! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Quasi-Monte Carlo helpers and engines corresponding to the first core slice
! of scipy.stats.qmc. Samples use SciPy's row-major mathematical convention:
! shape (n,d), with one point per row.

module scifort_qmc
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_ppf
    use scifort_linalg, only : cholesky_lower, symmetric_eigen_jacobi, linalg_status_success
    use scifort_random_variates, only : normal_rvs
    use scifort_random, only : rng_state, rng_uniform
    use scifort_qmc_sobol_data, only : sobol_maxdim, sobol_poly_value, sobol_vinit_value
    implicit none
    private

    integer, parameter, public :: qmc_status_success = 0
    integer, parameter, public :: qmc_status_invalid_shape = 1
    integer, parameter, public :: qmc_status_invalid_parameter = 2
    integer, parameter, public :: qmc_status_out_of_bounds = 3

    type, public :: qmc_halton
        integer :: d = 0 !! dimension of each generated point
        integer :: num_generated = 0 !! number of sequence points already consumed
        logical :: scramble = .false. !! whether Owen-style digit permutations are active
        integer, allocatable :: base(:) !! prime base associated with each coordinate
        integer, allocatable, private :: permutation_count(:)
        integer, allocatable, private :: permutations(:, :, :)
    contains
        procedure, public :: random => qmc_halton_random
        procedure, public :: fast_forward => qmc_halton_fast_forward
        procedure, public :: reset => qmc_halton_reset
    end type qmc_halton

    type, public :: qmc_latin_hypercube
        integer :: d = 0 !! dimension of each generated point
        integer :: num_generated = 0 !! cumulative number of points returned
        integer :: strength = 1 !! orthogonal-array strength, one or two
        logical :: scramble = .true. !! randomize positions within one-dimensional strata
    contains
        procedure, public :: random => qmc_latin_hypercube_random
        procedure, public :: reset => qmc_latin_hypercube_reset
    end type qmc_latin_hypercube

    type, public :: qmc_sobol
        integer :: d = 0 !! dimension of each generated point, at most 21201
        integer :: bits = 30 !! number of direction-number bits, from one through 64
        integer(int64) :: num_generated = 0_int64 !! number of sequence points already consumed
        logical :: scramble = .true. !! whether LMS plus digital-shift scrambling is active
        integer(int64), allocatable, private :: sv(:, :) !! Sobol direction numbers by dimension and bit
        integer(int64), allocatable, private :: shift(:) !! digital shift for each coordinate
        integer(int64), allocatable, private :: quasi(:) !! current integer-valued Sobol state
        real(dp), private :: scale_factor = 0.0_dp !! conversion factor from integer state to [0,1)
    contains
        procedure, public :: random => qmc_sobol_random
        procedure, public :: random_base2 => qmc_sobol_random_base2
        procedure, public :: fast_forward => qmc_sobol_fast_forward
        procedure, public :: reset => qmc_sobol_reset
    end type qmc_sobol

    type, public :: qmc_multinomial
        integer :: n_categories = 0 !! number of multinomial categories
        integer :: n_trials = 0 !! number of categorical draws per QMC sample
        real(dp), allocatable :: pvals(:) !! category probabilities
        real(dp), allocatable, private :: cumulative(:) !! cumulative probabilities for categorization
        type(qmc_sobol), private :: engine !! one-dimensional default Sobol engine
    contains
        procedure, public :: random => qmc_multinomial_random
        procedure, public :: reset => qmc_multinomial_reset
    end type qmc_multinomial

    type, public :: qmc_multivariate_normal
        integer :: d = 0 !! dimension of each multivariate-normal draw
        integer :: engine_d = 0 !! dimension consumed from the QMC engine
        logical :: inv_transform = .true. !! use inverse-normal transform instead of Box-Muller
        logical :: has_root = .false. !! whether a non-identity covariance root is stored
        real(dp), allocatable :: mean(:) !! mean vector
        real(dp), allocatable, private :: covariance_root(:, :) !! right-multiplication root of covariance
        type(qmc_sobol), private :: engine !! default Sobol engine
    contains
        procedure, public :: random => qmc_multivariate_normal_random
        procedure, public :: reset => qmc_multivariate_normal_reset
    end type qmc_multivariate_normal

    type, public :: qmc_poisson_disk
        integer :: d = 0 !! sampling dimension
        integer :: num_generated = 0 !! cumulative number of accepted points returned
        integer :: ncandidates = 30 !! candidates proposed around each active point
        logical :: surface = .false. !! sample candidate directions on a shell rather than in a ball
        real(dp) :: radius = 0.05_dp !! minimum accepted pairwise distance
        real(dp) :: radius_squared = 0.0025_dp !! squared minimum distance
        real(dp) :: radius_factor = 2.0_dp !! proposal-ball radius relative to minimum distance
        real(dp), allocatable :: lower(:) !! lower sampling bounds
        real(dp), allocatable :: upper(:) !! upper sampling bounds
        real(dp), allocatable, private :: pool(:, :) !! active Bridson points
        real(dp), allocatable, private :: accepted(:, :) !! all accepted points for neighbor checks
        integer, private :: pool_size = 0 !! active points currently stored
        integer, private :: accepted_size = 0 !! total accepted points currently stored
    contains
        procedure, public :: random => qmc_poisson_disk_random
        procedure, public :: fill_space => qmc_poisson_disk_fill_space
        procedure, public :: reset => qmc_poisson_disk_reset
    end type qmc_poisson_disk

    public :: qmc_discrepancy
    public :: qmc_geometric_discrepancy
    public :: qmc_halton_init
    public :: qmc_latin_hypercube_init
    public :: qmc_scale
    public :: qmc_sobol_init
    public :: qmc_multinomial_init
    public :: qmc_multivariate_normal_init
    public :: qmc_poisson_disk_init
    public :: qmc_to_integers
    public :: qmc_update_discrepancy
    public :: qmc_van_der_corput

contains

    subroutine qmc_scale(sample, lower, upper, scaled, status, reverse)
        real(dp), intent(in) :: sample(:, :) !! input sample with shape (n,d)
        real(dp), intent(in) :: lower(:) !! lower bound for each coordinate
        real(dp), intent(in) :: upper(:) !! upper bound for each coordinate
        real(dp), intent(out) :: scaled(:, :) !! transformed sample with shape (n,d)
        integer, intent(out) :: status !! qmc_status_* result
        logical, intent(in), optional :: reverse !! map bounds to unit cube when true

        integer :: j
        logical :: do_reverse

        status = qmc_status_success
        scaled = 0.0_dp
        do_reverse = .false.
        if (present(reverse)) do_reverse = reverse
        if (size(sample,2) /= size(lower) .or. size(lower) /= size(upper) .or. &
                any(shape(scaled) /= shape(sample))) then
            status = qmc_status_invalid_shape
            return
        end if
        if (.not. all(ieee_is_finite(sample)) .or. .not. all(ieee_is_finite(lower)) .or. &
                .not. all(ieee_is_finite(upper)) .or. any(upper <= lower)) then
            status = qmc_status_invalid_parameter
            return
        end if
        if (.not. do_reverse) then
            if (any(sample < 0.0_dp) .or. any(sample > 1.0_dp)) then
                status = qmc_status_out_of_bounds
                return
            end if
            do j = 1, size(sample,2)
                scaled(:,j) = sample(:,j) * (upper(j) - lower(j)) + lower(j)
            end do
        else
            do j = 1, size(sample,2)
                if (any(sample(:,j) < lower(j)) .or. any(sample(:,j) > upper(j))) then
                    status = qmc_status_out_of_bounds
                    return
                end if
                scaled(:,j) = (sample(:,j) - lower(j)) / (upper(j) - lower(j))
            end do
        end if
    end subroutine qmc_scale

    function qmc_discrepancy(sample, method, iterative, status) result(value)
        real(dp), intent(in) :: sample(:, :) !! unit-hypercube sample with shape (n,d)
        character(len=*), intent(in), optional :: method !! 'CD', 'WD', 'MD', or 'L2-star'
        logical, intent(in), optional :: iterative !! use n+1 denominators for update workflow
        integer, intent(out), optional :: status !! qmc_status_* result
        real(dp) :: value

        character(len=:), allocatable :: selected
        integer :: i
        integer :: j
        integer :: k
        integer :: nden
        logical :: iterative_local
        real(dp) :: a
        real(dp) :: b
        real(dp) :: disc1
        real(dp) :: disc2
        real(dp) :: prod
        real(dp) :: term

        if (present(status)) status = qmc_status_success
        if (.not. valid_unit_sample(sample) .or. size(sample,1) < 1 .or. size(sample,2) < 1) then
            value = quiet_nan(0.0_dp)
            if (present(status)) status = qmc_status_invalid_parameter
            return
        end if
        selected = 'cd'
        if (present(method)) selected = lower_ascii(trim(method))
        iterative_local = .false.
        if (present(iterative)) iterative_local = iterative
        nden = size(sample,1)
        if (iterative_local) nden = nden + 1

        select case (selected)
        case ('cd')
            disc1 = 0.0_dp
            do i = 1, size(sample,1)
                prod = 1.0_dp
                do k = 1, size(sample,2)
                    a = abs(sample(i,k) - 0.5_dp)
                    prod = prod * (1.0_dp + 0.5_dp * a - 0.5_dp * a * a)
                end do
                disc1 = disc1 + prod
            end do
            disc2 = 0.0_dp
            do i = 1, size(sample,1)
                do j = 1, size(sample,1)
                    prod = 1.0_dp
                    do k = 1, size(sample,2)
                        a = abs(sample(i,k) - 0.5_dp)
                        b = abs(sample(j,k) - 0.5_dp)
                        prod = prod * (1.0_dp + 0.5_dp * a + 0.5_dp * b - &
                            0.5_dp * abs(sample(i,k) - sample(j,k)))
                    end do
                    disc2 = disc2 + prod
                end do
            end do
            value = (13.0_dp / 12.0_dp) ** size(sample,2) - &
                2.0_dp * disc1 / real(nden,dp) + disc2 / real(nden,dp) ** 2
        case ('wd')
            disc2 = 0.0_dp
            do i = 1, size(sample,1)
                do j = 1, size(sample,1)
                    prod = 1.0_dp
                    do k = 1, size(sample,2)
                        term = abs(sample(i,k) - sample(j,k))
                        prod = prod * (1.5_dp - term + term * term)
                    end do
                    disc2 = disc2 + prod
                end do
            end do
            value = -(4.0_dp / 3.0_dp) ** size(sample,2) + &
                disc2 / real(nden,dp) ** 2
        case ('md')
            disc1 = 0.0_dp
            do i = 1, size(sample,1)
                prod = 1.0_dp
                do k = 1, size(sample,2)
                    a = abs(sample(i,k) - 0.5_dp)
                    prod = prod * (5.0_dp / 3.0_dp - 0.25_dp * a - 0.25_dp * a * a)
                end do
                disc1 = disc1 + prod
            end do
            disc2 = 0.0_dp
            do i = 1, size(sample,1)
                do j = 1, size(sample,1)
                    prod = 1.0_dp
                    do k = 1, size(sample,2)
                        a = abs(sample(i,k) - 0.5_dp)
                        b = abs(sample(j,k) - 0.5_dp)
                        term = abs(sample(i,k) - sample(j,k))
                        prod = prod * (15.0_dp / 8.0_dp - 0.25_dp * a - 0.25_dp * b - &
                            0.75_dp * term + 0.5_dp * term * term)
                    end do
                    disc2 = disc2 + prod
                end do
            end do
            value = (19.0_dp / 12.0_dp) ** size(sample,2) - &
                2.0_dp * disc1 / real(nden,dp) + disc2 / real(nden,dp) ** 2
        case ('l2-star', 'l2star')
            disc1 = 0.0_dp
            do i = 1, size(sample,1)
                prod = 1.0_dp
                do k = 1, size(sample,2)
                    prod = prod * (1.0_dp - sample(i,k) * sample(i,k))
                end do
                disc1 = disc1 + prod
            end do
            disc2 = 0.0_dp
            do i = 1, size(sample,1)
                do j = 1, size(sample,1)
                    prod = 1.0_dp
                    do k = 1, size(sample,2)
                        prod = prod * (1.0_dp - max(sample(i,k), sample(j,k)))
                    end do
                    disc2 = disc2 + prod
                end do
            end do
            term = 3.0_dp ** (-size(sample,2)) - 2.0_dp ** (1 - size(sample,2)) * &
                disc1 / real(nden,dp) + disc2 / real(nden,dp) ** 2
            value = sqrt(max(0.0_dp, term))
        case default
            value = quiet_nan(0.0_dp)
            if (present(status)) status = qmc_status_invalid_parameter
        end select
    end function qmc_discrepancy

    function qmc_update_discrepancy(x_new, sample, initial_disc, status) result(value)
        real(dp), intent(in) :: x_new(:) !! candidate point with d coordinates
        real(dp), intent(in) :: sample(:, :) !! existing unit-hypercube sample with shape (n,d)
        real(dp), intent(in) :: initial_disc !! centered discrepancy computed with iterative=true
        integer, intent(out), optional :: status !! qmc_status_* result
        real(dp) :: value

        integer :: i
        integer :: j
        integer :: n
        real(dp) :: disc1
        real(dp) :: disc2
        real(dp) :: disc3
        real(dp) :: prod
        real(dp), allocatable :: abs_center(:)

        if (present(status)) status = qmc_status_success
        if (size(x_new) /= size(sample,2) .or. .not. valid_unit_sample(sample) .or. &
                any(x_new < 0.0_dp) .or. any(x_new > 1.0_dp) .or. &
                .not. all(ieee_is_finite(x_new)) .or. .not. ieee_is_finite(initial_disc)) then
            value = quiet_nan(0.0_dp)
            if (present(status)) status = qmc_status_invalid_parameter
            return
        end if
        n = size(sample,1) + 1
        allocate(abs_center(size(x_new)))
        abs_center = abs(x_new - 0.5_dp)

        prod = 1.0_dp
        do j = 1, size(x_new)
            prod = prod * (1.0_dp + 0.5_dp * abs_center(j) - &
                0.5_dp * abs_center(j) * abs_center(j))
        end do
        disc1 = -2.0_dp * prod / real(n,dp)

        disc2 = 0.0_dp
        do i = 1, size(sample,1)
            prod = 1.0_dp
            do j = 1, size(x_new)
                prod = prod * (1.0_dp + 0.5_dp * abs_center(j) + &
                    0.5_dp * abs(sample(i,j) - 0.5_dp) - &
                    0.5_dp * abs(x_new(j) - sample(i,j)))
            end do
            disc2 = disc2 + prod
        end do
        disc2 = 2.0_dp * disc2 / real(n,dp) ** 2

        prod = 1.0_dp
        do j = 1, size(x_new)
            prod = prod * (1.0_dp + abs_center(j))
        end do
        disc3 = prod / real(n,dp) ** 2
        value = initial_disc + disc1 + disc2 + disc3
    end function qmc_update_discrepancy

    function qmc_geometric_discrepancy(sample, method, metric, status) result(value)
        real(dp), intent(in) :: sample(:, :) !! unit-hypercube sample with shape (n,d)
        character(len=*), intent(in), optional :: method !! 'mindist' or 'mst'
        character(len=*), intent(in), optional :: metric !! 'euclidean' or 'cityblock'
        integer, intent(out), optional :: status !! qmc_status_* result
        real(dp) :: value

        character(len=:), allocatable :: metric_local
        character(len=:), allocatable :: method_local
        integer :: i
        integer :: j
        integer :: k
        integer :: next
        integer :: n
        real(dp) :: dist
        real(dp) :: best
        real(dp) :: total
        real(dp), allocatable :: min_edge(:)
        logical :: found_positive
        logical, allocatable :: used(:)

        if (present(status)) status = qmc_status_success
        if (.not. valid_unit_sample(sample) .or. size(sample,1) < 2 .or. size(sample,2) < 1) then
            value = quiet_nan(0.0_dp)
            if (present(status)) status = qmc_status_invalid_parameter
            return
        end if
        method_local = 'mindist'
        if (present(method)) method_local = lower_ascii(trim(method))
        metric_local = 'euclidean'
        if (present(metric)) metric_local = lower_ascii(trim(metric))
        if (metric_local /= 'euclidean' .and. metric_local /= 'cityblock') then
            value = quiet_nan(0.0_dp)
            if (present(status)) status = qmc_status_invalid_parameter
            return
        end if

        if (method_local == 'mindist') then
            value = huge(1.0_dp)
            found_positive = .false.
            do i = 1, size(sample,1) - 1
                do j = i + 1, size(sample,1)
                    dist = point_distance(sample(i,:), sample(j,:), metric_local)
                    if (dist > 0.0_dp) then
                        value = min(value, dist)
                        found_positive = .true.
                    end if
                end do
            end do
            if (.not. found_positive) then
                value = quiet_nan(0.0_dp)
                if (present(status)) status = qmc_status_invalid_parameter
            end if
            return
        end if
        if (method_local /= 'mst') then
            value = quiet_nan(0.0_dp)
            if (present(status)) status = qmc_status_invalid_parameter
            return
        end if

        n = size(sample,1)
        allocate(min_edge(n), used(n))
        used = .false.
        min_edge = huge(1.0_dp)
        min_edge(1) = 0.0_dp
        total = 0.0_dp
        do k = 1, n
            next = 0
            best = huge(1.0_dp)
            do i = 1, n
                if (.not. used(i) .and. min_edge(i) < best) then
                    best = min_edge(i)
                    next = i
                end if
            end do
            if (next == 0) then
                value = quiet_nan(0.0_dp)
                if (present(status)) status = qmc_status_invalid_parameter
                return
            end if
            used(next) = .true.
            total = total + best
            do j = 1, n
                if (.not. used(j)) then
                    dist = point_distance(sample(next,:), sample(j,:), metric_local)
                    if (dist < min_edge(j)) min_edge(j) = dist
                end if
            end do
        end do
        value = total / real(n - 1,dp)
    end function qmc_geometric_discrepancy

    subroutine qmc_van_der_corput(n, base, sequence, status, start_index, scramble, state)
        integer, intent(in) :: n !! number of sequence values to generate
        integer, intent(in) :: base !! integer radix, at least two
        real(dp), intent(out) :: sequence(:) !! generated values in [0,1)
        integer, intent(out) :: status !! qmc_status_* result
        integer, intent(in), optional :: start_index !! zero-based initial sequence index
        logical, intent(in), optional :: scramble !! apply independent digit permutations when true
        type(rng_state), intent(inout), optional :: state !! explicit state used to create scrambling permutations

        integer :: count
        integer :: digit
        integer :: i
        integer :: index0
        integer :: j
        integer :: k
        integer :: q
        integer :: start
        integer, allocatable :: perms(:, :)
        logical :: do_scramble
        real(dp) :: place

        status = qmc_status_success
        sequence = 0.0_dp
        start = 0
        if (present(start_index)) start = start_index
        do_scramble = .false.
        if (present(scramble)) do_scramble = scramble
        if (n < 0 .or. base < 2 .or. start < 0 .or. size(sequence) /= n) then
            status = qmc_status_invalid_parameter
            return
        end if
        if (do_scramble .and. .not. present(state)) then
            status = qmc_status_invalid_parameter
            return
        end if

        if (do_scramble) then
            count = ceiling(54.0_dp / log(real(base,dp)) * log(2.0_dp)) - 1
            allocate(perms(count,0:base-1))
            do j = 1, count
                do k = 0, base - 1
                    perms(j,k) = k
                end do
                call shuffle_zero_based(state, perms(j,:))
            end do
        end if

        do i = 1, n
            index0 = start + i - 1
            q = index0
            place = 1.0_dp / real(base,dp)
            if (do_scramble) then
                do j = 1, count
                    digit = modulo(q, base)
                    sequence(i) = sequence(i) + real(perms(j,digit),dp) * place
                    place = place / real(base,dp)
                    q = q / base
                end do
            else
                do while (q > 0)
                    digit = modulo(q, base)
                    sequence(i) = sequence(i) + real(digit,dp) * place
                    place = place / real(base,dp)
                    q = q / base
                end do
            end if
        end do
    end subroutine qmc_van_der_corput

    subroutine qmc_halton_init(engine, d, status, scramble, state)
        type(qmc_halton), intent(out) :: engine !! initialized Halton engine
        integer, intent(in) :: d !! sequence dimension, at least one
        integer, intent(out) :: status !! qmc_status_* result
        logical, intent(in), optional :: scramble !! activate Owen-style scrambling
        type(rng_state), intent(inout), optional :: state !! explicit state used for scrambling permutations

        integer :: count
        integer :: i
        integer :: j
        integer :: k
        integer :: max_base
        integer :: max_count

        status = qmc_status_success
        if (d < 1) then
            status = qmc_status_invalid_parameter
            return
        end if
        engine%d = d
        engine%num_generated = 0
        engine%scramble = .false.
        if (present(scramble)) engine%scramble = scramble
        if (engine%scramble .and. .not. present(state)) then
            status = qmc_status_invalid_parameter
            return
        end if
        allocate(engine%base(d), engine%permutation_count(d))
        call first_primes(d, engine%base)
        engine%permutation_count = 0
        if (.not. engine%scramble) return

        max_base = maxval(engine%base)
        max_count = 0
        do i = 1, d
            count = ceiling(54.0_dp / log(real(engine%base(i),dp)) * log(2.0_dp)) - 1
            engine%permutation_count(i) = count
            max_count = max(max_count, count)
        end do
        allocate(engine%permutations(max_count,0:max_base-1,d))
        engine%permutations = -1
        do i = 1, d
            do j = 1, engine%permutation_count(i)
                do k = 0, engine%base(i) - 1
                    engine%permutations(j,k,i) = k
                end do
                call shuffle_zero_based(state, engine%permutations(j,0:engine%base(i)-1,i))
            end do
        end do
    end subroutine qmc_halton_init

    subroutine qmc_halton_random(self, n, sample, status)
        class(qmc_halton), intent(inout) :: self !! Halton engine advanced by n points
        integer, intent(in) :: n !! number of points to generate
        real(dp), intent(out) :: sample(:, :) !! generated sample with shape (n,d)
        integer, intent(out) :: status !! qmc_status_* result

        integer :: digit
        integer :: i
        integer :: j
        integer :: k
        integer :: q
        real(dp) :: place

        status = qmc_status_success
        sample = 0.0_dp
        if (self%d < 1 .or. n < 0 .or. size(sample,1) /= n .or. size(sample,2) /= self%d) then
            status = qmc_status_invalid_shape
            return
        end if
        do i = 1, n
            do j = 1, self%d
                q = self%num_generated + i - 1
                place = 1.0_dp / real(self%base(j),dp)
                if (self%scramble) then
                    do k = 1, self%permutation_count(j)
                        digit = modulo(q, self%base(j))
                        sample(i,j) = sample(i,j) + &
                            real(self%permutations(k,digit,j),dp) * place
                        place = place / real(self%base(j),dp)
                        q = q / self%base(j)
                    end do
                else
                    do while (q > 0)
                        digit = modulo(q, self%base(j))
                        sample(i,j) = sample(i,j) + real(digit,dp) * place
                        place = place / real(self%base(j),dp)
                        q = q / self%base(j)
                    end do
                end if
            end do
        end do
        self%num_generated = self%num_generated + n
    end subroutine qmc_halton_random

    subroutine qmc_halton_fast_forward(self, n, status)
        class(qmc_halton), intent(inout) :: self !! Halton engine whose index is advanced
        integer, intent(in) :: n !! number of sequence positions to skip
        integer, intent(out), optional :: status !! qmc_status_* result

        if (present(status)) status = qmc_status_success
        if (n < 0 .or. self%d < 1) then
            if (present(status)) status = qmc_status_invalid_parameter
            return
        end if
        self%num_generated = self%num_generated + n
    end subroutine qmc_halton_fast_forward

    subroutine qmc_halton_reset(self)
        class(qmc_halton), intent(inout) :: self !! Halton engine reset to its first point

        self%num_generated = 0
    end subroutine qmc_halton_reset

    subroutine qmc_sobol_init(engine, d, status, scramble, bits, state)
        type(qmc_sobol), intent(out) :: engine !! initialized Sobol engine
        integer, intent(in) :: d !! sequence dimension, from one through 21201
        integer, intent(out) :: status !! qmc_status_* result
        logical, intent(in), optional :: scramble !! activate LMS plus digital-shift scrambling
        integer, intent(in), optional :: bits !! number of sequence bits, from one through 64
        type(rng_state), intent(inout), optional :: state !! explicit RNG state used for scrambling

        integer :: selected_bits

        status = qmc_status_success
        selected_bits = 30
        if (present(bits)) selected_bits = bits
        if (d < 1 .or. d > sobol_maxdim .or. selected_bits < 1 .or. selected_bits > 64) then
            status = qmc_status_invalid_parameter
            return
        end if

        engine%d = d
        engine%bits = selected_bits
        engine%num_generated = 0_int64
        engine%scramble = .true.
        if (present(scramble)) engine%scramble = scramble
        if (engine%scramble .and. .not. present(state)) then
            status = qmc_status_invalid_parameter
            return
        end if

        allocate(engine%sv(d,selected_bits), engine%shift(d), engine%quasi(d))
        call initialize_sobol_directions(engine%sv, d, selected_bits)
        engine%shift = 0_int64
        if (engine%scramble) call scramble_sobol(state, engine%sv, engine%shift)
        engine%quasi = engine%shift
        engine%scale_factor = scale(1.0_dp, -selected_bits)
    end subroutine qmc_sobol_init

    subroutine qmc_sobol_random(self, n, sample, status)
        class(qmc_sobol), intent(inout) :: self !! Sobol engine advanced by n points
        integer, intent(in) :: n !! number of points to generate
        real(dp), intent(out) :: sample(:, :) !! generated sample with shape (n,d)
        integer, intent(out) :: status !! qmc_status_* result

        integer :: i
        integer :: j
        integer :: l
        integer(int64) :: previous_index

        status = qmc_status_success
        sample = 0.0_dp
        if (self%d < 1 .or. n < 0 .or. size(sample,1) /= n .or. size(sample,2) /= self%d) then
            status = qmc_status_invalid_shape
            return
        end if
        if (int(n,int64) > huge(0_int64) - self%num_generated) then
            status = qmc_status_invalid_parameter
            return
        end if
        if (self%bits < 63) then
            if (self%num_generated + int(n,int64) > ishft(1_int64,self%bits)) then
                status = qmc_status_invalid_parameter
                return
            end if
        end if

        do i = 1, n
            if (self%num_generated == 0_int64) then
                do j = 1, self%d
                    sample(i,j) = uint64_to_unit(self%quasi(j), self%scale_factor)
                end do
                self%num_generated = 1_int64
            else
                previous_index = self%num_generated - 1_int64
                l = low_zero_bit(previous_index)
                if (l > self%bits) then
                    status = qmc_status_invalid_parameter
                    return
                end if
                do j = 1, self%d
                    self%quasi(j) = ieor(self%quasi(j), self%sv(j,l))
                    sample(i,j) = uint64_to_unit(self%quasi(j), self%scale_factor)
                end do
                self%num_generated = self%num_generated + 1_int64
            end if
        end do
    end subroutine qmc_sobol_random

    subroutine qmc_sobol_random_base2(self, m, sample, status)
        class(qmc_sobol), intent(inout) :: self !! Sobol engine advanced by 2**m points
        integer, intent(in) :: m !! base-two logarithm of number of requested points
        real(dp), intent(out) :: sample(:, :) !! generated sample with shape (2**m,d)
        integer, intent(out) :: status !! qmc_status_* result

        integer :: n
        integer(int64) :: total

        status = qmc_status_success
        if (m < 0 .or. m > 30) then
            status = qmc_status_invalid_parameter
            sample = 0.0_dp
            return
        end if
        n = ishft(1,m)
        if (size(sample,1) /= n .or. size(sample,2) /= self%d) then
            status = qmc_status_invalid_shape
            sample = 0.0_dp
            return
        end if
        total = self%num_generated + int(n,int64)
        if (total <= 0_int64 .or. popcnt(total) /= 1) then
            status = qmc_status_invalid_parameter
            sample = 0.0_dp
            return
        end if
        call self%random(n, sample, status)
    end subroutine qmc_sobol_random_base2

    subroutine qmc_sobol_fast_forward(self, n, status)
        class(qmc_sobol), intent(inout) :: self !! Sobol engine whose index is advanced
        integer, intent(in) :: n !! number of sequence positions to skip
        integer, intent(out), optional :: status !! qmc_status_* result

        integer :: i
        integer :: j
        integer :: l
        integer(int64) :: previous_index

        if (present(status)) status = qmc_status_success
        if (n < 0 .or. self%d < 1) then
            if (present(status)) status = qmc_status_invalid_parameter
            return
        end if
        if (int(n,int64) > huge(0_int64) - self%num_generated) then
            if (present(status)) status = qmc_status_invalid_parameter
            return
        end if
        if (self%bits < 63) then
            if (self%num_generated + int(n,int64) > ishft(1_int64,self%bits)) then
                if (present(status)) status = qmc_status_invalid_parameter
                return
            end if
        end if

        do i = 1, n
            if (self%num_generated == 0_int64) then
                self%num_generated = 1_int64
            else
                previous_index = self%num_generated - 1_int64
                l = low_zero_bit(previous_index)
                if (l > self%bits) then
                    if (present(status)) status = qmc_status_invalid_parameter
                    return
                end if
                do j = 1, self%d
                    self%quasi(j) = ieor(self%quasi(j), self%sv(j,l))
                end do
                self%num_generated = self%num_generated + 1_int64
            end if
        end do
    end subroutine qmc_sobol_fast_forward

    subroutine qmc_sobol_reset(self)
        class(qmc_sobol), intent(inout) :: self !! Sobol engine reset to its first point

        self%num_generated = 0_int64
        if (allocated(self%quasi) .and. allocated(self%shift)) self%quasi = self%shift
    end subroutine qmc_sobol_reset

    subroutine initialize_sobol_directions(sv, dim, bits)
        integer, intent(in) :: dim !! number of dimensions
        integer, intent(in) :: bits !! number of bits
        integer(int64), intent(out) :: sv(dim,bits) !! initialized direction-number matrix

        integer :: d
        integer :: j
        integer :: k
        integer :: m
        integer(int64) :: newv
        integer(int64) :: p
        integer(int64) :: pow2

        sv = 0_int64
        sv(1,:) = 1_int64
        do d = 2, dim
            p = int(sobol_poly_value(d),int64)
            m = bit_size(p) - leadz(p) - 1
            do j = 1, min(m,bits)
                sv(d,j) = int(sobol_vinit_value(d,j),int64)
            end do
            do j = m + 1, bits
                newv = sv(d,j-m)
                pow2 = 1_int64
                do k = 1, m
                    pow2 = ishft(pow2,1)
                    if (btest(p,m-k)) then
                        newv = ieor(newv, pow2 * sv(d,j-k))
                    end if
                end do
                sv(d,j) = newv
            end do
        end do
        do j = 1, bits
            do d = 1, dim
                sv(d,j) = ishft(sv(d,j), bits-j)
            end do
        end do
    end subroutine initialize_sobol_directions

    subroutine scramble_sobol(state, sv, shift_bits)
        type(rng_state), intent(inout) :: state !! explicit state advanced by scrambling draws
        integer(int64), intent(inout) :: sv(:, :) !! direction numbers scrambled in place
        integer(int64), intent(out) :: shift_bits(:) !! digital shift generated for each dimension

        integer :: bits
        integer :: d
        integer :: j
        integer :: k
        integer :: p
        integer :: parity
        integer(int64) :: original
        integer(int64) :: scrambled
        logical, allocatable :: ltm(:, :)

        bits = size(sv,2)
        allocate(ltm(bits,bits))
        do d = 1, size(sv,1)
            shift_bits(d) = 0_int64
            do k = 0, bits - 1
                if (rng_uniform(state) >= 0.5_dp) shift_bits(d) = ibset(shift_bits(d),k)
            end do

            ltm = .false.
            do p = 1, bits
                ltm(p,p) = .true.
                do k = 1, p - 1
                    ltm(p,k) = rng_uniform(state) >= 0.5_dp
                end do
            end do
            do j = 1, bits
                original = sv(d,j)
                scrambled = 0_int64
                do p = 1, bits
                    parity = 0
                    do k = 1, p
                        if (ltm(p,k) .and. btest(original,bits-k)) parity = 1 - parity
                    end do
                    if (parity == 1) scrambled = ibset(scrambled,bits-p)
                end do
                sv(d,j) = scrambled
            end do
        end do
    end subroutine scramble_sobol

    pure function low_zero_bit(x) result(position)
        integer(int64), intent(in) :: x !! nonnegative integer whose rightmost zero bit is sought
        integer :: position

        integer(int64) :: work

        work = x
        position = 1
        do while (btest(work,0))
            work = ishft(work,-1)
            position = position + 1
        end do
    end function low_zero_bit

    pure function uint64_to_unit(x, factor) result(value)
        integer(int64), intent(in) :: x !! unsigned 64-bit pattern stored in signed int64
        real(dp), intent(in) :: factor !! scale factor 2**(-bits)
        real(dp) :: value

        real(dp), parameter :: two63 = 9223372036854775808.0_dp

        if (btest(x,63)) then
            value = (real(ibclr(x,63),dp) + two63) * factor
        else
            value = real(x,dp) * factor
        end if
    end function uint64_to_unit

    subroutine qmc_multinomial_init(dist, pvals, n_trials, state, status, scramble, bits)
        type(qmc_multinomial), intent(out) :: dist !! initialized multinomial QMC sampler
        real(dp), intent(in) :: pvals(:) !! nonnegative category probabilities summing to one
        integer, intent(in) :: n_trials !! categorical trials per returned sample
        type(rng_state), intent(inout), optional :: state !! explicit state used for Sobol scrambling
        integer, intent(out) :: status !! qmc_status_* result
        logical, intent(in), optional :: scramble !! whether to scramble the default Sobol engine
        integer, intent(in), optional :: bits !! Sobol direction-number bits

        integer :: i
        logical :: do_scramble
        real(dp) :: total

        status = qmc_status_success
        if (size(pvals) < 1 .or. n_trials < 0 .or. .not. all(ieee_is_finite(pvals)) .or. any(pvals < 0.0_dp)) then
            status = qmc_status_invalid_parameter
            return
        end if
        total = sum(pvals)
        if (abs(total - 1.0_dp) > 1.0e-8_dp + 1.0e-5_dp) then
            status = qmc_status_invalid_parameter
            return
        end if

        dist%n_categories = size(pvals)
        dist%n_trials = n_trials
        allocate(dist%pvals(size(pvals)), dist%cumulative(size(pvals)))
        dist%pvals = pvals
        dist%cumulative(1) = pvals(1)
        do i = 2, size(pvals)
            dist%cumulative(i) = dist%cumulative(i-1) + pvals(i)
        end do
        dist%cumulative(size(pvals)) = 1.0_dp

        do_scramble = .true.
        if (present(scramble)) do_scramble = scramble
        if (present(bits)) then
            if (present(state)) then
                call qmc_sobol_init(dist%engine, 1, status, scramble=do_scramble, bits=bits, state=state)
            else
                call qmc_sobol_init(dist%engine, 1, status, scramble=do_scramble, bits=bits)
            end if
        else
            if (present(state)) then
                call qmc_sobol_init(dist%engine, 1, status, scramble=do_scramble, state=state)
            else
                call qmc_sobol_init(dist%engine, 1, status, scramble=do_scramble)
            end if
        end if
    end subroutine qmc_multinomial_init

    subroutine qmc_multinomial_random(self, n, sample, status)
        class(qmc_multinomial), intent(inout) :: self !! sampler advanced by n multinomial draws
        integer, intent(in) :: n !! number of QMC multinomial samples
        integer, intent(out) :: sample(:, :) !! category counts with shape (n,k)
        integer, intent(out) :: status !! qmc_status_* result

        integer :: category
        integer :: i
        integer :: j
        real(dp), allocatable :: draws(:, :)

        status = qmc_status_success
        sample = 0
        if (n < 0 .or. self%n_categories < 1 .or. size(sample,1) /= n .or. &
                size(sample,2) /= self%n_categories) then
            status = qmc_status_invalid_shape
            return
        end if
        allocate(draws(self%n_trials,1))
        do i = 1, n
            call self%engine%random(self%n_trials, draws, status)
            if (status /= qmc_status_success) return
            do j = 1, self%n_trials
                category = 1
                do while (category < self%n_categories .and. draws(j,1) > self%cumulative(category))
                    category = category + 1
                end do
                sample(i,category) = sample(i,category) + 1
            end do
        end do
    end subroutine qmc_multinomial_random

    subroutine qmc_multinomial_reset(self)
        class(qmc_multinomial), intent(inout) :: self !! sampler reset to its initial QMC point

        call self%engine%reset()
    end subroutine qmc_multinomial_reset

    subroutine qmc_multivariate_normal_init(dist, mean, status, state, cov, cov_root, &
            inv_transform, scramble, bits)
        type(qmc_multivariate_normal), intent(out) :: dist !! initialized multivariate-normal QMC sampler
        real(dp), intent(in) :: mean(:) !! mean vector
        integer, intent(out) :: status !! qmc_status_* result
        type(rng_state), intent(inout), optional :: state !! explicit state used for Sobol scrambling
        real(dp), intent(in), optional :: cov(:, :) !! symmetric positive-semidefinite covariance matrix
        real(dp), intent(in), optional :: cov_root(:, :) !! right-multiplication covariance root
        logical, intent(in), optional :: inv_transform !! inverse-normal transform when true
        logical, intent(in), optional :: scramble !! whether to scramble the default Sobol engine
        integer, intent(in), optional :: bits !! Sobol direction-number bits

        integer :: d
        integer :: i
        integer :: info
        integer :: j
        logical :: do_scramble
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: lower(:, :)
        real(dp) :: symmetry_tol

        status = qmc_status_success
        d = size(mean)
        if (d < 1 .or. .not. all(ieee_is_finite(mean)) .or. (present(cov) .and. present(cov_root))) then
            status = qmc_status_invalid_parameter
            return
        end if

        dist%d = d
        dist%inv_transform = .true.
        if (present(inv_transform)) dist%inv_transform = inv_transform
        dist%engine_d = d
        if (.not. dist%inv_transform .and. modulo(d,2) /= 0) dist%engine_d = d + 1
        allocate(dist%mean(d))
        dist%mean = mean
        dist%has_root = .false.

        if (present(cov)) then
            if (size(cov,1) /= d .or. size(cov,2) /= d .or. .not. all(ieee_is_finite(cov))) then
                status = qmc_status_invalid_shape
                return
            end if
            symmetry_tol = maxval(abs(cov - transpose(cov)) - &
                (1.0e-8_dp + 1.0e-5_dp * abs(transpose(cov))))
            if (symmetry_tol > 0.0_dp) then
                status = qmc_status_invalid_parameter
                return
            end if
            allocate(dist%covariance_root(d,d), lower(d,d))
            call cholesky_lower(cov, lower, info)
            if (info == linalg_status_success) then
                dist%covariance_root = transpose(lower)
            else
                allocate(eigenvalues(d), eigenvectors(d,d))
                call symmetric_eigen_jacobi(cov, eigenvalues, eigenvectors, info)
                if (info /= linalg_status_success .or. minval(eigenvalues) < -1.0e-8_dp) then
                    status = qmc_status_invalid_parameter
                    return
                end if
                do i = 1, d
                    do j = 1, d
                        dist%covariance_root(i,j) = sqrt(max(0.0_dp,eigenvalues(i))) * eigenvectors(j,i)
                    end do
                end do
            end if
            dist%has_root = .true.
        else if (present(cov_root)) then
            if (size(cov_root,1) /= d .or. size(cov_root,2) /= d .or. .not. all(ieee_is_finite(cov_root))) then
                status = qmc_status_invalid_shape
                return
            end if
            allocate(dist%covariance_root(d,d))
            dist%covariance_root = cov_root
            dist%has_root = .true.
        end if

        do_scramble = .true.
        if (present(scramble)) do_scramble = scramble
        if (present(bits)) then
            if (present(state)) then
                call qmc_sobol_init(dist%engine, dist%engine_d, status, scramble=do_scramble, bits=bits, state=state)
            else
                call qmc_sobol_init(dist%engine, dist%engine_d, status, scramble=do_scramble, bits=bits)
            end if
        else
            if (present(state)) then
                call qmc_sobol_init(dist%engine, dist%engine_d, status, scramble=do_scramble, state=state)
            else
                call qmc_sobol_init(dist%engine, dist%engine_d, status, scramble=do_scramble)
            end if
        end if
    end subroutine qmc_multivariate_normal_init

    subroutine qmc_multivariate_normal_random(self, n, sample, status)
        class(qmc_multivariate_normal), intent(inout) :: self !! sampler advanced by n QMC points
        integer, intent(in) :: n !! number of multivariate-normal samples
        real(dp), intent(out) :: sample(:, :) !! generated values with shape (n,d)
        integer, intent(out) :: status !! qmc_status_* result

        integer :: i
        integer :: j
        integer :: pair
        real(dp), allocatable :: base(:, :)
        real(dp), allocatable :: unit(:, :)
        real(dp) :: angle
        real(dp) :: p
        real(dp) :: radius
        real(dp), parameter :: clip_scale = 1.0_dp - 1.0e-10_dp
        real(dp), parameter :: two_pi = 6.2831853071795864769252867665590058_dp

        status = qmc_status_success
        sample = 0.0_dp
        if (n < 0 .or. self%d < 1 .or. size(sample,1) /= n .or. size(sample,2) /= self%d) then
            status = qmc_status_invalid_shape
            return
        end if
        allocate(unit(n,self%engine_d), base(n,self%d))
        call self%engine%random(n, unit, status)
        if (status /= qmc_status_success) return

        if (self%inv_transform) then
            do j = 1, self%d
                do i = 1, n
                    p = 0.5_dp + clip_scale * (unit(i,j) - 0.5_dp)
                    base(i,j) = normal_ppf(p, 0.0_dp, 1.0_dp)
                end do
            end do
        else
            base = 0.0_dp
            pair = 0
            do j = 1, self%engine_d, 2
                pair = pair + 1
                do i = 1, n
                    radius = sqrt(-2.0_dp * log(unit(i,j)))
                    angle = two_pi * unit(i,j+1)
                    if (2*pair-1 <= self%d) base(i,2*pair-1) = radius * cos(angle)
                    if (2*pair <= self%d) base(i,2*pair) = radius * sin(angle)
                end do
            end do
        end if

        if (self%has_root) then
            do i = 1, n
                do j = 1, self%d
                    sample(i,j) = self%mean(j) + dot_product(base(i,:), self%covariance_root(:,j))
                end do
            end do
        else
            do j = 1, self%d
                sample(:,j) = base(:,j) + self%mean(j)
            end do
        end if
    end subroutine qmc_multivariate_normal_random

    subroutine qmc_multivariate_normal_reset(self)
        class(qmc_multivariate_normal), intent(inout) :: self !! sampler reset to its initial QMC point

        call self%engine%reset()
    end subroutine qmc_multivariate_normal_reset

    subroutine qmc_poisson_disk_init(engine, d, status, radius, hypersphere, ncandidates, lower, upper)
        type(qmc_poisson_disk), intent(out) :: engine !! initialized Poisson-disk engine
        integer, intent(in) :: d !! sampling dimension
        integer, intent(out) :: status !! qmc_status_* result
        real(dp), intent(in), optional :: radius !! minimum point separation, default 0.05
        character(len=*), intent(in), optional :: hypersphere !! 'volume' or 'surface'
        integer, intent(in), optional :: ncandidates !! candidates generated per active center
        real(dp), intent(in), optional :: lower(:) !! lower bounds, default zero
        real(dp), intent(in), optional :: upper(:) !! upper bounds, default one

        character(len=:), allocatable :: method

        status = qmc_status_success
        if (d < 1) then
            status = qmc_status_invalid_parameter
            return
        end if
        engine%d = d
        engine%radius = 0.05_dp
        if (present(radius)) engine%radius = radius
        engine%ncandidates = 30
        if (present(ncandidates)) engine%ncandidates = ncandidates
        if (.not. ieee_is_finite(engine%radius) .or. engine%radius <= 0.0_dp .or. engine%ncandidates < 1) then
            status = qmc_status_invalid_parameter
            return
        end if
        engine%surface = .false.
        if (present(hypersphere)) then
            method = lower_ascii(trim(hypersphere))
            select case (method)
            case ('volume')
                engine%surface = .false.
            case ('surface')
                engine%surface = .true.
            case default
                status = qmc_status_invalid_parameter
                return
            end select
        end if
        if (engine%surface) then
            engine%radius_factor = 1.001_dp
        else
            engine%radius_factor = 2.0_dp
        end if
        engine%radius_squared = engine%radius * engine%radius

        allocate(engine%lower(d), engine%upper(d))
        engine%lower = 0.0_dp
        engine%upper = 1.0_dp
        if (present(lower)) then
            if (size(lower) /= d) then
                status = qmc_status_invalid_shape
                return
            end if
            engine%lower = lower
        end if
        if (present(upper)) then
            if (size(upper) /= d) then
                status = qmc_status_invalid_shape
                return
            end if
            engine%upper = upper
        end if
        if (.not. all(ieee_is_finite(engine%lower)) .or. .not. all(ieee_is_finite(engine%upper)) .or. &
                any(engine%upper <= engine%lower)) then
            status = qmc_status_invalid_parameter
            return
        end if
        allocate(engine%pool(16,d), engine%accepted(16,d))
        engine%pool = 0.0_dp
        engine%accepted = 0.0_dp
        engine%pool_size = 0
        engine%accepted_size = 0
        engine%num_generated = 0
    end subroutine qmc_poisson_disk_init

    subroutine qmc_poisson_disk_random(self, state, n, sample, status, n_drawn)
        class(qmc_poisson_disk), intent(inout) :: self !! Poisson-disk engine advanced by accepted points
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by proposals
        integer, intent(in) :: n !! maximum number of points requested
        real(dp), intent(out) :: sample(:, :) !! output buffer with shape (n,d); unused rows are NaN
        integer, intent(out) :: status !! qmc_status_* result
        integer, intent(out), optional :: n_drawn !! actual number of accepted points returned

        integer :: candidate_index
        integer :: center_index
        integer :: i
        integer :: produced
        real(dp), allocatable :: candidate(:)
        real(dp), allocatable :: center(:)

        status = qmc_status_success
        sample = quiet_nan(0.0_dp)
        if (present(n_drawn)) n_drawn = 0
        if (n < 0 .or. self%d < 1 .or. size(sample,1) /= n .or. size(sample,2) /= self%d) then
            status = qmc_status_invalid_shape
            return
        end if
        if (n == 0) return
        allocate(candidate(self%d), center(self%d))
        produced = 0

        if (self%accepted_size == 0) then
            do i = 1, self%d
                candidate(i) = self%lower(i) + rng_uniform(state) * (self%upper(i) - self%lower(i))
            end do
            call poisson_append_accepted(self, candidate)
            call poisson_append_pool(self, candidate)
            produced = 1
            sample(1,:) = candidate
            self%num_generated = self%num_generated + 1
        end if

        do while (self%pool_size > 0 .and. produced < n)
            center_index = 1 + int(floor(rng_uniform(state) * real(self%pool_size,dp)))
            center_index = min(self%pool_size,max(1,center_index))
            center = self%pool(center_index,:)
            call poisson_remove_pool(self, center_index)

            do candidate_index = 1, self%ncandidates
                call poisson_candidate(self, state, center, candidate)
                if (.not. poisson_in_bounds(self, candidate)) cycle
                if (poisson_has_neighbor(self, candidate)) cycle
                call poisson_append_accepted(self, candidate)
                call poisson_append_pool(self, candidate)
                produced = produced + 1
                sample(produced,:) = candidate
                self%num_generated = self%num_generated + 1
                if (produced >= n) exit
            end do
        end do
        if (present(n_drawn)) n_drawn = produced
    end subroutine qmc_poisson_disk_random

    subroutine qmc_poisson_disk_fill_space(self, state, sample, status)
        class(qmc_poisson_disk), intent(inout) :: self !! engine exhausted until no active centers remain
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by proposals
        real(dp), allocatable, intent(out) :: sample(:, :) !! all newly returned points with shape (n,d)
        integer, intent(out) :: status !! qmc_status_* result

        integer, parameter :: block_size = 256
        integer :: got
        integer :: total
        real(dp), allocatable :: block(:, :)
        real(dp), allocatable :: work(:, :)

        status = qmc_status_success
        total = 0
        allocate(sample(0,self%d), block(block_size,self%d))
        do
            call self%random(state, block_size, block, status, got)
            if (status /= qmc_status_success) return
            if (got > 0) then
                allocate(work(total+got,self%d))
                if (total > 0) work(1:total,:) = sample
                work(total+1:total+got,:) = block(1:got,:)
                call move_alloc(work,sample)
                total = total + got
            end if
            if (self%pool_size == 0) exit
        end do
    end subroutine qmc_poisson_disk_fill_space

    subroutine qmc_poisson_disk_reset(self)
        class(qmc_poisson_disk), intent(inout) :: self !! engine reset to an empty sampling state

        self%num_generated = 0
        self%pool_size = 0
        self%accepted_size = 0
        if (allocated(self%pool)) self%pool = 0.0_dp
        if (allocated(self%accepted)) self%accepted = 0.0_dp
    end subroutine qmc_poisson_disk_reset

    subroutine poisson_candidate(self, state, center, candidate)
        class(qmc_poisson_disk), intent(in) :: self !! Poisson-disk geometry and proposal settings
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by proposal generation
        real(dp), intent(in) :: center(:) !! active center point
        real(dp), intent(out) :: candidate(:) !! proposed point

        integer :: i
        real(dp) :: norm
        real(dp) :: proposal_radius
        real(dp), allocatable :: direction(:)

        allocate(direction(self%d))
        do
            do i = 1, self%d
                direction(i) = normal_rvs(state, 0.0_dp, 1.0_dp)
            end do
            norm = sqrt(sum(direction * direction))
            if (norm > 0.0_dp .and. ieee_is_finite(norm)) exit
        end do
        direction = direction / norm
        if (self%surface) then
            proposal_radius = self%radius * self%radius_factor
        else
            proposal_radius = self%radius * self%radius_factor * rng_uniform(state) ** (1.0_dp / real(self%d,dp))
        end if
        candidate = center + proposal_radius * direction
    end subroutine poisson_candidate

    pure function poisson_in_bounds(self, point) result(valid)
        class(qmc_poisson_disk), intent(in) :: self !! Poisson-disk bounds
        real(dp), intent(in) :: point(:) !! candidate point
        logical :: valid

        valid = all(point >= self%lower) .and. all(point <= self%upper)
    end function poisson_in_bounds

    function poisson_has_neighbor(self, point) result(found)
        class(qmc_poisson_disk), intent(in) :: self !! accepted-point collection
        real(dp), intent(in) :: point(:) !! candidate point
        logical :: found

        integer :: i

        found = .false.
        do i = 1, self%accepted_size
            if (sum((point - self%accepted(i,:)) ** 2) < self%radius_squared) then
                found = .true.
                return
            end if
        end do
    end function poisson_has_neighbor

    subroutine poisson_append_pool(self, point)
        class(qmc_poisson_disk), intent(inout) :: self !! active-point collection grown as needed
        real(dp), intent(in) :: point(:) !! point appended to active pool

        real(dp), allocatable :: work(:, :)
        integer :: capacity

        capacity = size(self%pool,1)
        if (self%pool_size >= capacity) then
            allocate(work(max(2*capacity,1),self%d))
            work = 0.0_dp
            if (self%pool_size > 0) work(1:self%pool_size,:) = self%pool(1:self%pool_size,:)
            call move_alloc(work,self%pool)
        end if
        self%pool_size = self%pool_size + 1
        self%pool(self%pool_size,:) = point
    end subroutine poisson_append_pool

    subroutine poisson_append_accepted(self, point)
        class(qmc_poisson_disk), intent(inout) :: self !! accepted-point collection grown as needed
        real(dp), intent(in) :: point(:) !! point appended to accepted set

        real(dp), allocatable :: work(:, :)
        integer :: capacity

        capacity = size(self%accepted,1)
        if (self%accepted_size >= capacity) then
            allocate(work(max(2*capacity,1),self%d))
            work = 0.0_dp
            if (self%accepted_size > 0) work(1:self%accepted_size,:) = self%accepted(1:self%accepted_size,:)
            call move_alloc(work,self%accepted)
        end if
        self%accepted_size = self%accepted_size + 1
        self%accepted(self%accepted_size,:) = point
    end subroutine poisson_append_accepted

    subroutine poisson_remove_pool(self, index)
        class(qmc_poisson_disk), intent(inout) :: self !! active pool modified in place
        integer, intent(in) :: index !! one-based active point removed

        if (index < self%pool_size) self%pool(index,:) = self%pool(self%pool_size,:)
        self%pool_size = self%pool_size - 1
    end subroutine poisson_remove_pool

    subroutine qmc_latin_hypercube_init(engine, d, status, scramble, strength)
        type(qmc_latin_hypercube), intent(out) :: engine !! initialized Latin-hypercube engine
        integer, intent(in) :: d !! design dimension, at least one
        integer, intent(out) :: status !! qmc_status_* result
        logical, intent(in), optional :: scramble !! randomize positions inside strata
        integer, intent(in), optional :: strength !! one for plain LHS or two for OA-LHS

        status = qmc_status_success
        if (d < 1) then
            status = qmc_status_invalid_parameter
            return
        end if
        engine%d = d
        engine%num_generated = 0
        engine%scramble = .true.
        if (present(scramble)) engine%scramble = scramble
        engine%strength = 1
        if (present(strength)) engine%strength = strength
        if (engine%strength /= 1 .and. engine%strength /= 2) then
            status = qmc_status_invalid_parameter
        end if
    end subroutine qmc_latin_hypercube_init

    subroutine qmc_latin_hypercube_random(self, state, n, sample, status)
        class(qmc_latin_hypercube), intent(inout) :: self !! design engine whose count is advanced
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by randomization
        integer, intent(in) :: n !! number of points, p**2 for strength two
        real(dp), intent(out) :: sample(:, :) !! generated sample with shape (n,d)
        integer, intent(out) :: status !! qmc_status_* result

        if (n < 1 .or. self%d < 1 .or. size(sample,1) /= n .or. size(sample,2) /= self%d) then
            sample = 0.0_dp
            status = qmc_status_invalid_shape
            return
        end if
        if (self%strength == 1) then
            call plain_lhs(state, n, self%d, self%scramble, sample)
            status = qmc_status_success
        else if (self%strength == 2) then
            call oa_lhs(state, n, self%d, self%scramble, sample, status)
        else
            sample = 0.0_dp
            status = qmc_status_invalid_parameter
        end if
        if (status == qmc_status_success) self%num_generated = self%num_generated + n
    end subroutine qmc_latin_hypercube_random

    subroutine qmc_latin_hypercube_reset(self)
        class(qmc_latin_hypercube), intent(inout) :: self !! design engine count reset to zero

        self%num_generated = 0
    end subroutine qmc_latin_hypercube_reset

    subroutine qmc_to_integers(sample, lower, upper, values, status, endpoint)
        real(dp), intent(in) :: sample(:, :) !! unit-hypercube sample with shape (n,d)
        integer(int64), intent(in) :: lower(:) !! inclusive lower integer bound per coordinate
        integer(int64), intent(in) :: upper(:) !! exclusive upper bound, or inclusive when endpoint=true
        integer(int64), intent(out) :: values(:, :) !! integer-valued mapped sample
        integer, intent(out) :: status !! qmc_status_* result
        logical, intent(in), optional :: endpoint !! include the upper bound when true

        integer :: i
        integer :: j
        integer(int64) :: width
        logical :: include_upper

        status = qmc_status_success
        values = 0_int64
        include_upper = .false.
        if (present(endpoint)) include_upper = endpoint
        if (size(sample,2) /= size(lower) .or. size(lower) /= size(upper) .or. &
                any(shape(values) /= shape(sample))) then
            status = qmc_status_invalid_shape
            return
        end if
        if (.not. valid_unit_sample(sample)) then
            status = qmc_status_out_of_bounds
            return
        end if
        do j = 1, size(lower)
            width = upper(j) - lower(j)
            if (include_upper) width = width + 1_int64
            if (width <= 0_int64) then
                status = qmc_status_invalid_parameter
                return
            end if
            do i = 1, size(sample,1)
                values(i,j) = lower(j) + int(floor(real(width,dp) * sample(i,j)), int64)
                if (include_upper) values(i,j) = min(values(i,j), upper(j))
            end do
        end do
    end subroutine qmc_to_integers

    subroutine plain_lhs(state, n, d, scramble, sample)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by cell randomization
        integer, intent(in) :: n !! number of strata per coordinate
        integer, intent(in) :: d !! number of coordinates
        logical, intent(in) :: scramble !! use uniform offsets rather than cell centers
        real(dp), intent(out) :: sample(n,d) !! Latin-hypercube sample

        integer :: i
        integer :: j
        integer, allocatable :: perm(:)
        real(dp) :: offset

        allocate(perm(n))
        do j = 1, d
            do i = 1, n
                perm(i) = i
            end do
            call shuffle_one_based(state, perm)
            do i = 1, n
                offset = 0.5_dp
                if (scramble) offset = rng_uniform(state)
                sample(i,j) = (real(perm(i),dp) - offset) / real(n,dp)
            end do
        end do
    end subroutine plain_lhs

    subroutine oa_lhs(state, n, d, scramble, sample, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by OA/LHS scrambling
        integer, intent(in) :: n !! required square of a prime p
        integer, intent(in) :: d !! dimension, no greater than p+1
        logical, intent(in) :: scramble !! randomize within one-dimensional strata
        real(dp), intent(out) :: sample(n,d) !! strength-two orthogonal-array LHS
        integer, intent(out) :: status !! qmc_status_* result

        integer :: col
        integer :: idx
        integer :: p
        integer :: row
        integer :: symbol
        integer :: u
        integer :: v
        integer, allocatable :: oa(:, :)
        integer, allocatable :: perm(:)
        integer, allocatable :: rows(:)
        real(dp), allocatable :: lhs(:, :)

        status = qmc_status_success
        sample = 0.0_dp
        p = nint(sqrt(real(n,dp)))
        if (p * p /= n .or. .not. is_prime(p) .or. d > p + 1) then
            status = qmc_status_invalid_parameter
            return
        end if
        allocate(oa(n,p+1), perm(0:p-1), rows(p), lhs(p,1))
        row = 0
        do u = 0, p - 1
            do v = 0, p - 1
                row = row + 1
                oa(row,1) = u
                oa(row,2) = v
                do col = 3, p + 1
                    oa(row,col) = modulo(u + (col - 2) * v, p)
                end do
            end do
        end do

        do col = 1, d
            do symbol = 0, p - 1
                perm(symbol) = symbol
            end do
            call shuffle_zero_based(state, perm)
            do row = 1, n
                oa(row,col) = perm(oa(row,col))
            end do
        end do

        do col = 1, d
            do symbol = 0, p - 1
                idx = 0
                do row = 1, n
                    if (oa(row,col) == symbol) then
                        idx = idx + 1
                        rows(idx) = row
                    end if
                end do
                call plain_lhs(state, p, 1, scramble, lhs)
                do idx = 1, p
                    sample(rows(idx),col) = (real(symbol,dp) + lhs(idx,1)) / real(p,dp)
                end do
            end do
        end do
    end subroutine oa_lhs

    subroutine first_primes(n, primes)
        integer, intent(in) :: n !! number of primes requested
        integer, intent(out) :: primes(n) !! first n primes in increasing order

        integer :: candidate
        integer :: found

        candidate = 2
        found = 0
        do while (found < n)
            if (is_prime(candidate)) then
                found = found + 1
                primes(found) = candidate
            end if
            candidate = candidate + 1
        end do
    end subroutine first_primes

    pure function is_prime(n) result(prime)
        integer, intent(in) :: n !! integer tested for primality
        logical :: prime

        integer :: divisor

        if (n < 2) then
            prime = .false.
            return
        end if
        if (n == 2) then
            prime = .true.
            return
        end if
        if (modulo(n,2) == 0) then
            prime = .false.
            return
        end if
        divisor = 3
        do while (divisor * divisor <= n)
            if (modulo(n,divisor) == 0) then
                prime = .false.
                return
            end if
            divisor = divisor + 2
        end do
        prime = .true.
    end function is_prime

    subroutine shuffle_one_based(state, values)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by shuffle draws
        integer, intent(inout) :: values(:) !! one-based list randomized in place

        integer :: i
        integer :: j
        integer :: tmp

        do i = size(values), 2, -1
            j = 1 + int(floor(rng_uniform(state) * real(i,dp)))
            j = min(i, max(1,j))
            tmp = values(i)
            values(i) = values(j)
            values(j) = tmp
        end do
    end subroutine shuffle_one_based

    subroutine shuffle_zero_based(state, values)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by shuffle draws
        integer, intent(inout) :: values(0:) !! zero-based permutation randomized in place

        integer :: i
        integer :: j
        integer :: tmp

        do i = ubound(values,1), 1, -1
            j = int(floor(rng_uniform(state) * real(i + 1,dp)))
            j = min(i, max(0,j))
            tmp = values(i)
            values(i) = values(j)
            values(j) = tmp
        end do
    end subroutine shuffle_zero_based

    pure function point_distance(x, y, metric) result(value)
        real(dp), intent(in) :: x(:) !! first point
        real(dp), intent(in) :: y(:) !! second point
        character(len=*), intent(in) :: metric !! supported distance metric name
        real(dp) :: value

        if (metric == 'cityblock') then
            value = sum(abs(x - y))
        else
            value = sqrt(sum((x - y) ** 2))
        end if
    end function point_distance

    pure function valid_unit_sample(sample) result(valid)
        real(dp), intent(in) :: sample(:, :) !! candidate unit-hypercube points
        logical :: valid

        valid = size(sample,1) >= 1 .and. size(sample,2) >= 1 .and. &
            all(ieee_is_finite(sample)) .and. all(sample >= 0.0_dp) .and. all(sample <= 1.0_dp)
    end function valid_unit_sample

    pure function lower_ascii(text) result(lower)
        character(len=*), intent(in) :: text !! ASCII text to convert to lower case
        character(len=len(text)) :: lower

        integer :: c
        integer :: i

        lower = text
        do i = 1, len(text)
            c = iachar(text(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) lower(i:i) = achar(c + 32)
        end do
    end function lower_ascii

end module scifort_qmc
