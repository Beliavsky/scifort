! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Adapted portions: SciPy 1.17.0 scipy/stats/_qmvnt.py and _qmvnt_cy.pyx (box transform).
! Retained BSD terms: THIRD_PARTY_LICENSES.md; details: CODE_PROVENANCE.md.
!
! Box probabilities for multivariate normal and Student t laws.
!
! The conditioning transformation and covariance permutation follow the
! Genz/Bretz strategy used by SciPy's scipy.stats._qmvnt implementation.
! SciFort uses randomized, tent-transformed Halton points instead of SciPy's
! FFT-based component-by-component lattice construction so this module stays
! self contained and does not add an FFT dependency.

module scifort_multivariate_integration
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_gamma, only : gammaincinv
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_cdf, normal_ppf
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_student_t, only : t_cdf
    implicit none
    private

    integer, parameter, public :: qmv_status_success = 0
    integer, parameter, public :: qmv_status_invalid_shape = 1
    integer, parameter, public :: qmv_status_invalid_parameter = 2
    integer, parameter, public :: qmv_status_not_psd = 3

    integer, parameter :: default_batches = 10
    integer, parameter :: minimum_points_per_batch = 64
    real(dp), parameter :: cholesky_tol = 1.0e-10_dp
    real(dp), parameter :: sqrt_two_pi = 2.506628274631000502415765284811045_dp

    public :: qmvn_box
    public :: qmvt_box

contains

    subroutine qmvn_box(covar, low, high, probability, estimated_error, n_samples, status, &
            state, maxpts, abseps)
        real(dp), intent(in) :: covar(:, :) !! positive-semidefinite covariance matrix
        real(dp), intent(in) :: low(:) !! lower integration limits
        real(dp), intent(in) :: high(:) !! upper integration limits
        real(dp), intent(out) :: probability !! estimated probability of the box
        real(dp), intent(out) :: estimated_error !! three-standard-error batch estimate
        integer, intent(out) :: n_samples !! number of low-discrepancy points evaluated
        integer, intent(out) :: status !! zero on success; nonzero for invalid input
        type(rng_state), intent(inout), optional :: state !! optional RNG for randomized shifts
        integer, intent(in), optional :: maxpts !! approximate maximum number of integration points
        real(dp), intent(in), optional :: abseps !! requested absolute error target

        integer :: limit
        real(dp) :: target
        type(rng_state) :: local_state

        call validate_box_input(covar, low, high, status)
        if (status /= qmv_status_success) then
            probability = quiet_nan(0.0_dp)
            estimated_error = probability
            n_samples = 0
            return
        end if

        limit = default_limit(size(low), 1000000, maxpts)
        target = 1.0e-5_dp
        if (present(abseps)) target = abseps
        if (.not. ieee_is_finite(target) .or. target <= 0.0_dp) then
            probability = quiet_nan(0.0_dp)
            estimated_error = probability
            n_samples = 0
            status = qmv_status_invalid_parameter
            return
        end if

        if (present(state)) then
            call qmvn_auto(covar, low, high, state, limit, target, probability, &
                estimated_error, n_samples, status)
        else
            call rng_seed(local_state, 170141183)
            call qmvn_auto(covar, low, high, local_state, limit, target, probability, &
                estimated_error, n_samples, status)
        end if
    end subroutine qmvn_box

    subroutine qmvt_box(covar, low, high, df, probability, estimated_error, n_samples, status, &
            state, maxpts, abseps)
        real(dp), intent(in) :: covar(:, :) !! positive-semidefinite shape matrix
        real(dp), intent(in) :: low(:) !! lower integration limits
        real(dp), intent(in) :: high(:) !! upper integration limits
        real(dp), intent(in) :: df !! degrees of freedom, > 0; +infinity gives normal
        real(dp), intent(out) :: probability !! estimated probability of the box
        real(dp), intent(out) :: estimated_error !! three-standard-error batch estimate
        integer, intent(out) :: n_samples !! number of low-discrepancy points evaluated
        integer, intent(out) :: status !! zero on success; nonzero for invalid input
        type(rng_state), intent(inout), optional :: state !! optional RNG for randomized shifts
        integer, intent(in), optional :: maxpts !! approximate maximum number of integration points
        real(dp), intent(in), optional :: abseps !! requested absolute error target

        integer :: limit
        real(dp) :: target
        type(rng_state) :: local_state

        call validate_box_input(covar, low, high, status)
        if (status /= qmv_status_success) then
            probability = quiet_nan(0.0_dp)
            estimated_error = probability
            n_samples = 0
            return
        end if
        if (ieee_is_nan(df) .or. df <= 0.0_dp) then
            probability = quiet_nan(0.0_dp)
            estimated_error = probability
            n_samples = 0
            status = qmv_status_invalid_parameter
            return
        end if

        if (.not. ieee_is_finite(df)) then
            if (present(state)) then
                call qmvn_box(covar, low, high, probability, estimated_error, n_samples, status, &
                    state=state, maxpts=maxpts, abseps=abseps)
            else
                call qmvn_box(covar, low, high, probability, estimated_error, n_samples, status, &
                    maxpts=maxpts, abseps=abseps)
            end if
            return
        end if

        limit = default_limit(size(low), 1000, maxpts)
        target = 1.0e-5_dp
        if (present(abseps)) target = abseps
        if (.not. ieee_is_finite(target) .or. target <= 0.0_dp) then
            probability = quiet_nan(0.0_dp)
            estimated_error = probability
            n_samples = 0
            status = qmv_status_invalid_parameter
            return
        end if

        if (present(state)) then
            call qmvt_auto(covar, low, high, df, state, limit, target, probability, &
                estimated_error, n_samples, status)
        else
            call rng_seed(local_state, 104729)
            call qmvt_auto(covar, low, high, df, local_state, limit, target, probability, &
                estimated_error, n_samples, status)
        end if
    end subroutine qmvt_box

    subroutine qmvn_auto(covar, low, high, state, limit, target, probability, &
            estimated_error, n_samples, status)
        real(dp), intent(in) :: covar(:, :) !! covariance matrix
        real(dp), intent(in) :: low(:) !! lower limits
        real(dp), intent(in) :: high(:) !! upper limits
        type(rng_state), intent(inout) :: state !! RNG used for randomized shifts
        integer, intent(in) :: limit !! integration point limit
        real(dp), intent(in) :: target !! absolute error target
        real(dp), intent(out) :: probability !! probability estimate
        real(dp), intent(out) :: estimated_error !! error estimate
        integer, intent(out) :: n_samples !! points evaluated
        integer, intent(out) :: status !! integration status

        integer :: m
        integer :: n
        integer :: ni
        integer :: remaining
        real(dp), allocatable :: cho(:, :)
        real(dp), allocatable :: lo(:)
        real(dp), allocatable :: hi(:)
        real(dp) :: ei
        real(dp) :: pi
        real(dp) :: weight

        n = size(low)
        if (n == 1) then
            call univariate_normal_box(covar(1, 1), low(1), high(1), probability, status)
            estimated_error = epsilon(1.0_dp)
            n_samples = 0
            return
        end if

        allocate(cho(n, n), lo(n), hi(n))
        call permuted_cholesky(covar, low, high, cho, lo, hi, status)
        if (status /= qmv_status_success) then
            probability = quiet_nan(0.0_dp)
            estimated_error = probability
            n_samples = 0
            return
        end if

        probability = 0.0_dp
        estimated_error = 1.0_dp
        n_samples = 0
        m = max(minimum_points_per_batch, (n * 1000) / default_batches)
        do while (n_samples < limit .and. estimated_error > target)
            remaining = limit - n_samples
            if (remaining < default_batches) exit
            m = min(m, max(1, remaining / default_batches))
            call qmvn_qmc(cho, lo, hi, state, m, default_batches, pi, ei, ni)
            if (n_samples == 0) then
                probability = pi
                estimated_error = ei
            else
                if (estimated_error <= tiny(1.0_dp)) then
                    weight = 0.0_dp
                else
                    weight = 1.0_dp / (1.0_dp + (ei / estimated_error) ** 2)
                end if
                probability = probability + weight * (pi - probability)
                estimated_error = sqrt(max(0.0_dp, weight)) * ei
            end if
            n_samples = n_samples + ni
            if (m >= max(1, (limit - n_samples) / default_batches)) exit
            m = max(m + 1, nint(sqrt(2.0_dp) * real(m, dp)))
        end do
        probability = max(0.0_dp, min(1.0_dp, probability))
        status = qmv_status_success
    end subroutine qmvn_auto

    subroutine qmvt_auto(covar, low, high, df, state, limit, target, probability, &
            estimated_error, n_samples, status)
        real(dp), intent(in) :: covar(:, :) !! shape matrix
        real(dp), intent(in) :: low(:) !! lower limits
        real(dp), intent(in) :: high(:) !! upper limits
        real(dp), intent(in) :: df !! degrees of freedom
        type(rng_state), intent(inout) :: state !! RNG used for randomized shifts
        integer, intent(in) :: limit !! integration point limit
        real(dp), intent(in) :: target !! absolute error target
        real(dp), intent(out) :: probability !! probability estimate
        real(dp), intent(out) :: estimated_error !! error estimate
        integer, intent(out) :: n_samples !! points evaluated
        integer, intent(out) :: status !! integration status

        integer :: m
        integer :: n
        integer :: ni
        integer :: remaining
        real(dp), allocatable :: cho(:, :)
        real(dp), allocatable :: lo(:)
        real(dp), allocatable :: hi(:)
        real(dp) :: ei
        real(dp) :: pi
        real(dp) :: scale
        real(dp) :: weight

        n = size(low)
        if (n == 1) then
            call univariate_t_box(covar(1, 1), low(1), high(1), df, probability, status)
            estimated_error = epsilon(1.0_dp)
            n_samples = 0
            return
        end if

        scale = max(1.0_dp, sqrt(df))
        allocate(cho(n, n), lo(n), hi(n))
        call permuted_cholesky(covar, low / scale, high / scale, cho, lo, hi, status)
        if (status /= qmv_status_success) then
            probability = quiet_nan(0.0_dp)
            estimated_error = probability
            n_samples = 0
            return
        end if

        probability = 0.0_dp
        estimated_error = 1.0_dp
        n_samples = 0
        m = max(minimum_points_per_batch, (n * 1000) / default_batches)
        do while (n_samples < limit .and. estimated_error > target)
            remaining = limit - n_samples
            if (remaining < default_batches) exit
            m = min(m, max(1, remaining / default_batches))
            call qmvt_qmc(cho, lo, hi, df, state, m, default_batches, pi, ei, ni)
            if (n_samples == 0) then
                probability = pi
                estimated_error = ei
            else
                if (estimated_error <= tiny(1.0_dp)) then
                    weight = 0.0_dp
                else
                    weight = 1.0_dp / (1.0_dp + (ei / estimated_error) ** 2)
                end if
                probability = probability + weight * (pi - probability)
                estimated_error = sqrt(max(0.0_dp, weight)) * ei
            end if
            n_samples = n_samples + ni
            if (m >= max(1, (limit - n_samples) / default_batches)) exit
            m = max(m + 1, nint(sqrt(2.0_dp) * real(m, dp)))
        end do
        probability = max(0.0_dp, min(1.0_dp, probability))
        status = qmv_status_success
    end subroutine qmvt_auto

    subroutine qmvn_qmc(cho, lo, hi, state, n_qmc, n_batches, probability, &
            estimated_error, n_samples)
        real(dp), intent(in) :: cho(:, :) !! permuted/scaled lower Cholesky factor
        real(dp), intent(in) :: lo(:) !! permuted/scaled lower limits
        real(dp), intent(in) :: hi(:) !! permuted/scaled upper limits
        type(rng_state), intent(inout) :: state !! RNG used for shifts
        integer, intent(in) :: n_qmc !! low-discrepancy points per batch
        integer, intent(in) :: n_batches !! number of independently shifted batches
        real(dp), intent(out) :: probability !! mean batch probability
        real(dp), intent(out) :: estimated_error !! three-standard-error estimate
        integer, intent(out) :: n_samples !! total integrand evaluations

        integer, allocatable :: bases(:)
        real(dp), allocatable :: batch(:)
        real(dp), allocatable :: shifts(:)
        real(dp), allocatable :: y(:)
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        real(dp) :: c
        real(dp) :: ct
        real(dp) :: d
        real(dp) :: dc
        real(dp) :: p
        real(dp) :: pv
        real(dp) :: s
        real(dp) :: x
        real(dp) :: z

        n = size(lo)
        allocate(batch(n_batches), shifts(max(1, n - 1)), y(max(1, n - 1)))
        allocate(bases(max(1, n - 1)))
        do i = 1, max(1, n - 1)
            bases(i) = nth_prime(i)
        end do

        do j = 1, n_batches
            do i = 1, size(shifts)
                shifts(i) = rng_uniform(state)
            end do
            batch(j) = 0.0_dp
            do k = 1, n_qmc
                call first_normal_interval(cho(1, 1), lo(1), hi(1), c, d)
                dc = max(0.0_dp, d - c)
                pv = dc
                if (pv <= 0.0_dp) then
                    cycle
                end if
                y = 0.0_dp
                do i = 2, n
                    z = modulo(radical_inverse(k, bases(i - 1)) + shifts(i - 1), 1.0_dp)
                    x = abs(2.0_dp * z - 1.0_dp)
                    p = c + x * dc
                    y(i - 1) = finite_normal_ppf(p)
                    s = dot_product(cho(i, 1:i - 1), y(1:i - 1))
                    ct = cho(i, i)
                    if (abs(ct) <= cholesky_tol) then
                        if (lo(i) <= s .and. s <= hi(i)) then
                            c = 0.0_dp
                            d = 1.0_dp
                        else
                            c = 0.0_dp
                            d = 0.0_dp
                        end if
                    else
                        c = normal_cdf((lo(i) - s) / ct)
                        d = normal_cdf((hi(i) - s) / ct)
                    end if
                    dc = max(0.0_dp, d - c)
                    pv = pv * dc
                    if (pv <= 0.0_dp) exit
                end do
                batch(j) = batch(j) + pv
            end do
            batch(j) = batch(j) / real(n_qmc, dp)
        end do

        probability = sum(batch) / real(n_batches, dp)
        estimated_error = three_standard_error(batch, probability)
        n_samples = n_qmc * n_batches
    end subroutine qmvn_qmc

    subroutine qmvt_qmc(cho, lo, hi, df, state, n_qmc, n_batches, probability, &
            estimated_error, n_samples)
        real(dp), intent(in) :: cho(:, :) !! permuted/scaled lower Cholesky factor
        real(dp), intent(in) :: lo(:) !! permuted/scaled lower limits
        real(dp), intent(in) :: hi(:) !! permuted/scaled upper limits
        real(dp), intent(in) :: df !! degrees of freedom
        type(rng_state), intent(inout) :: state !! RNG used for shifts
        integer, intent(in) :: n_qmc !! low-discrepancy points per batch
        integer, intent(in) :: n_batches !! number of independently shifted batches
        real(dp), intent(out) :: probability !! mean batch probability
        real(dp), intent(out) :: estimated_error !! three-standard-error estimate
        integer, intent(out) :: n_samples !! total integrand evaluations

        integer, allocatable :: bases(:)
        real(dp), allocatable :: batch(:)
        real(dp), allocatable :: c(:)
        real(dp), allocatable :: dc(:)
        real(dp), allocatable :: shifts(:)
        real(dp), allocatable :: s(:)
        integer :: i
        integer :: ip
        integer :: j
        integer :: k
        integer :: n
        real(dp) :: ck
        real(dp) :: dk
        real(dp) :: hi_scaled
        real(dp) :: lo_scaled
        real(dp) :: p
        real(dp) :: pv
        real(dp) :: r
        real(dp) :: x
        real(dp) :: y
        real(dp) :: z

        n = size(lo)
        allocate(batch(n_batches), shifts(n), bases(n), c(n_qmc), dc(n_qmc), s(n))
        do i = 1, n
            bases(i) = nth_prime(i)
        end do

        do j = 1, n_batches
            do i = 1, n
                shifts(i) = rng_uniform(state)
            end do
            batch(j) = 0.0_dp
            do k = 1, n_qmc
                s = 0.0_dp
                pv = 1.0_dp
                c(k) = 1.0_dp
                dc(k) = 0.0_dp

                z = modulo(radical_inverse(k, bases(1)) + shifts(1), 1.0_dp)
                x = abs(2.0_dp * z - 1.0_dp)
                r = sqrt(2.0_dp * gammaincinv(0.5_dp * df, x))

                lo_scaled = lo(1) * r
                hi_scaled = hi(1) * r
                ck = clipped_normal_cdf(lo_scaled)
                dk = clipped_normal_cdf(hi_scaled)
                c(k) = ck
                dc(k) = max(0.0_dp, dk - ck)
                pv = dc(k)

                do i = 2, n
                    if (pv <= 0.0_dp) exit
                    z = modulo(radical_inverse(k, bases(i)) + shifts(i), 1.0_dp)
                    x = abs(2.0_dp * z - 1.0_dp)
                    p = c(k) + x * dc(k)
                    y = finite_normal_ppf(p)
                    do ip = i, n
                        s(ip) = s(ip) + cho(ip, i - 1) * y
                    end do

                    lo_scaled = lo(i) * r - s(i)
                    hi_scaled = hi(i) * r - s(i)
                    ck = clipped_normal_cdf(lo_scaled)
                    dk = clipped_normal_cdf(hi_scaled)
                    c(k) = ck
                    dc(k) = max(0.0_dp, dk - ck)
                    pv = pv * dc(k)
                end do
                batch(j) = batch(j) + pv
            end do
            batch(j) = batch(j) / real(n_qmc, dp)
        end do

        probability = sum(batch) / real(n_batches, dp)
        estimated_error = three_standard_error(batch, probability)
        n_samples = n_qmc * n_batches
    end subroutine qmvt_qmc

    subroutine permuted_cholesky(covar, low, high, cho, new_low, new_high, status)
        real(dp), intent(in) :: covar(:, :) !! PSD covariance/shape matrix; lower triangle is authoritative
        real(dp), intent(in) :: low(:) !! lower limits
        real(dp), intent(in) :: high(:) !! upper limits
        real(dp), intent(out) :: cho(:, :) !! scaled and permuted Cholesky-like factor
        real(dp), intent(out) :: new_low(:) !! correspondingly scaled/permuted lower limits
        real(dp), intent(out) :: new_high(:) !! correspondingly scaled/permuted upper limits
        integer, intent(out) :: status !! zero on success

        integer :: i
        integer :: im
        integer :: k
        integer :: n
        real(dp), allocatable :: dc(:)
        real(dp), allocatable :: y(:)
        real(dp) :: ci
        real(dp) :: ck
        real(dp) :: de
        real(dp) :: dem
        real(dp) :: epk
        real(dp) :: hi_i
        real(dp) :: hi_m
        real(dp) :: lo_i
        real(dp) :: lo_m
        real(dp) :: s
        real(dp) :: tmp

        n = size(low)
        if (size(covar, 1) /= n .or. size(covar, 2) /= n .or. size(high) /= n .or. &
                size(cho, 1) /= n .or. size(cho, 2) /= n .or. size(new_low) /= n .or. &
                size(new_high) /= n) then
            status = qmv_status_invalid_shape
            return
        end if

        cho = 0.0_dp
        do i = 1, n
            cho(i, i) = covar(i, i)
            do k = 1, i - 1
                cho(i, k) = covar(i, k)
                cho(k, i) = covar(i, k)
            end do
        end do
        new_low = low
        new_high = high
        allocate(dc(n), y(n))
        y = 0.0_dp

        do i = 1, n
            if (cho(i, i) < -cholesky_tol) then
                status = qmv_status_not_psd
                return
            end if
            dc(i) = sqrt(max(cho(i, i), 0.0_dp))
            if (dc(i) == 0.0_dp) dc(i) = 1.0_dp
            new_low(i) = new_low(i) / dc(i)
            new_high(i) = new_high(i) / dc(i)
        end do
        do i = 1, n
            do k = 1, n
                cho(i, k) = cho(i, k) / (dc(i) * dc(k))
            end do
        end do

        do k = 1, n
            epk = real(k, dp) * cholesky_tol
            im = k
            ck = 0.0_dp
            dem = 1.0_dp
            lo_m = 0.0_dp
            hi_m = 0.0_dp
            do i = k, n
                if (cho(i, i) > cholesky_tol) then
                    ci = sqrt(cho(i, i))
                    if (k > 1) then
                        s = dot_product(cho(i, 1:k - 1), y(1:k - 1))
                    else
                        s = 0.0_dp
                    end if
                    lo_i = (new_low(i) - s) / ci
                    hi_i = (new_high(i) - s) / ci
                    de = normal_cdf(hi_i) - normal_cdf(lo_i)
                    if (de <= dem) then
                        ck = ci
                        dem = de
                        lo_m = lo_i
                        hi_m = hi_i
                        im = i
                    end if
                end if
            end do

            if (im > k) then
                tmp = cho(im, im)
                cho(im, im) = cho(k, k)
                cho(k, k) = tmp
                if (k > 1) call swap_rows_prefix(cho, im, k, k - 1)
                if (im < n) call swap_column_tail(cho, im, k, im + 1, n)
                if (im > k + 1) call swap_cross_segment(cho, im, k)
                tmp = new_low(k)
                new_low(k) = new_low(im)
                new_low(im) = tmp
                tmp = new_high(k)
                new_high(k) = new_high(im)
                new_high(im) = tmp
            end if

            if (ck > epk) then
                cho(k, k) = ck
                if (k < n) cho(k, k + 1:n) = 0.0_dp
                do i = k + 1, n
                    cho(i, k) = cho(i, k) / ck
                    if (i >= k + 1) then
                        cho(i, k + 1:i) = cho(i, k + 1:i) - &
                            cho(i, k) * cho(k + 1:i, k)
                    end if
                end do
                if (abs(dem) > cholesky_tol) then
                    y(k) = (exp(-0.5_dp * lo_m * lo_m) - exp(-0.5_dp * hi_m * hi_m)) / &
                        (sqrt_two_pi * dem)
                else
                    y(k) = 0.5_dp * (lo_m + hi_m)
                    if (lo_m < -10.0_dp) then
                        y(k) = hi_m
                    else if (hi_m > 10.0_dp) then
                        y(k) = lo_m
                    end if
                end if
                cho(k, 1:k) = cho(k, 1:k) / ck
                new_low(k) = new_low(k) / ck
                new_high(k) = new_high(k) / ck
            else
                cho(k:n, k) = 0.0_dp
                y(k) = 0.5_dp * (new_low(k) + new_high(k))
            end if
        end do
        status = qmv_status_success
    end subroutine permuted_cholesky

    subroutine validate_box_input(covar, low, high, status)
        real(dp), intent(in) :: covar(:, :) !! covariance or shape matrix
        real(dp), intent(in) :: low(:) !! lower bounds
        real(dp), intent(in) :: high(:) !! upper bounds
        integer, intent(out) :: status !! zero if dimensions and values are acceptable

        integer :: n

        n = size(low)
        if (n < 1 .or. size(high) /= n .or. size(covar, 1) /= n .or. size(covar, 2) /= n) then
            status = qmv_status_invalid_shape
            return
        end if
        if (any(.not. ieee_is_finite(covar)) .or. any(ieee_is_nan(low)) .or. any(ieee_is_nan(high))) then
            status = qmv_status_invalid_parameter
            return
        end if
        status = qmv_status_success
    end subroutine validate_box_input

    subroutine univariate_normal_box(variance, low, high, probability, status)
        real(dp), intent(in) :: variance !! nonnegative variance
        real(dp), intent(in) :: low !! lower limit
        real(dp), intent(in) :: high !! upper limit
        real(dp), intent(out) :: probability !! exact one-dimensional probability
        integer, intent(out) :: status !! zero on success

        real(dp) :: sd

        if (variance < 0.0_dp) then
            probability = quiet_nan(0.0_dp)
            status = qmv_status_not_psd
            return
        end if
        if (variance == 0.0_dp) then
            if (low <= 0.0_dp .and. 0.0_dp <= high) then
                probability = 1.0_dp
            else
                probability = 0.0_dp
            end if
        else
            sd = sqrt(variance)
            probability = normal_cdf(high / sd) - normal_cdf(low / sd)
        end if
        probability = max(0.0_dp, min(1.0_dp, probability))
        status = qmv_status_success
    end subroutine univariate_normal_box

    subroutine univariate_t_box(shape, low, high, df, probability, status)
        real(dp), intent(in) :: shape !! nonnegative one-dimensional shape
        real(dp), intent(in) :: low !! lower limit
        real(dp), intent(in) :: high !! upper limit
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(out) :: probability !! exact one-dimensional probability
        integer, intent(out) :: status !! zero on success

        real(dp) :: scale

        if (shape < 0.0_dp) then
            probability = quiet_nan(0.0_dp)
            status = qmv_status_not_psd
            return
        end if
        if (shape == 0.0_dp) then
            if (low <= 0.0_dp .and. 0.0_dp <= high) then
                probability = 1.0_dp
            else
                probability = 0.0_dp
            end if
        else
            scale = sqrt(shape)
            probability = t_cdf(high, df, scale=scale) - t_cdf(low, df, scale=scale)
        end if
        probability = max(0.0_dp, min(1.0_dp, probability))
        status = qmv_status_success
    end subroutine univariate_t_box

    pure subroutine first_normal_interval(ct, low, high, c, d)
        real(dp), intent(in) :: ct !! conditional standard deviation
        real(dp), intent(in) :: low !! lower bound
        real(dp), intent(in) :: high !! upper bound
        real(dp), intent(out) :: c !! lower normal CDF
        real(dp), intent(out) :: d !! upper normal CDF

        if (abs(ct) <= cholesky_tol) then
            if (low <= 0.0_dp) then
                c = 0.0_dp
            else
                c = 1.0_dp
            end if
            if (high < 0.0_dp) then
                d = 0.0_dp
            else
                d = 1.0_dp
            end if
        else
            c = normal_cdf(low / ct)
            d = normal_cdf(high / ct)
        end if
    end subroutine first_normal_interval

    pure function finite_normal_ppf(p) result(x)
        real(dp), intent(in) :: p !! normal probability, nominally in [0,1]
        real(dp) :: x
        real(dp) :: q

        q = max(0.0_dp, min(1.0_dp, p))
        x = normal_ppf(q)
        if (.not. ieee_is_finite(x)) then
            if (q <= 0.5_dp) then
                x = -1.0e100_dp
            else
                x = 1.0e100_dp
            end if
        end if
    end function finite_normal_ppf

    pure function clipped_normal_cdf(x) result(p)
        real(dp), intent(in) :: x !! standardized normal argument
        real(dp) :: p

        if (x < -9.0_dp) then
            p = 0.0_dp
        else if (x < 9.0_dp) then
            p = normal_cdf(x)
        else
            p = 1.0_dp
        end if
    end function clipped_normal_cdf

    pure function radical_inverse(index, base) result(value)
        integer, intent(in) :: index !! positive low-discrepancy point index
        integer, intent(in) :: base !! prime Halton base
        real(dp) :: value

        integer :: k
        real(dp) :: factor

        k = index
        factor = 1.0_dp / real(base, dp)
        value = 0.0_dp
        do while (k > 0)
            value = value + factor * real(mod(k, base), dp)
            k = k / base
            factor = factor / real(base, dp)
        end do
    end function radical_inverse

    pure integer function nth_prime(n) result(prime)
        integer, intent(in) :: n !! one-based prime index

        integer :: candidate
        integer :: count

        candidate = 1
        count = 0
        do while (count < n)
            candidate = candidate + 1
            if (is_prime(candidate)) count = count + 1
        end do
        prime = candidate
    end function nth_prime

    pure logical function is_prime(n) result(prime)
        integer, intent(in) :: n !! positive integer candidate

        integer :: divisor

        if (n < 2) then
            prime = .false.
            return
        end if
        if (n == 2) then
            prime = .true.
            return
        end if
        if (mod(n, 2) == 0) then
            prime = .false.
            return
        end if
        divisor = 3
        do while (divisor * divisor <= n)
            if (mod(n, divisor) == 0) then
                prime = .false.
                return
            end if
            divisor = divisor + 2
        end do
        prime = .true.
    end function is_prime

    pure function three_standard_error(values, mean_value) result(error)
        real(dp), intent(in) :: values(:) !! independently shifted batch estimates
        real(dp), intent(in) :: mean_value !! arithmetic mean of batch estimates
        real(dp) :: error

        integer :: n
        real(dp) :: ss

        n = size(values)
        if (n <= 1) then
            error = 1.0_dp
            return
        end if
        ss = sum((values - mean_value) ** 2)
        error = 3.0_dp * sqrt(max(0.0_dp, ss / real(n * (n - 1), dp)))
    end function three_standard_error

    pure integer function default_limit(n, multiplier, requested) result(limit)
        integer, intent(in) :: n !! distribution dimension
        integer, intent(in) :: multiplier !! SciPy-compatible default multiplier
        integer, intent(in), optional :: requested !! caller-supplied maximum

        integer :: default_value

        if (n > huge(1) / max(1, multiplier)) then
            default_value = huge(1)
        else
            default_value = max(1, n * multiplier)
        end if
        limit = default_value
        if (present(requested)) then
            if (requested > 0) limit = requested
        end if
        limit = max(default_batches, limit)
    end function default_limit

    pure subroutine swap_rows_prefix(a, row1, row2, last_col)
        real(dp), intent(inout) :: a(:, :) !! matrix whose row prefixes are swapped
        integer, intent(in) :: row1 !! first row
        integer, intent(in) :: row2 !! second row
        integer, intent(in) :: last_col !! last prefix column

        integer :: j
        real(dp) :: tmp

        do j = 1, last_col
            tmp = a(row1, j)
            a(row1, j) = a(row2, j)
            a(row2, j) = tmp
        end do
    end subroutine swap_rows_prefix

    pure subroutine swap_column_tail(a, col1, col2, first_row, last_row)
        real(dp), intent(inout) :: a(:, :) !! matrix whose column tails are swapped
        integer, intent(in) :: col1 !! first column
        integer, intent(in) :: col2 !! second column
        integer, intent(in) :: first_row !! first row to swap
        integer, intent(in) :: last_row !! final row to swap

        integer :: i
        real(dp) :: tmp

        do i = first_row, last_row
            tmp = a(i, col1)
            a(i, col1) = a(i, col2)
            a(i, col2) = tmp
        end do
    end subroutine swap_column_tail

    pure subroutine swap_cross_segment(a, im, k)
        real(dp), intent(inout) :: a(:, :) !! lower-triangle work matrix
        integer, intent(in) :: im !! pivot row/column moved to k
        integer, intent(in) :: k !! current factorization index

        integer :: j
        real(dp) :: tmp

        do j = k + 1, im - 1
            tmp = a(j, k)
            a(j, k) = a(im, j)
            a(im, j) = tmp
        end do
    end subroutine swap_cross_segment

end module scifort_multivariate_integration
