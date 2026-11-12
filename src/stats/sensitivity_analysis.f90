! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2001-2002 Enthought, Inc. 2003-2026, SciPy Developers.
! Copyright (c) 2026 SciFort contributors
!
! Sobol sensitivity indices matching scipy.stats.sobol_indices 1.17 semantics.
! The Saltelli 2010 estimator, A/B sampling layout, mean normalization, and
! constant-output handling follow SciPy's scipy/stats/_sensitivity_analysis.py.
! See CODE_PROVENANCE.md and THIRD_PARTY_LICENSES.md.

module scifort_sensitivity_analysis
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_qmc, only : qmc_sobol, qmc_sobol_init, qmc_status_success
    use scifort_random, only : rng_state
    implicit none
    private

    integer, parameter, public :: sobol_status_success = 0
    integer, parameter, public :: sobol_status_invalid_shape = 1
    integer, parameter, public :: sobol_status_invalid_parameter = 2
    integer, parameter, public :: sobol_status_model_failure = 3
    integer, parameter, public :: sobol_status_qmc_failure = 4

    abstract interface
        pure function sobol_ppf_callback(q, params) result(x)
            import :: dp
            real(dp), intent(in) :: q !! probability in the unit interval
            real(dp), intent(in) :: params(:) !! distribution-specific parameters
            real(dp) :: x
        end function sobol_ppf_callback

        subroutine sobol_model_callback(x, y, status)
            import :: dp
            real(dp), intent(in) :: x(:, :) !! model inputs with shape (d,n)
            real(dp), allocatable, intent(out) :: y(:, :) !! outputs with shape (s,n)
            integer, intent(out) :: status !! zero on successful model evaluation
        end subroutine sobol_model_callback
    end interface

    type, public :: sobol_distribution
        procedure(sobol_ppf_callback), pointer, nopass :: ppf => null() !! inverse CDF callback
        real(dp), allocatable :: params(:) !! parameters supplied to the inverse CDF callback
    end type sobol_distribution

    type, public :: sobol_result
        real(dp), allocatable :: first_order(:, :) !! first-order indices with shape (s,d)
        real(dp), allocatable :: total_order(:, :) !! total-order indices with shape (s,d)
        integer :: status = sobol_status_success !! sobol_status_* result
    end type sobol_result

    interface sobol_indices
        module procedure sobol_indices_values
        module procedure sobol_indices_model
    end interface sobol_indices

    public :: sobol_indices
    public :: sobol_uniform_distribution

contains

    function sobol_indices_values(f_a, f_b, f_ab, status) result(res)
        real(dp), intent(in) :: f_a(:, :) !! outputs at A with shape (s,n)
        real(dp), intent(in) :: f_b(:, :) !! outputs at B with shape (s,n)
        real(dp), intent(in) :: f_ab(:, :, :) !! outputs at AB pages with shape (d,s,n)
        integer, intent(out), optional :: status !! sobol_status_* result
        type(sobol_result) :: res

        integer :: d
        integer :: n
        integer :: s
        integer :: stat
        integer :: j
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:, :)
        real(dp), allocatable :: ab(:, :, :)
        real(dp), allocatable :: center(:)

        s = size(f_a, 1)
        n = size(f_a, 2)
        d = size(f_ab, 1)
        stat = sobol_status_success

        if (s < 1 .or. n < 1 .or. d < 1 .or. size(f_b,1) /= s .or. &
                size(f_b,2) /= n .or. size(f_ab,2) /= s .or. size(f_ab,3) /= n) then
            stat = sobol_status_invalid_shape
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if
        if (.not. is_power_of_two(n) .or. .not. all(ieee_is_finite(f_a)) .or. &
                .not. all(ieee_is_finite(f_b)) .or. .not. all(ieee_is_finite(f_ab))) then
            stat = sobol_status_invalid_parameter
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if

        allocate(a(s,n), b(s,n), ab(d,s,n), center(s))
        a = f_a
        b = f_b
        ab = f_ab
        do j = 1, s
            center(j) = (sum(a(j,:)) + sum(b(j,:))) / real(2*n, dp)
            a(j,:) = a(j,:) - center(j)
            b(j,:) = b(j,:) - center(j)
            ab(:,j,:) = ab(:,j,:) - center(j)
        end do

        call saltelli_2010(a, b, ab, res%first_order, res%total_order, stat)
        res%status = stat
        if (present(status)) status = stat
    end function sobol_indices_values

    function sobol_indices_model(func, n, dists, state, status) result(res)
        procedure(sobol_model_callback) :: func !! vectorized model evaluated on shape (d,n)
        integer, intent(in) :: n !! number of A and B samples; must be a positive power of two
        type(sobol_distribution), intent(in) :: dists(:) !! independent input distributions
        type(rng_state), intent(inout) :: state !! explicit RNG state used to scramble Sobol points
        integer, intent(out), optional :: status !! sobol_status_* result
        type(sobol_result) :: res

        type(qmc_sobol) :: engine
        integer :: d
        integer :: i
        integer :: j
        integer :: qstat
        integer :: s
        integer :: stat
        real(dp), allocatable :: sample(:, :)
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:, :)
        real(dp), allocatable :: ab_flat(:, :)
        real(dp), allocatable :: f_a(:, :)
        real(dp), allocatable :: f_b(:, :)
        real(dp), allocatable :: f_flat(:, :)
        real(dp), allocatable :: f_ab(:, :, :)

        d = size(dists)
        stat = sobol_status_success
        if (d < 1 .or. n < 1 .or. .not. is_power_of_two(n)) then
            stat = sobol_status_invalid_parameter
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if
        do i = 1, d
            if (.not. associated(dists(i)%ppf)) then
                stat = sobol_status_invalid_parameter
                call invalidate_result(res, stat)
                if (present(status)) status = stat
                return
            end if
        end do

        allocate(sample(n,2*d), a(d,n), b(d,n))
        call qmc_sobol_init(engine, 2*d, qstat, scramble=.true., bits=64, state=state)
        if (qstat /= qmc_status_success) then
            stat = sobol_status_qmc_failure
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if
        call engine%random(n, sample, qstat)
        if (qstat /= qmc_status_success) then
            stat = sobol_status_qmc_failure
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if

        do i = 1, d
            do j = 1, n
                a(i,j) = dists(i)%ppf(sample(j,i), dists(i)%params)
                b(i,j) = dists(i)%ppf(sample(j,d+i), dists(i)%params)
            end do
        end do
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) then
            stat = sobol_status_invalid_parameter
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if

        call func(a, f_a, stat)
        if (stat /= 0 .or. .not. allocated(f_a)) then
            stat = sobol_status_model_failure
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if
        if (size(f_a,2) /= n .or. size(f_a,1) < 1 .or. .not. all(ieee_is_finite(f_a))) then
            stat = sobol_status_model_failure
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if
        s = size(f_a,1)

        call func(b, f_b, stat)
        if (stat /= 0 .or. .not. allocated(f_b)) then
            stat = sobol_status_model_failure
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if
        if (size(f_b,1) /= s .or. size(f_b,2) /= n .or. .not. all(ieee_is_finite(f_b))) then
            stat = sobol_status_model_failure
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if

        allocate(ab_flat(d,n*d))
        do i = 1, d
            ab_flat(:,(i-1)*n+1:i*n) = a
            ab_flat(i,(i-1)*n+1:i*n) = b(i,:)
        end do
        call func(ab_flat, f_flat, stat)
        if (stat /= 0 .or. .not. allocated(f_flat)) then
            stat = sobol_status_model_failure
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if
        if (size(f_flat,1) /= s .or. size(f_flat,2) /= n*d .or. &
                .not. all(ieee_is_finite(f_flat))) then
            stat = sobol_status_model_failure
            call invalidate_result(res, stat)
            if (present(status)) status = stat
            return
        end if

        allocate(f_ab(d,s,n))
        do i = 1, d
            f_ab(i,:,:) = f_flat(:,(i-1)*n+1:i*n)
        end do
        res = sobol_indices_values(f_a, f_b, f_ab, stat)
        if (present(status)) status = stat
    end function sobol_indices_model

    function sobol_uniform_distribution(lower, upper) result(dist)
        real(dp), intent(in) :: lower !! finite lower endpoint
        real(dp), intent(in) :: upper !! finite upper endpoint, greater than lower
        type(sobol_distribution) :: dist

        allocate(dist%params(2))
        dist%params = [lower, upper]
        dist%ppf => uniform_ppf
    end function sobol_uniform_distribution

    subroutine saltelli_2010(f_a, f_b, f_ab, first_order, total_order, status)
        real(dp), intent(in) :: f_a(:, :) !! centered outputs at A, shape (s,n)
        real(dp), intent(in) :: f_b(:, :) !! centered outputs at B, shape (s,n)
        real(dp), intent(in) :: f_ab(:, :, :) !! centered outputs at AB pages, shape (d,s,n)
        real(dp), allocatable, intent(out) :: first_order(:, :) !! first-order indices, shape (s,d)
        real(dp), allocatable, intent(out) :: total_order(:, :) !! total-order indices, shape (s,d)
        integer, intent(out) :: status !! sobol_status_* result

        integer :: d
        integer :: i
        integer :: j
        integer :: n
        integer :: s
        real(dp) :: mean_value
        real(dp) :: var

        s = size(f_a,1)
        n = size(f_a,2)
        d = size(f_ab,1)
        allocate(first_order(s,d), total_order(s,d))
        first_order = 0.0_dp
        total_order = 0.0_dp
        status = sobol_status_success

        do j = 1, s
            mean_value = (sum(f_a(j,:)) + sum(f_b(j,:))) / real(2*n, dp)
            var = (sum((f_a(j,:) - mean_value)**2) + &
                sum((f_b(j,:) - mean_value)**2)) / real(2*n, dp)
            if (.not. ieee_is_finite(var) .or. var <= 0.0_dp) cycle
            do i = 1, d
                first_order(j,i) = sum(f_b(j,:) * (f_ab(i,j,:) - f_a(j,:))) / &
                    real(n,dp) / var
                total_order(j,i) = 0.5_dp * sum((f_a(j,:) - f_ab(i,j,:))**2) / &
                    real(n,dp) / var
                if (.not. ieee_is_finite(first_order(j,i))) first_order(j,i) = 0.0_dp
                if (.not. ieee_is_finite(total_order(j,i))) total_order(j,i) = 0.0_dp
            end do
        end do
    end subroutine saltelli_2010

    pure function uniform_ppf(q, params) result(x)
        real(dp), intent(in) :: q !! probability in the unit interval
        real(dp), intent(in) :: params(:) !! [lower, upper] endpoints
        real(dp) :: x

        if (size(params) /= 2 .or. .not. ieee_is_finite(q) .or. q < 0.0_dp .or. q > 1.0_dp .or. &
                .not. all(ieee_is_finite(params)) .or. params(2) <= params(1)) then
            x = quiet_nan(0.0_dp)
            return
        end if
        x = params(1) + q * (params(2) - params(1))
    end function uniform_ppf

    pure logical function is_power_of_two(n) result(ok)
        integer, intent(in) :: n !! candidate positive integer

        ok = n > 0
        if (ok) ok = iand(n, n - 1) == 0
    end function is_power_of_two

    subroutine invalidate_result(res, status)
        type(sobol_result), intent(out) :: res !! result object reset to empty arrays
        integer, intent(in) :: status !! failure code to store

        allocate(res%first_order(0,0), res%total_order(0,0))
        res%status = status
    end subroutine invalidate_result

end module scifort_sensitivity_analysis
