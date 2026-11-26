! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_wishart
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_chi2, only : chi2_ppf
    use scifort_kinds, only : dp
    use scifort_matrix_distribution_helpers, only : log_multivariate_gamma, matrix_status_success, &
        multivariate_digamma_sum, solve_spd_matrix, spd_factor
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_ppf
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: wishart_status_success = 0
    integer, parameter, public :: wishart_status_invalid_shape = 1
    integer, parameter, public :: wishart_status_invalid_parameter = 2
    integer, parameter, public :: wishart_status_linalg_failure = 3

    public :: wishart_entropy
    public :: wishart_logpdf
    public :: wishart_mean
    public :: wishart_mode
    public :: wishart_pdf
    public :: wishart_rvs
    public :: wishart_rvs_array
    public :: wishart_var

contains

    function wishart_logpdf(x, df, scale) result(y)
        !! Evaluate the log density at one positive-definite matrix.
        !! @dimension p = size(scale,1)
        !! @require p >= 1
        !! Finite symmetric inputs are required; factorization uses lower triangles.
        !! The trace solve also reads the full x matrix; do not omit its upper triangle.
        real(dp), intent(in) :: x(:, :) !! shape (p,p); finite symmetric positive-definite evaluation matrix
        real(dp), intent(in) :: df !! finite degrees of freedom, > p-1
        real(dp), intent(in) :: scale(:, :) !! shape (p,p); finite symmetric positive-definite scale matrix
        real(dp) :: y

        integer :: info
        integer :: p
        real(dp), allocatable :: cx(:, :)
        real(dp), allocatable :: cs(:, :)
        real(dp), allocatable :: solved(:, :)
        real(dp) :: logdet_scale
        real(dp) :: logdet_x
        real(dp) :: trace_term

        p = size(scale,1)
        if (.not. valid_shapes(x, scale) .or. .not. valid_df(df, p)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(cx(p,p), cs(p,p), solved(p,p))
        call spd_factor(scale, cs, logdet_scale, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call spd_factor(x, cx, logdet_x, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        call solve_spd_matrix(cs, x, solved, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        trace_term = trace_matrix(solved)
        y = 0.5_dp * (df - real(p,dp) - 1.0_dp) * logdet_x - 0.5_dp * trace_term - &
            0.5_dp * df * real(p,dp) * log(2.0_dp) - 0.5_dp * df * logdet_scale - &
            log_multivariate_gamma(0.5_dp * df, p)
    end function wishart_logpdf

    function wishart_pdf(x, df, scale) result(y)
        !! Evaluate the density at one positive-definite matrix.
        !! @dimension p = size(scale,1)
        !! @require p >= 1
        !! Finite symmetric inputs are required; factorization uses lower triangles.
        !! The trace solve also reads the full x matrix; do not omit its upper triangle.
        real(dp), intent(in) :: x(:, :) !! shape (p,p); finite symmetric positive-definite evaluation matrix
        real(dp), intent(in) :: df !! finite degrees of freedom, > p-1
        real(dp), intent(in) :: scale(:, :) !! shape (p,p); finite symmetric positive-definite scale matrix
        real(dp) :: y

        y = exp(wishart_logpdf(x, df, scale))
    end function wishart_pdf

    function wishart_mean(df, scale) result(mu)
        !! Return the entrywise mean matrix.
        !! @dimension p = size(scale,1)
        !! @require p >= 1
        !! Scale validation factors its lower triangle; the formula reads all entries.
        real(dp), intent(in) :: df !! finite degrees of freedom, > p-1
        real(dp), intent(in) :: scale(:, :) !! shape (p,p); finite symmetric positive-definite scale matrix
        real(dp) :: mu(size(scale,1), size(scale,2)) !! shape (p,p); entrywise means

        if (.not. valid_scale_and_df(scale, df)) then
            mu = quiet_nan(0.0_dp)
        else
            mu = df * scale
        end if
    end function wishart_mean

    function wishart_mode(df, scale) result(mode)
        !! Return the mode matrix.
        !! @dimension p = size(scale,1)
        !! @require p >= 1
        !! The mode is NaN for df < p+1; at df = p+1 it is the zero matrix.
        !! Scale validation factors its lower triangle; the formula reads all entries.
        real(dp), intent(in) :: df !! finite degrees of freedom, > p-1; mode requires df >= p+1
        real(dp), intent(in) :: scale(:, :) !! shape (p,p); finite symmetric positive-definite scale matrix
        real(dp) :: mode(size(scale,1), size(scale,2)) !! shape (p,p); mode matrix

        integer :: p

        p = size(scale,1)
        if (.not. valid_scale_and_df(scale, df) .or. df < real(p + 1,dp)) then
            mode = quiet_nan(0.0_dp)
        else
            mode = (df - real(p + 1,dp)) * scale
        end if
    end function wishart_mode

    function wishart_var(df, scale) result(v)
        !! Return the entrywise variances, not a covariance tensor.
        !! @dimension p = size(scale,1)
        !! @require p >= 1
        !! Scale validation factors its lower triangle; the formula reads all entries.
        real(dp), intent(in) :: df !! finite degrees of freedom, > p-1
        real(dp), intent(in) :: scale(:, :) !! shape (p,p); finite symmetric positive-definite scale matrix
        real(dp) :: v(size(scale,1), size(scale,2)) !! shape (p,p); entrywise variances

        integer :: i
        integer :: j

        if (.not. valid_scale_and_df(scale, df)) then
            v = quiet_nan(0.0_dp)
            return
        end if
        do j = 1, size(scale,2)
            do i = 1, size(scale,1)
                v(i,j) = df * (scale(i,j)**2 + scale(i,i) * scale(j,j))
            end do
        end do
    end function wishart_var

    function wishart_entropy(df, scale) result(y)
        !! Return the differential entropy.
        !! @dimension p = size(scale,1)
        !! @require p >= 1
        real(dp), intent(in) :: df !! finite degrees of freedom, > p-1
        real(dp), intent(in) :: scale(:, :) !! shape (p,p); finite symmetric positive-definite scale matrix
        real(dp) :: y

        integer :: info
        integer :: p
        real(dp), allocatable :: c(:, :)
        real(dp) :: logdet_scale

        p = size(scale,1)
        if (.not. valid_scale_and_df(scale, df)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(c(p,p))
        call spd_factor(scale, c, logdet_scale, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        y = 0.5_dp * real(p + 1,dp) * logdet_scale + &
            0.5_dp * real(p * (p + 1),dp) * log(2.0_dp) + &
            log_multivariate_gamma(0.5_dp * df, p) - &
            0.5_dp * (df - real(p + 1,dp)) * multivariate_digamma_sum(0.5_dp * df, p) + &
            0.5_dp * df * real(p,dp)
    end function wishart_entropy

    subroutine wishart_rvs(state, df, scale, sample, status)
        !! Draw one random positive-definite matrix.
        !! p = size(scale,1), p >= 1; scale and each draw have shape (p,p).
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(in) :: df !! finite degrees of freedom, > p-1
        real(dp), intent(in) :: scale(:, :) !! shape (p,p); finite SPD scale; lower triangle defines symmetry
        real(dp), intent(out) :: sample(:, :) !! shape (p,p); random SPD matrix
        integer, intent(out), optional :: status !! wishart_status_* result

        integer :: i
        integer :: info
        integer :: j
        integer :: p
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: c(:, :)
        real(dp), allocatable :: ca(:, :)
        real(dp) :: dummy

        p = size(scale,1)
        sample = 0.0_dp
        if (size(scale,2) /= p .or. size(sample,1) /= p .or. size(sample,2) /= p .or. p < 1) then
            if (present(status)) status = wishart_status_invalid_shape
            return
        end if
        if (.not. valid_df(df,p) .or. .not. all(ieee_is_finite(scale))) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = wishart_status_invalid_parameter
            return
        end if
        allocate(a(p,p), c(p,p), ca(p,p))
        call spd_factor(scale, c, dummy, info)
        if (info /= matrix_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = wishart_status_linalg_failure
            return
        end if
        a = 0.0_dp
        do i = 1, p
            a(i,i) = sqrt(chi2_ppf(rng_uniform(state), df - real(i,dp) + 1.0_dp))
            do j = 1, i - 1
                a(i,j) = normal_ppf(rng_uniform(state))
            end do
        end do
        ca = matmul(c, a)
        sample = matmul(ca, transpose(ca))
        if (present(status)) status = wishart_status_success
    end subroutine wishart_rvs

    subroutine wishart_rvs_array(state, df, scale, samples, status)
        !! Draw random matrices into the third dimension.
        !! p = size(scale,1), p >= 1; scale and each draw have shape (p,p).
        !! b = size(samples,3), b >= 0; empty batches are permitted.
        !! With b = 0 only the leading sample extents are checked; df and scale
        !! squareness, finiteness and positive-definiteness are not validated.
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(in) :: df !! finite degrees of freedom, > p-1
        real(dp), intent(in) :: scale(:, :) !! shape (p,p); finite SPD scale; lower triangle defines symmetry
        real(dp), intent(out) :: samples(:, :, :) !! shape (p,p,b); random SPD matrices; b may be zero
        integer, intent(out), optional :: status !! wishart_status_* result

        integer :: info
        integer :: k

        if (size(samples,1) /= size(scale,1) .or. size(samples,2) /= size(scale,1)) then
            samples = quiet_nan(0.0_dp)
            if (present(status)) status = wishart_status_invalid_shape
            return
        end if
        do k = 1, size(samples,3)
            call wishart_rvs(state, df, scale, samples(:,:,k), info)
            if (info /= wishart_status_success) then
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = wishart_status_success
    end subroutine wishart_rvs_array

    pure logical function valid_shapes(x, scale) result(valid)
        !! Check finite equal nonempty square shapes; invalid inputs return false.
        !! Candidate matrix extents are unrestricted, including zero extents.
        real(dp), intent(in) :: x(:, :) !! candidate matrix; arbitrary extents to be checked
        real(dp), intent(in) :: scale(:, :) !! candidate scale; arbitrary extents to be checked

        valid = size(scale,1) >= 1 .and. size(scale,2) == size(scale,1) .and. &
            size(x,1) == size(scale,1) .and. size(x,2) == size(scale,1) .and. &
            all(ieee_is_finite(x)) .and. all(ieee_is_finite(scale))
    end function valid_shapes

    pure logical function valid_df(df, p) result(valid)
        !! Check finite degrees of freedom for a positive matrix dimension.
        real(dp), intent(in) :: df !! degrees of freedom
        integer, intent(in) :: p !! matrix dimension

        valid = ieee_is_finite(df) .and. p >= 1 .and. df > real(p - 1,dp)
    end function valid_df

    pure logical function valid_scale_and_df(scale, df) result(valid)
        !! Check a finite nonempty square scale, degrees of freedom and SPD factor.
        !! Candidate scale extents are unrestricted; invalid inputs return false.
        real(dp), intent(in) :: scale(:, :) !! candidate scale; arbitrary extents to be checked
        real(dp), intent(in) :: df !! degrees of freedom

        integer :: info
        integer :: p
        real(dp), allocatable :: c(:, :)
        real(dp) :: dummy

        p = size(scale,1)
        valid = p >= 1 .and. size(scale,2) == p .and. all(ieee_is_finite(scale)) .and. valid_df(df,p)
        if (.not. valid) return
        allocate(c(p,p))
        call spd_factor(scale, c, dummy, info)
        valid = info == matrix_status_success
    end function valid_scale_and_df

    pure function trace_matrix(a) result(y)
        !! Sum the diagonal through min(r,c), where r = size(a,1), c = size(a,2).
        !! Rectangular matrices and empty extents are permitted; an empty trace is zero.
        real(dp), intent(in) :: a(:, :) !! shape (r,c); matrix whose available diagonal is summed
        real(dp) :: y

        integer :: i

        y = 0.0_dp
        do i = 1, min(size(a,1),size(a,2))
            y = y + a(i,i)
        end do
    end function trace_matrix

end module scifort_wishart
