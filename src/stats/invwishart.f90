! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_invwishart
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_matrix_distribution_helpers, only : log_multivariate_gamma, matrix_status_success, &
        multivariate_digamma_sum, solve_lower_matrix, solve_spd_matrix, spd_factor
    use scifort_math, only : quiet_nan
    use scifort_random, only : rng_state
    use scifort_wishart, only : wishart_rvs, wishart_status_success
    implicit none
    private

    integer, parameter, public :: invwishart_status_success = 0
    integer, parameter, public :: invwishart_status_invalid_shape = 1
    integer, parameter, public :: invwishart_status_invalid_parameter = 2
    integer, parameter, public :: invwishart_status_linalg_failure = 3

    public :: invwishart_entropy
    public :: invwishart_logpdf
    public :: invwishart_mean
    public :: invwishart_mode
    public :: invwishart_pdf
    public :: invwishart_rvs
    public :: invwishart_rvs_array
    public :: invwishart_var

contains

    function invwishart_logpdf(x, df, scale) result(y)
        real(dp), intent(in) :: x(:, :) !! SPD matrix quantile
        real(dp), intent(in) :: df !! degrees of freedom; must exceed dimension-1
        real(dp), intent(in) :: scale(:, :) !! SPD scale matrix; lower triangle defines symmetry
        real(dp) :: y

        integer :: info
        integer :: p
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: cx(:, :)
        real(dp), allocatable :: cs(:, :)
        real(dp) :: logdet_scale
        real(dp) :: logdet_x
        real(dp) :: trace_term

        p = size(scale,1)
        if (.not. valid_shapes(x, scale) .or. .not. valid_df(df, p)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        allocate(a(p,p), cx(p,p), cs(p,p))
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
        ! tr(scale * x^{-1}) = || L_x^{-1} L_scale ||_F^2.
        call solve_lower_matrix(cx, cs, a, info)
        if (info /= matrix_status_success) then
            y = quiet_nan(0.0_dp)
            return
        end if
        trace_term = sum(a * a)
        y = 0.5_dp * df * logdet_scale - 0.5_dp * trace_term - &
            0.5_dp * df * real(p,dp) * log(2.0_dp) - &
            0.5_dp * (df + real(p + 1,dp)) * logdet_x - &
            log_multivariate_gamma(0.5_dp * df, p)
    end function invwishart_logpdf

    function invwishart_pdf(x, df, scale) result(y)
        real(dp), intent(in) :: x(:, :) !! SPD matrix quantile
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in) :: scale(:, :) !! SPD scale matrix
        real(dp) :: y

        y = exp(invwishart_logpdf(x, df, scale))
    end function invwishart_pdf

    function invwishart_mean(df, scale) result(mu)
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in) :: scale(:, :) !! SPD scale matrix
        real(dp) :: mu(size(scale,1), size(scale,2))

        integer :: p

        p = size(scale,1)
        if (.not. valid_scale_and_df(scale, df) .or. df <= real(p + 1,dp)) then
            mu = quiet_nan(0.0_dp)
        else
            mu = scale / (df - real(p + 1,dp))
        end if
    end function invwishart_mean

    function invwishart_mode(df, scale) result(mode)
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in) :: scale(:, :) !! SPD scale matrix
        real(dp) :: mode(size(scale,1), size(scale,2))

        integer :: p

        p = size(scale,1)
        if (.not. valid_scale_and_df(scale, df)) then
            mode = quiet_nan(0.0_dp)
        else
            mode = scale / (df + real(p + 1,dp))
        end if
    end function invwishart_mode

    function invwishart_var(df, scale) result(v)
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in) :: scale(:, :) !! SPD scale matrix
        real(dp) :: v(size(scale,1), size(scale,2))

        integer :: i
        integer :: j
        integer :: p
        real(dp) :: denom
        real(dp) :: d

        p = size(scale,1)
        if (.not. valid_scale_and_df(scale, df) .or. df <= real(p + 3,dp)) then
            v = quiet_nan(0.0_dp)
            return
        end if
        d = df - real(p,dp)
        denom = d * (d - 1.0_dp)**2 * (d - 3.0_dp)
        do j = 1, p
            do i = 1, p
                v(i,j) = ((d + 1.0_dp) * scale(i,j)**2 + &
                    (d - 1.0_dp) * scale(i,i) * scale(j,j)) / denom
            end do
        end do
    end function invwishart_var

    function invwishart_entropy(df, scale) result(y)
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in) :: scale(:, :) !! SPD scale matrix
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
        y = log_multivariate_gamma(0.5_dp * df, p) + &
            0.5_dp * df * real(p,dp) + &
            0.5_dp * real(p + 1,dp) * (logdet_scale - log(2.0_dp)) - &
            0.5_dp * (df + real(p + 1,dp)) * multivariate_digamma_sum(0.5_dp * df, p)
    end function invwishart_entropy

    subroutine invwishart_rvs(state, df, scale, sample, status)
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in) :: scale(:, :) !! SPD scale matrix
        real(dp), intent(out) :: sample(:, :) !! random SPD matrix
        integer, intent(out), optional :: status !! invwishart_status_* result

        integer :: i
        integer :: info
        integer :: p
        real(dp), allocatable :: cscale(:, :)
        real(dp), allocatable :: cw(:, :)
        real(dp), allocatable :: ident(:, :)
        real(dp), allocatable :: invscale(:, :)
        real(dp), allocatable :: w(:, :)
        real(dp) :: dummy

        p = size(scale,1)
        sample = 0.0_dp
        if (p < 1 .or. size(scale,2) /= p .or. size(sample,1) /= p .or. size(sample,2) /= p) then
            if (present(status)) status = invwishart_status_invalid_shape
            return
        end if
        if (.not. valid_df(df,p) .or. .not. all(ieee_is_finite(scale))) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = invwishart_status_invalid_parameter
            return
        end if
        allocate(cscale(p,p), cw(p,p), ident(p,p), invscale(p,p), w(p,p))
        call spd_factor(scale, cscale, dummy, info)
        if (info /= matrix_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = invwishart_status_linalg_failure
            return
        end if
        ident = 0.0_dp
        do i = 1, p
            ident(i,i) = 1.0_dp
        end do
        call solve_spd_matrix(cscale, ident, invscale, info)
        if (info /= matrix_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = invwishart_status_linalg_failure
            return
        end if
        call wishart_rvs(state, df, invscale, w, info)
        if (info /= wishart_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = invwishart_status_linalg_failure
            return
        end if
        call spd_factor(w, cw, dummy, info)
        if (info /= matrix_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = invwishart_status_linalg_failure
            return
        end if
        call solve_spd_matrix(cw, ident, sample, info)
        if (info /= matrix_status_success) then
            sample = quiet_nan(0.0_dp)
            if (present(status)) status = invwishart_status_linalg_failure
            return
        end if
        if (present(status)) status = invwishart_status_success
    end subroutine invwishart_rvs

    subroutine invwishart_rvs_array(state, df, scale, samples, status)
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in) :: scale(:, :) !! SPD scale matrix
        real(dp), intent(out) :: samples(:, :, :) !! random matrices indexed by third dimension
        integer, intent(out), optional :: status !! invwishart_status_* result

        integer :: info
        integer :: k

        if (size(samples,1) /= size(scale,1) .or. size(samples,2) /= size(scale,1)) then
            samples = quiet_nan(0.0_dp)
            if (present(status)) status = invwishart_status_invalid_shape
            return
        end if
        do k = 1, size(samples,3)
            call invwishart_rvs(state, df, scale, samples(:,:,k), info)
            if (info /= invwishart_status_success) then
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = invwishart_status_success
    end subroutine invwishart_rvs_array

    pure logical function valid_shapes(x, scale) result(valid)
        real(dp), intent(in) :: x(:, :) !! matrix quantile
        real(dp), intent(in) :: scale(:, :) !! scale matrix

        valid = size(scale,1) >= 1 .and. size(scale,2) == size(scale,1) .and. &
            size(x,1) == size(scale,1) .and. size(x,2) == size(scale,1) .and. &
            all(ieee_is_finite(x)) .and. all(ieee_is_finite(scale))
    end function valid_shapes

    pure logical function valid_df(df, p) result(valid)
        real(dp), intent(in) :: df !! degrees of freedom
        integer, intent(in) :: p !! matrix dimension

        valid = ieee_is_finite(df) .and. p >= 1 .and. df > real(p - 1,dp)
    end function valid_df

    pure logical function valid_scale_and_df(scale, df) result(valid)
        real(dp), intent(in) :: scale(:, :) !! scale matrix
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

end module scifort_invwishart
