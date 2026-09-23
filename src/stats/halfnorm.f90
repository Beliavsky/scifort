! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Half-normal distribution with lower support endpoint loc and scale > 0.
! With z = (x - loc) / scale >= 0, f(x) = sqrt(2/pi)*exp(-z**2/2)/scale.
! This matches scipy.stats.halfnorm.

module scifort_halfnorm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_log_two, scifort_sqrt_two
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_isf, normal_logpdf, normal_logsf, normal_pdf
    implicit none
    private

    public :: halfnorm_cdf
    public :: halfnorm_isf
    public :: halfnorm_logcdf
    public :: halfnorm_logpdf
    public :: halfnorm_logsf
    public :: halfnorm_pdf
    public :: halfnorm_ppf
    public :: halfnorm_sf

contains

    pure elemental function halfnorm_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu) then
            y = 0.0_dp
        else
            y = 2.0_dp * normal_pdf(x, mu, sigma)
        end if
    end function halfnorm_pdf

    pure elemental function halfnorm_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu) then
            y = negative_infinity(x)
        else
            y = scifort_log_two + normal_logpdf(x, mu, sigma)
        end if
    end function halfnorm_logpdf

    pure elemental function halfnorm_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else
            y = erf(z / scifort_sqrt_two)
        end if
    end function halfnorm_cdf

    pure elemental function halfnorm_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 1.0_dp
        else
            y = erfc(z / scifort_sqrt_two)
        end if
    end function halfnorm_sf

    pure elemental function halfnorm_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: cdf_value
        real(dp) :: mu
        real(dp) :: sf_value
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = negative_infinity(x)
            return
        end if
        sf_value = erfc(z / scifort_sqrt_two)
        if (sf_value < 0.5_dp) then
            y = log1p_safe(-sf_value)
        else
            cdf_value = erf(z / scifort_sqrt_two)
            y = log(cdf_value)
        end if
    end function halfnorm_logcdf

    pure elemental function halfnorm_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else
            y = scifort_log_two + normal_logsf(x, mu, sigma)
        end if
    end function halfnorm_logsf

    pure elemental function halfnorm_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            x = normal_isf(0.5_dp * (1.0_dp - p), mu, sigma)
        end if
    end function halfnorm_ppf

    pure elemental function halfnorm_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = mu
        else
            x = normal_isf(0.5_dp * p, mu, sigma)
        end if
    end function halfnorm_isf

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when loc is absent
        real(dp), intent(out) :: sigma !! scale, 1 when scale is absent

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_halfnorm
