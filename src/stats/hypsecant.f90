! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Hyperbolic-secant distribution. In standardized coordinates
! z = (x - loc) / scale, f(z) = sech(z) / pi.
! This matches scipy.stats.hypsecant.

module scifort_hypsecant
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_log_two, scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: hypsecant_cdf
    public :: hypsecant_isf
    public :: hypsecant_logcdf
    public :: hypsecant_logpdf
    public :: hypsecant_logsf
    public :: hypsecant_pdf
    public :: hypsecant_ppf
    public :: hypsecant_sf

contains

    pure elemental function hypsecant_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            y = exp(-log_cosh(z)) / (scifort_pi * sigma)
        end if
    end function hypsecant_pdf

    pure elemental function hypsecant_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            y = -log(scifort_pi) - log_cosh(z) - log(sigma)
        end if
    end function hypsecant_logpdf

    pure elemental function hypsecant_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = (2.0_dp / scifort_pi) * atan(exp(z))
            else
                y = 1.0_dp - (2.0_dp / scifort_pi) * atan(exp(-z))
            end if
        end if
    end function hypsecant_cdf

    pure elemental function hypsecant_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z >= 0.0_dp) then
                y = (2.0_dp / scifort_pi) * atan(exp(-z))
            else
                y = 1.0_dp - (2.0_dp / scifort_pi) * atan(exp(z))
            end if
        end if
    end function hypsecant_sf

    pure elemental function hypsecant_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: tail
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < -18.0_dp) then
                y = log(2.0_dp / scifort_pi) + z
            else if (z <= 0.0_dp) then
                y = log((2.0_dp / scifort_pi) * atan(exp(z)))
            else
                tail = (2.0_dp / scifort_pi) * atan(exp(-z))
                y = log1p_safe(-tail)
            end if
        end if
    end function hypsecant_logcdf

    pure elemental function hypsecant_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: tail
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z > 18.0_dp) then
                y = log(2.0_dp / scifort_pi) - z
            else if (z >= 0.0_dp) then
                y = log((2.0_dp / scifort_pi) * atan(exp(-z)))
            else
                tail = (2.0_dp / scifort_pi) * atan(exp(z))
                y = log1p_safe(-tail)
            end if
        end if
    end function hypsecant_logsf

    pure elemental function hypsecant_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else if (p <= 0.5_dp) then
            z = log(tan(0.5_dp * scifort_pi * p))
            x = mu + sigma * z
        else
            z = -log(tan(0.5_dp * scifort_pi * (1.0_dp - p)))
            x = mu + sigma * z
        end if
    end function hypsecant_ppf

    pure elemental function hypsecant_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else if (p <= 0.5_dp) then
            z = -log(tan(0.5_dp * scifort_pi * p))
            x = mu + sigma * z
        else
            z = log(tan(0.5_dp * scifort_pi * (1.0_dp - p)))
            x = mu + sigma * z
        end if
    end function hypsecant_isf

    pure elemental function log_cosh(x) result(y)
        real(dp), intent(in) :: x !! real argument of log(cosh(x))
        real(dp) :: y

        real(dp) :: a

        a = abs(x)
        if (a > 0.5_dp * log(huge(1.0_dp))) then
            y = a - scifort_log_two
        else
            y = a + log1p_safe(exp(-2.0_dp * a)) - scifort_log_two
        end if
    end function log_cosh

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

end module scifort_hypsecant
