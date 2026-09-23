! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Lomax (Pareto type II) distribution with shape c > 0, lower endpoint loc,
! and scale > 0. With z = (x - loc) / scale >= 0,
! f(x) = (c / scale) * (1 + z)**(-c - 1). This matches scipy.stats.lomax.

module scifort_lomax
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: lomax_cdf
    public :: lomax_isf
    public :: lomax_logcdf
    public :: lomax_logpdf
    public :: lomax_logsf
    public :: lomax_pdf
    public :: lomax_ppf
    public :: lomax_sf

contains

    pure elemental function lomax_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp) then
                y = 0.0_dp
            else
                y = exp(log(c) - log(sigma) - (c + 1.0_dp) * log1p_safe(z))
            end if
        end if
    end function lomax_pdf

    pure elemental function lomax_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp) then
                y = negative_infinity(x)
            else
                y = log(c) - log(sigma) - (c + 1.0_dp) * log1p_safe(z)
            end if
        end if
    end function lomax_logpdf

    pure elemental function lomax_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logsf

        logsf = lomax_logsf(x, c, loc, scale)
        if (ieee_is_nan(logsf)) then
            y = logsf
        else if (logsf == 0.0_dp) then
            y = 0.0_dp
        else if (logsf == negative_infinity(logsf)) then
            y = 1.0_dp
        else
            y = -expm1_safe(logsf)
        end if
    end function lomax_cdf

    pure elemental function lomax_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logsf

        logsf = lomax_logsf(x, c, loc, scale)
        if (ieee_is_nan(logsf)) then
            y = logsf
        else if (logsf == 0.0_dp) then
            y = 1.0_dp
        else if (logsf == negative_infinity(logsf)) then
            y = 0.0_dp
        else
            y = exp(logsf)
        end if
    end function lomax_sf

    pure elemental function lomax_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logsf

        logsf = lomax_logsf(x, c, loc, scale)
        if (ieee_is_nan(logsf)) then
            y = logsf
        else if (logsf == 0.0_dp) then
            y = negative_infinity(x)
        else if (logsf == negative_infinity(logsf)) then
            y = 0.0_dp
        else if (-logsf < scifort_log_two) then
            y = log(-expm1_safe(logsf))
        else
            y = log1p_safe(-exp(logsf))
        end if
    end function lomax_logcdf

    pure elemental function lomax_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp) then
                y = 0.0_dp
            else
                y = -c * log1p_safe(z)
            end if
        end if
    end function lomax_logsf

    pure elemental function lomax_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: exponent
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            exponent = -log1p_safe(-p) / c
            if (exponent > log(huge(1.0_dp))) then
                x = positive_infinity(p)
            else
                z = expm1_safe(exponent)
                if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
                    x = positive_infinity(p)
                else
                    x = mu + sigma * z
                end if
            end if
        end if
    end function lomax_ppf

    pure elemental function lomax_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: exponent
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            exponent = -log(p) / c
            if (exponent > log(huge(1.0_dp))) then
                x = positive_infinity(p)
            else
                z = expm1_safe(exponent)
                if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
                    x = positive_infinity(p)
                else
                    x = mu + sigma * z
                end if
            end if
        end if
    end function lomax_isf

    pure elemental logical function valid_shape(c) result(valid)
        real(dp), intent(in) :: c !! shape parameter to check

        valid = ieee_is_finite(c) .and. c > 0.0_dp
    end function valid_shape

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

end module scifort_lomax
