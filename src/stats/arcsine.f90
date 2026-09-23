! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Arcsine distribution. In standardized coordinates z = (x - loc) / scale,
! support is [0, 1] and f(z) = 1 / (pi*sqrt(z*(1-z))).
! This matches scipy.stats.arcsine.

module scifort_arcsine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_log_two, scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: arcsine_cdf
    public :: arcsine_isf
    public :: arcsine_logcdf
    public :: arcsine_logpdf
    public :: arcsine_logsf
    public :: arcsine_pdf
    public :: arcsine_ppf
    public :: arcsine_sf

contains

    pure elemental function arcsine_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
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
        if (z < 0.0_dp .or. z > 1.0_dp) then
            y = 0.0_dp
        else if (z == 0.0_dp .or. z == 1.0_dp) then
            y = positive_infinity(x)
        else
            y = 1.0_dp / (scifort_pi * sigma * sqrt(z * (1.0_dp - z)))
        end if
    end function arcsine_pdf

    pure elemental function arcsine_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
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
        if (z < 0.0_dp .or. z > 1.0_dp) then
            y = negative_infinity(x)
        else if (z == 0.0_dp .or. z == 1.0_dp) then
            y = positive_infinity(x)
        else
            y = -log(scifort_pi) - log(sigma) - &
                0.5_dp * (log(z) + log1p_safe(-z))
        end if
    end function arcsine_logpdf

    pure elemental function arcsine_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
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
        else if (z >= 1.0_dp) then
            y = 1.0_dp
        else
            y = (2.0_dp / scifort_pi) * asin(sqrt(z))
        end if
    end function arcsine_cdf

    pure elemental function arcsine_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
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
        else if (z >= 1.0_dp) then
            y = 0.0_dp
        else
            y = (2.0_dp / scifort_pi) * asin(sqrt(1.0_dp - z))
        end if
    end function arcsine_sf

    pure elemental function arcsine_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
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
            y = negative_infinity(x)
        else if (z >= 1.0_dp) then
            y = 0.0_dp
        else
            y = scifort_log_two - log(scifort_pi) + log(asin(sqrt(z)))
        end if
    end function arcsine_logcdf

    pure elemental function arcsine_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
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
        else if (z >= 1.0_dp) then
            y = negative_infinity(x)
        else
            y = scifort_log_two - log(scifort_pi) + log(asin(sqrt(1.0_dp - z)))
        end if
    end function arcsine_logsf

    pure elemental function arcsine_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: s

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = mu + sigma
        else
            s = sin(0.5_dp * scifort_pi * p)
            x = mu + sigma * s * s
        end if
    end function arcsine_ppf

    pure elemental function arcsine_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: s

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu + sigma
        else if (p == 1.0_dp) then
            x = mu
        else
            s = sin(0.5_dp * scifort_pi * p)
            x = mu + sigma * (1.0_dp - s * s)
        end if
    end function arcsine_isf

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

end module scifort_arcsine
