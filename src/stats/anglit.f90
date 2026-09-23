! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Anglit distribution. In standardized coordinates z = (x - loc) / scale,
! support is [-pi/4, pi/4] and f(z) = cos(2*z).
! This matches scipy.stats.anglit.

module scifort_anglit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan, &
        valid_loc_scale
    implicit none
    private

    public :: anglit_cdf
    public :: anglit_isf
    public :: anglit_logcdf
    public :: anglit_logpdf
    public :: anglit_logsf
    public :: anglit_pdf
    public :: anglit_ppf
    public :: anglit_sf

contains

    pure elemental function anglit_pdf(x, loc, scale) result(y)
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
            if (abs(z) >= 0.25_dp * scifort_pi) then
                y = 0.0_dp
            else
                y = cos(2.0_dp * z) / sigma
            end if
        end if
    end function anglit_pdf

    pure elemental function anglit_logpdf(x, loc, scale) result(y)
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
            if (abs(z) >= 0.25_dp * scifort_pi) then
                y = negative_infinity(x)
            else
                y = log(cos(2.0_dp * z)) - log(sigma)
            end if
        end if
    end function anglit_logpdf

    pure elemental function anglit_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: s
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= -0.25_dp * scifort_pi) then
                y = 0.0_dp
            else if (z >= 0.25_dp * scifort_pi) then
                y = 1.0_dp
            else if (z <= 0.0_dp) then
                s = sin(z + 0.25_dp * scifort_pi)
                y = s * s
            else
                s = sin(0.25_dp * scifort_pi - z)
                y = 1.0_dp - s * s
            end if
        end if
    end function anglit_cdf

    pure elemental function anglit_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else
            y = anglit_cdf(2.0_dp * mu - x, mu, sigma)
        end if
    end function anglit_sf

    pure elemental function anglit_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
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
        if (z <= -0.25_dp * scifort_pi) then
            y = negative_infinity(x)
        else if (z >= 0.25_dp * scifort_pi) then
            y = 0.0_dp
        else if (z <= 0.0_dp) then
            cdf_value = anglit_cdf(x, mu, sigma)
            y = log(cdf_value)
        else
            sf_value = anglit_sf(x, mu, sigma)
            y = log1p_safe(-sf_value)
        end if
    end function anglit_logcdf

    pure elemental function anglit_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else
            y = anglit_logcdf(2.0_dp * mu - x, mu, sigma)
        end if
    end function anglit_logsf

    pure elemental function anglit_ppf(p, loc, scale) result(x)
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
        else if (p <= 0.0_dp) then
            x = mu - 0.25_dp * scifort_pi * sigma
        else if (p >= 1.0_dp) then
            x = mu + 0.25_dp * scifort_pi * sigma
        else if (p <= 0.5_dp) then
            z = asin(sqrt(p)) - 0.25_dp * scifort_pi
            x = mu + sigma * z
        else
            z = asin(sqrt(1.0_dp - p)) - 0.25_dp * scifort_pi
            x = mu - sigma * z
        end if
    end function anglit_ppf

    pure elemental function anglit_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else
            x = 2.0_dp * mu - anglit_ppf(p, mu, sigma)
        end if
    end function anglit_isf

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

end module scifort_anglit
