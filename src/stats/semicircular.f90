! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Semicircular distribution. In standardized coordinates z = (x - loc) / scale,
! support is [-1, 1] and f(z) = 2*sqrt(1-z**2)/pi.
! This matches scipy.stats.semicircular.

module scifort_semicircular
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: semicircular_cdf
    public :: semicircular_isf
    public :: semicircular_logcdf
    public :: semicircular_logpdf
    public :: semicircular_logsf
    public :: semicircular_pdf
    public :: semicircular_ppf
    public :: semicircular_sf

contains

    pure elemental function semicircular_pdf(x, loc, scale) result(y)
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
            if (abs(z) >= 1.0_dp) then
                y = 0.0_dp
            else
                y = (2.0_dp / (scifort_pi * sigma)) * &
                    sqrt((1.0_dp - abs(z)) * (1.0_dp + abs(z)))
            end if
        end if
    end function semicircular_pdf

    pure elemental function semicircular_logpdf(x, loc, scale) result(y)
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
            if (abs(z) >= 1.0_dp) then
                y = negative_infinity(x)
            else
                y = log(2.0_dp / scifort_pi) - log(sigma) + &
                    0.5_dp * (log1p_safe(-abs(z)) + log1p_safe(abs(z)))
            end if
        end if
    end function semicircular_logpdf

    pure elemental function semicircular_cdf(x, loc, scale) result(y)
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
            y = standard_cdf(z)
        end if
    end function semicircular_cdf

    pure elemental function semicircular_sf(x, loc, scale) result(y)
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
            y = standard_cdf(-z)
        end if
    end function semicircular_sf

    pure elemental function semicircular_logcdf(x, loc, scale) result(y)
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
        if (z <= -1.0_dp) then
            y = negative_infinity(x)
        else if (z >= 1.0_dp) then
            y = 0.0_dp
        else if (z <= 0.0_dp) then
            cdf_value = standard_cdf(z)
            y = log(cdf_value)
        else
            sf_value = standard_cdf(-z)
            y = log1p_safe(-sf_value)
        end if
    end function semicircular_logcdf

    pure elemental function semicircular_logsf(x, loc, scale) result(y)
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
            y = semicircular_logcdf(2.0_dp * mu - x, mu, sigma)
        end if
    end function semicircular_logsf

    pure elemental function semicircular_ppf(p, loc, scale) result(x)
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
            x = mu - sigma
        else if (p >= 1.0_dp) then
            x = mu + sigma
        else if (p < 0.5_dp) then
            z = lower_quantile(p)
            x = mu + sigma * z
        else
            z = lower_quantile(1.0_dp - p)
            x = mu - sigma * z
        end if
    end function semicircular_ppf

    pure elemental function semicircular_isf(p, loc, scale) result(x)
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
            x = 2.0_dp * mu - semicircular_ppf(p, mu, sigma)
        end if
    end function semicircular_isf

    pure elemental function standard_cdf(z) result(y)
        real(dp), intent(in) :: z !! standardized variate
        real(dp) :: y

        real(dp) :: theta

        if (z <= -1.0_dp) then
            y = 0.0_dp
        else if (z >= 1.0_dp) then
            y = 1.0_dp
        else if (z < 0.0_dp) then
            theta = acos(-z)
            y = theta_minus_sin_cos(theta) / scifort_pi
        else
            theta = acos(z)
            y = 1.0_dp - theta_minus_sin_cos(theta) / scifort_pi
        end if
    end function standard_cdf

    pure elemental function theta_minus_sin_cos(theta) result(y)
        real(dp), intent(in) :: theta !! angle in [0, pi/2]
        real(dp) :: y

        real(dp) :: t2

        if (abs(theta) < 1.0e-3_dp) then
            t2 = theta * theta
            y = (2.0_dp * theta * t2 / 3.0_dp) * (1.0_dp - t2 / 5.0_dp + &
                2.0_dp * t2 * t2 / 105.0_dp - t2 * t2 * t2 / 945.0_dp)
        else
            y = theta - sin(theta) * cos(theta)
        end if
    end function theta_minus_sin_cos

    pure elemental function lower_quantile(p) result(z)
        real(dp), intent(in) :: p !! lower-tail probability strictly between 0 and 0.5
        real(dp) :: z

        integer :: iteration
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid

        lo = -1.0_dp
        hi = 0.0_dp
        do iteration = 1, 80
            mid = 0.5_dp * (lo + hi)
            if (standard_cdf(mid) < p) then
                lo = mid
            else
                hi = mid
            end if
        end do
        z = 0.5_dp * (lo + hi)
    end function lower_quantile

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

end module scifort_semicircular
