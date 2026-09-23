! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Cosine distribution. In standardized coordinates z = (x - loc) / scale,
! support is [-pi, pi] and f(z) = (1 + cos(z)) / (2*pi).
! This matches scipy.stats.cosine.

module scifort_cosine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: cosine_cdf
    public :: cosine_isf
    public :: cosine_logcdf
    public :: cosine_logpdf
    public :: cosine_logsf
    public :: cosine_pdf
    public :: cosine_ppf
    public :: cosine_sf

contains

    pure elemental function cosine_pdf(x, loc, scale) result(y)
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
            if (abs(z) >= scifort_pi) then
                y = 0.0_dp
            else
                y = cos(0.5_dp * z)**2 / (scifort_pi * sigma)
            end if
        end if
    end function cosine_pdf

    pure elemental function cosine_logpdf(x, loc, scale) result(y)
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
            if (abs(z) >= scifort_pi) then
                y = negative_infinity(x)
            else
                y = 2.0_dp * log(abs(cos(0.5_dp * z))) - &
                    log(scifort_pi) - log(sigma)
            end if
        end if
    end function cosine_logpdf

    pure elemental function cosine_cdf(x, loc, scale) result(y)
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
    end function cosine_cdf

    pure elemental function cosine_sf(x, loc, scale) result(y)
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
    end function cosine_sf

    pure elemental function cosine_logcdf(x, loc, scale) result(y)
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
        if (z <= -scifort_pi) then
            y = negative_infinity(x)
        else if (z >= scifort_pi) then
            y = 0.0_dp
        else if (z <= 0.0_dp) then
            cdf_value = standard_cdf(z)
            y = log(cdf_value)
        else
            sf_value = standard_cdf(-z)
            y = log1p_safe(-sf_value)
        end if
    end function cosine_logcdf

    pure elemental function cosine_logsf(x, loc, scale) result(y)
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
            y = cosine_logcdf(2.0_dp * mu - x, mu, sigma)
        end if
    end function cosine_logsf

    pure elemental function cosine_ppf(p, loc, scale) result(x)
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
            x = mu - scifort_pi * sigma
        else if (p >= 1.0_dp) then
            x = mu + scifort_pi * sigma
        else if (p < 0.5_dp) then
            z = lower_quantile(p)
            x = mu + sigma * z
        else
            z = lower_quantile(1.0_dp - p)
            x = mu - sigma * z
        end if
    end function cosine_ppf

    pure elemental function cosine_isf(p, loc, scale) result(x)
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
            x = 2.0_dp * mu - cosine_ppf(p, mu, sigma)
        end if
    end function cosine_isf

    pure elemental function standard_cdf(z) result(y)
        real(dp), intent(in) :: z !! standardized variate
        real(dp) :: y

        real(dp) :: t

        if (z <= -scifort_pi) then
            y = 0.0_dp
        else if (z >= scifort_pi) then
            y = 1.0_dp
        else if (z < 0.0_dp) then
            t = scifort_pi + z
            y = x_minus_sin(t) / (2.0_dp * scifort_pi)
        else
            y = 0.5_dp + (z + sin(z)) / (2.0_dp * scifort_pi)
        end if
    end function standard_cdf

    pure elemental function x_minus_sin(x) result(y)
        real(dp), intent(in) :: x !! nonnegative argument
        real(dp) :: y

        real(dp) :: x2

        if (abs(x) < 1.0e-3_dp) then
            x2 = x * x
            y = (x * x2 / 6.0_dp) * (1.0_dp - x2 / 20.0_dp + &
                x2 * x2 / 840.0_dp - x2 * x2 * x2 / 60480.0_dp)
        else
            y = x - sin(x)
        end if
    end function x_minus_sin

    pure elemental function lower_quantile(p) result(z)
        real(dp), intent(in) :: p !! lower-tail probability strictly between 0 and 0.5
        real(dp) :: z

        integer :: iteration
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid

        lo = -scifort_pi
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

end module scifort_cosine
