! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Wrapped Cauchy distribution matching scipy.stats.wrapcauchy.
! Standardized support is 0 <= z <= 2*pi with 0 < c < 1.
module scifort_wrapcauchy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: wrapcauchy_cdf, wrapcauchy_isf, wrapcauchy_logcdf, wrapcauchy_logpdf
    public :: wrapcauchy_logsf, wrapcauchy_pdf, wrapcauchy_ppf, wrapcauchy_sf

contains

    pure elemental function wrapcauchy_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y
        real(dp) :: h, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp .or. z > 2.0_dp * scifort_pi .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                h = hypot(1.0_dp - c, 2.0_dp * sqrt(c) * sin(0.5_dp * z))
                y = log(1.0_dp - c) + log(1.0_dp + c) - &
                    log(2.0_dp * scifort_pi) - 2.0_dp * log(h) - log(sigma)
            end if
        end if
    end function wrapcauchy_logpdf

    pure elemental function wrapcauchy_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y, ly
        ly = wrapcauchy_logpdf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function wrapcauchy_pdf

    pure elemental function wrapcauchy_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (z >= 2.0_dp * scifort_pi .or. .not. ieee_is_finite(z)) then
                y = 1.0_dp
            else
                y = standard_cdf(z, c)
            end if
        end if
    end function wrapcauchy_cdf

    pure elemental function wrapcauchy_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 1.0_dp
            else if (z >= 2.0_dp * scifort_pi .or. .not. ieee_is_finite(z)) then
                y = 0.0_dp
            else
                y = standard_cdf(2.0_dp * scifort_pi - z, c)
            end if
        end if
    end function wrapcauchy_sf

    pure elemental function wrapcauchy_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y, p
        p = wrapcauchy_cdf(x, c, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(p)
        end if
    end function wrapcauchy_logcdf

    pure elemental function wrapcauchy_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y, p
        p = wrapcauchy_sf(x, c, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(p)
        end if
    end function wrapcauchy_logsf

    pure elemental function wrapcauchy_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, theta, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = mu + sigma * (2.0_dp * scifort_pi)
        else
            theta = scifort_pi * p
            z = 2.0_dp * atan2((1.0_dp - c) * sin(theta), &
                (1.0_dp + c) * cos(theta))
            x = mu + sigma * z
        end if
    end function wrapcauchy_ppf

    pure elemental function wrapcauchy_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, theta, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = mu + sigma * (2.0_dp * scifort_pi)
        else
            theta = scifort_pi * p
            z = 2.0_dp * atan2((1.0_dp - c) * sin(theta), &
                (1.0_dp + c) * cos(theta))
            x = mu + sigma * (2.0_dp * scifort_pi - z)
        end if
    end function wrapcauchy_isf

    pure elemental function standard_cdf(z, c) result(p)
        real(dp), intent(in) :: z !! standardized angle strictly between 0 and 2*pi
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp) :: p
        p = atan2((1.0_dp + c) * sin(0.5_dp * z), &
            (1.0_dp - c) * cos(0.5_dp * z)) / scifort_pi
    end function standard_cdf

    pure elemental logical function valid_shape(c) result(ok)
        real(dp), intent(in) :: c !! candidate wrapped-Cauchy concentration parameter
        ok = ieee_is_finite(c) .and. c > 0.0_dp .and. c < 1.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! resolved location, 0 when absent
        real(dp), intent(out) :: sigma !! resolved scale, 1 when absent
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_wrapcauchy
