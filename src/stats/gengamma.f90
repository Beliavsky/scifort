! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized gamma distribution matching scipy.stats.gengamma.
module scifort_gengamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv, log_gammainc, log_gammaincc
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, &
        valid_loc_scale
    implicit none
    private

    public :: gengamma_cdf, gengamma_isf, gengamma_logcdf, gengamma_logpdf
    public :: gengamma_logsf, gengamma_pdf, gengamma_ppf, gengamma_sf

contains

    pure elemental function gengamma_logpdf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, logz, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else if (x == mu) then
            if (c < 0.0_dp .or. c * a > 1.0_dp) then
                y = negative_infinity(x)
            else if (c * a < 1.0_dp) then
                y = positive_infinity(x)
            else
                y = log(abs(c)) - log_gamma(a) - log(sigma)
            end if
        else
            z = (x - mu) / sigma
            logz = log(z)
            t = power_from_log(c * logz)
            if (.not. ieee_is_finite(t)) then
                y = negative_infinity(x)
            else
                y = log(abs(c)) + (c * a - 1.0_dp) * logz - t - &
                    log_gamma(a) - log(sigma)
            end if
        end if
    end function gengamma_logpdf

    pure elemental function gengamma_pdf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = gengamma_logpdf(x, a, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function gengamma_pdf

    pure elemental function gengamma_logcdf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            t = power_from_log(c * log(z))
            if (c > 0.0_dp) then
                y = log_gammainc(a, t)
            else
                y = log_gammaincc(a, t)
            end if
        end if
    end function gengamma_logcdf

    pure elemental function gengamma_logsf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            t = power_from_log(c * log(z))
            if (c > 0.0_dp) then
                y = log_gammaincc(a, t)
            else
                y = log_gammainc(a, t)
            end if
        end if
    end function gengamma_logsf

    pure elemental function gengamma_cdf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            t = power_from_log(c * log(z))
            if (c > 0.0_dp) then
                y = gammainc(a, t)
            else
                y = gammaincc(a, t)
            end if
        end if
    end function gengamma_cdf

    pure elemental function gengamma_sf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            t = power_from_log(c * log(z))
            if (c > 0.0_dp) then
                y = gammaincc(a, t)
            else
                y = gammainc(a, t)
            end if
        end if
    end function gengamma_sf

    pure elemental function gengamma_ppf(p, a, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, c) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            if (c > 0.0_dp) then
                t = gammaincinv(a, p)
            else
                t = gammainccinv(a, p)
            end if
            z = root_power(t, c)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function gengamma_ppf

    pure elemental function gengamma_isf(p, a, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, c) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            if (c > 0.0_dp) then
                t = gammainccinv(a, p)
            else
                t = gammaincinv(a, p)
            end if
            z = root_power(t, c)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function gengamma_isf

    pure elemental function power_from_log(logt) result(t)
        real(dp), intent(in) :: logt !! logarithm of a positive transformed value
        real(dp) :: t
        if (logt > log(huge(1.0_dp))) then
            t = positive_infinity(logt)
        else if (logt < log(tiny(1.0_dp))) then
            t = 0.0_dp
        else
            t = exp(logt)
        end if
    end function power_from_log

    pure elemental function root_power(t, c) result(z)
        real(dp), intent(in) :: t !! nonnegative incomplete-gamma inverse value
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp) :: z, logz
        if (t == 0.0_dp) then
            if (c > 0.0_dp) then
                z = 0.0_dp
            else
                z = positive_infinity(t)
            end if
        else if (.not. ieee_is_finite(t)) then
            if (c > 0.0_dp) then
                z = positive_infinity(t)
            else
                z = 0.0_dp
            end if
        else
            logz = log(t) / c
            z = power_from_log(logz)
        end if
    end function root_power

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive scale parameter
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z) .or. z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive

    pure elemental logical function valid_shape(a, c) result(ok)
        real(dp), intent(in) :: a !! candidate positive gamma-shape parameter
        real(dp), intent(in) :: c !! candidate nonzero power parameter
        ok = ieee_is_finite(a) .and. a > 0.0_dp .and. ieee_is_finite(c) .and. c /= 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_gengamma
