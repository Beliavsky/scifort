! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Exponentiated Weibull distribution matching scipy.stats.exponweib.
! Standard CDF is [1-exp(-z**c)]**a for z>0, a>0, c>0.
module scifort_exponweib
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private
    public :: exponweib_pdf, exponweib_logpdf, exponweib_cdf, exponweib_sf
    public :: exponweib_logcdf, exponweib_logsf, exponweib_ppf, exponweib_isf
contains
    pure elemental function exponweib_logpdf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive finite exponentiation shape parameter
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lz, lr, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z == 0.0_dp) then
                y = endpoint_logpdf(a, c, sigma, x)
            else
                lz = log(z)
                t = safe_power_from_log(c * lz, x)
                if (.not. ieee_is_finite(t)) then
                    y = negative_infinity(x)
                else
                    lr = log(-expm1_safe(-t))
                    y = log(a) + log(c) + (a - 1.0_dp) * lr - t + &
                        (c - 1.0_dp) * lz - log(sigma)
                end if
            end if
        end if
    end function exponweib_logpdf

    pure elemental function exponweib_pdf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive finite exponentiation shape parameter
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = exponweib_logpdf(x, a, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function exponweib_pdf

    pure elemental function exponweib_logcdf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive finite exponentiation shape parameter
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lz, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            lz = log(z)
            t = safe_power_from_log(c * lz, x)
            if (.not. ieee_is_finite(t)) then
                y = 0.0_dp
            else
                y = a * log(-expm1_safe(-t))
            end if
        end if
    end function exponweib_logcdf

    pure elemental function exponweib_cdf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability
        real(dp), intent(in) :: a !! positive finite exponentiation shape parameter
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = exponweib_logcdf(x, a, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function exponweib_cdf

    pure elemental function exponweib_logsf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive finite exponentiation shape parameter
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = exponweib_logcdf(x, a, c, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 0.0_dp
        else if (lc == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(-expm1_safe(lc))
        end if
    end function exponweib_logsf

    pure elemental function exponweib_sf(x, a, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability
        real(dp), intent(in) :: a !! positive finite exponentiation shape parameter
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = exponweib_logcdf(x, a, c, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 1.0_dp
        else
            y = -expm1_safe(lc)
        end if
    end function exponweib_sf

    pure elemental function exponweib_ppf(p, a, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive finite exponentiation shape parameter
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: logq, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            logq = log(p) / a
            t = -log(-expm1_safe(logq))
            z = exp(log(t) / c)
            x = affine(mu, sigma, z, p)
        end if
    end function exponweib_ppf

    pure elemental function exponweib_isf(p, a, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive finite exponentiation shape parameter
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: logq, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            logq = log1p_safe(-p) / a
            t = -log(-expm1_safe(logq))
            z = exp(log(t) / c)
            x = affine(mu, sigma, z, p)
        end if
    end function exponweib_isf

    pure elemental function endpoint_logpdf(a, c, sigma, seed) result(y)
        real(dp), intent(in) :: a !! positive exponentiation shape parameter
        real(dp), intent(in) :: c !! positive Weibull shape parameter
        real(dp), intent(in) :: sigma !! positive scale parameter
        real(dp), intent(in) :: seed !! seed used to form infinities
        real(dp) :: y
        real(dp) :: ac
        ac = a * c
        if (ac < 1.0_dp) then
            y = positive_infinity(seed)
        else if (ac > 1.0_dp) then
            y = negative_infinity(seed)
        else
            y = log(a) + log(c) - log(sigma)
        end if
    end function endpoint_logpdf

    pure elemental function safe_power_from_log(log_value, seed) result(value)
        real(dp), intent(in) :: log_value !! logarithm of nonnegative power value
        real(dp), intent(in) :: seed !! seed used to construct positive infinity
        real(dp) :: value
        if (log_value > log(huge(1.0_dp))) then
            value = positive_infinity(seed)
        else if (log_value < log(tiny(1.0_dp))) then
            value = 0.0_dp
        else
            value = exp(log_value)
        end if
    end function safe_power_from_log

    pure elemental function valid_shape(a) result(ok)
        real(dp), intent(in) :: a !! candidate shape parameter
        logical :: ok
        ok = ieee_is_finite(a) .and. a > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! requested location
        real(dp), intent(in), optional :: scale !! requested scale
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure elemental function affine(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! standardized quantile
        real(dp), intent(in) :: seed !! seed used to construct positive infinity
        real(dp) :: x
        if (.not. ieee_is_finite(z)) then
            x = positive_infinity(seed)
        else if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine
end module scifort_exponweib
