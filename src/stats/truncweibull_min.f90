! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Doubly truncated Weibull-minimum distribution matching scipy.stats.truncweibull_min.
module scifort_truncweibull_min
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: truncweibull_min_cdf, truncweibull_min_isf, truncweibull_min_logcdf
    public :: truncweibull_min_logpdf, truncweibull_min_logsf, truncweibull_min_pdf
    public :: truncweibull_min_ppf, truncweibull_min_sf

contains

    pure elemental function truncweibull_min_logpdf(x, c, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! finite standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! finite standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, zc
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < a .or. z > b .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else if (z == 0.0_dp) then
                if (c < 1.0_dp) then
                    y = positive_infinity(x)
                else if (c == 1.0_dp) then
                    y = -log_normalizer(c, a, b) - log(sigma)
                else
                    y = negative_infinity(x)
                end if
            else
                zc = positive_power(z, c)
                y = log(c) + (c - 1.0_dp) * log(z) - zc - &
                    log_normalizer(c, a, b) - log(sigma)
            end if
        end if
    end function truncweibull_min_logpdf

    pure elemental function truncweibull_min_pdf(x, c, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! finite standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! finite standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = truncweibull_min_logpdf(x, c, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else if (ly >= log(huge(1.0_dp))) then
            y = huge(1.0_dp)
        else
            y = exp(ly)
        end if
    end function truncweibull_min_pdf

    pure elemental function truncweibull_min_logcdf(x, c, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! finite standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! finite standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ac, mu, sigma, z, zc
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= a) then
                y = negative_infinity(x)
            else if (z >= b) then
                y = 0.0_dp
            else
                ac = positive_power(a, c)
                zc = positive_power(z, c)
                y = -ac + log_one_minus_exp(-(zc - ac)) - log_normalizer(c, a, b)
            end if
        end if
    end function truncweibull_min_logcdf

    pure elemental function truncweibull_min_cdf(x, c, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! finite standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! finite standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = truncweibull_min_logcdf(x, c, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function truncweibull_min_cdf

    pure elemental function truncweibull_min_logsf(x, c, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! finite standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! finite standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, bc, mu, sigma, z, zc
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= a) then
                y = 0.0_dp
            else if (z >= b) then
                y = negative_infinity(x)
            else
                zc = positive_power(z, c)
                bc = positive_power(b, c)
                y = -zc + log_one_minus_exp(-(bc - zc)) - log_normalizer(c, a, b)
            end if
        end if
    end function truncweibull_min_logsf

    pure elemental function truncweibull_min_sf(x, c, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! finite standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! finite standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = truncweibull_min_logsf(x, c, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function truncweibull_min_sf

    pure elemental function truncweibull_min_ppf(p, c, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! finite standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! finite standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, ac, bc, loge, mu, sigma, z, zc
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, a, b) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu + sigma * a
        else if (p == 1.0_dp) then
            x = mu + sigma * b
        else
            ac = positive_power(a, c)
            bc = positive_power(b, c)
            loge = logaddexp(log1p_safe(-p) - ac, log(p) - bc)
            zc = -loge
            z = exp(log(zc) / c)
            x = mu + sigma * z
        end if
    end function truncweibull_min_ppf

    pure elemental function truncweibull_min_isf(p, c, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! finite standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! finite standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, ac, bc, loge, mu, sigma, z, zc
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, a, b) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu + sigma * a
        else if (p == 0.0_dp) then
            x = mu + sigma * b
        else
            ac = positive_power(a, c)
            bc = positive_power(b, c)
            loge = logaddexp(log(p) - ac, log1p_safe(-p) - bc)
            zc = -loge
            z = exp(log(zc) / c)
            x = mu + sigma * z
        end if
    end function truncweibull_min_isf

    pure elemental function log_normalizer(c, a, b) result(y)
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point
        real(dp) :: y, ac, bc
        ac = positive_power(a, c)
        bc = positive_power(b, c)
        y = -ac + log_one_minus_exp(-(bc - ac))
    end function log_normalizer

    pure elemental function positive_power(x, c) result(y)
        real(dp), intent(in) :: x !! nonnegative finite base
        real(dp), intent(in) :: c !! positive finite exponent
        real(dp) :: y, e
        if (x == 0.0_dp) then
            y = 0.0_dp
        else
            e = c * log(x)
            if (e > log(huge(1.0_dp))) then
                y = huge(1.0_dp)
            else if (e < log(tiny(1.0_dp))) then
                y = 0.0_dp
            else
                y = exp(e)
            end if
        end if
    end function positive_power

    pure elemental function log_one_minus_exp(t) result(y)
        real(dp), intent(in) :: t !! nonpositive exponent in log(1-exp(t))
        real(dp) :: y
        if (t == 0.0_dp) then
            y = negative_infinity(t)
        else if (t < -0.69314718055994530942_dp) then
            y = log(1.0_dp - exp(t))
        else
            y = log(-expm1_safe(t))
        end if
    end function log_one_minus_exp

    pure elemental function logaddexp(a, b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y, m
        m = max(a, b)
        if (m == negative_infinity(m)) then
            y = m
        else
            y = m + log(exp(a - m) + exp(b - m))
        end if
    end function logaddexp

    pure elemental logical function valid_shapes(c, a, b) result(ok)
        real(dp), intent(in) :: c !! candidate Weibull shape parameter
        real(dp), intent(in) :: a !! candidate lower standardized truncation point
        real(dp), intent(in) :: b !! candidate upper standardized truncation point
        ok = ieee_is_finite(c) .and. ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
            c > 0.0_dp .and. a >= 0.0_dp .and. b > a
    end function valid_shapes

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

end module scifort_truncweibull_min
