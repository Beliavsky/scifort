! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Skew-Cauchy distribution with -1 < a < 1, matching scipy.stats.skewcauchy.

module scifort_skewcauchy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: skewcauchy_cdf
    public :: skewcauchy_isf
    public :: skewcauchy_logcdf
    public :: skewcauchy_logpdf
    public :: skewcauchy_logsf
    public :: skewcauchy_pdf
    public :: skewcauchy_ppf
    public :: skewcauchy_sf

contains

    pure elemental function skewcauchy_pdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy

        logy = skewcauchy_logpdf(x, a, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function skewcauchy_pdf

    pure elemental function skewcauchy_logpdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: h, mu, q, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            h = side_scale(z, a)
            q = z / h
            y = -log(scifort_pi) - log(sigma) - log1p_square(q)
        end if
    end function skewcauchy_logpdf

    pure elemental function skewcauchy_cdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: h, mu, sigma, tail, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 0.0_dp
        else if (x == positive_infinity(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                h = 1.0_dp - a
                y = h * atan(h, -z) / scifort_pi
            else
                h = 1.0_dp + a
                tail = h * atan(h, z) / scifort_pi
                y = 1.0_dp - tail
            end if
        end if
    end function skewcauchy_cdf

    pure elemental function skewcauchy_sf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: h, lower, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 1.0_dp
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z >= 0.0_dp) then
                h = 1.0_dp + a
                y = h * atan(h, z) / scifort_pi
            else
                h = 1.0_dp - a
                lower = h * atan(h, -z) / scifort_pi
                y = 1.0_dp - lower
            end if
        end if
    end function skewcauchy_sf

    pure elemental function skewcauchy_logcdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: p

        p = skewcauchy_cdf(x, a, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(x)
        else if (p < 0.5_dp) then
            y = log(p)
        else
            y = log1p_safe(-skewcauchy_sf(x, a, loc, scale))
        end if
    end function skewcauchy_logcdf

    pure elemental function skewcauchy_logsf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: p

        p = skewcauchy_sf(x, a, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(x)
        else if (p < 0.5_dp) then
            y = log(p)
        else
            y = log1p_safe(-skewcauchy_cdf(x, a, loc, scale))
        end if
    end function skewcauchy_logsf

    pure elemental function skewcauchy_ppf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: cut, h, mu, sigma, theta, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            cut = 0.5_dp * (1.0_dp - a)
            if (p < cut) then
                h = 1.0_dp - a
                theta = scifort_pi * (p - cut) / h
            else
                h = 1.0_dp + a
                theta = scifort_pi * (p - cut) / h
            end if
            z = h * tan(theta)
            x = affine(mu, sigma, z, p)
        end if
    end function skewcauchy_ppf

    pure elemental function skewcauchy_isf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            z = skewcauchy_ppf(p, -a)
            if (ieee_is_nan(z)) then
                x = z
            else
                x = affine(mu, sigma, -z, p)
            end if
        end if
    end function skewcauchy_isf

    pure elemental function side_scale(z, a) result(h)
        real(dp), intent(in) :: z !! standardized point
        real(dp), intent(in) :: a !! skewness parameter
        real(dp) :: h

        if (z < 0.0_dp) then
            h = 1.0_dp - a
        else if (z > 0.0_dp) then
            h = 1.0_dp + a
        else
            h = 1.0_dp
        end if
    end function side_scale

    pure elemental function log1p_square(q) result(y)
        real(dp), intent(in) :: q !! real argument in log(1 + q**2)
        real(dp) :: y
        real(dp) :: aq, invq

        aq = abs(q)
        if (aq <= 1.0_dp) then
            y = log1p_safe(q * q)
        else
            invq = 1.0_dp / aq
            y = 2.0_dp * log(aq) + log1p_safe(invq * invq)
        end if
    end function log1p_square

    pure elemental logical function valid_shape(a) result(ok)
        real(dp), intent(in) :: a !! candidate skewness parameter
        ok = ieee_is_finite(a) .and. abs(a) < 1.0_dp
    end function valid_shape

    pure elemental function affine(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! standardized quantile
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: x

        if (.not. ieee_is_finite(z)) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else if (abs(z) > (huge(1.0_dp) - abs(mu)) / sigma) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else
            x = mu + sigma * z
        end if
    end function affine

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

end module scifort_skewcauchy
