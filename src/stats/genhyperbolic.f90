! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized hyperbolic distribution matching scipy.stats.genhyperbolic.
module scifort_genhyperbolic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_k, only : besselk_log_derivative_x, log_besselk_order_derivative, &
        log_besselk_scaled
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: genhyperbolic_cdf, genhyperbolic_isf, genhyperbolic_logcdf, genhyperbolic_logpdf
    public :: genhyperbolic_logsf, genhyperbolic_pdf, genhyperbolic_ppf, genhyperbolic_sf
    public :: genhyperbolic_logpdf_derivatives

    real(dp), parameter :: gl_x(8) = [ &
        0.0950125098376374401853193354250_dp, 0.281603550779258913230460501460_dp, &
        0.458016777657227386342419442984_dp, 0.617876244402643748446671764049_dp, &
        0.755404408355003033895101194847_dp, 0.865631202387831743880467897712_dp, &
        0.944575023073232576077988415535_dp, 0.989400934991649932596154173450_dp]
    real(dp), parameter :: gl_w(8) = [ &
        0.189450610455068496285396723208_dp, 0.182603415044923588866763667969_dp, &
        0.169156519395002538189312079030_dp, 0.149595988816576732081501730547_dp, &
        0.124628971255533872052476282192_dp, 0.0951585116824927848099251076022_dp, &
        0.0622535239386478928628438369944_dp, 0.0271524594117540948517805724560_dp]

contains

    pure elemental function genhyperbolic_logpdf(x, p, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape constrained by a and p
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, r, q
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            r = hypot(1.0_dp, z)
            q = p - 0.5_dp
            y = log_normalizer(p, a, b) + log_besselk_scaled(q, a * r) + &
                stable_exponent(z, a, b) + q * log(r) - log(sigma)
        end if
    end function genhyperbolic_logpdf

    pure elemental function genhyperbolic_pdf(x, p, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape constrained by a and p
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = genhyperbolic_logpdf(x, p, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function genhyperbolic_pdf

    pure elemental function genhyperbolic_cdf(x, p, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape constrained by a and p
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then; y = 0.0_dp; else; y = 1.0_dp; end if
        else
            z = (x - mu) / sigma
            y = genhyperbolic_tail_y(asinh(z), p, a, b, .false.)
        end if
    end function genhyperbolic_cdf

    pure elemental function genhyperbolic_sf(x, p, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape constrained by a and p
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then; y = 1.0_dp; else; y = 0.0_dp; end if
        else
            z = (x - mu) / sigma
            y = genhyperbolic_tail_y(asinh(z), p, a, b, .true.)
        end if
    end function genhyperbolic_sf

    pure elemental function genhyperbolic_logcdf(x, p, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape constrained by a and p
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then; y = negative_infinity(x); else; y = 0.0_dp; end if
        else
            z = (x - mu) / sigma
            y = genhyperbolic_logtail_y(asinh(z), p, a, b, .false.)
        end if
    end function genhyperbolic_logcdf

    pure elemental function genhyperbolic_logsf(x, p, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape constrained by a and p
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then; y = 0.0_dp; else; y = negative_infinity(x); end if
        else
            z = (x - mu) / sigma
            y = genhyperbolic_logtail_y(asinh(z), p, a, b, .true.)
        end if
    end function genhyperbolic_logsf

    pure elemental function genhyperbolic_ppf(probability, p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape constrained by a and p
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, bound, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, a, b) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            x = negative_infinity(probability)
        else if (probability >= 1.0_dp) then
            x = positive_infinity(probability)
        else
            bound = integration_bound(p, a, b)
            lo = -bound; hi = bound
            if (probability <= 0.5_dp) then
                do iter = 1, 90
                    mid = 0.5_dp * (lo + hi)
                    if (genhyperbolic_logtail_y(mid, p, a, b, .false.) < log(probability)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            else
                do iter = 1, 90
                    mid = 0.5_dp * (lo + hi)
                    if (genhyperbolic_logtail_y(mid, p, a, b, .true.) > log1p_safe(-probability)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            end if
            x = mu + sigma * sinh(0.5_dp * (lo + hi))
        end if
    end function genhyperbolic_ppf

    pure elemental function genhyperbolic_isf(probability, p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape constrained by a and p
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, bound, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, a, b) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            x = negative_infinity(probability)
        else if (probability <= 0.0_dp) then
            x = positive_infinity(probability)
        else
            bound = integration_bound(p, a, b)
            lo = -bound; hi = bound
            if (probability <= 0.5_dp) then
                do iter = 1, 90
                    mid = 0.5_dp * (lo + hi)
                    if (genhyperbolic_logtail_y(mid, p, a, b, .true.) > log(probability)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            else
                do iter = 1, 90
                    mid = 0.5_dp * (lo + hi)
                    if (genhyperbolic_logtail_y(mid, p, a, b, .false.) < log1p_safe(-probability)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            end if
            x = mu + sigma * sinh(0.5_dp * (lo + hi))
        end if
    end function genhyperbolic_isf

    pure subroutine genhyperbolic_logpdf_derivatives(z, p, a, b, logf, dp_shape, da, db, dz)
        real(dp), intent(in) :: z !! standardized observation
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape constrained by a and p
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: dp_shape !! derivative with respect to p
        real(dp), intent(out) :: da !! derivative with respect to a
        real(dp), intent(out) :: db !! derivative with respect to b
        real(dp), intent(out) :: dz !! derivative with respect to standardized observation
        real(dp) :: d2, d, q, r, dlogk_norm, dlogk_obs
        logf = genhyperbolic_logpdf(z, p, a, b)
        d2 = (a - b) * (a + b)
        if (d2 <= 0.0_dp) then
            dp_shape = quiet_nan(z); da = dp_shape; db = dp_shape; dz = dp_shape
            return
        end if
        d = sqrt(d2)
        q = p - 0.5_dp
        r = hypot(1.0_dp, z)
        dlogk_norm = besselk_log_derivative_x(p, d)
        dlogk_obs = besselk_log_derivative_x(q, a * r)
        dp_shape = log(d) - log(a) - log_besselk_order_derivative(p, d) + &
            log_besselk_order_derivative(q, a * r) + log(r)
        da = p * a / d2 - q / a - dlogk_norm * a / d + dlogk_obs * r
        db = -p * b / d2 + dlogk_norm * b / d + z
        dz = b + dlogk_obs * a * z / r + q * z / (r * r)
    end subroutine genhyperbolic_logpdf_derivatives

    pure function genhyperbolic_tail_y(y, p, a, b, upper) result(probability)
        real(dp), intent(in) :: y !! asinh of standardized variate
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        logical, intent(in) :: upper !! true for upper tail and false for lower tail
        real(dp) :: probability, log_probability
        log_probability = genhyperbolic_logtail_y(y, p, a, b, upper)
        if (log_probability == negative_infinity(log_probability)) then
            probability = 0.0_dp
        else
            probability = exp(log_probability)
        end if
        probability = min(1.0_dp, max(0.0_dp, probability))
    end function genhyperbolic_tail_y

    pure function genhyperbolic_logtail_y(y, p, a, b, upper) result(log_probability)
        real(dp), intent(in) :: y !! asinh of standardized variate
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        logical, intent(in) :: upper !! true for upper tail and false for lower tail
        real(dp) :: log_probability, bound, left, right, slope
        bound = integration_bound(p, a, b)
        slope = log_density_slope(y, p, a, b)
        if (upper) then
            if (slope <= 0.0_dp) then
                if (y >= bound) then
                    log_probability = negative_infinity(y)
                    return
                end if
                left = y
                right = tail_cutoff(y, p, a, b, 1, bound)
            else
                left = y
                right = bound
            end if
        else
            if (slope >= 0.0_dp) then
                if (y <= -bound) then
                    log_probability = negative_infinity(y)
                    return
                end if
                left = tail_cutoff(y, p, a, b, -1, bound)
                right = y
            else
                left = -bound
                right = y
            end if
        end if
        log_probability = integrate_log_y_density(left, right, p, a, b)
        log_probability = min(0.0_dp, log_probability)
    end function genhyperbolic_logtail_y

    pure function tail_cutoff(y, p, a, b, direction, bound) result(cutoff)
        real(dp), intent(in) :: y !! tail endpoint in asinh coordinate
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        integer, intent(in) :: direction !! -1 for lower tail or +1 for upper tail
        real(dp), intent(in) :: bound !! absolute integration limit
        real(dp) :: cutoff, target, step, near, far, mid
        integer :: iter
        target = transformed_log_density(y, p, a, b) - 45.0_dp
        step = 0.5_dp; near = y; far = y + real(direction, dp) * step
        do iter = 1, 30
            if (abs(far) >= bound) then
                far = sign(bound, far)
                exit
            end if
            if (transformed_log_density(far, p, a, b) <= target) exit
            step = 2.0_dp * step
            far = y + real(direction, dp) * step
        end do
        if (transformed_log_density(far, p, a, b) > target) then
            cutoff = far
            return
        end if
        do iter = 1, 50
            mid = 0.5_dp * (near + far)
            if (transformed_log_density(mid, p, a, b) > target) then
                near = mid
            else
                far = mid
            end if
        end do
        cutoff = far
    end function tail_cutoff

    pure elemental function log_density_slope(y, p, a, b) result(value)
        real(dp), intent(in) :: y !! asinh-coordinate value
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        real(dp) :: value, c, s, q
        c = cosh(y); s = sinh(y); q = p - 0.5_dp
        value = a * s * besselk_log_derivative_x(q, a * c) + b * c + (p + 0.5_dp) * tanh(y)
    end function log_density_slope

    pure function integrate_log_y_density(left, right, p, a, b) result(log_value)
        real(dp), intent(in) :: left !! lower integration bound in asinh coordinate
        real(dp), intent(in) :: right !! upper integration bound in asinh coordinate
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        real(dp) :: log_value, l, r, mid, half, y1, y2, log_weight, segment_width
        integer :: nseg, j, i
        if (right <= left) then
            log_value = negative_infinity(left)
            return
        end if
        segment_width = min(0.4_dp, 4.0_dp / sqrt(max(1.0_dp, hypot(a, b))))
        nseg = max(1, ceiling((right - left) / segment_width))
        log_value = negative_infinity(left)
        do j = 1, nseg
            l = left + (right-left) * real(j-1,dp) / real(nseg,dp)
            r = left + (right-left) * real(j,dp) / real(nseg,dp)
            mid = 0.5_dp * (l+r); half = 0.5_dp * (r-l)
            do i = 1, 8
                log_weight = log(half * gl_w(i))
                y1 = mid - half*gl_x(i); y2 = mid + half*gl_x(i)
                log_value = logaddexp_pair(log_value, log_weight + transformed_log_density(y1, p, a, b))
                log_value = logaddexp_pair(log_value, log_weight + transformed_log_density(y2, p, a, b))
            end do
        end do
    end function integrate_log_y_density

    pure elemental function transformed_log_density(y, p, a, b) result(v)
        real(dp), intent(in) :: y !! asinh-coordinate integration variable
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        real(dp) :: v, c, q
        c = cosh(y); q = p - 0.5_dp
        v = log_normalizer(p, a, b) + log_besselk_scaled(q, a*c) + &
            stable_hyperbolic_exponent(y, a, b) + (p + 0.5_dp) * log(c)
    end function transformed_log_density

    pure elemental function log_normalizer(p, a, b) result(v)
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        real(dp) :: v, d2, d, nu, q
        q = p - 0.5_dp
        d2 = max(0.0_dp, (a - b) * (a + b))
        if (d2 == 0.0_dp) then
            nu = -p
            v = -(nu - 1.0_dp) * log(2.0_dp) - log_gamma(nu) - &
                0.5_dp * log(2.0_dp * scifort_pi) - q * log(a)
        else
            d = sqrt(d2)
            v = p * log(d) - 0.5_dp * log(2.0_dp * scifort_pi) - q * log(a) - &
                log_besselk_scaled(p, d) + d
        end if
    end function log_normalizer

    pure elemental function stable_exponent(z, a, b) result(v)
        real(dp), intent(in) :: z !! standardized observation
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        real(dp) :: v, r, az
        r = hypot(1.0_dp, z); az = abs(z)
        if (z >= 0.0_dp) then
            v = -(a - b) * z - a / (r + az)
        else
            v = -(a + b) * az - a / (r + az)
        end if
    end function stable_exponent

    pure elemental function stable_hyperbolic_exponent(y, a, b) result(v)
        real(dp), intent(in) :: y !! asinh-coordinate value
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        real(dp) :: v, ep, em
        ep = exp(y); em = exp(-y)
        v = -0.5_dp * ((a - b) * ep + (a + b) * em)
    end function stable_hyperbolic_exponent

    pure elemental function logaddexp_pair(a, b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y, m
        if (.not. ieee_is_finite(a)) then
            y = b
        else if (.not. ieee_is_finite(b)) then
            y = a
        else
            m = max(a, b)
            y = m + log1p_safe(exp(-abs(a - b)))
        end if
    end function logaddexp_pair

    pure elemental function integration_bound(p, a, b) result(bound)
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        real(dp) :: bound, gap
        gap = min(a - b, a + b)
        if (gap > 0.0_dp) then
            bound = max(14.0_dp, log(2.0_dp/min(gap,1.0_dp)) + 10.0_dp)
            bound = min(100.0_dp, bound)
        else
            bound = max(30.0_dp, 50.0_dp / max(abs(p), 0.08_dp))
            bound = min(650.0_dp, bound)
        end if
    end function integration_bound

    pure elemental logical function valid_shapes(p, a, b)
        real(dp), intent(in) :: p !! real tail shape parameter
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! skew shape
        valid_shapes = ieee_is_finite(p) .and. ieee_is_finite(a) .and. ieee_is_finite(b) .and. a > 0.0_dp .and. &
            ((p >= 0.0_dp .and. abs(b) < a) .or. (p < 0.0_dp .and. abs(b) <= a))
    end function valid_shapes

    pure elemental logical function valid_probability(q)
        real(dp), intent(in) :: q !! probability to validate
        valid_probability = ieee_is_finite(q) .and. q >= 0.0_dp .and. q <= 1.0_dp
    end function valid_probability

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location parameter
        real(dp), intent(in), optional :: scale !! optional positive scale parameter
        real(dp), intent(out) :: mu !! resolved location parameter
        real(dp), intent(out) :: sigma !! resolved scale parameter
        mu = 0.0_dp; sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_genhyperbolic
