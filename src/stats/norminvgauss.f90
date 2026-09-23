! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Normal inverse Gaussian distribution matching scipy.stats.norminvgauss.
module scifort_norminvgauss
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_k, only : besselk_log_derivative_x, log_besselk
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: norminvgauss_cdf, norminvgauss_isf, norminvgauss_logcdf, norminvgauss_logpdf
    public :: norminvgauss_logsf, norminvgauss_pdf, norminvgauss_ppf, norminvgauss_sf
    public :: norminvgauss_logpdf_derivatives

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

    pure elemental function norminvgauss_logpdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, r, delta
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            r = hypot(1.0_dp, z)
            delta = sqrt((a - b) * (a + b))
            y = log(a) + log_besselk(1.0_dp, a * r) - log(scifort_pi) - log(r) + &
                delta + b * z - log(sigma)
        end if
    end function norminvgauss_logpdf

    pure elemental function norminvgauss_pdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = norminvgauss_logpdf(x, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function norminvgauss_pdf

    pure elemental function norminvgauss_cdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then; y = 0.0_dp; else; y = 1.0_dp; end if
        else
            z = (x - mu) / sigma
            y = norminvgauss_tail_y(asinh(z), a, b, .false.)
        end if
    end function norminvgauss_cdf

    pure elemental function norminvgauss_sf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then; y = 1.0_dp; else; y = 0.0_dp; end if
        else
            z = (x - mu) / sigma
            y = norminvgauss_tail_y(asinh(z), a, b, .true.)
        end if
    end function norminvgauss_sf

    pure elemental function norminvgauss_logcdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = negative_infinity(x)
            else
                y = 0.0_dp
            end if
        else
            z = (x - mu) / sigma
            y = norminvgauss_logtail_y(asinh(z), a, b, .false.)
        end if
    end function norminvgauss_logcdf

    pure elemental function norminvgauss_logsf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = 0.0_dp
            else
                y = negative_infinity(x)
            end if
        else
            z = (x - mu) / sigma
            y = norminvgauss_logtail_y(asinh(z), a, b, .true.)
        end if
    end function norminvgauss_logsf

    pure elemental function norminvgauss_ppf(q, a, b, loc, scale) result(x)
        real(dp), intent(in) :: q !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, bound, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma) .or. .not. valid_probability(q)) then
            x = quiet_nan(q)
        else if (q <= 0.0_dp) then
            x = negative_infinity(q)
        else if (q >= 1.0_dp) then
            x = positive_infinity(q)
        else
            bound = integration_bound(a, b)
            lo = -bound
            hi = bound
            if (q <= 0.5_dp) then
                do iter = 1, 80
                    mid = 0.5_dp * (lo + hi)
                    if (norminvgauss_logtail_y(mid, a, b, .false.) < log(q)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            else
                do iter = 1, 80
                    mid = 0.5_dp * (lo + hi)
                    if (norminvgauss_logtail_y(mid, a, b, .true.) > log1p_safe(-q)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            end if
            x = mu + sigma * sinh(0.5_dp * (lo + hi))
        end if
    end function norminvgauss_ppf

    pure elemental function norminvgauss_isf(q, a, b, loc, scale) result(x)
        real(dp), intent(in) :: q !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, bound, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma) .or. .not. valid_probability(q)) then
            x = quiet_nan(q)
        else if (q >= 1.0_dp) then
            x = negative_infinity(q)
        else if (q <= 0.0_dp) then
            x = positive_infinity(q)
        else
            bound = integration_bound(a, b)
            lo = -bound
            hi = bound
            if (q <= 0.5_dp) then
                do iter = 1, 80
                    mid = 0.5_dp * (lo + hi)
                    if (norminvgauss_logtail_y(mid, a, b, .true.) > log(q)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            else
                do iter = 1, 80
                    mid = 0.5_dp * (lo + hi)
                    if (norminvgauss_logtail_y(mid, a, b, .false.) < log1p_safe(-q)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            end if
            x = mu + sigma * sinh(0.5_dp * (lo + hi))
        end if
    end function norminvgauss_isf

    pure subroutine norminvgauss_logpdf_derivatives(z, a, b, logf, da, db, dz)
        real(dp), intent(in) :: z !! standardized observation
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: da !! derivative with respect to a
        real(dp), intent(out) :: db !! derivative with respect to b
        real(dp), intent(out) :: dz !! derivative with respect to z
        real(dp) :: r, delta, dlogk
        if (.not. valid_shapes(a, b) .or. .not. ieee_is_finite(z)) then
            logf = quiet_nan(z); da = logf; db = logf; dz = logf; return
        end if
        r = hypot(1.0_dp, z)
        delta = sqrt((a - b) * (a + b))
        dlogk = besselk_log_derivative_x(1.0_dp, a * r)
        logf = log(a) + log_besselk(1.0_dp, a * r) - log(scifort_pi) - log(r) + delta + b*z
        da = 1.0_dp / a + r * dlogk + a / delta
        db = z - b / delta
        dz = a * z / r * dlogk - z / (r * r) + b
    end subroutine norminvgauss_logpdf_derivatives

    pure function norminvgauss_tail_y(y, a, b, upper) result(probability)
        real(dp), intent(in) :: y !! asinh of standardized variate
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        logical, intent(in) :: upper !! true for upper tail and false for lower tail
        real(dp) :: probability, log_probability
        log_probability = norminvgauss_logtail_y(y, a, b, upper)
        if (log_probability == negative_infinity(log_probability)) then
            probability = 0.0_dp
        else
            probability = exp(log_probability)
        end if
        probability = min(1.0_dp, max(0.0_dp, probability))
    end function norminvgauss_tail_y

    pure function norminvgauss_logtail_y(y, a, b, upper) result(log_probability)
        real(dp), intent(in) :: y !! asinh of standardized variate
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        logical, intent(in) :: upper !! true for upper tail and false for lower tail
        real(dp) :: log_probability, bound, left, right, slope
        bound = integration_bound(a, b)
        slope = norminvgauss_log_density_slope(y, a, b)
        if (upper) then
            if (slope <= 0.0_dp) then
                if (y >= 100.0_dp) then
                    log_probability = negative_infinity(y)
                    return
                end if
                left = y
                right = norminvgauss_tail_cutoff(y, a, b, 1)
            else
                if (y <= -bound) then
                    log_probability = 0.0_dp
                    return
                end if
                left = y
                right = bound
            end if
        else
            if (slope >= 0.0_dp) then
                if (y <= -100.0_dp) then
                    log_probability = negative_infinity(y)
                    return
                end if
                left = norminvgauss_tail_cutoff(y, a, b, -1)
                right = y
            else
                if (y >= bound) then
                    log_probability = 0.0_dp
                    return
                end if
                left = -bound
                right = y
            end if
        end if
        log_probability = integrate_log_y_density(left, right, a, b)
        log_probability = min(0.0_dp, log_probability)
    end function norminvgauss_logtail_y

    pure function norminvgauss_tail_cutoff(y, a, b, direction) result(cutoff)
        real(dp), intent(in) :: y !! tail endpoint in asinh coordinate
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        integer, intent(in) :: direction !! -1 for lower tail or +1 for upper tail
        real(dp) :: cutoff, delta, target, step, near, far, mid
        integer :: iter
        delta = sqrt((a - b) * (a + b))
        target = transformed_log_density(y, a, b, delta) - 45.0_dp
        step = 0.5_dp
        near = y
        far = y + real(direction, dp) * step
        do iter = 1, 20
            if (abs(far) >= 100.0_dp) then
                far = sign(100.0_dp, far)
                exit
            end if
            if (transformed_log_density(far, a, b, delta) <= target) exit
            step = 2.0_dp * step
            far = y + real(direction, dp) * step
        end do
        if (transformed_log_density(far, a, b, delta) > target) then
            cutoff = far
            return
        end if
        do iter = 1, 50
            mid = 0.5_dp * (near + far)
            if (transformed_log_density(mid, a, b, delta) > target) then
                near = mid
            else
                far = mid
            end if
        end do
        cutoff = far
    end function norminvgauss_tail_cutoff

    pure elemental function norminvgauss_log_density_slope(y, a, b) result(value)
        real(dp), intent(in) :: y !! asinh-coordinate value
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp) :: value, c, s
        c = cosh(y)
        s = sinh(y)
        value = a * s * besselk_log_derivative_x(1.0_dp, a * c) + b * c
    end function norminvgauss_log_density_slope

    pure function integrate_log_y_density(left, right, a, b) result(log_value)
        real(dp), intent(in) :: left !! lower integration bound in asinh coordinate
        real(dp), intent(in) :: right !! upper integration bound in asinh coordinate
        real(dp), intent(in) :: a !! positive tail shape with abs(b) < a
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp) :: log_value, delta, l, r, mid, half, y1, y2, log_weight
        real(dp) :: segment_width
        integer :: nseg, j, i
        if (right <= left) then
            log_value = negative_infinity(left)
            return
        end if
        delta = sqrt((a - b) * (a + b))
        segment_width = min(0.4_dp, 4.0_dp / sqrt(max(1.0_dp, hypot(a, b))))
        nseg = max(1, ceiling((right - left) / segment_width))
        log_value = negative_infinity(left)
        do j = 1, nseg
            l = left + (right-left) * real(j-1,dp) / real(nseg,dp)
            r = left + (right-left) * real(j,dp) / real(nseg,dp)
            mid = 0.5_dp * (l+r)
            half = 0.5_dp * (r-l)
            do i = 1, 8
                log_weight = log(half * gl_w(i))
                y1 = mid - half*gl_x(i)
                y2 = mid + half*gl_x(i)
                log_value = logaddexp_pair(log_value, &
                    log_weight + transformed_log_density(y1, a, b, delta))
                log_value = logaddexp_pair(log_value, &
                    log_weight + transformed_log_density(y2, a, b, delta))
            end do
        end do
    end function integrate_log_y_density

    pure elemental function transformed_log_density(y, a, b, delta) result(v)
        real(dp), intent(in) :: y !! asinh-coordinate integration variable
        real(dp), intent(in) :: a !! positive tail shape
        real(dp), intent(in) :: b !! skew shape
        real(dp), intent(in) :: delta !! sqrt(a^2-b^2)
        real(dp) :: v, c
        c = cosh(y)
        v = log(a) - log(scifort_pi) + log_besselk(1.0_dp, a*c) + delta + b*sinh(y)
    end function transformed_log_density

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

    pure elemental function integration_bound(a, b) result(bound)
        real(dp), intent(in) :: a !! positive tail shape
        real(dp), intent(in) :: b !! skew shape
        real(dp) :: bound, gap
        gap = max(tiny(1.0_dp), min(a-b, a+b))
        bound = max(12.0_dp, log(2.0_dp/min(gap,1.0_dp)) + 8.0_dp)
        bound = min(80.0_dp, bound)
    end function integration_bound

    pure elemental logical function valid_shapes(a, b)
        real(dp), intent(in) :: a !! positive tail shape
        real(dp), intent(in) :: b !! skew shape
        valid_shapes = ieee_is_finite(a) .and. ieee_is_finite(b) .and. a > 0.0_dp .and. abs(b) < a
    end function valid_shapes

    pure elemental logical function valid_probability(q)
        real(dp), intent(in) :: q !! probability to validate
        valid_probability = q >= 0.0_dp .and. q <= 1.0_dp
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

end module scifort_norminvgauss
