! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Gauss hypergeometric distribution matching scipy.stats.gausshyper.
!
! The density on [0,1] is proportional to
! x**(a-1) (1-x)**(b-1) (1+z*x)**(-c), with a,b > 0 and z > -1.
! Normalization and tails are evaluated directly with a double-exponential
! endpoint transform.  This avoids depending on a separate 2F1 implementation
! and remains accurate when the beta factors are singular at either endpoint.
module scifort_gausshyper
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: gausshyper_cdf, gausshyper_isf, gausshyper_logcdf, gausshyper_logpdf
    public :: gausshyper_logsf, gausshyper_pdf, gausshyper_ppf, gausshyper_sf
    public :: gausshyper_logpdf_derivatives

    real(dp), parameter :: quadrature_h_target = 0.03125_dp
    real(dp), parameter :: endpoint_log_cutoff = 55.0_dp

contains

    pure elemental function gausshyper_logpdf(x, a, b, c, z, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, q, lognorm

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b, c, z) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        q = (x - mu) / sigma
        if (q < 0.0_dp .or. q > 1.0_dp .or. .not. ieee_is_finite(q)) then
            y = negative_infinity(x)
            return
        end if
        lognorm = log_segment_integral(0.0_dp, 1.0_dp, a, b, c, z)
        if (q <= 0.0_dp) then
            if (a < 1.0_dp) then
                y = positive_infinity(x)
            else if (a > 1.0_dp) then
                y = negative_infinity(x)
            else
                y = -lognorm - log(sigma)
            end if
        else if (q >= 1.0_dp) then
            if (b < 1.0_dp) then
                y = positive_infinity(x)
            else if (b > 1.0_dp) then
                y = negative_infinity(x)
            else
                y = -c * log1p_safe(z) - lognorm - log(sigma)
            end if
        else
            y = (a - 1.0_dp) * log(q) + (b - 1.0_dp) * log1p_safe(-q) - &
                c * log1p_safe(z * q) - lognorm - log(sigma)
        end if
    end function gausshyper_logpdf

    pure elemental function gausshyper_pdf(x, a, b, c, z, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = gausshyper_logpdf(x, a, b, c, z, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function gausshyper_pdf

    pure elemental function gausshyper_logcdf(x, a, b, c, z, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, q, lognorm, logupper

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b, c, z) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        q = (x - mu) / sigma
        if (q <= 0.0_dp) then
            y = negative_infinity(x)
        else if (q >= 1.0_dp) then
            y = 0.0_dp
        else
            lognorm = log_segment_integral(0.0_dp, 1.0_dp, a, b, c, z)
            y = log_segment_integral(0.0_dp, q, a, b, c, z) - lognorm
            if (y > -log(2.0_dp)) then
                logupper = log_segment_integral(q, 1.0_dp, a, b, c, z) - lognorm
                y = log1p_safe(-exp(logupper))
            end if
            y = min(0.0_dp, y)
        end if
    end function gausshyper_logcdf

    pure elemental function gausshyper_logsf(x, a, b, c, z, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, q, lognorm, loglower

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b, c, z) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        q = (x - mu) / sigma
        if (q <= 0.0_dp) then
            y = 0.0_dp
        else if (q >= 1.0_dp) then
            y = negative_infinity(x)
        else
            lognorm = log_segment_integral(0.0_dp, 1.0_dp, a, b, c, z)
            y = log_segment_integral(q, 1.0_dp, a, b, c, z) - lognorm
            if (y > -log(2.0_dp)) then
                loglower = log_segment_integral(0.0_dp, q, a, b, c, z) - lognorm
                y = log1p_safe(-exp(loglower))
            end if
            y = min(0.0_dp, y)
        end if
    end function gausshyper_logsf

    pure elemental function gausshyper_cdf(x, a, b, c, z, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lcdf, lsf

        lcdf = gausshyper_logcdf(x, a, b, c, z, loc, scale)
        if (ieee_is_nan(lcdf)) then
            y = lcdf
        else if (lcdf < -log(2.0_dp)) then
            y = exp(lcdf)
        else
            lsf = gausshyper_logsf(x, a, b, c, z, loc, scale)
            y = -expm1_safe(lsf)
        end if
    end function gausshyper_cdf

    pure elemental function gausshyper_sf(x, a, b, c, z, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lcdf, lsf

        lsf = gausshyper_logsf(x, a, b, c, z, loc, scale)
        if (ieee_is_nan(lsf)) then
            y = lsf
        else if (lsf < -log(2.0_dp)) then
            y = exp(lsf)
        else
            lcdf = gausshyper_logcdf(x, a, b, c, z, loc, scale)
            y = -expm1_safe(lcdf)
        end if
    end function gausshyper_sf

    pure elemental function gausshyper_ppf(probability, a, b, c, z, loc, scale) result(x)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, q

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b, c, z) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            x = mu
        else if (probability >= 1.0_dp) then
            x = mu + sigma
        else
            q = standard_quantile(probability, a, b, c, z, .false.)
            x = mu + sigma * q
        end if
    end function gausshyper_ppf

    pure elemental function gausshyper_isf(probability, a, b, c, z, loc, scale) result(x)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, q

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b, c, z) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            x = mu
        else if (probability <= 0.0_dp) then
            x = mu + sigma
        else
            q = standard_quantile(probability, a, b, c, z, .true.)
            x = mu + sigma * q
        end if
    end function gausshyper_isf

    pure subroutine gausshyper_logpdf_derivatives(q, a, b, c, z, logf, da, db, dc, dz_shape, dq)
        real(dp), intent(in) :: q !! standardized observation strictly inside (0,1)
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: da !! derivative with respect to a
        real(dp), intent(out) :: db !! derivative with respect to b
        real(dp), intent(out) :: dc !! derivative with respect to c
        real(dp), intent(out) :: dz_shape !! derivative with respect to z
        real(dp), intent(out) :: dq !! derivative with respect to standardized observation
        real(dp) :: e_logx, e_log1mx, e_logtilt, e_ratio, lognorm

        if (.not. valid_shapes(a, b, c, z) .or. q <= 0.0_dp .or. q >= 1.0_dp .or. &
            .not. ieee_is_finite(q)) then
            logf = quiet_nan(q); da = logf; db = logf; dc = logf; dz_shape = logf; dq = logf
            return
        end if
        call normalization_moments(a, b, c, z, lognorm, e_logx, e_log1mx, e_logtilt, e_ratio)
        logf = (a - 1.0_dp) * log(q) + (b - 1.0_dp) * log1p_safe(-q) - &
            c * log1p_safe(z * q) - lognorm
        da = log(q) - e_logx
        db = log1p_safe(-q) - e_log1mx
        dc = -log1p_safe(z * q) + e_logtilt
        dz_shape = -c * q / (1.0_dp + z * q) + c * e_ratio
        dq = (a - 1.0_dp) / q - (b - 1.0_dp) / (1.0_dp - q) - &
            c * z / (1.0_dp + z * q)
    end subroutine gausshyper_logpdf_derivatives

    pure function standard_quantile(probability, a, b, c, z, upper) result(q)
        real(dp), intent(in) :: probability !! requested tail probability in (0,1)
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        logical, intent(in) :: upper !! true when probability is an upper tail
        real(dp) :: q, lo, hi, mid, lognorm, logtarget, value
        integer :: iter

        lognorm = log_segment_integral(0.0_dp, 1.0_dp, a, b, c, z)
        logtarget = log(probability)
        lo = 0.0_dp; hi = 1.0_dp
        do iter = 1, 60
            mid = 0.5_dp * (lo + hi)
            if (upper) then
                value = log_segment_integral(mid, 1.0_dp, a, b, c, z) - lognorm
                if (value > logtarget) then
                    lo = mid
                else
                    hi = mid
                end if
            else
                value = log_segment_integral(0.0_dp, mid, a, b, c, z) - lognorm
                if (value < logtarget) then
                    lo = mid
                else
                    hi = mid
                end if
            end if
        end do
        q = 0.5_dp * (lo + hi)
    end function standard_quantile

    pure function log_segment_integral(lo, hi, a, b, c, z) result(log_value)
        real(dp), intent(in) :: lo !! lower integration endpoint in [0,1)
        real(dp), intent(in) :: hi !! upper integration endpoint in (0,1]
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp) :: log_value, tmax, h, t, li
        integer :: j, n

        if (hi <= lo) then
            log_value = negative_infinity(lo)
            return
        end if
        tmax = segment_bound(lo, hi, a, b)
        n = max(2, ceiling(2.0_dp * tmax / quadrature_h_target))
        if (mod(n, 2) /= 0) n = n + 1
        h = 2.0_dp * tmax / real(n, dp)
        log_value = negative_infinity(lo)
        do j = 0, n
            t = -tmax + h * real(j, dp)
            li = transformed_log_integrand(t, lo, hi, a, b, c, z)
            if (j == 0 .or. j == n) li = li - log(2.0_dp)
            log_value = logaddexp_pair(log_value, li)
        end do
        log_value = log_value + log(h)
    end function log_segment_integral

    pure subroutine normalization_moments(a, b, c, z, lognorm, e_logx, e_log1mx, e_logtilt, e_ratio)
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(out) :: lognorm !! log normalization integral
        real(dp), intent(out) :: e_logx !! normalized expectation of log(X)
        real(dp), intent(out) :: e_log1mx !! normalized expectation of log(1-X)
        real(dp), intent(out) :: e_logtilt !! normalized expectation of log(1+zX)
        real(dp), intent(out) :: e_ratio !! normalized expectation of X/(1+zX)
        real(dp) :: tmax, h, t, li, weight, x, lx, l1
        integer :: j, n

        lognorm = log_segment_integral(0.0_dp, 1.0_dp, a, b, c, z)
        tmax = segment_bound(0.0_dp, 1.0_dp, a, b)
        n = max(2, ceiling(2.0_dp * tmax / quadrature_h_target))
        if (mod(n, 2) /= 0) n = n + 1
        h = 2.0_dp * tmax / real(n, dp)
        e_logx = 0.0_dp; e_log1mx = 0.0_dp; e_logtilt = 0.0_dp; e_ratio = 0.0_dp
        do j = 0, n
            t = -tmax + h * real(j, dp)
            call transformed_point(t, 0.0_dp, 1.0_dp, x, lx, l1)
            li = transformed_log_integrand(t, 0.0_dp, 1.0_dp, a, b, c, z)
            if (j == 0 .or. j == n) li = li - log(2.0_dp)
            weight = exp(li + log(h) - lognorm)
            e_logx = e_logx + weight * lx
            e_log1mx = e_log1mx + weight * l1
            e_logtilt = e_logtilt + weight * log1p_safe(z * x)
            e_ratio = e_ratio + weight * x / (1.0_dp + z * x)
        end do
    end subroutine normalization_moments

    pure function transformed_log_integrand(t, lo, hi, a, b, c, z) result(value)
        real(dp), intent(in) :: t !! double-exponential integration coordinate
        real(dp), intent(in) :: lo !! lower x endpoint
        real(dp), intent(in) :: hi !! upper x endpoint
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp) :: value, x, lx, l1, lr, lrc, logjac

        call transformed_point_full(t, lo, hi, x, lx, l1, lr, lrc)
        logjac = log(hi - lo) + log(scifort_pi) + log(cosh(t)) + lr + lrc
        value = (a - 1.0_dp) * lx + (b - 1.0_dp) * l1 - &
            c * log1p_safe(z * x) + logjac
    end function transformed_log_integrand

    pure subroutine transformed_point(t, lo, hi, x, lx, l1)
        real(dp), intent(in) :: t !! double-exponential integration coordinate
        real(dp), intent(in) :: lo !! lower x endpoint
        real(dp), intent(in) :: hi !! upper x endpoint
        real(dp), intent(out) :: x !! transformed x coordinate
        real(dp), intent(out) :: lx !! log(x) formed stably near zero
        real(dp), intent(out) :: l1 !! log(1-x) formed stably near one
        real(dp) :: lr, lrc
        call transformed_point_full(t, lo, hi, x, lx, l1, lr, lrc)
    end subroutine transformed_point

    pure subroutine transformed_point_full(t, lo, hi, x, lx, l1, lr, lrc)
        real(dp), intent(in) :: t !! double-exponential integration coordinate
        real(dp), intent(in) :: lo !! lower x endpoint
        real(dp), intent(in) :: hi !! upper x endpoint
        real(dp), intent(out) :: x !! transformed x coordinate
        real(dp), intent(out) :: lx !! log(x) formed stably near zero
        real(dp), intent(out) :: l1 !! log(1-x) formed stably near one
        real(dp), intent(out) :: lr !! log of logistic transform r
        real(dp), intent(out) :: lrc !! log(1-r)
        real(dp) :: e, ld, r, v, span

        v = scifort_pi * sinh(t)
        if (v >= 0.0_dp) then
            e = exp(-v)
            ld = log1p_safe(e)
            lr = -ld
            lrc = -v - ld
            if (lrc > -36.0_dp) then
                r = 1.0_dp - exp(lrc)
            else
                r = 1.0_dp
            end if
        else
            e = exp(v)
            ld = log1p_safe(e)
            lr = v - ld
            lrc = -ld
            if (lr > log(tiny(1.0_dp))) then
                r = exp(lr)
            else
                r = 0.0_dp
            end if
        end if
        span = hi - lo
        if (lo == 0.0_dp) then
            lx = log(hi) + lr
        else
            lx = log(lo + span * r)
        end if
        if (hi == 1.0_dp) then
            l1 = log1p_safe(-lo) + lrc
        else
            l1 = log1p_safe(-(lo + span * r))
        end if
        if (r <= 0.0_dp) then
            x = lo
        else if (r >= 1.0_dp) then
            x = hi
        else
            x = lo + span * r
        end if
    end subroutine transformed_point_full

    pure elemental function segment_bound(lo, hi, a, b) result(bound)
        real(dp), intent(in) :: lo !! lower integration endpoint
        real(dp), intent(in) :: hi !! upper integration endpoint
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp) :: bound, edge_rate

        edge_rate = 1.0_dp
        if (lo <= 0.0_dp) edge_rate = min(edge_rate, a)
        if (hi >= 1.0_dp) edge_rate = min(edge_rate, b)
        bound = 4.0_dp
        if (edge_rate < 1.0_dp) then
            bound = max(bound, asinh(endpoint_log_cutoff / (scifort_pi * edge_rate)))
        end if
    end function segment_bound

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

    pure elemental logical function valid_shapes(a, b, c, z)
        real(dp), intent(in) :: a !! first beta shape to validate
        real(dp), intent(in) :: b !! second beta shape to validate
        real(dp), intent(in) :: c !! tilt exponent to validate
        real(dp), intent(in) :: z !! tilt parameter to validate
        valid_shapes = ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
            ieee_is_finite(c) .and. ieee_is_finite(z) .and. a > 0.0_dp .and. &
            b > 0.0_dp .and. z > -1.0_dp
    end function valid_shapes

    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability to validate
        valid_probability = ieee_is_finite(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location parameter
        real(dp), intent(in), optional :: scale !! optional positive scale parameter
        real(dp), intent(out) :: mu !! resolved location parameter
        real(dp), intent(out) :: sigma !! resolved positive scale parameter
        mu = 0.0_dp; sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_gausshyper
