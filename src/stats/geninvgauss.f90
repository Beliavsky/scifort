! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized inverse Gaussian distribution matching scipy.stats.geninvgauss.
module scifort_geninvgauss
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_k, only : besselk_log_derivative_x, log_besselk, &
        log_besselk_order_derivative
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: geninvgauss_cdf, geninvgauss_isf, geninvgauss_logcdf, geninvgauss_logpdf
    public :: geninvgauss_logsf, geninvgauss_pdf, geninvgauss_ppf, geninvgauss_sf
    public :: geninvgauss_logpdf_derivatives

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

    pure elemental function geninvgauss_logpdf(x, p, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                y = (p - 1.0_dp) * log(z) - 0.5_dp * b * (z + 1.0_dp / z) - &
                    log(2.0_dp) - log_besselk(p, b) - log(sigma)
            end if
        end if
    end function geninvgauss_logpdf

    pure elemental function geninvgauss_pdf(x, p, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = geninvgauss_logpdf(x, p, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function geninvgauss_pdf

    pure elemental function geninvgauss_cdf(x, p, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, yy
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (.not. ieee_is_finite(z)) then
                y = 1.0_dp
            else
                yy = log(z)
                y = geninvgauss_tail_y(yy, p, b, .false.)
            end if
        end if
    end function geninvgauss_cdf

    pure elemental function geninvgauss_sf(x, p, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, yy
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 1.0_dp
            else if (.not. ieee_is_finite(z)) then
                y = 0.0_dp
            else
                yy = log(z)
                y = geninvgauss_tail_y(yy, p, b, .true.)
            end if
        end if
    end function geninvgauss_sf

    pure elemental function geninvgauss_logcdf(x, p, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = negative_infinity(x)
            else if (.not. ieee_is_finite(z)) then
                y = 0.0_dp
            else
                y = geninvgauss_logtail_y(log(z), p, b, .false.)
            end if
        end if
    end function geninvgauss_logcdf

    pure elemental function geninvgauss_logsf(x, p, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (.not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                y = geninvgauss_logtail_y(log(z), p, b, .true.)
            end if
        end if
    end function geninvgauss_logsf

    pure elemental function geninvgauss_ppf(q, p, b, loc, scale) result(x)
        real(dp), intent(in) :: q !! lower-tail probability in [0,1]
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, bound, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, b) .or. .not. valid_loc_scale(mu, sigma) .or. .not. valid_probability(q)) then
            x = quiet_nan(q)
        else if (q <= 0.0_dp) then
            x = mu
        else if (q >= 1.0_dp) then
            x = positive_infinity(q)
        else
            bound = integration_bound(p, b)
            lo = -bound
            hi = bound
            if (q <= 0.5_dp) then
                do iter = 1, 90
                    mid = 0.5_dp * (lo + hi)
                    if (geninvgauss_logtail_y(mid, p, b, .false.) < log(q)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            else
                do iter = 1, 90
                    mid = 0.5_dp * (lo + hi)
                    if (geninvgauss_logtail_y(mid, p, b, .true.) > log1p_safe(-q)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            end if
            x = mu + sigma * exp(0.5_dp * (lo + hi))
        end if
    end function geninvgauss_ppf

    pure elemental function geninvgauss_isf(q, p, b, loc, scale) result(x)
        real(dp), intent(in) :: q !! upper-tail probability in [0,1]
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, bound, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(p, b) .or. .not. valid_loc_scale(mu, sigma) .or. .not. valid_probability(q)) then
            x = quiet_nan(q)
        else if (q >= 1.0_dp) then
            x = mu
        else if (q <= 0.0_dp) then
            x = positive_infinity(q)
        else
            bound = integration_bound(p, b)
            lo = -bound
            hi = bound
            if (q <= 0.5_dp) then
                do iter = 1, 90
                    mid = 0.5_dp * (lo + hi)
                    if (geninvgauss_logtail_y(mid, p, b, .true.) > log(q)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            else
                do iter = 1, 90
                    mid = 0.5_dp * (lo + hi)
                    if (geninvgauss_logtail_y(mid, p, b, .false.) < log1p_safe(-q)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
            end if
            x = mu + sigma * exp(0.5_dp * (lo + hi))
        end if
    end function geninvgauss_isf

    pure subroutine geninvgauss_logpdf_derivatives(z, p, b, logf, dp_shape, db_shape, dz)
        real(dp), intent(in) :: z !! positive standardized observation
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: dp_shape !! derivative with respect to p
        real(dp), intent(out) :: db_shape !! derivative with respect to b
        real(dp), intent(out) :: dz !! derivative with respect to z
        if (.not. valid_shapes(p, b) .or. z <= 0.0_dp .or. .not. ieee_is_finite(z)) then
            logf = quiet_nan(z)
            dp_shape = logf; db_shape = logf; dz = logf
            return
        end if
        logf = (p - 1.0_dp) * log(z) - 0.5_dp * b * (z + 1.0_dp / z) - &
            log(2.0_dp) - log_besselk(p, b)
        dp_shape = log(z) - log_besselk_order_derivative(p, b)
        db_shape = -0.5_dp * (z + 1.0_dp / z) - besselk_log_derivative_x(p, b)
        dz = (p - 1.0_dp) / z - 0.5_dp * b * (1.0_dp - 1.0_dp / (z * z))
    end subroutine geninvgauss_logpdf_derivatives

    pure function geninvgauss_tail_y(y, p, b, upper) result(probability)
        real(dp), intent(in) :: y !! log of the positive standardized variate
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        logical, intent(in) :: upper !! true for upper tail and false for lower tail
        real(dp) :: probability, log_probability
        log_probability = geninvgauss_logtail_y(y, p, b, upper)
        if (log_probability == negative_infinity(log_probability)) then
            probability = 0.0_dp
        else
            probability = exp(log_probability)
        end if
        probability = min(1.0_dp, max(0.0_dp, probability))
    end function geninvgauss_tail_y

    pure function geninvgauss_logtail_y(y, p, b, upper) result(log_probability)
        real(dp), intent(in) :: y !! log of the positive standardized variate
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        logical, intent(in) :: upper !! true for upper tail and false for lower tail
        real(dp) :: log_probability, bound, mode_y, a, c
        bound = integration_bound(p, b)
        mode_y = asinh(p / b)
        if (upper) then
            if (y >= mode_y) then
                if (y >= 100.0_dp) then
                    log_probability = negative_infinity(y)
                    return
                end if
                a = y
                c = geninvgauss_tail_cutoff(y, p, b, 1)
            else
                if (y <= -bound) then
                    log_probability = 0.0_dp
                    return
                end if
                a = y
                c = bound
            end if
        else
            if (y <= mode_y) then
                if (y <= -100.0_dp) then
                    log_probability = negative_infinity(y)
                    return
                end if
                a = geninvgauss_tail_cutoff(y, p, b, -1)
                c = y
            else
                if (y >= bound) then
                    log_probability = 0.0_dp
                    return
                end if
                a = -bound
                c = y
            end if
        end if
        log_probability = integrate_log_y_density(a, c, p, b)
        log_probability = min(0.0_dp, log_probability)
    end function geninvgauss_logtail_y

    pure function geninvgauss_tail_cutoff(y, p, b, direction) result(cutoff)
        real(dp), intent(in) :: y !! tail endpoint in log coordinate
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        integer, intent(in) :: direction !! -1 for lower tail or +1 for upper tail
        real(dp) :: cutoff, target, step, near, far, mid
        integer :: iter
        target = geninvgauss_log_kernel(y, p, b) - 45.0_dp
        step = 0.5_dp
        near = y
        far = y + real(direction, dp) * step
        do iter = 1, 20
            if (abs(far) >= 100.0_dp) then
                far = sign(100.0_dp, far)
                exit
            end if
            if (geninvgauss_log_kernel(far, p, b) <= target) exit
            step = 2.0_dp * step
            far = y + real(direction, dp) * step
        end do
        if (geninvgauss_log_kernel(far, p, b) > target) then
            cutoff = far
            return
        end if
        do iter = 1, 50
            mid = 0.5_dp * (near + far)
            if (geninvgauss_log_kernel(mid, p, b) > target) then
                near = mid
            else
                far = mid
            end if
        end do
        cutoff = far
    end function geninvgauss_tail_cutoff

    pure elemental function geninvgauss_log_kernel(y, p, b) result(value)
        real(dp), intent(in) :: y !! log-coordinate value
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp) :: value
        value = p * y - b * cosh(y)
    end function geninvgauss_log_kernel

    pure function integrate_log_y_density(a, bnd, p, b) result(log_value)
        real(dp), intent(in) :: a !! lower integration bound in log coordinate
        real(dp), intent(in) :: bnd !! upper integration bound in log coordinate
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp) :: log_value, norm, left, right, mid, half, y1, y2, log_weight
        real(dp) :: segment_width
        integer :: nseg, j, i
        if (bnd <= a) then
            log_value = negative_infinity(a)
            return
        end if
        norm = log(2.0_dp) + log_besselk(p, b)
        segment_width = min(0.5_dp, 4.0_dp / sqrt(max(1.0_dp, hypot(p, b))))
        nseg = max(1, ceiling((bnd - a) / segment_width))
        log_value = negative_infinity(a)
        do j = 1, nseg
            left = a + (bnd - a) * real(j - 1, dp) / real(nseg, dp)
            right = a + (bnd - a) * real(j, dp) / real(nseg, dp)
            mid = 0.5_dp * (left + right)
            half = 0.5_dp * (right - left)
            do i = 1, 8
                log_weight = log(half * gl_w(i))
                y1 = mid - half * gl_x(i)
                y2 = mid + half * gl_x(i)
                log_value = logaddexp_pair(log_value, &
                    log_weight + geninvgauss_log_kernel(y1, p, b) - norm)
                log_value = logaddexp_pair(log_value, &
                    log_weight + geninvgauss_log_kernel(y2, p, b) - norm)
            end do
        end do
    end function integrate_log_y_density

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

    pure elemental function integration_bound(p, b) result(bound)
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        real(dp) :: bound, mode_y
        mode_y = asinh(p / b)
        bound = abs(mode_y) + max(12.0_dp, log(2.0_dp / min(b, 1.0_dp)) + 8.0_dp)
        bound = min(100.0_dp, bound)
    end function integration_bound

    pure elemental logical function valid_shapes(p, b)
        real(dp), intent(in) :: p !! real shape parameter
        real(dp), intent(in) :: b !! strictly positive shape parameter
        valid_shapes = ieee_is_finite(p) .and. ieee_is_finite(b) .and. b > 0.0_dp
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

end module scifort_geninvgauss
