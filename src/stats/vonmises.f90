! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Circular von Mises distribution matching scipy.stats.vonmises.
module scifort_vonmises
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan, valid_loc_scale
    use scifort_special_elementary, only : i0e, i1e
    implicit none
    private

    public :: vonmises_cdf, vonmises_isf, vonmises_logcdf, vonmises_logpdf
    public :: vonmises_logsf, vonmises_pdf, vonmises_ppf, vonmises_sf
    public :: vonmises_i1_i0_ratio, vonmises_standard_cdf, vonmises_standard_logpdf
    public :: vonmises_line_cdf_standard, vonmises_line_ppf_standard

contains

    pure elemental function vonmises_standard_logpdf(z, kappa) result(y)
        real(dp), intent(in) :: z !! standardized angular variate
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp) :: y

        if (.not. valid_shape(kappa) .or. ieee_is_nan(z)) then
            y = quiet_nan(z)
        else
            y = kappa * (cos(z) - 1.0_dp) - log(2.0_dp * scifort_pi) - log(i0e(kappa))
        end if
    end function vonmises_standard_logpdf

    pure elemental function vonmises_logpdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(kappa) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            y = vonmises_standard_logpdf(z, kappa) - log(sigma)
        end if
    end function vonmises_logpdf

    pure elemental function vonmises_pdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y, ly

        ly = vonmises_logpdf(x, kappa, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else
            y = exp(ly)
        end if
    end function vonmises_pdf

    pure elemental function vonmises_standard_cdf(z, kappa) result(p)
        real(dp), intent(in) :: z !! standardized angular variate on the real line
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp) :: p
        real(dp) :: period, reduced, turns

        if (.not. valid_shape(kappa) .or. ieee_is_nan(z)) then
            p = quiet_nan(z)
        else if (.not. ieee_is_finite(z)) then
            if (z > 0.0_dp) then
                p = huge(1.0_dp)
            else
                p = -huge(1.0_dp)
            end if
        else
            period = 2.0_dp * scifort_pi
            turns = floor((z + scifort_pi) / period)
            reduced = z - turns * period
            if (reduced >= scifort_pi) then
                reduced = reduced - period
                turns = turns + 1.0_dp
            else if (reduced < -scifort_pi) then
                reduced = reduced + period
                turns = turns - 1.0_dp
            end if
            p = turns + vonmises_line_cdf_standard(reduced, kappa)
        end if
    end function vonmises_standard_cdf

    pure elemental function vonmises_cdf(x, kappa, loc, scale) result(p)
        real(dp), intent(in) :: x !! point x of the lower-tail circular CDF
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: p
        real(dp) :: mu, sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(kappa) .or. .not. valid_loc_scale(mu, sigma)) then
            p = quiet_nan(x)
        else
            p = vonmises_standard_cdf((x - mu) / sigma, kappa)
        end if
    end function vonmises_cdf

    pure elemental function vonmises_sf(x, kappa, loc, scale) result(p)
        real(dp), intent(in) :: x !! point x of one minus the circular CDF
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: p

        p = 1.0_dp - vonmises_cdf(x, kappa, loc, scale)
    end function vonmises_sf

    pure elemental function vonmises_logcdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log of the circular CDF
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y, p

        p = vonmises_cdf(x, kappa, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else if (p < 0.0_dp) then
            y = quiet_nan(p)
        else
            y = log(p)
        end if
    end function vonmises_logcdf

    pure elemental function vonmises_logsf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log of one minus the circular CDF
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y, p

        p = vonmises_sf(x, kappa, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else if (p < 0.0_dp) then
            y = quiet_nan(p)
        else
            y = log(p)
        end if
    end function vonmises_logsf

    pure elemental function vonmises_ppf(p, kappa, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(kappa) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            z = vonmises_line_ppf_standard(p, kappa)
            x = mu + sigma * z
        end if
    end function vonmises_ppf

    pure elemental function vonmises_isf(p, kappa, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: x

        x = vonmises_ppf(1.0_dp - p, kappa, loc, scale)
    end function vonmises_isf

    pure elemental function vonmises_i1_i0_ratio(kappa) result(r)
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp) :: r

        if (.not. valid_shape(kappa)) then
            r = quiet_nan(kappa)
        else if (kappa == 0.0_dp) then
            r = 0.0_dp
        else
            r = i1e(kappa) / i0e(kappa)
        end if
    end function vonmises_i1_i0_ratio

    pure elemental function vonmises_line_cdf_standard(z, kappa) result(p)
        real(dp), intent(in) :: z !! standardized variate in [-pi,pi]
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp) :: p
        real(dp) :: area, az, fa, fb, fm, whole

        if (z <= -scifort_pi) then
            p = 0.0_dp
        else if (z >= scifort_pi) then
            p = 1.0_dp
        else if (z == 0.0_dp) then
            p = 0.5_dp
        else if (kappa == 0.0_dp) then
            p = 0.5_dp + z / (2.0_dp * scifort_pi)
        else
            az = abs(z)
            fa = standard_pdf(0.0_dp, kappa)
            fb = standard_pdf(az, kappa)
            fm = standard_pdf(0.5_dp * az, kappa)
            whole = az * (fa + 4.0_dp * fm + fb) / 6.0_dp
            area = adaptive_simpson(0.0_dp, az, fa, fb, fm, whole, &
                2.0e-14_dp, 28, kappa)
            if (z > 0.0_dp) then
                p = 0.5_dp + area
            else
                p = 0.5_dp - area
            end if
            p = min(1.0_dp, max(0.0_dp, p))
        end if
    end function vonmises_line_cdf_standard

    pure elemental function vonmises_line_ppf_standard(p, kappa) result(z)
        real(dp), intent(in) :: p !! probability in [0,1]
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp) :: z
        real(dp) :: lo, hi, mid
        integer :: iter

        if (p == 0.0_dp) then
            z = -scifort_pi
        else if (p == 1.0_dp) then
            z = scifort_pi
        else if (p == 0.5_dp) then
            z = 0.0_dp
        else
            lo = -scifort_pi
            hi = scifort_pi
            do iter = 1, 75
                mid = 0.5_dp * (lo + hi)
                if (vonmises_line_cdf_standard(mid, kappa) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            z = 0.5_dp * (lo + hi)
        end if
    end function vonmises_line_ppf_standard

    pure recursive function adaptive_simpson(a, b, fa, fb, fm, whole, tol, depth, kappa) result(value)
        real(dp), intent(in) :: a !! left integration endpoint
        real(dp), intent(in) :: b !! right integration endpoint
        real(dp), intent(in) :: fa !! density at a
        real(dp), intent(in) :: fb !! density at b
        real(dp), intent(in) :: fm !! density at the midpoint
        real(dp), intent(in) :: whole !! Simpson estimate on [a,b]
        real(dp), intent(in) :: tol !! absolute error target
        integer, intent(in) :: depth !! remaining recursion depth
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp) :: value
        real(dp) :: m, lm, rm, fl, fr, left, right, delta

        m = 0.5_dp * (a + b)
        lm = 0.5_dp * (a + m)
        rm = 0.5_dp * (m + b)
        fl = standard_pdf(lm, kappa)
        fr = standard_pdf(rm, kappa)
        left = (m - a) * (fa + 4.0_dp * fl + fm) / 6.0_dp
        right = (b - m) * (fm + 4.0_dp * fr + fb) / 6.0_dp
        delta = left + right - whole
        if (depth <= 0 .or. abs(delta) <= 15.0_dp * tol) then
            value = left + right + delta / 15.0_dp
        else
            value = adaptive_simpson(a, m, fa, fm, fl, left, 0.5_dp * tol, depth - 1, kappa) + &
                adaptive_simpson(m, b, fm, fb, fr, right, 0.5_dp * tol, depth - 1, kappa)
        end if
    end function adaptive_simpson

    pure elemental function standard_pdf(z, kappa) result(y)
        real(dp), intent(in) :: z !! standardized angular variate
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp) :: y

        y = exp(vonmises_standard_logpdf(z, kappa))
    end function standard_pdf

    pure elemental logical function valid_shape(kappa) result(ok)
        real(dp), intent(in) :: kappa !! candidate concentration parameter
        ok = ieee_is_finite(kappa) .and. kappa >= 0.0_dp
    end function valid_shape

    pure elemental subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location parameter
        real(dp), intent(in), optional :: scale !! optional positive scale parameter
        real(dp), intent(out) :: mu !! resolved location parameter
        real(dp), intent(out) :: sigma !! resolved scale parameter

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_vonmises
