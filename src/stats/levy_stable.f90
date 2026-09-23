! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2001-2002 Enthought, Inc. 2003-2026, SciPy Developers.
! Copyright (c) 2026 SciFort contributors
!
! Levy-stable distribution matching scipy.stats.levy_stable in its default S1
! parameterization.  The S1-to-S0 conversion, difficult-input rounding, special
! cases, and Nolan/Zolotarev piecewise formulas are adapted from SciPy 1.17.0
! scipy.stats._levy_stable and its levyst.c helper.  SciFort uses its own
! Gauss-Legendre/adaptive quadrature implementation around those formulas.
! See THIRD_PARTY_LICENSES.md and CODE_PROVENANCE.md.
module scifort_levy_stable
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_pdf
    implicit none
    private

    public :: levy_stable_cdf, levy_stable_isf, levy_stable_logcdf
    public :: levy_stable_logpdf, levy_stable_logsf, levy_stable_pdf
    public :: levy_stable_ppf, levy_stable_sf

    real(dp), parameter :: half_pi = 0.5_dp * scifort_pi
    real(dp), parameter :: inv_pi = 1.0_dp / scifort_pi
    real(dp), parameter :: sqrt_two = 1.414213562373095048801688724209698079_dp
    real(dp), parameter :: sqrt_two_pi = 2.506628274631000502415765284811045253_dp
    real(dp), parameter :: x_tol_near_zeta = 0.005_dp
    real(dp), parameter :: alpha_tol_near_one = 0.005_dp
    real(dp), parameter :: quad_rel_tol = 2.0e-13_dp
    integer, parameter :: max_quad_depth = 11

    ! Positive abscissas and weights of the 16-point Gauss-Legendre rule.
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

    pure elemental function levy_stable_pdf(x, alpha, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(alpha, beta) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = standardized_x(x, alpha, beta, mu, sigma)
            y = standard_pdf_s1(z, alpha, beta) / sigma
        end if
    end function levy_stable_pdf

    pure elemental function levy_stable_logpdf(x, alpha, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, f

        f = levy_stable_pdf(x, alpha, beta, loc, scale)
        if (ieee_is_nan(f)) then
            y = f
        else if (f <= 0.0_dp) then
            y = negative_infinity(f)
        else
            y = log(f)
        end if
    end function levy_stable_logpdf

    pure elemental function levy_stable_cdf(x, alpha, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(alpha, beta) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = 0.0_dp
            else
                y = 1.0_dp
            end if
        else
            z = standardized_x(x, alpha, beta, mu, sigma)
            y = standard_cdf_s1(z, alpha, beta)
        end if
    end function levy_stable_cdf

    pure elemental function levy_stable_sf(x, alpha, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(alpha, beta) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = 1.0_dp
            else
                y = 0.0_dp
            end if
        else
            z = standardized_x(x, alpha, beta, mu, sigma)
            ! Stable laws obey X(alpha,beta) == -X(alpha,-beta) in S1.
            y = standard_cdf_s1(-z, alpha, -beta)
        end if
    end function levy_stable_sf

    pure elemental function levy_stable_logcdf(x, alpha, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p, q

        p = levy_stable_cdf(x, alpha, beta, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else if (p < 0.5_dp) then
            y = log(p)
        else
            q = levy_stable_sf(x, alpha, beta, loc, scale)
            y = log1p_safe(-q)
        end if
    end function levy_stable_logcdf

    pure elemental function levy_stable_logsf(x, alpha, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p, q

        q = levy_stable_sf(x, alpha, beta, loc, scale)
        if (ieee_is_nan(q)) then
            y = q
        else if (q <= 0.0_dp) then
            y = negative_infinity(q)
        else if (q < 0.5_dp) then
            y = log(q)
        else
            p = levy_stable_cdf(x, alpha, beta, loc, scale)
            y = log1p_safe(-p)
        end if
    end function levy_stable_logsf

    pure elemental function levy_stable_ppf(probability, alpha, beta, loc, scale) result(x)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(alpha, beta) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            x = negative_infinity(probability)
        else if (probability >= 1.0_dp) then
            x = positive_infinity(probability)
        else
            if (probability <= 0.5_dp) then
                z = standard_lower_quantile(probability, alpha, beta)
            else
                z = -standard_lower_quantile(1.0_dp - probability, alpha, -beta)
            end if
            x = unstandardized_x(z, alpha, beta, mu, sigma)
        end if
    end function levy_stable_ppf

    pure elemental function levy_stable_isf(probability, alpha, beta, loc, scale) result(x)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(alpha, beta) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            x = negative_infinity(probability)
        else if (probability <= 0.0_dp) then
            x = positive_infinity(probability)
        else
            if (probability <= 0.5_dp) then
                z = -standard_lower_quantile(probability, alpha, -beta)
            else
                z = standard_lower_quantile(1.0_dp - probability, alpha, beta)
            end if
            x = unstandardized_x(z, alpha, beta, mu, sigma)
        end if
    end function levy_stable_isf

    pure function standard_pdf_s1(x, alpha, beta) result(value)
        real(dp), intent(in) :: x !! standardized S1 variate
        real(dp), intent(in) :: alpha !! stability parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness parameter in [-1,1]
        real(dp) :: value, zeta, x0

        if (alpha == 2.0_dp) then
            value = normal_pdf(x / sqrt_two) / sqrt_two
            return
        else if (alpha == 1.0_dp .and. beta == 0.0_dp) then
            value = inv_pi / (1.0_dp + x * x)
            return
        end if
        zeta = -beta * tan(half_pi * alpha)
        if (alpha == 1.0_dp) then
            x0 = x
        else
            x0 = x + zeta
        end if
        value = standard_pdf_s0(x0, alpha, beta)
    end function standard_pdf_s1

    pure function standard_cdf_s1(x, alpha, beta) result(value)
        real(dp), intent(in) :: x !! standardized S1 variate
        real(dp), intent(in) :: alpha !! stability parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness parameter in [-1,1]
        real(dp) :: value, zeta, x0

        if (alpha == 2.0_dp) then
            value = normal_cdf(x / sqrt_two)
            return
        else if (alpha == 1.0_dp .and. beta == 0.0_dp) then
            value = 0.5_dp + atan(x) * inv_pi
            return
        end if
        zeta = -beta * tan(half_pi * alpha)
        if (alpha == 1.0_dp) then
            x0 = x
        else
            x0 = x + zeta
        end if
        value = standard_cdf_s0(x0, alpha, beta)
    end function standard_cdf_s1

    pure recursive function standard_pdf_s0(x0_in, alpha_in, beta_in) result(value)
        real(dp), intent(in) :: x0_in !! standardized S0/Zolotarev-M variate
        real(dp), intent(in) :: alpha_in !! input stability parameter in (0,2]
        real(dp), intent(in) :: beta_in !! input skewness parameter in [-1,1]
        real(dp) :: value
        real(dp) :: x0, alpha, beta, zeta, xi, c2, integral_value, xp

        alpha = alpha_in
        beta = beta_in
        zeta = -beta * tan(half_pi * alpha)
        x0 = x0_in
        call round_difficult_input(x0, alpha, zeta)

        if (alpha == 2.0_dp) then
            value = normal_pdf(x0 / sqrt_two) / sqrt_two
            return
        else if (alpha == 0.5_dp .and. beta == 1.0_dp) then
            xp = x0 + 1.0_dp
            if (xp <= 0.0_dp) then
                value = 0.0_dp
            else
                value = exp(-0.5_dp / xp) / (sqrt_two_pi * xp**1.5_dp)
            end if
            return
        else if (alpha == 1.0_dp .and. beta == 0.0_dp) then
            value = inv_pi / (1.0_dp + x0 * x0)
            return
        end if

        call nolan_constants(alpha, beta, x0, zeta, xi, c2=c2)
        call round_x_near_zeta(x0, alpha, zeta)
        if (x0 == zeta) then
            value = exp(log_gamma(1.0_dp + 1.0_dp / alpha) + log(max(tiny(1.0_dp), cos(xi))) &
                - log(scifort_pi) - 0.5_dp * log1p_safe(zeta * zeta) / alpha)
            return
        else if (x0 < zeta) then
            value = standard_pdf_s0(-x0, alpha, -beta)
            return
        end if

        if (abs((-xi) - half_pi) <= 2.0e-14_dp) then
            value = 0.0_dp
            return
        end if
        integral_value = nolan_integral(alpha, beta, x0, xi, .true.)
        value = max(0.0_dp, c2 * integral_value)
    end function standard_pdf_s0

    pure recursive function standard_cdf_s0(x0_in, alpha_in, beta_in) result(value)
        real(dp), intent(in) :: x0_in !! standardized S0/Zolotarev-M variate
        real(dp), intent(in) :: alpha_in !! input stability parameter in (0,2]
        real(dp), intent(in) :: beta_in !! input skewness parameter in [-1,1]
        real(dp) :: value
        real(dp) :: x0, alpha, beta, zeta, xi, c1, c3, integral_value, xp

        alpha = alpha_in
        beta = beta_in
        zeta = -beta * tan(half_pi * alpha)
        x0 = x0_in
        call round_difficult_input(x0, alpha, zeta)

        if (alpha == 2.0_dp) then
            value = normal_cdf(x0 / sqrt_two)
            return
        else if (alpha == 0.5_dp .and. beta == 1.0_dp) then
            xp = x0 + 1.0_dp
            if (xp <= 0.0_dp) then
                value = 0.0_dp
            else
                value = erfc(sqrt(0.5_dp / xp))
            end if
            return
        else if (alpha == 1.0_dp .and. beta == 0.0_dp) then
            value = 0.5_dp + atan(x0) * inv_pi
            return
        end if

        call nolan_constants(alpha, beta, x0, zeta, xi, c1=c1, c3=c3)
        call round_x_near_zeta(x0, alpha, zeta)
        if ((alpha == 1.0_dp .and. beta < 0.0_dp) .or. x0 < zeta) then
            value = 1.0_dp - standard_cdf_s0(-x0, alpha, -beta)
            value = min(1.0_dp, max(0.0_dp, value))
            return
        else if (x0 == zeta) then
            value = 0.5_dp - xi * inv_pi
            return
        end if

        if (abs((-xi) - half_pi) <= 2.0e-14_dp) then
            value = min(1.0_dp, max(0.0_dp, c1))
            return
        end if
        integral_value = nolan_integral(alpha, beta, x0, xi, .false.)
        value = min(1.0_dp, max(0.0_dp, c1 + c3 * integral_value))
    end function standard_cdf_s0

    pure subroutine round_difficult_input(x0, alpha, zeta)
        real(dp), intent(inout) :: x0 !! S0 variate, rounded near zeta when needed
        real(dp), intent(inout) :: alpha !! stability parameter, rounded near one
        real(dp), intent(in) :: zeta !! zeta computed before alpha rounding

        if (abs(alpha - 1.0_dp) < alpha_tol_near_one) alpha = 1.0_dp
        call round_x_near_zeta(x0, alpha, zeta)
    end subroutine round_difficult_input

    pure subroutine round_x_near_zeta(x0, alpha, zeta)
        real(dp), intent(inout) :: x0 !! S0 variate to condition near zeta
        real(dp), intent(in) :: alpha !! positive stability parameter
        real(dp), intent(in) :: zeta !! Nolan zeta location
        real(dp) :: tolerance_scale

        tolerance_scale = alpha**(1.0_dp / alpha)
        if (abs(x0 - zeta) < x_tol_near_zeta * tolerance_scale) x0 = zeta
    end subroutine round_x_near_zeta

    pure subroutine nolan_constants(alpha, beta, x0, zeta, xi, c1, c2, c3)
        real(dp), intent(in) :: alpha !! rounded stability parameter
        real(dp), intent(in) :: beta !! skewness parameter
        real(dp), intent(in) :: x0 !! S0 variate
        real(dp), intent(out) :: zeta !! Nolan zeta
        real(dp), intent(out) :: xi !! Nolan xi angle
        real(dp), intent(out), optional :: c1 !! CDF additive constant
        real(dp), intent(out), optional :: c2 !! PDF multiplicative constant
        real(dp), intent(out), optional :: c3 !! CDF integral multiplier

        zeta = -beta * tan(half_pi * alpha)
        if (alpha == 1.0_dp) then
            xi = half_pi
            if (present(c1)) c1 = 0.0_dp
            if (present(c2)) c2 = 0.5_dp / abs(beta)
            if (present(c3)) c3 = inv_pi
        else
            xi = atan(-zeta) / alpha
            if (alpha < 1.0_dp) then
                if (present(c1)) c1 = 0.5_dp - xi * inv_pi
                if (present(c3)) c3 = inv_pi
            else
                if (present(c1)) c1 = 1.0_dp
                if (present(c3)) c3 = -inv_pi
            end if
            if (present(c2)) c2 = alpha * inv_pi / (abs(alpha - 1.0_dp) * (x0 - zeta))
        end if
    end subroutine nolan_constants

    pure function nolan_integral(alpha, beta, x0, xi, want_pdf) result(value)
        real(dp), intent(in) :: alpha !! rounded stability parameter
        real(dp), intent(in) :: beta !! skewness parameter after any reflection
        real(dp), intent(in) :: x0 !! S0 variate on Nolan's x0 > zeta branch
        real(dp), intent(in) :: xi !! Nolan xi angle
        logical, intent(in) :: want_pdf !! true integrates g*exp(-g), false exp(-g)
        real(dp) :: value
        real(dp) :: points(16), levels(10), left, right, root, tol
        integer :: npoints, i

        left = -xi
        right = half_pi
        if (right <= left) then
            value = 0.0_dp
            return
        end if

        levels = [1.0e-8_dp, 1.0e-5_dp, 1.0e-3_dp, 1.0e-2_dp, 1.0e-1_dp, &
            1.0_dp, 5.0_dp, 10.0_dp, 50.0_dp, 100.0_dp]
        npoints = 2
        points(1) = left
        points(2) = right
        if (left < 0.0_dp .and. 0.0_dp < right) call add_point(points, npoints, 0.0_dp)
        do i = 1, size(levels)
            root = nolan_level_root(alpha, beta, x0, xi, log(levels(i)))
            if (root > left .and. root < right) call add_point(points, npoints, root)
        end do
        call sort_points(points, npoints)

        value = 0.0_dp
        tol = quad_rel_tol / real(max(1, npoints - 1), dp)
        do i = 1, npoints - 1
            if (points(i + 1) > points(i)) then
                value = value + adaptive_nolan(points(i), points(i + 1), alpha, beta, x0, xi, &
                    want_pdf, tol, 0)
            end if
        end do
    end function nolan_integral

    pure function nolan_level_root(alpha, beta, x0, xi, target_log_g) result(root)
        real(dp), intent(in) :: alpha !! rounded stability parameter
        real(dp), intent(in) :: beta !! skewness parameter after any reflection
        real(dp), intent(in) :: x0 !! S0 variate on Nolan's x0 > zeta branch
        real(dp), intent(in) :: xi !! Nolan xi angle
        real(dp), intent(in) :: target_log_g !! logarithm of desired g level
        real(dp) :: root, lo, hi, mid, lg
        logical :: increasing
        integer :: iter

        lo = -xi
        hi = half_pi
        increasing = alpha <= 1.0_dp
        do iter = 1, 90
            mid = 0.5_dp * (lo + hi)
            lg = nolan_log_g(mid, alpha, beta, x0)
            if (increasing) then
                if (lg < target_log_g) then
                    lo = mid
                else
                    hi = mid
                end if
            else
                if (lg > target_log_g) then
                    lo = mid
                else
                    hi = mid
                end if
            end if
        end do
        root = 0.5_dp * (lo + hi)
    end function nolan_level_root

    pure recursive function adaptive_nolan(a, b, alpha, beta, x0, xi, want_pdf, tol, depth) result(value)
        real(dp), intent(in) :: a !! left angular integration endpoint
        real(dp), intent(in) :: b !! right angular integration endpoint
        real(dp), intent(in) :: alpha !! rounded stability parameter
        real(dp), intent(in) :: beta !! skewness parameter after reflection
        real(dp), intent(in) :: x0 !! S0 variate
        real(dp), intent(in) :: xi !! Nolan xi angle (retained for provenance/debug symmetry)
        logical, intent(in) :: want_pdf !! true for PDF integrand
        real(dp), intent(in) :: tol !! local relative/absolute quadrature target
        integer, intent(in) :: depth !! current adaptive recursion depth
        real(dp) :: value, whole, left_value, right_value, split, refined, error_est

        if (b <= a) then
            value = 0.0_dp
            return
        end if
        split = 0.5_dp * (a + b)
        whole = gl16_nolan(a, b, alpha, beta, x0, want_pdf)
        left_value = gl16_nolan(a, split, alpha, beta, x0, want_pdf)
        right_value = gl16_nolan(split, b, alpha, beta, x0, want_pdf)
        refined = left_value + right_value
        error_est = abs(refined - whole)
        if (depth >= max_quad_depth .or. error_est <= tol * (1.0_dp + abs(refined))) then
            value = refined
        else
            value = adaptive_nolan(a, split, alpha, beta, x0, xi, want_pdf, 0.5_dp * tol, depth + 1) + &
                adaptive_nolan(split, b, alpha, beta, x0, xi, want_pdf, 0.5_dp * tol, depth + 1)
        end if
    end function adaptive_nolan

    pure function gl16_nolan(a, b, alpha, beta, x0, want_pdf) result(value)
        real(dp), intent(in) :: a !! left angular integration endpoint
        real(dp), intent(in) :: b !! right angular integration endpoint
        real(dp), intent(in) :: alpha !! rounded stability parameter
        real(dp), intent(in) :: beta !! skewness parameter after reflection
        real(dp), intent(in) :: x0 !! S0 variate
        logical, intent(in) :: want_pdf !! true for PDF integrand
        real(dp) :: value, mid, half, t1, t2
        integer :: i

        mid = 0.5_dp * (a + b)
        half = 0.5_dp * (b - a)
        value = 0.0_dp
        do i = 1, 8
            t1 = mid - half * gl_x(i)
            t2 = mid + half * gl_x(i)
            value = value + gl_w(i) * nolan_integrand(t1, alpha, beta, x0, want_pdf)
            value = value + gl_w(i) * nolan_integrand(t2, alpha, beta, x0, want_pdf)
        end do
        value = half * value
    end function gl16_nolan

    pure function nolan_integrand(theta, alpha, beta, x0, want_pdf) result(value)
        real(dp), intent(in) :: theta !! Nolan angular integration coordinate
        real(dp), intent(in) :: alpha !! rounded stability parameter
        real(dp), intent(in) :: beta !! skewness parameter after reflection
        real(dp), intent(in) :: x0 !! S0 variate
        logical, intent(in) :: want_pdf !! true returns g*exp(-g), false exp(-g)
        real(dp) :: value, lg, g

        lg = nolan_log_g(theta, alpha, beta, x0)
        if (want_pdf) then
            if (lg > 6.0_dp .or. lg < log(tiny(1.0_dp))) then
                value = 0.0_dp
            else
                g = exp(lg)
                value = g * exp(-g)
            end if
        else
            if (lg > 6.0_dp) then
                value = 0.0_dp
            else if (lg < -36.0_dp) then
                value = 1.0_dp
            else
                g = exp(lg)
                value = exp(-g)
            end if
        end if
    end function nolan_integrand

    pure function nolan_log_g(theta, alpha, beta, x0) result(value)
        real(dp), intent(in) :: theta !! interior Nolan angular coordinate
        real(dp), intent(in) :: alpha !! rounded stability parameter
        real(dp), intent(in) :: beta !! skewness parameter after reflection
        real(dp), intent(in) :: x0 !! S0 variate on x0 > zeta branch
        real(dp) :: value
        real(dp) :: zeta, aangle, d, s1, ctheta, c2theta, exponent
        real(dp) :: two_beta_over_pi, first_factor

        if (alpha == 1.0_dp) then
            two_beta_over_pi = 2.0_dp * beta * inv_pi
            first_factor = 1.0_dp + theta * two_beta_over_pi
            ctheta = cos(theta)
            if (first_factor <= 0.0_dp .or. ctheta <= 0.0_dp) then
                value = negative_infinity(theta)
                return
            end if
            value = log(first_factor) + (half_pi / beta + theta) * tan(theta) &
                - x0 / two_beta_over_pi - log(ctheta)
            return
        end if

        zeta = -beta * tan(half_pi * alpha)
        aangle = atan(-zeta)
        d = x0 - zeta
        s1 = sin(aangle + alpha * theta)
        ctheta = cos(theta)
        c2theta = cos(aangle + (alpha - 1.0_dp) * theta)
        if (d <= 0.0_dp .or. s1 <= 0.0_dp .or. ctheta <= 0.0_dp .or. c2theta <= 0.0_dp) then
            if (alpha < 1.0_dp) then
                value = negative_infinity(theta)
            else
                value = positive_infinity(theta)
            end if
            return
        end if
        exponent = alpha / (alpha - 1.0_dp)
        value = -0.5_dp * log1p_safe(zeta * zeta) / (alpha - 1.0_dp) + &
            exponent * (log(ctheta) + log(d) - log(s1)) + log(c2theta) - log(ctheta)
    end function nolan_log_g

    pure function standard_lower_quantile(probability, alpha, beta) result(x)
        real(dp), intent(in) :: probability !! lower-tail probability in (0,0.5]
        real(dp), intent(in) :: alpha !! stability parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness parameter in [-1,1]
        real(dp) :: x, lo, hi, mid, p
        integer :: iter

        lo = -1.0_dp
        hi = 1.0_dp
        p = standard_cdf_s1(lo, alpha, beta)
        do iter = 1, 256
            if (p <= probability) exit
            hi = lo
            lo = 2.0_dp * lo
            if (.not. ieee_is_finite(lo)) then
                x = negative_infinity(probability)
                return
            end if
            p = standard_cdf_s1(lo, alpha, beta)
        end do
        p = standard_cdf_s1(hi, alpha, beta)
        do iter = 1, 256
            if (p >= probability) exit
            lo = hi
            hi = 2.0_dp * hi
            if (.not. ieee_is_finite(hi)) then
                x = positive_infinity(probability)
                return
            end if
            p = standard_cdf_s1(hi, alpha, beta)
        end do

        do iter = 1, 72
            mid = lo + 0.5_dp * (hi - lo)
            if (standard_cdf_s1(mid, alpha, beta) < probability) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x = lo + 0.5_dp * (hi - lo)
    end function standard_lower_quantile

    pure elemental function standardized_x(x, alpha, beta, loc, scale) result(z)
        real(dp), intent(in) :: x !! physical variate
        real(dp), intent(in) :: alpha !! stability parameter
        real(dp), intent(in) :: beta !! skewness parameter
        real(dp), intent(in) :: loc !! S1 location
        real(dp), intent(in) :: scale !! positive scale
        real(dp) :: z, shift

        shift = 0.0_dp
        if (alpha == 1.0_dp) shift = 2.0_dp * beta * scale * log(scale) * inv_pi
        z = (x - loc - shift) / scale
    end function standardized_x

    pure elemental function unstandardized_x(z, alpha, beta, loc, scale) result(x)
        real(dp), intent(in) :: z !! standardized S1 variate
        real(dp), intent(in) :: alpha !! stability parameter
        real(dp), intent(in) :: beta !! skewness parameter
        real(dp), intent(in) :: loc !! S1 location
        real(dp), intent(in) :: scale !! positive scale
        real(dp) :: x, shift

        shift = 0.0_dp
        if (alpha == 1.0_dp) shift = 2.0_dp * beta * scale * log(scale) * inv_pi
        x = loc + shift + scale * z
    end function unstandardized_x

    pure subroutine add_point(points, npoints, x)
        real(dp), intent(inout) :: points(:) !! point buffer to append to
        integer, intent(inout) :: npoints !! current number of buffered points
        real(dp), intent(in) :: x !! point to append if capacity permits

        if (npoints < size(points)) then
            npoints = npoints + 1
            points(npoints) = x
        end if
    end subroutine add_point

    pure subroutine sort_points(points, npoints)
        real(dp), intent(inout) :: points(:) !! point buffer whose active prefix is sorted
        integer, intent(in) :: npoints !! active number of points
        real(dp) :: key
        integer :: i, j

        do i = 2, npoints
            key = points(i)
            j = i - 1
            do while (j >= 1)
                if (points(j) <= key) exit
                points(j + 1) = points(j)
                j = j - 1
            end do
            points(j + 1) = key
        end do
    end subroutine sort_points

    pure elemental logical function valid_shapes(alpha, beta) result(ok)
        real(dp), intent(in) :: alpha !! candidate stability parameter
        real(dp), intent(in) :: beta !! candidate skewness parameter
        ok = ieee_is_finite(alpha) .and. ieee_is_finite(beta) .and. &
            alpha > 0.0_dp .and. alpha <= 2.0_dp .and. beta >= -1.0_dp .and. beta <= 1.0_dp
    end function valid_shapes

    pure elemental logical function valid_probability(p) result(ok)
        real(dp), intent(in) :: p !! probability to validate
        ok = .not. ieee_is_nan(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional S1 location parameter
        real(dp), intent(in), optional :: scale !! optional positive scale parameter
        real(dp), intent(out) :: mu !! resolved location parameter
        real(dp), intent(out) :: sigma !! resolved scale parameter

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_levy_stable
