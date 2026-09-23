! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Landau distribution matching scipy.stats.landau.
!
! The standard density is evaluated by three complementary representations:
! a left-tail asymptotic expansion, exponentially damped characteristic-function
! quadrature in the central region, and the defining Landau integral after a
! scale change on the heavy right tail.  The left-tail coefficients are the
! classical Kolbig-Schorr asymptotic coefficients as tabulated by Yoshimura
! (2024); no source code from SciPy, Boost, or Yoshimura's implementation is
! copied or adapted.
module scifort_landau
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi, scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: landau_cdf, landau_isf, landau_logcdf, landau_logpdf
    public :: landau_logsf, landau_pdf, landau_ppf, landau_sf
    public :: landau_logpdf_derivatives

    real(dp), parameter :: half_pi = 0.5_dp * scifort_pi
    real(dp), parameter :: inv_pi = 1.0_dp / scifort_pi
    real(dp), parameter :: alpha = 2.0_dp / scifort_pi
    real(dp), parameter :: log_half_pi = 0.451582705289454864726195229894882143572_dp
    real(dp), parameter :: log_two_over_pi = -0.451582705289454864726195229894882143572_dp
    real(dp), parameter :: left_switch = -2.7_dp
    real(dp), parameter :: right_switch = 1.0_dp
    real(dp), parameter :: log_quad_lower = -40.0_dp
    real(dp), parameter :: log_quad_upper = 3.688879454113936302852455697600717343752_dp ! log(40)
    real(dp), parameter :: central_segment_width = 0.5_dp
    real(dp), parameter :: right_segment_width = 1.0_dp
    real(dp), parameter :: extreme_right = 1.0e100_dp

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

    ! PDF expansion coefficients a_n^- in sigma^{-n}.
    real(dp), parameter :: left_pdf_coef(0:8) = [ &
        1.0_dp, &
        1.0_dp/24.0_dp, &
        -23.0_dp/1152.0_dp, &
        11237.0_dp/414720.0_dp, &
        -2482411.0_dp/39813120.0_dp, &
        272785979.0_dp/1337720832.0_dp, &
        -4175309343349.0_dp/4815794995200.0_dp, &
        525035501918789.0_dp/115579079884800.0_dp, &
        -628141988536245979.0_dp/22191183337881600.0_dp ]

    ! CDF expansion coefficients A_n^- in sigma^{-n}.
    real(dp), parameter :: left_cdf_coef(0:8) = [ &
        1.0_dp, &
        -11.0_dp/24.0_dp, &
        769.0_dp/1152.0_dp, &
        -680863.0_dp/414720.0_dp, &
        226287557.0_dp/39813120.0_dp, &
        -169709463197.0_dp/6688604160.0_dp, &
        667874164916771.0_dp/4815794995200.0_dp, &
        -103663334225097487.0_dp/115579079884800.0_dp, &
        21235294185086305043.0_dp/3170169048268800.0_dp ]

contains

    pure elemental function landau_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, dz

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            call standard_logpdf_derivative(z, y, dz)
            y = y - log(sigma)
        end if
    end function landau_logpdf

    pure elemental function landau_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = landau_logpdf(x, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function landau_pdf

    pure elemental function landau_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, cdf_value, sf_value

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
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
            if (z <= left_switch) then
                y = left_logcdf(z)
            else if (z < right_switch) then
                call central_values(z, cdf_value=cdf_value)
                if (cdf_value <= 0.0_dp) then
                    y = negative_infinity(x)
                else
                    y = log(cdf_value)
                end if
            else
                call right_values(z, sf_value=sf_value)
                if (sf_value >= 1.0_dp) then
                    y = negative_infinity(x)
                else
                    y = log1p_safe(-sf_value)
                end if
            end if
        end if
    end function landau_logcdf

    pure elemental function landau_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, lcdf, cdf_value, sf_value

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
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
            if (z <= left_switch) then
                lcdf = left_logcdf(z)
                if (lcdf < log(tiny(1.0_dp))) then
                    y = -exp(lcdf)
                else
                    y = log1p_safe(-exp(lcdf))
                end if
            else if (z < right_switch) then
                call central_values(z, cdf_value=cdf_value)
                if (cdf_value <= 0.0_dp) then
                    y = 0.0_dp
                else if (cdf_value >= 1.0_dp) then
                    y = negative_infinity(x)
                else
                    y = log1p_safe(-cdf_value)
                end if
            else
                if (z >= extreme_right) then
                    y = log_two_over_pi - log(z)
                else
                    call right_values(z, sf_value=sf_value)
                    if (sf_value <= 0.0_dp) then
                        y = negative_infinity(x)
                    else
                        y = log(sf_value)
                    end if
                end if
            end if
        end if
    end function landau_logsf

    pure elemental function landau_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lcdf, lsf

        lcdf = landau_logcdf(x, loc, scale)
        if (ieee_is_nan(lcdf)) then
            y = lcdf
        else if (lcdf < -log(2.0_dp)) then
            y = exp(lcdf)
        else
            lsf = landau_logsf(x, loc, scale)
            y = -expm1_safe(lsf)
        end if
    end function landau_cdf

    pure elemental function landau_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lcdf, lsf

        lsf = landau_logsf(x, loc, scale)
        if (ieee_is_nan(lsf)) then
            y = lsf
        else if (lsf < -log(2.0_dp)) then
            y = exp(lsf)
        else
            lcdf = landau_logcdf(x, loc, scale)
            y = -expm1_safe(lcdf)
        end if
    end function landau_sf

    pure elemental function landau_ppf(probability, loc, scale) result(x)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma) .or. .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            x = negative_infinity(probability)
        else if (probability >= 1.0_dp) then
            x = positive_infinity(probability)
        else if (probability <= 0.5_dp) then
            z = standard_quantile(probability, .false.)
            x = mu + sigma * z
        else
            z = standard_quantile(1.0_dp - probability, .true.)
            x = mu + sigma * z
        end if
    end function landau_ppf

    pure elemental function landau_isf(probability, loc, scale) result(x)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma) .or. .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            x = negative_infinity(probability)
        else if (probability <= 0.0_dp) then
            x = positive_infinity(probability)
        else if (probability <= 0.5_dp) then
            z = standard_quantile(probability, .true.)
            x = mu + sigma * z
        else
            z = standard_quantile(1.0_dp - probability, .false.)
            x = mu + sigma * z
        end if
    end function landau_isf

    pure subroutine landau_logpdf_derivatives(z, logf, dz)
        real(dp), intent(in) :: z !! standardized observation
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: dz !! derivative of standardized log density with respect to z

        if (.not. ieee_is_finite(z)) then
            logf = negative_infinity(z)
            dz = quiet_nan(z)
        else
            call standard_logpdf_derivative(z, logf, dz)
        end if
    end subroutine landau_logpdf_derivatives

    pure subroutine standard_logpdf_derivative(z, logf, derivative)
        real(dp), intent(in) :: z !! standardized observation
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: derivative !! derivative of log density with respect to z
        real(dp) :: pdf_value, slope

        if (z <= left_switch) then
            call left_logpdf_derivative(z, logf, derivative)
        else if (z >= extreme_right) then
            logf = log_two_over_pi - 2.0_dp * log(z)
            derivative = -2.0_dp / z
        else if (z < right_switch) then
            call central_values(z, pdf_value=pdf_value, slope=slope)
            pdf_value = max(pdf_value, tiny(1.0_dp))
            logf = log(pdf_value)
            derivative = slope / pdf_value
        else
            call right_values(z, pdf_value=pdf_value, slope=slope)
            pdf_value = max(pdf_value, tiny(1.0_dp))
            logf = log(pdf_value)
            derivative = slope / pdf_value
        end if
    end subroutine standard_logpdf_derivative

    pure subroutine left_logpdf_derivative(z, logf, derivative)
        real(dp), intent(in) :: z !! standardized left-tail observation
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: derivative !! derivative of log density with respect to z
        real(dp) :: logsigma, sigma, invsigma, poly, dpoly, dlogp_dlogsigma

        logsigma = -half_pi * z - log_half_pi - 1.0_dp
        if (logsigma >= log(huge(1.0_dp))) then
            logf = negative_infinity(z)
            derivative = positive_infinity(z)
            return
        end if
        sigma = exp(logsigma)
        invsigma = exp(-logsigma)
        call polynomial_and_derivative(invsigma, left_pdf_coef, poly, dpoly)
        if (poly <= 0.0_dp) then
            logf = negative_infinity(z)
            derivative = positive_infinity(z)
            return
        end if
        logf = log_half_pi - scifort_log_sqrt_two_pi + 0.5_dp * logsigma - sigma + log(poly)
        dlogp_dlogsigma = 0.5_dp - sigma - invsigma * dpoly / poly
        derivative = -half_pi * dlogp_dlogsigma
    end subroutine left_logpdf_derivative

    pure elemental function left_logcdf(z) result(value)
        real(dp), intent(in) :: z !! standardized left-tail observation
        real(dp) :: value, logsigma, sigma, invsigma, poly, dummy

        logsigma = -half_pi * z - log_half_pi - 1.0_dp
        if (logsigma >= log(huge(1.0_dp))) then
            value = negative_infinity(z)
            return
        end if
        sigma = exp(logsigma)
        invsigma = exp(-logsigma)
        call polynomial_and_derivative(invsigma, left_cdf_coef, poly, dummy)
        if (poly <= 0.0_dp) then
            value = negative_infinity(z)
        else
            value = -scifort_log_sqrt_two_pi - 0.5_dp * logsigma - sigma + log(poly)
        end if
    end function left_logcdf

    pure subroutine polynomial_and_derivative(x, coefficients, value, derivative)
        real(dp), intent(in) :: x !! polynomial argument
        real(dp), intent(in) :: coefficients(0:) !! coefficients in ascending powers
        real(dp), intent(out) :: value !! polynomial value
        real(dp), intent(out) :: derivative !! derivative with respect to x
        integer :: k

        value = coefficients(ubound(coefficients,1))
        derivative = 0.0_dp
        do k = ubound(coefficients,1)-1, 0, -1
            derivative = derivative * x + value
            value = value * x + coefficients(k)
        end do
    end subroutine polynomial_and_derivative

    pure subroutine central_values(z, pdf_value, slope, cdf_value)
        real(dp), intent(in) :: z !! standardized observation in the central integration region
        real(dp), intent(out), optional :: pdf_value !! standardized density
        real(dp), intent(out), optional :: slope !! derivative of standardized density
        real(dp), intent(out), optional :: cdf_value !! lower-tail probability
        real(dp) :: pdf_sum, slope_sum, cdf_sum, l, r, mid, half, s1, s2
        real(dp) :: t1, t2, phase1, phase2, e1, e2, weight
        integer :: nseg, j, i

        nseg = ceiling((log_quad_upper-log_quad_lower) / central_segment_width)
        pdf_sum = 0.0_dp
        slope_sum = 0.0_dp
        cdf_sum = 0.0_dp
        do j = 1, nseg
            l = log_quad_lower + (log_quad_upper-log_quad_lower) * real(j-1,dp) / real(nseg,dp)
            r = log_quad_lower + (log_quad_upper-log_quad_lower) * real(j,dp) / real(nseg,dp)
            mid = 0.5_dp * (l+r)
            half = 0.5_dp * (r-l)
            do i = 1, 8
                weight = half * gl_w(i)
                s1 = mid - half * gl_x(i)
                s2 = mid + half * gl_x(i)
                t1 = exp(s1)
                t2 = exp(s2)
                phase1 = t1 * (z + alpha * s1)
                phase2 = t2 * (z + alpha * s2)
                e1 = exp(-t1)
                e2 = exp(-t2)
                pdf_sum = pdf_sum + weight * (t1*e1*cos(phase1) + t2*e2*cos(phase2))
                slope_sum = slope_sum - weight * (t1*t1*e1*sin(phase1) + t2*t2*e2*sin(phase2))
                cdf_sum = cdf_sum + weight * (e1*sin(phase1) + e2*sin(phase2))
            end do
        end do
        if (present(pdf_value)) pdf_value = max(0.0_dp, inv_pi * pdf_sum)
        if (present(slope)) slope = inv_pi * slope_sum
        if (present(cdf_value)) cdf_value = min(1.0_dp, max(0.0_dp, 0.5_dp + inv_pi * cdf_sum))
    end subroutine central_values

    pure subroutine right_values(z, pdf_value, slope, sf_value)
        real(dp), intent(in) :: z !! standardized observation in the heavy right-tail region
        real(dp), intent(out), optional :: pdf_value !! standardized density
        real(dp), intent(out), optional :: slope !! derivative of standardized density
        real(dp), intent(out), optional :: sf_value !! upper-tail probability
        real(dp) :: lambda, pdf_sum, slope_sum, sf_sum, l, r, mid, half, s1, s2
        real(dp) :: u1, u2, t1, t2, q1, q2, base1, base2, weight
        integer :: nseg, j, i

        if (z >= extreme_right) then
            if (present(pdf_value)) pdf_value = exp(log_two_over_pi - 2.0_dp*log(z))
            if (present(slope)) slope = -2.0_dp * exp(log_two_over_pi - 3.0_dp*log(z))
            if (present(sf_value)) sf_value = exp(log_two_over_pi - log(z))
            return
        end if
        lambda = log_half_pi + half_pi * z
        nseg = ceiling((log_quad_upper-log_quad_lower) / right_segment_width)
        pdf_sum = 0.0_dp
        slope_sum = 0.0_dp
        sf_sum = 0.0_dp
        do j = 1, nseg
            l = log_quad_lower + (log_quad_upper-log_quad_lower) * real(j-1,dp) / real(nseg,dp)
            r = log_quad_lower + (log_quad_upper-log_quad_lower) * real(j,dp) / real(nseg,dp)
            mid = 0.5_dp * (l+r)
            half = 0.5_dp * (r-l)
            do i = 1, 8
                weight = half * gl_w(i)
                s1 = mid - half * gl_x(i)
                s2 = mid + half * gl_x(i)
                u1 = exp(s1)
                u2 = exp(s2)
                t1 = u1 / lambda
                t2 = u2 / lambda
                q1 = -u1 - t1 * log(t1)
                q2 = -u2 - t2 * log(t2)
                base1 = exp(q1) * sin(scifort_pi*t1)
                base2 = exp(q2) * sin(scifort_pi*t2)
                pdf_sum = pdf_sum + weight * (u1*base1 + u2*base2)
                slope_sum = slope_sum + weight * (u1*u1*base1 + u2*u2*base2)
                sf_sum = sf_sum + weight * (base1 + base2)
            end do
        end do
        if (present(pdf_value)) pdf_value = max(0.0_dp, pdf_sum / (2.0_dp*lambda))
        if (present(slope)) slope = -scifort_pi * slope_sum / (4.0_dp*lambda*lambda)
        if (present(sf_value)) sf_value = min(1.0_dp, max(0.0_dp, inv_pi*sf_sum))
    end subroutine right_values

    pure function standard_quantile(probability, upper) result(z)
        real(dp), intent(in) :: probability !! requested tail probability in (0,0.5]
        logical, intent(in) :: upper !! true for upper-tail inversion
        real(dp) :: z, lo, hi, mid, logtarget, value
        integer :: iter

        logtarget = log(probability)
        lo = -4.0_dp
        hi = 4.0_dp
        if (upper) then
            do iter = 1, 1100
                value = landau_logsf(hi)
                if (value <= logtarget) exit
                if (hi >= 0.25_dp*huge(1.0_dp)) then
                    z = positive_infinity(probability)
                    return
                end if
                hi = 2.0_dp * hi + 1.0_dp
            end do
            do iter = 1, 80
                mid = 0.5_dp * (lo+hi)
                if (landau_logsf(mid) > logtarget) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        else
            do iter = 1, 1100
                value = landau_logcdf(lo)
                if (value <= logtarget) exit
                if (lo <= -0.25_dp*huge(1.0_dp)) then
                    z = negative_infinity(probability)
                    return
                end if
                lo = 2.0_dp * lo - 1.0_dp
            end do
            do iter = 1, 80
                mid = 0.5_dp * (lo+hi)
                if (landau_logcdf(mid) < logtarget) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        end if
        z = 0.5_dp * (lo+hi)
    end function standard_quantile

    pure elemental function valid_probability(p) result(ok)
        real(dp), intent(in) :: p !! candidate probability
        logical :: ok
        ok = .not. ieee_is_nan(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

    pure elemental subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location parameter
        real(dp), intent(in), optional :: scale !! optional scale parameter
        real(dp), intent(out) :: mu !! resolved location parameter
        real(dp), intent(out) :: sigma !! resolved scale parameter
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_landau
