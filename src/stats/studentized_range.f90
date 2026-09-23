! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Studentized-range distribution matching scipy.stats.studentized_range.
! The finite-df implementation uses the defining scale-mixture representation:
! conditionally on S=sqrt(chi-square(df)/df), q*S is the range of k standard
! normals.  The normal-range integrals are the formulas in Batista et al.
! (2017), while the df -> infinity limit is the Lund & Lund (1983) formula.
! SciPy 1.17.0 is used only as an API and numerical-parity reference.
module scifort_studentized_range
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_pdf, normal_sf
    implicit none
    private

    public :: studentized_range_cdf, studentized_range_isf
    public :: studentized_range_logcdf, studentized_range_logpdf
    public :: studentized_range_logsf, studentized_range_pdf
    public :: studentized_range_ppf, studentized_range_sf

    real(dp), parameter :: asymptotic_df = 100000.0_dp
    real(dp), parameter :: log_two = 0.693147180559945309417232121458176568_dp
    real(dp), parameter :: sqrt_pi = 1.772453850905516027298167483341145183_dp
    real(dp), parameter :: z_cutoff = 11.0_dp
    real(dp), parameter :: z_panel_width = 0.75_dp
    real(dp), parameter :: outer_panel_width = 0.5_dp

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

    pure elemental function studentized_range_pdf(x, k, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, q

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(k, df) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            q = (x - mu) / sigma
            if (q < 0.0_dp .or. q == positive_infinity(q)) then
                y = 0.0_dp
            else if (q == 0.0_dp) then
                y = pdf_at_zero(k, df) / sigma
            else if (.not. ieee_is_finite(q)) then
                y = 0.0_dp
            else if (df >= asymptotic_df) then
                y = range_pdf(q, k) / sigma
            else
                y = finite_df_value(q, k, df, .true.) / sigma
            end if
        end if
    end function studentized_range_pdf

    pure elemental function studentized_range_logpdf(x, k, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p

        p = studentized_range_pdf(x, k, df, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else if (.not. ieee_is_finite(p)) then
            y = positive_infinity(p)
        else
            y = log(p)
        end if
    end function studentized_range_logpdf

    pure elemental function studentized_range_cdf(x, k, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, q

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(k, df) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            q = (x - mu) / sigma
            if (q <= 0.0_dp) then
                y = 0.0_dp
            else if (.not. ieee_is_finite(q)) then
                y = 1.0_dp
            else if (df >= asymptotic_df) then
                y = range_cdf(q, k)
            else
                y = finite_df_value(q, k, df, .false.)
            end if
            y = min(1.0_dp, max(0.0_dp, y))
        end if
    end function studentized_range_cdf

    pure elemental function studentized_range_sf(x, k, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, cdf_value

        cdf_value = studentized_range_cdf(x, k, df, loc, scale)
        if (ieee_is_nan(cdf_value)) then
            y = cdf_value
        else
            y = max(0.0_dp, 1.0_dp - cdf_value)
        end if
    end function studentized_range_sf

    pure elemental function studentized_range_logcdf(x, k, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p

        p = studentized_range_cdf(x, k, df, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function studentized_range_logcdf

    pure elemental function studentized_range_logsf(x, k, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p

        p = studentized_range_cdf(x, k, df, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p >= 1.0_dp) then
            y = negative_infinity(p)
        else
            y = log1p_safe(-p)
        end if
    end function studentized_range_logsf

    pure elemental function studentized_range_ppf(probability, k, df, loc, scale) result(x)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, q

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(k, df) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability == 0.0_dp) then
            x = mu
        else if (probability == 1.0_dp) then
            x = positive_infinity(probability)
        else
            q = standard_quantile(probability, k, df, .false.)
            x = mu + sigma * q
        end if
    end function studentized_range_ppf

    pure elemental function studentized_range_isf(probability, k, df, loc, scale) result(x)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, q

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(k, df) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability == 1.0_dp) then
            x = mu
        else if (probability == 0.0_dp) then
            x = positive_infinity(probability)
        else
            q = standard_quantile(probability, k, df, .true.)
            x = mu + sigma * q
        end if
    end function studentized_range_isf

    pure function standard_quantile(probability, k, df, upper) result(q)
        real(dp), intent(in) :: probability !! probability strictly between zero and one
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        logical, intent(in) :: upper !! true to invert the survival function
        real(dp) :: q, lo, hi, mid, value, target
        integer :: iter

        target = probability
        lo = 0.0_dp
        hi = 1.0_dp
        do iter = 1, 256
            if (upper) then
                value = standard_sf(hi, k, df)
                if (value <= target) exit
            else
                value = standard_cdf(hi, k, df)
                if (value >= target) exit
            end if
            hi = 2.0_dp * hi
            if (.not. ieee_is_finite(hi)) then
                q = positive_infinity(probability)
                return
            end if
        end do
        do iter = 1, 64
            mid = 0.5_dp * (lo + hi)
            if (upper) then
                if (standard_sf(mid, k, df) > target) then
                    lo = mid
                else
                    hi = mid
                end if
            else
                if (standard_cdf(mid, k, df) < target) then
                    lo = mid
                else
                    hi = mid
                end if
            end if
        end do
        q = 0.5_dp * (lo + hi)
    end function standard_quantile

    pure function standard_cdf(q, k, df) result(p)
        real(dp), intent(in) :: q !! standardized nonnegative variate
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp) :: p
        if (q <= 0.0_dp) then
            p = 0.0_dp
        else if (df >= asymptotic_df) then
            p = range_cdf(q, k)
        else
            p = finite_df_value(q, k, df, .false.)
        end if
        p = min(1.0_dp, max(0.0_dp, p))
    end function standard_cdf

    pure function standard_sf(q, k, df) result(p)
        real(dp), intent(in) :: q !! standardized nonnegative variate
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp) :: p
        p = max(0.0_dp, 1.0_dp - standard_cdf(q, k, df))
    end function standard_sf

    pure function finite_df_value(q, k, df, want_pdf) result(value)
        real(dp), intent(in) :: q !! standardized positive variate
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! finite degrees of freedom in (0,100000)
        logical, intent(in) :: want_pdf !! true for density, false for CDF
        real(dp) :: value
        real(dp) :: coeff, left_y, right_y, tlo, thi
        real(dp) :: left, right, mid, half, t1, t2
        integer :: nseg, j, i

        coeff = df + k - 1.0_dp
        left_y = -max(8.0_dp, 48.0_dp / coeff)
        right_y = max(4.0_dp, 0.5_dp * log(max(1.0_dp, 100.0_dp / df)))
        tlo = asinh(left_y)
        thi = asinh(right_y)
        nseg = max(1, ceiling((thi - tlo) / outer_panel_width))
        value = 0.0_dp
        do j = 1, nseg
            left = tlo + (thi - tlo) * real(j - 1, dp) / real(nseg, dp)
            right = tlo + (thi - tlo) * real(j, dp) / real(nseg, dp)
            mid = 0.5_dp * (left + right)
            half = 0.5_dp * (right - left)
            do i = 1, 8
                t1 = mid - half * gl_x(i)
                t2 = mid + half * gl_x(i)
                value = value + half * gl_w(i) * outer_integrand(t1, q, k, df, want_pdf)
                value = value + half * gl_w(i) * outer_integrand(t2, q, k, df, want_pdf)
            end do
        end do
    end function finite_df_value

    pure function outer_integrand(t, q, k, df, want_pdf) result(value)
        real(dp), intent(in) :: t !! double-exponential log-scale coordinate
        real(dp), intent(in) :: q !! standardized positive variate
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! finite degrees of freedom in (0,100000)
        logical, intent(in) :: want_pdf !! true for density, false for CDF
        real(dp) :: value, y, s, r, log_weight, conditional

        y = sinh(t)
        if (y > 350.0_dp) then
            value = 0.0_dp
            return
        end if
        s = exp(y)
        if (s == 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (.not. ieee_is_finite(s)) then
            value = 0.0_dp
            return
        end if
        r = q * s
        log_weight = scaled_chi_log_constant(df) + df * y - &
            0.5_dp * df * s * s + log_cosh_safe(t)
        if (want_pdf) then
            if (log_weight + y < log(tiny(1.0_dp))) then
                value = 0.0_dp
                return
            end if
            conditional = range_pdf(r, k)
            value = exp(log_weight + y) * conditional
        else
            if (log_weight < log(tiny(1.0_dp))) then
                value = 0.0_dp
                return
            end if
            conditional = range_cdf(r, k)
            value = exp(log_weight) * conditional
        end if
    end function outer_integrand

    pure elemental function scaled_chi_log_constant(df) result(value)
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp) :: value, a
        a = 0.5_dp * df
        value = log_two + a * log(a) - log_gamma(a)
    end function scaled_chi_log_constant

    pure function range_cdf(r, k) result(value)
        real(dp), intent(in) :: r !! nonnegative range threshold
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp) :: value
        real(dp) :: left, right, mid, half, z1, z2
        integer :: nseg, j, i

        if (r <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        nseg = max(1, ceiling((2.0_dp * z_cutoff) / z_panel_width))
        value = 0.0_dp
        do j = 1, nseg
            left = -z_cutoff + 2.0_dp * z_cutoff * real(j - 1, dp) / real(nseg, dp)
            right = -z_cutoff + 2.0_dp * z_cutoff * real(j, dp) / real(nseg, dp)
            mid = 0.5_dp * (left + right)
            half = 0.5_dp * (right - left)
            do i = 1, 8
                z1 = mid - half * gl_x(i)
                z2 = mid + half * gl_x(i)
                value = value + half * gl_w(i) * range_cdf_integrand(z1, r, k)
                value = value + half * gl_w(i) * range_cdf_integrand(z2, r, k)
            end do
        end do
        value = min(1.0_dp, max(0.0_dp, value))
    end function range_cdf

    pure function range_pdf(r, k) result(value)
        real(dp), intent(in) :: r !! nonnegative range threshold
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp) :: value
        real(dp) :: center, left_bound, right_bound, left, right, mid, half, z1, z2
        integer :: nseg, j, i

        if (r < 0.0_dp) then
            value = 0.0_dp
            return
        else if (r == 0.0_dp) then
            if (k < 2.0_dp) then
                value = positive_infinity(r)
            else if (k == 2.0_dp) then
                value = 1.0_dp / sqrt_pi
            else
                value = 0.0_dp
            end if
            return
        end if
        center = -0.5_dp * r
        left_bound = center - z_cutoff
        right_bound = center + z_cutoff
        nseg = max(1, ceiling((right_bound - left_bound) / z_panel_width))
        value = 0.0_dp
        do j = 1, nseg
            left = left_bound + (right_bound - left_bound) * real(j - 1, dp) / real(nseg, dp)
            right = left_bound + (right_bound - left_bound) * real(j, dp) / real(nseg, dp)
            mid = 0.5_dp * (left + right)
            half = 0.5_dp * (right - left)
            do i = 1, 8
                z1 = mid - half * gl_x(i)
                z2 = mid + half * gl_x(i)
                value = value + half * gl_w(i) * range_pdf_integrand(z1, r, k)
                value = value + half * gl_w(i) * range_pdf_integrand(z2, r, k)
            end do
        end do
        value = max(0.0_dp, value)
    end function range_pdf

    pure elemental function range_cdf_integrand(z, r, k) result(value)
        real(dp), intent(in) :: z !! standard-normal integration coordinate
        real(dp), intent(in) :: r !! positive range threshold
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp) :: value, d
        d = normal_interval_probability(z, z + r)
        if (d <= 0.0_dp) then
            value = 0.0_dp
        else
            value = k * normal_pdf(z) * exp((k - 1.0_dp) * log(d))
        end if
    end function range_cdf_integrand

    pure elemental function range_pdf_integrand(z, r, k) result(value)
        real(dp), intent(in) :: z !! standard-normal integration coordinate
        real(dp), intent(in) :: r !! positive range threshold
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp) :: value, d, exponent
        d = normal_interval_probability(z, z + r)
        if (d <= 0.0_dp) then
            value = 0.0_dp
        else
            exponent = (k - 2.0_dp) * log(d)
            value = k * (k - 1.0_dp) * normal_pdf(z) * normal_pdf(z + r) * exp(exponent)
        end if
    end function range_pdf_integrand

    pure elemental function normal_interval_probability(a, b) result(p)
        real(dp), intent(in) :: a !! lower standard-normal endpoint
        real(dp), intent(in) :: b !! upper standard-normal endpoint, >= a
        real(dp) :: p
        if (a >= 0.0_dp) then
            p = normal_sf(a) - normal_sf(b)
        else
            p = normal_cdf(b) - normal_cdf(a)
        end if
        p = max(0.0_dp, p)
    end function normal_interval_probability

    pure elemental function pdf_at_zero(k, df) result(value)
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp) :: value
        if (k < 2.0_dp) then
            value = positive_infinity(k)
        else if (k > 2.0_dp) then
            value = 0.0_dp
        else if (df >= asymptotic_df) then
            value = 1.0_dp / sqrt_pi
        else
            value = exp(0.5_dp * (log(2.0_dp) - log(df)) + &
                log_gamma(0.5_dp * (df + 1.0_dp)) - log_gamma(0.5_dp * df)) / sqrt_pi
        end if
    end function pdf_at_zero

    pure elemental function log_cosh_safe(x) result(y)
        real(dp), intent(in) :: x !! real argument
        real(dp) :: y, ax
        ax = abs(x)
        y = ax + log1p_safe(exp(-2.0_dp * ax)) - log_two
    end function log_cosh_safe

    pure elemental logical function valid_shapes(k, df) result(ok)
        real(dp), intent(in) :: k !! candidate number-of-means shape
        real(dp), intent(in) :: df !! candidate degrees of freedom
        ok = ieee_is_finite(k) .and. ieee_is_finite(df) .and. k > 1.0_dp .and. df > 0.0_dp
    end function valid_shapes

    pure elemental logical function valid_probability(p) result(ok)
        real(dp), intent(in) :: p !! probability to validate
        ok = .not. ieee_is_nan(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location parameter
        real(dp), intent(in), optional :: scale !! optional positive scale parameter
        real(dp), intent(out) :: mu !! resolved location parameter
        real(dp), intent(out) :: sigma !! resolved scale parameter
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_studentized_range
