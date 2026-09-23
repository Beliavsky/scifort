! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Independent implementation of the Birnbaum-Tingey/Smirnov finite-sample law.
! See CODE_PROVENANCE.md.

module scifort_ksone
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: ksone_cdf, ksone_isf, ksone_logcdf, ksone_logpdf
    public :: ksone_logsf, ksone_pdf, ksone_ppf, ksone_sf
    public :: ksone_logpdf_derivative

contains

    pure elemental function ksone_logsf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z, s0, s1, s2, ls
        integer :: ni

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else if (z >= 1.0_dp) then
            y = negative_infinity(z)
        else
            ni = int(n)
            call standard_sums(z, ni, s0, s1, s2, ls)
            if (s0 <= 0.0_dp) then
                y = negative_infinity(z)
            else
                y = min(0.0_dp, ls + log(s0))
            end if
        end if
    end function ksone_logsf

    pure elemental function ksone_sf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y, ly
        ly = ksone_logsf(x, n, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function ksone_sf

    pure elemental function ksone_logcdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y, lq
        real(dp) :: mu, sigma, z
        integer :: ni

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = negative_infinity(z)
        else if (z >= 1.0_dp) then
            y = 0.0_dp
        else if (z <= 1.0_dp / n) then
            ni = int(n)
            y = log(z) + real(ni - 1, dp) * log1p_safe(z)
        else
            lq = ksone_logsf(x, n, loc, scale)
            if (lq < -0.69314718055994530942_dp) then
                y = log1p_safe(-exp(lq))
            else
                y = log(-expm1_safe(lq))
            end if
        end if
    end function ksone_logcdf

    pure elemental function ksone_cdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y, ly
        ly = ksone_logcdf(x, n, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function ksone_cdf

    pure elemental function ksone_logpdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z, s0, s1, s2, ls
        integer :: ni

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (z < 0.0_dp .or. z > 1.0_dp) then
            y = negative_infinity(z)
        else if (int(n) == 1) then
            y = -log(sigma)
        else if (z <= 0.0_dp .or. z >= 1.0_dp) then
            y = negative_infinity(z)
        else if (at_knot(z, int(n))) then
            y = log(standard_pdf_knot(z, int(n))) - log(sigma)
        else
            ni = int(n)
            call standard_sums(z, ni, s0, s1, s2, ls)
            if (s1 >= 0.0_dp) then
                y = negative_infinity(z)
            else
                y = ls + log(-s1) - log(sigma)
            end if
        end if
    end function ksone_logpdf

    pure elemental function ksone_pdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y, ly
        ly = ksone_logpdf(x, n, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function ksone_pdf

    pure elemental function ksone_ppf(p, n, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = mu
        else if (p >= 1.0_dp) then
            x = mu + sigma
        else if (p > 0.5_dp) then
            x = ksone_isf(1.0_dp - p, n, mu, sigma)
        else
            lo = 0.0_dp
            hi = 1.0_dp
            do iter = 1, 100
                mid = 0.5_dp * (lo + hi)
                if (ksone_cdf(mu + sigma * mid, n, mu, sigma) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function ksone_ppf

    pure elemental function ksone_isf(p, n, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p >= 1.0_dp) then
            x = mu
        else if (p <= 0.0_dp) then
            x = mu + sigma
        else
            lo = 0.0_dp
            hi = 1.0_dp
            do iter = 1, 100
                mid = 0.5_dp * (lo + hi)
                if (ksone_sf(mu + sigma * mid, n, mu, sigma) > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function ksone_isf

    pure elemental function ksone_logpdf_derivative(x, n) result(y)
        real(dp), intent(in) :: x !! standardized interior point
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp) :: y
        real(dp) :: s0, s1, s2, ls
        integer :: ni

        if (.not. valid_shape(n) .or. x <= 0.0_dp .or. x >= 1.0_dp) then
            y = quiet_nan(x)
            return
        end if
        ni = int(n)
        if (ni == 1) then
            y = 0.0_dp
        else if (at_knot(x, ni)) then
            y = quiet_nan(x)
        else
            call standard_sums(x, ni, s0, s1, s2, ls)
            if (s1 >= 0.0_dp) then
                y = quiet_nan(x)
            else
                y = s2 / s1
            end if
        end if
    end function ksone_logpdf_derivative

    pure subroutine standard_sums(x, n, s0, s1, s2, log_scale)
        real(dp), intent(in) :: x !! standardized interior point, not a knot
        integer, intent(in) :: n !! positive integer sample size
        real(dp), intent(out) :: s0 !! scaled survival sum
        real(dp), intent(out) :: s1 !! scaled first derivative of survival
        real(dp), intent(out) :: s2 !! scaled second derivative of survival
        real(dp), intent(out) :: log_scale !! common logarithmic scale
        integer :: j, m
        real(dp) :: a, b, logt, maxlog, g, gp, w
        real(dp) :: logs(0:n), gs(0:n), gps(0:n)

        m = min(n, int(floor(real(n, dp) * (1.0_dp - x))))
        maxlog = -huge(1.0_dp)
        logs = -huge(1.0_dp)
        gs = 0.0_dp
        gps = 0.0_dp
        do j = 0, m
            a = x + real(j, dp) / real(n, dp)
            b = 1.0_dp - x - real(j, dp) / real(n, dp)
            if (b <= 0.0_dp) cycle
            if (j == 0) then
                logt = real(n, dp) * log(b)
                g = -real(n, dp) / b
                gp = -real(n, dp) / (b * b)
            else
                logt = log(x) + log_binomial(n, j) + real(j - 1, dp) * log(a) + &
                    real(n - j, dp) * log(b)
                g = 1.0_dp / x + real(j - 1, dp) / a - real(n - j, dp) / b
                gp = -1.0_dp / (x * x) - real(j - 1, dp) / (a * a) - &
                    real(n - j, dp) / (b * b)
            end if
            logs(j) = logt
            gs(j) = g
            gps(j) = gp
            maxlog = max(maxlog, logt)
        end do

        s0 = 0.0_dp
        s1 = 0.0_dp
        s2 = 0.0_dp
        if (maxlog <= -huge(1.0_dp) / 2.0_dp) then
            log_scale = maxlog
            return
        end if
        do j = 0, m
            if (logs(j) <= -huge(1.0_dp) / 2.0_dp) cycle
            w = exp(logs(j) - maxlog)
            s0 = s0 + w
            s1 = s1 + w * gs(j)
            s2 = s2 + w * (gs(j) * gs(j) + gps(j))
        end do
        log_scale = maxlog
    end subroutine standard_sums

    pure function standard_pdf_knot(x, n) result(f)
        real(dp), intent(in) :: x !! standardized knot
        integer, intent(in) :: n !! positive sample size
        real(dp) :: f
        real(dp) :: epsx
        epsx = max(64.0_dp * epsilon(1.0_dp), 8.0_dp * epsilon(1.0_dp) * abs(x))
        f = 0.5_dp * (standard_pdf_regular(max(0.0_dp, x - epsx), n) + &
            standard_pdf_regular(min(1.0_dp, x + epsx), n))
    end function standard_pdf_knot

    pure function standard_pdf_regular(x, n) result(f)
        real(dp), intent(in) :: x !! standardized non-knot point
        integer, intent(in) :: n !! positive sample size
        real(dp) :: f, s0, s1, s2, ls
        if (n == 1) then
            f = 1.0_dp
        else if (x <= 0.0_dp .or. x >= 1.0_dp) then
            f = 0.0_dp
        else
            call standard_sums(x, n, s0, s1, s2, ls)
            if (s1 < 0.0_dp) then
                f = exp(ls) * (-s1)
            else
                f = 0.0_dp
            end if
        end if
    end function standard_pdf_regular

    pure elemental function at_knot(x, n) result(ok)
        real(dp), intent(in) :: x !! standardized point
        integer, intent(in) :: n !! positive sample size
        logical :: ok
        real(dp) :: t
        t = real(n, dp) * (1.0_dp - x)
        ok = abs(t - anint(t)) <= 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(t))
    end function at_knot

    pure elemental function valid_shape(n) result(ok)
        real(dp), intent(in) :: n !! candidate sample size
        logical :: ok
        ok = ieee_is_finite(n) .and. n >= 1.0_dp .and. n == floor(n)
    end function valid_shape

    pure elemental function log_binomial(n, k) result(y)
        integer, intent(in) :: n !! nonnegative total count
        integer, intent(in) :: k !! count in [0,n]
        real(dp) :: y
        y = log_gamma(real(n + 1, dp)) - log_gamma(real(k + 1, dp)) - &
            log_gamma(real(n - k + 1, dp))
    end function log_binomial

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! requested location, absent for zero
        real(dp), intent(in), optional :: scale !! requested scale, absent for one
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_ksone
