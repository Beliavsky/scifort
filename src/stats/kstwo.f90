! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2001-2002 Enthought, Inc. 2003-2026, SciPy Developers.
! Copyright (c) 2026 SciFort contributors
!
! The regime selection, Durbin/Marsaglia-Tsang-Wang matrix algorithm, and
! Pelz-Good expansion in this file are adapted from scipy.stats._ksstats in
! SciPy 1.17.0. See THIRD_PARTY_LICENSES.md and CODE_PROVENANCE.md.

module scifort_kstwo
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_ksone, only : ksone_pdf, ksone_sf
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    real(dp), parameter :: log_two_pi = 1.83787706640934548356_dp
    real(dp), parameter :: sqrt_two_pi = 2.50662827463100050242_dp
    real(dp), parameter :: sqrt_three = 1.73205080756887729353_dp
    real(dp), parameter :: pi2 = scifort_pi * scifort_pi
    real(dp), parameter :: pi4 = pi2 * pi2
    real(dp), parameter :: pi6 = pi4 * pi2
    real(dp), parameter :: ep128 = 2.0_dp**128
    real(dp), parameter :: em128 = 2.0_dp**(-128)

    public :: kstwo_cdf, kstwo_isf, kstwo_logcdf, kstwo_logpdf
    public :: kstwo_logsf, kstwo_pdf, kstwo_ppf, kstwo_sf
    public :: kstwo_logpdf_derivative

contains

    pure elemental function kstwo_cdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            y = standard_cdf(z, int(n))
        end if
    end function kstwo_cdf

    pure elemental function kstwo_sf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            y = standard_sf(z, int(n))
        end if
    end function kstwo_sf

    pure elemental function kstwo_logcdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y, p, q

        p = kstwo_cdf(x, n, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else if (p < 0.5_dp) then
            y = log(p)
        else
            q = kstwo_sf(x, n, loc, scale)
            y = log1p_safe(-q)
        end if
    end function kstwo_logcdf

    pure elemental function kstwo_logsf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y, p, q

        q = kstwo_sf(x, n, loc, scale)
        if (ieee_is_nan(q)) then
            y = q
        else if (q <= 0.0_dp) then
            y = negative_infinity(q)
        else if (q < 0.5_dp) then
            y = log(q)
        else
            p = kstwo_cdf(x, n, loc, scale)
            y = log1p_safe(-p)
        end if
    end function kstwo_logsf

    pure elemental function kstwo_pdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            y = standard_pdf(z, int(n)) / sigma
        end if
    end function kstwo_pdf

    pure elemental function kstwo_logpdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: y, f

        f = kstwo_pdf(x, n, loc, scale)
        if (ieee_is_nan(f)) then
            y = f
        else if (f <= 0.0_dp) then
            y = negative_infinity(f)
        else
            y = log(f)
        end if
    end function kstwo_logpdf

    pure elemental function kstwo_ppf(p, n, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: ni, iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
            return
        end if
        ni = int(n)
        if (p <= 0.0_dp) then
            x = mu + sigma * (0.5_dp / real(ni, dp))
        else if (p >= 1.0_dp) then
            x = mu + sigma
        else if (p > 0.5_dp) then
            x = kstwo_isf(1.0_dp - p, n, mu, sigma)
        else
            lo = 0.5_dp / real(ni, dp)
            hi = 1.0_dp
            do iter = 1, 80
                mid = 0.5_dp * (lo + hi)
                if (standard_cdf(mid, ni) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function kstwo_ppf

    pure elemental function kstwo_isf(p, n, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: ni, iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
            return
        end if
        ni = int(n)
        if (p >= 1.0_dp) then
            x = mu + sigma * (0.5_dp / real(ni, dp))
        else if (p <= 0.0_dp) then
            x = mu + sigma
        else
            lo = 0.5_dp / real(ni, dp)
            hi = 1.0_dp
            do iter = 1, 80
                mid = 0.5_dp * (lo + hi)
                if (standard_sf(mid, ni) > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function kstwo_isf

    pure elemental function kstwo_logpdf_derivative(z, n) result(g)
        real(dp), intent(in) :: z !! standardized interior point
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp) :: g
        real(dp) :: h, fm2, fm1, fp1, fp2, f0, dfdz
        integer :: ni

        if (.not. valid_shape(n) .or. z <= 0.5_dp / n .or. z >= 1.0_dp) then
            g = quiet_nan(z)
            return
        end if
        ni = int(n)
        h = derivative_step(z, ni)
        f0 = standard_pdf(z, ni)
        if (h <= 0.0_dp .or. f0 <= 0.0_dp) then
            g = quiet_nan(z)
            return
        end if
        fm2 = standard_pdf(z - 2.0_dp * h, ni)
        fm1 = standard_pdf(z - h, ni)
        fp1 = standard_pdf(z + h, ni)
        fp2 = standard_pdf(z + 2.0_dp * h, ni)
        dfdz = (fm2 - 8.0_dp * fm1 + 8.0_dp * fp1 - fp2) / (12.0_dp * h)
        g = dfdz / f0
    end function kstwo_logpdf_derivative

    pure elemental function standard_cdf(x, n) result(p)
        real(dp), intent(in) :: x !! standardized KS statistic
        integer, intent(in) :: n !! positive integer sample size
        real(dp) :: p
        real(dp) :: t, nx2, q

        if (x >= 1.0_dp) then
            p = 1.0_dp
            return
        else if (x <= 0.0_dp) then
            p = 0.0_dp
            return
        end if
        t = real(n, dp) * x
        if (t <= 1.0_dp) then
            if (t <= 0.5_dp) then
                p = 0.0_dp
            else if (n <= 140) then
                p = small_t_cdf(n, t)
            else if (n <= 140) then
                p = small_t_cdf(n, t)
            else
                p = exp(log_nfactorial_div_n_pow_n(n) + real(n, dp) * log(2.0_dp * t - 1.0_dp))
            end if
            p = clip_probability(p)
            return
        end if
        if (t >= real(n - 1, dp)) then
            q = 2.0_dp * (1.0_dp - x)**n
            p = clip_probability(1.0_dp - q)
            return
        end if
        if (x >= 0.5_dp) then
            q = 2.0_dp * ksone_sf(x, real(n, dp))
            p = clip_probability(1.0_dp - q)
            return
        end if

        nx2 = t * x
        if (n <= 140) then
            if (nx2 <= 4.0_dp) then
                p = dmtw_cdf(n, x)
            else
                q = 2.0_dp * ksone_sf(x, real(n, dp))
                p = 1.0_dp - q
            end if
        else
            if (nx2 >= 18.0_dp) then
                p = 1.0_dp
            else if (n <= 100000 .and. real(n, dp) * x**1.5_dp <= 1.4_dp) then
                p = dmtw_cdf(n, x)
            else
                p = pelz_good_cdf(n, x)
            end if
        end if
        p = clip_probability(p)
    end function standard_cdf

    pure elemental function standard_sf(x, n) result(q)
        real(dp), intent(in) :: x !! standardized KS statistic
        integer, intent(in) :: n !! positive integer sample size
        real(dp) :: q
        real(dp) :: t, nx2, p

        if (x >= 1.0_dp) then
            q = 0.0_dp
            return
        else if (x <= 0.0_dp) then
            q = 1.0_dp
            return
        end if
        t = real(n, dp) * x
        if (t <= 1.0_dp) then
            if (t <= 0.5_dp) then
                q = 1.0_dp
            else if (n <= 140) then
                p = small_t_cdf(n, t)
                q = 1.0_dp - p
            else
                p = exp(log_nfactorial_div_n_pow_n(n) + real(n, dp) * log(2.0_dp * t - 1.0_dp))
                q = 1.0_dp - p
            end if
            q = clip_probability(q)
            return
        end if
        if (t >= real(n - 1, dp)) then
            q = clip_probability(2.0_dp * (1.0_dp - x)**n)
            return
        end if
        if (x >= 0.5_dp) then
            q = clip_probability(2.0_dp * ksone_sf(x, real(n, dp)))
            return
        end if

        nx2 = t * x
        if (n <= 140) then
            if (nx2 <= 4.0_dp) then
                q = 1.0_dp - dmtw_cdf(n, x)
            else
                q = 2.0_dp * ksone_sf(x, real(n, dp))
            end if
        else if (nx2 >= 370.0_dp) then
            q = 0.0_dp
        else if (nx2 >= 2.2_dp) then
            q = 2.0_dp * ksone_sf(x, real(n, dp))
        else
            q = 1.0_dp - standard_cdf(x, n)
        end if
        q = clip_probability(q)
    end function standard_sf

    pure elemental function standard_pdf(x, n) result(f)
        real(dp), intent(in) :: x !! standardized KS statistic
        integer, intent(in) :: n !! positive integer sample size
        real(dp) :: f
        real(dp) :: t, h, fm2, fm1, fp1, fp2

        if (x <= 0.5_dp / real(n, dp) .or. x >= 1.0_dp) then
            f = 0.0_dp
            return
        end if
        t = real(n, dp) * x
        if (t <= 1.0_dp) then
            if (n <= 140) then
                f = small_t_pdf(n, t)
            else
                f = exp(log_nfactorial_div_n_pow_n(n) + real(n - 1, dp) * log(2.0_dp * t - 1.0_dp)) * &
                    2.0_dp * real(n, dp)**2
            end if
        else if (t >= real(n - 1, dp)) then
            f = 2.0_dp * real(n, dp) * (1.0_dp - x)**(n - 1)
        else if (x >= 0.5_dp) then
            f = 2.0_dp * ksone_pdf(x, real(n, dp))
        else
            h = cdf_derivative_step(x, n)
            if (h <= 0.0_dp) then
                f = 0.0_dp
            else
                fm2 = standard_cdf(x - 2.0_dp * h, n)
                fm1 = standard_cdf(x - h, n)
                fp1 = standard_cdf(x + h, n)
                fp2 = standard_cdf(x + 2.0_dp * h, n)
                f = (fm2 - 8.0_dp * fm1 + 8.0_dp * fp1 - fp2) / (12.0_dp * h)
                f = max(0.0_dp, f)
            end if
        end if
    end function standard_pdf

    pure elemental function cdf_derivative_step(x, n) result(h)
        real(dp), intent(in) :: x !! standardized interior point
        integer, intent(in) :: n !! positive integer sample size
        real(dp) :: h
        real(dp) :: left, right, knot, rn
        integer :: j

        rn = real(n, dp)
        h = max(32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x)), x / 65536.0_dp)
        j = floor(rn * x)
        left = max(0.5_dp / rn, real(j, dp) / rn)
        right = min(1.0_dp, real(j + 1, dp) / rn)
        if (abs(rn * x - real(nint(rn * x), dp)) <= 64.0_dp * epsilon(1.0_dp) * rn) then
            knot = real(nint(rn * x), dp) / rn
            left = max(0.5_dp / rn, knot - 1.0_dp / rn)
            right = min(1.0_dp, knot + 1.0_dp / rn)
        end if
        h = min(h, 0.24_dp * max(0.0_dp, x - left))
        h = min(h, 0.24_dp * max(0.0_dp, right - x))
        h = min(h, 0.24_dp * max(0.0_dp, 0.5_dp - x))
        if (h <= 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x))) then
            h = min(x / 65536.0_dp, 0.1_dp / rn)
        end if
    end function cdf_derivative_step

    pure elemental function derivative_step(x, n) result(h)
        real(dp), intent(in) :: x !! standardized interior point
        integer, intent(in) :: n !! positive integer sample size
        real(dp) :: h
        real(dp) :: lo, hi

        lo = 0.5_dp / real(n, dp)
        hi = 1.0_dp
        h = max(1.0e-5_dp * max(1.0_dp, abs(x)), 256.0_dp * epsilon(1.0_dp))
        h = min(h, 0.20_dp * (x - lo))
        h = min(h, 0.20_dp * (hi - x))
    end function derivative_step

    pure function dmtw_cdf(n, d) result(p)
        integer, intent(in) :: n !! positive integer sample size
        real(dp), intent(in) :: d !! standardized KS statistic
        real(dp) :: p
        real(dp), allocatable :: hmat(:, :), hpwr(:, :), v(:), w(:), tmp(:, :)
        real(dp) :: nd, h, fac, tt
        integer :: k, m, i, j, nn, expnt, hexpnt

        nd = real(n, dp) * d
        if (d >= 1.0_dp) then
            p = 1.0_dp
            return
        else if (nd <= 0.5_dp) then
            p = 0.0_dp
            return
        end if
        k = ceiling(nd)
        h = real(k, dp) - nd
        m = 2 * k - 1
        allocate(hmat(m, m), hpwr(m, m), tmp(m, m), v(m), w(m))
        hmat = 0.0_dp
        fac = 1.0_dp
        do j = 1, m
            w(j) = fac
            fac = fac / real(j, dp)
            v(j) = (1.0_dp - h**j) * fac
        end do
        tt = max(2.0_dp * h - 1.0_dp, 0.0_dp)**m - 2.0_dp * h**m
        v(m) = (1.0_dp + tt) * fac
        do i = 1, m - 1
            hmat(i:m, i + 1) = w(1:m - i + 1)
        end do
        hmat(:, 1) = v
        do j = 1, m
            hmat(m, j) = v(m - j + 1)
        end do

        hpwr = 0.0_dp
        do i = 1, m
            hpwr(i, i) = 1.0_dp
        end do
        nn = n
        expnt = 0
        hexpnt = 0
        do while (nn > 0)
            if (mod(nn, 2) == 1) then
                tmp = matmul(hpwr, hmat)
                hpwr = tmp
                expnt = expnt + hexpnt
            end if
            tmp = matmul(hmat, hmat)
            hmat = tmp
            hexpnt = 2 * hexpnt
            if (abs(hmat(k, k)) > ep128) then
                hmat = hmat / ep128
                hexpnt = hexpnt + 128
            end if
            nn = nn / 2
        end do

        p = hpwr(k, k)
        do i = 1, n
            p = real(i, dp) * p / real(n, dp)
            if (abs(p) < em128 .and. abs(p) > 0.0_dp) then
                p = p * ep128
                expnt = expnt - 128
            end if
        end do
        if (expnt /= 0 .and. abs(p) > 0.0_dp) p = scale(p, expnt)
        p = clip_probability(p)
    end function dmtw_cdf

    pure elemental function pelz_good_cdf(n, x) result(p)
        integer, intent(in) :: n !! positive integer sample size
        real(dp), intent(in) :: x !! standardized KS statistic
        real(dp) :: p
        real(dp) :: z, z2, z3, z4, z6, qlog, q, qpower
        real(dp) :: k1a, k1b, k2a, k2b, k2c, k3a, k3b, k3c, k3d
        real(dp) :: kval(4), coeff(4), msq, m4, m6, ks, ks2, kpi, qp, sqrt3z
        real(dp) :: k2extra, k3extra
        integer :: k, maxk, m

        if (x <= 0.0_dp) then
            p = 0.0_dp
            return
        else if (x >= 1.0_dp) then
            p = 1.0_dp
            return
        end if
        z = sqrt(real(n, dp)) * x
        z2 = z * z
        z3 = z2 * z
        z4 = z2 * z2
        z6 = z3 * z3
        qlog = -pi2 / (8.0_dp * z2)
        if (qlog < log(tiny(1.0_dp))) then
            p = 0.0_dp
            return
        end if
        q = exp(qlog)
        k1a = -z2
        k1b = pi2 / 4.0_dp
        k2a = 6.0_dp * z6 + 2.0_dp * z4
        k2b = (2.0_dp * z4 - 5.0_dp * z2) * pi2 / 4.0_dp
        k2c = pi4 * (1.0_dp - 2.0_dp * z2) / 16.0_dp
        k3d = pi6 * (5.0_dp - 30.0_dp * z2) / 64.0_dp
        k3c = pi4 * (-60.0_dp * z2 + 212.0_dp * z4) / 16.0_dp
        k3b = pi2 * (135.0_dp * z4 - 96.0_dp * z6) / 4.0_dp
        k3a = -30.0_dp * z6 - 90.0_dp * z4 * z4

        kval = 0.0_dp
        maxk = ceiling(16.0_dp * z / scifort_pi)
        do k = maxk, 1, -1
            m = 2 * k - 1
            msq = real(m * m, dp)
            m4 = msq * msq
            m6 = m4 * msq
            qpower = q**(8 * k)
            coeff = [1.0_dp, k1a + k1b * msq, k2a + k2b * msq + k2c * m4, &
                k3a + k3b * msq + k3c * m4 + k3d * m6]
            kval = kval * qpower + coeff
        end do
        kval = kval * q * sqrt_two_pi
        kval(1) = kval(1) / z
        kval(2) = kval(2) / (6.0_dp * z4)
        kval(3) = kval(3) / (72.0_dp * z4 * z3)
        kval(4) = kval(4) / (6480.0_dp * z4 * z6)

        q = exp(-pi2 / (2.0_dp * z2))
        sqrt3z = sqrt_three * z
        k2extra = 0.0_dp
        k3extra = 0.0_dp
        do k = maxk, 1, -1
            ks = real(k, dp)
            ks2 = ks * ks
            kpi = scifort_pi * ks
            qp = q**(k * k)
            k2extra = k2extra + ks2 * qp
            k3extra = k3extra + (sqrt3z + kpi) * (sqrt3z - kpi) * ks2 * qp
        end do
        kval(3) = kval(3) + k2extra * pi2 * sqrt_two_pi / (-36.0_dp * z3)
        kval(4) = kval(4) + k3extra * pi2 * sqrt_two_pi / (216.0_dp * z6)
        kval(2) = kval(2) / sqrt(real(n, dp))
        kval(3) = kval(3) / real(n, dp)
        kval(4) = kval(4) / (real(n, dp) * sqrt(real(n, dp)))
        p = clip_probability(sum(kval))
    end function pelz_good_cdf

    pure elemental function small_t_cdf(n, t) result(p)
        integer, intent(in) :: n !! positive integer sample size no larger than 140
        real(dp), intent(in) :: t !! product n*x in (1/2,1]
        real(dp) :: p
        integer :: i

        p = 1.0_dp
        do i = 1, n
            p = p * real(i, dp) / real(n, dp) * (2.0_dp * t - 1.0_dp)
        end do
    end function small_t_cdf

    pure elemental function small_t_pdf(n, t) result(f)
        integer, intent(in) :: n !! positive integer sample size no larger than 140
        real(dp), intent(in) :: t !! product n*x in (1/2,1]
        real(dp) :: f
        integer :: i

        f = 1.0_dp
        do i = 1, n - 1
            f = f * real(i, dp) / real(n, dp) * (2.0_dp * t - 1.0_dp)
        end do
        f = f * 2.0_dp * real(n, dp)**2
    end function small_t_pdf

    pure elemental function log_nfactorial_div_n_pow_n(n) result(y)
        integer, intent(in) :: n !! positive integer sample size
        real(dp) :: y
        real(dp) :: rn, r2, poly

        rn = 1.0_dp / real(n, dp)
        r2 = rn * rn
        poly = -0.029550653594771242_dp
        poly = poly * r2 + 0.0064102564102564103_dp
        poly = poly * r2 - 0.0019175269175269175_dp
        poly = poly * r2 + 0.00084175084175084175_dp
        poly = poly * r2 - 0.00059523809523809524_dp
        poly = poly * r2 + 0.00079365079365079365_dp
        poly = poly * r2 - 0.0027777777777777778_dp
        poly = poly * r2 + 0.083333333333333333_dp
        y = 0.5_dp * log(real(n, dp)) - real(n, dp) + 0.5_dp * log_two_pi + rn * poly
    end function log_nfactorial_div_n_pow_n

    pure elemental function clip_probability(p) result(q)
        real(dp), intent(in) :: p !! raw probability
        real(dp) :: q
        q = min(1.0_dp, max(0.0_dp, p))
    end function clip_probability

    pure elemental logical function valid_shape(n)
        real(dp), intent(in) :: n !! candidate positive integer sample size
        integer :: ni

        if (.not. ieee_is_finite(n) .or. n < 1.0_dp .or. n > real(huge(1), dp)) then
            valid_shape = .false.
        else
            ni = nint(n)
            valid_shape = abs(n - real(ni, dp)) <= 0.25_dp * spacing(n)
        end if
    end function valid_shape

    pure elemental subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location
        real(dp), intent(in), optional :: scale !! optional positive scale
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_kstwo
