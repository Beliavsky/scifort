! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Skew-normal distribution matching scipy.stats.skewnorm.
module scifort_skewnorm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_logcdf, normal_logpdf
    implicit none
    private

    public :: skewnorm_cdf, skewnorm_isf, skewnorm_logcdf, skewnorm_logpdf
    public :: skewnorm_logsf, skewnorm_pdf, skewnorm_ppf, skewnorm_sf

    real(dp), parameter :: gl_x(32) = [ &
        -9.972638618494815699e-01_dp, -9.856115115452683817e-01_dp, &
        -9.647622555875063899e-01_dp, -9.349060759377396668e-01_dp, &
        -8.963211557660520912e-01_dp, -8.493676137325699704e-01_dp, &
        -7.944837959679423856e-01_dp, -7.321821187402897113e-01_dp, &
        -6.630442669302152314e-01_dp, -5.877157572407623043e-01_dp, &
        -5.068999089322293594e-01_dp, -4.213512761306353327e-01_dp, &
        -3.318686022821276671e-01_dp, -2.392873622521370647e-01_dp, &
        -1.444719615827964876e-01_dp, -4.830766568773832426e-02_dp, &
         4.830766568773832426e-02_dp,  1.444719615827964876e-01_dp, &
         2.392873622521370647e-01_dp,  3.318686022821276671e-01_dp, &
         4.213512761306353327e-01_dp,  5.068999089322293594e-01_dp, &
         5.877157572407623043e-01_dp,  6.630442669302152314e-01_dp, &
         7.321821187402897113e-01_dp,  7.944837959679423856e-01_dp, &
         8.493676137325699704e-01_dp,  8.963211557660520912e-01_dp, &
         9.349060759377396668e-01_dp,  9.647622555875063899e-01_dp, &
         9.856115115452683817e-01_dp,  9.972638618494815699e-01_dp ]
    real(dp), parameter :: gl_w(32) = [ &
        7.018610009470505756e-03_dp, 1.627439473090574323e-02_dp, &
        2.539206530926202410e-02_dp, 3.427386291302176452e-02_dp, &
        4.283589802222683568e-02_dp, 5.099805926237609144e-02_dp, &
        5.868409347853556501e-02_dp, 6.582222277636168295e-02_dp, &
        7.234579410884833806e-02_dp, 7.819389578707022781e-02_dp, &
        8.331192422694670696e-02_dp, 8.765209300440378326e-02_dp, &
        9.117387869576377979e-02_dp, 9.384439908080451087e-02_dp, &
        9.563872007927470831e-02_dp, 9.654008851472765940e-02_dp, &
        9.654008851472765940e-02_dp, 9.563872007927470831e-02_dp, &
        9.384439908080451087e-02_dp, 9.117387869576377979e-02_dp, &
        8.765209300440378326e-02_dp, 8.331192422694670696e-02_dp, &
        7.819389578707022781e-02_dp, 7.234579410884833806e-02_dp, &
        6.582222277636168295e-02_dp, 5.868409347853556501e-02_dp, &
        5.099805926237609144e-02_dp, 4.283589802222683568e-02_dp, &
        3.427386291302176452e-02_dp, 2.539206530926202410e-02_dp, &
        1.627439473090574323e-02_dp, 7.018610009470505756e-03_dp ]

contains

    pure elemental function skewnorm_logpdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! finite skewness shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            y = log(2.0_dp) + normal_logpdf(z) + normal_logcdf(a * z) - log(sigma)
        end if
    end function skewnorm_logpdf

    pure elemental function skewnorm_pdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! finite skewness shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = skewnorm_logpdf(x, a, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function skewnorm_pdf

    pure elemental function skewnorm_cdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! finite skewness shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 0.0_dp
        else if (x == positive_infinity(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            y = standard_cdf(z, a)
        end if
    end function skewnorm_cdf

    pure elemental function skewnorm_sf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! finite skewness shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 1.0_dp
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            y = standard_cdf(-z, -a)
        end if
    end function skewnorm_sf

    pure elemental function skewnorm_logcdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! finite skewness shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp .and. a > 0.0_dp .and. ieee_is_finite(z)) then
                y = tail_logcdf(z, a)
            else
                p = skewnorm_cdf(x, a, mu, sigma)
                if (p == 0.0_dp) then
                    y = negative_infinity(p)
                else
                    y = log(p)
                end if
            end if
        end if
    end function skewnorm_logcdf

    pure elemental function skewnorm_logsf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! finite skewness shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z > 0.0_dp .and. a < 0.0_dp .and. ieee_is_finite(z)) then
                y = tail_logcdf(-z, -a)
            else
                p = skewnorm_sf(x, a, mu, sigma)
                if (p == 0.0_dp) then
                    y = negative_infinity(p)
                else
                    y = log(p)
                end if
            end if
        end if
    end function skewnorm_logsf

    pure elemental function skewnorm_ppf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! finite skewness shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid, fmid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(a) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            lo = -40.0_dp
            hi = 40.0_dp
            do iter = 1, 140
                mid = 0.5_dp * (lo + hi)
                fmid = standard_cdf(mid, a)
                if (fmid < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function skewnorm_ppf

    pure elemental function skewnorm_isf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! finite skewness shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid, fmid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(a) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            lo = -40.0_dp
            hi = 40.0_dp
            do iter = 1, 140
                mid = 0.5_dp * (lo + hi)
                fmid = standard_cdf(-mid, -a)
                if (fmid > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function skewnorm_isf


    pure elemental function standard_cdf(z, a) result(p)
        real(dp), intent(in) :: z !! standardized evaluation point
        real(dp), intent(in) :: a !! finite skewness parameter
        real(dp) :: p
        if (z < 0.0_dp .and. a > 0.0_dp) then
            p = exp(tail_logcdf(z, a))
        else
            p = normal_cdf(z) - 2.0_dp * owens_t(z, a)
            p = min(1.0_dp, max(0.0_dp, p))
        end if
    end function standard_cdf

    pure elemental function tail_logcdf(z, a) result(logp)
        real(dp), intent(in) :: z !! negative standardized point
        real(dp), intent(in) :: a !! positive skewness parameter
        real(dp) :: logp
        real(dp) :: u, t, logterm, maxlog, sumv, center
        real(dp) :: logs(64)
        integer :: i, j, k
        maxlog = negative_infinity(z)
        k = 0
        do j = 1, 2
            center = 0.25_dp + 0.5_dp * real(j - 1, dp)
            do i = 1, size(gl_x)
                k = k + 1
                u = center + 0.25_dp * gl_x(i)
                t = z - u / (1.0_dp - u)
                logterm = log(2.0_dp) + normal_logpdf(t) + normal_logcdf(a * t) - &
                    2.0_dp * log(1.0_dp - u) + log(0.25_dp * gl_w(i))
                logs(k) = logterm
                maxlog = max(maxlog, logterm)
            end do
        end do
        sumv = 0.0_dp
        do i = 1, size(logs)
            sumv = sumv + exp(logs(i) - maxlog)
        end do
        logp = maxlog + log(sumv)
    end function tail_logcdf

    pure elemental function owens_t(h, a) result(t)
        real(dp), intent(in) :: h !! first Owen-T argument
        real(dp), intent(in) :: a !! second Owen-T argument
        real(dp) :: t
        real(dp) :: half, center, theta, ct, sum
        integer :: i

        if (a == 0.0_dp) then
            t = 0.0_dp
            return
        end if
        half = 0.5_dp * atan(a)
        center = half
        sum = 0.0_dp
        do i = 1, size(gl_x)
            theta = center + half * gl_x(i)
            ct = cos(theta)
            if (ct /= 0.0_dp) then
                sum = sum + gl_w(i) * exp(-0.5_dp * (h / ct)**2)
            end if
        end do
        t = half * sum / (2.0_dp * scifort_pi)
    end function owens_t

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

end module scifort_skewnorm
