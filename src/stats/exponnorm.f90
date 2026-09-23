! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Exponentially modified normal distribution, matching scipy.stats.exponnorm.

module scifort_exponnorm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_logcdf, normal_logpdf, normal_logsf
    implicit none
    private

    public :: exponnorm_cdf, exponnorm_isf, exponnorm_logcdf, exponnorm_logpdf
    public :: exponnorm_logsf, exponnorm_pdf, exponnorm_ppf, exponnorm_sf

contains

    pure elemental function exponnorm_pdf(x, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: k !! positive finite exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: y, ly
        ly = exponnorm_logpdf(x, k, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function exponnorm_pdf

    pure elemental function exponnorm_logpdf(x, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: k !! positive finite exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: y, mu, sigma, z, invk
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(k) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            invk = 1.0_dp / k
            y = invk * (0.5_dp * invk - z) + normal_logcdf(z - invk) - &
                log(k) - log(sigma)
        end if
    end function exponnorm_logpdf

    pure elemental function exponnorm_logcdf(x, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: k !! positive finite exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: y, mu, sigma, z, invk, a, b
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(k) .or. .not. valid_loc_scale(mu, sigma)) then
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
            invk = 1.0_dp / k
            a = normal_logcdf(z)
            b = invk * (0.5_dp * invk - z) + normal_logcdf(z - invk)
            y = logdiffexp_pair(a, b)
        end if
    end function exponnorm_logcdf

    pure elemental function exponnorm_logsf(x, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: k !! positive finite exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: y, mu, sigma, z, invk, a, b
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(k) .or. .not. valid_loc_scale(mu, sigma)) then
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
            invk = 1.0_dp / k
            a = normal_logsf(z)
            b = invk * (0.5_dp * invk - z) + normal_logcdf(z - invk)
            y = min(0.0_dp, logaddexp_pair(a, b))
        end if
    end function exponnorm_logsf

    pure elemental function exponnorm_cdf(x, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: k !! positive finite exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: y, ly
        ly = exponnorm_logcdf(x, k, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function exponnorm_cdf

    pure elemental function exponnorm_sf(x, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: k !! positive finite exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: y, ly
        ly = exponnorm_logsf(x, k, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function exponnorm_sf

    pure elemental function exponnorm_ppf(p, k, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: k !! positive finite exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(k) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            z = standard_quantile(p, k, .true.)
            x = affine(mu, sigma, z, p)
        end if
    end function exponnorm_ppf

    pure elemental function exponnorm_isf(p, k, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: k !! positive finite exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(k) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            z = standard_quantile(p, k, .false.)
            x = affine(mu, sigma, z, p)
        end if
    end function exponnorm_isf

    pure elemental function standard_quantile(p, k, lower_tail) result(z)
        real(dp), intent(in) :: p !! requested tail probability strictly between zero and one
        real(dp), intent(in) :: k !! positive exponential-to-normal scale ratio
        logical, intent(in) :: lower_tail !! true for CDF inversion, false for SF inversion
        real(dp) :: z, lo, hi, mid, value
        integer :: iter
        lo = -1.0_dp
        hi = max(1.0_dp, k + 1.0_dp)
        do
            if (lower_tail) then
                value = standard_cdf(lo, k)
                if (value <= p) exit
            else
                value = standard_sf(lo, k)
                if (value >= p) exit
            end if
            if (abs(lo) > 0.25_dp * huge(lo)) exit
            lo = 2.0_dp * lo
        end do
        do
            if (lower_tail) then
                value = standard_cdf(hi, k)
                if (value >= p) exit
            else
                value = standard_sf(hi, k)
                if (value <= p) exit
            end if
            if (hi > 0.25_dp * huge(hi)) exit
            hi = 2.0_dp * hi
        end do
        do iter = 1, 100
            mid = 0.5_dp * (lo + hi)
            if (lower_tail) then
                value = standard_cdf(mid, k)
                if (value < p) then
                    lo = mid
                else
                    hi = mid
                end if
            else
                value = standard_sf(mid, k)
                if (value > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end if
        end do
        z = 0.5_dp * (lo + hi)
    end function standard_quantile

    pure elemental function standard_cdf(z, k) result(p)
        real(dp), intent(in) :: z !! standardized variate
        real(dp), intent(in) :: k !! positive shape parameter
        real(dp) :: p, lp
        lp = standard_logcdf(z, k)
        if (lp == negative_infinity(lp)) then
            p = 0.0_dp
        else
            p = exp(lp)
        end if
    end function standard_cdf

    pure elemental function standard_sf(z, k) result(p)
        real(dp), intent(in) :: z !! standardized variate
        real(dp), intent(in) :: k !! positive shape parameter
        real(dp) :: p, lp
        lp = standard_logsf(z, k)
        if (lp == negative_infinity(lp)) then
            p = 0.0_dp
        else
            p = exp(lp)
        end if
    end function standard_sf

    pure elemental function standard_logcdf(z, k) result(y)
        real(dp), intent(in) :: z !! standardized variate
        real(dp), intent(in) :: k !! positive shape parameter
        real(dp) :: y, invk, a, b
        invk = 1.0_dp / k
        a = normal_logcdf(z)
        b = invk * (0.5_dp * invk - z) + normal_logcdf(z - invk)
        y = logdiffexp_pair(a, b)
    end function standard_logcdf

    pure elemental function standard_logsf(z, k) result(y)
        real(dp), intent(in) :: z !! standardized variate
        real(dp), intent(in) :: k !! positive shape parameter
        real(dp) :: y, invk, a, b
        invk = 1.0_dp / k
        a = normal_logsf(z)
        b = invk * (0.5_dp * invk - z) + normal_logcdf(z - invk)
        y = min(0.0_dp, logaddexp_pair(a, b))
    end function standard_logsf

    pure elemental function logaddexp_pair(a, b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y, m
        m = max(a, b)
        if (m == negative_infinity(m)) then
            y = m
        else
            y = m + log(exp(a - m) + exp(b - m))
        end if
    end function logaddexp_pair

    pure elemental function logdiffexp_pair(a, b) result(y)
        real(dp), intent(in) :: a !! logarithm of minuend
        real(dp), intent(in) :: b !! logarithm of subtrahend no larger than a
        real(dp) :: y, d
        if (b == negative_infinity(b)) then
            y = a
        else if (b >= a) then
            if (b - a <= 64.0_dp * epsilon(1.0_dp)) then
                y = negative_infinity(a)
            else
                y = quiet_nan(a)
            end if
        else
            d = b - a
            y = a + log(-expm1_safe(d))
        end if
    end function logdiffexp_pair

    pure elemental logical function valid_shape(k) result(ok)
        real(dp), intent(in) :: k !! candidate shape parameter
        ok = ieee_is_finite(k) .and. k > 0.0_dp
    end function valid_shape

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

    pure elemental function affine(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! standardized quantile
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z)) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else if (abs(z) > (huge(1.0_dp) - abs(mu)) / sigma) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else
            x = mu + sigma * z
        end if
    end function affine

end module scifort_exponnorm
