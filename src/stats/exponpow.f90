! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Exponential-power distribution matching scipy.stats.exponpow.
! Standard density is b*z**(b-1)*exp(1+z**b-exp(z**b)), z>=0, b>0.
module scifort_exponpow
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private
    public :: exponpow_pdf, exponpow_logpdf, exponpow_cdf, exponpow_sf
    public :: exponpow_logcdf, exponpow_logsf, exponpow_ppf, exponpow_isf
contains
    pure elemental function exponpow_logpdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: b !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lz, mu, sigma, xb, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z == 0.0_dp) then
                if (b < 1.0_dp) then
                    y = positive_infinity(x)
                else if (b == 1.0_dp) then
                    y = -log(sigma)
                else
                    y = negative_infinity(x)
                end if
            else
                lz = log(z)
                xb = safe_power_from_log(b * lz, x)
                if (.not. ieee_is_finite(xb) .or. xb > log(huge(1.0_dp))) then
                    y = negative_infinity(x)
                else
                    y = 1.0_dp + log(b) + (b - 1.0_dp) * lz + xb - exp(xb) - log(sigma)
                end if
            end if
        end if
    end function exponpow_logpdf

    pure elemental function exponpow_pdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: b !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = exponpow_logpdf(x, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function exponpow_pdf

    pure elemental function exponpow_logsf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: b !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lz, mu, sigma, xb, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            lz = log(z)
            xb = safe_power_from_log(b * lz, x)
            if (.not. ieee_is_finite(xb) .or. xb > log(huge(1.0_dp))) then
                y = negative_infinity(x)
            else
                y = -expm1_safe(xb)
            end if
        end if
    end function exponpow_logsf

    pure elemental function exponpow_sf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability
        real(dp), intent(in) :: b !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = exponpow_logsf(x, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function exponpow_sf

    pure elemental function exponpow_cdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability
        real(dp), intent(in) :: b !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ls
        ls = exponpow_logsf(x, b, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == negative_infinity(ls)) then
            y = 1.0_dp
        else
            y = -expm1_safe(ls)
        end if
    end function exponpow_cdf

    pure elemental function exponpow_logcdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: b !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ls
        ls = exponpow_logsf(x, b, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == 0.0_dp) then
            y = negative_infinity(x)
        else if (ls == negative_infinity(ls)) then
            y = 0.0_dp
        else
            y = log(-expm1_safe(ls))
        end if
    end function exponpow_logcdf

    pure elemental function exponpow_ppf(p, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: b !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            t = log1p_safe(-log1p_safe(-p))
            z = exp(log(t) / b)
            x = affine(mu, sigma, z, p)
        end if
    end function exponpow_ppf

    pure elemental function exponpow_isf(p, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: b !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            t = log1p_safe(-log(p))
            z = exp(log(t) / b)
            x = affine(mu, sigma, z, p)
        end if
    end function exponpow_isf

    pure elemental function safe_power_from_log(log_value, seed) result(value)
        real(dp), intent(in) :: log_value !! logarithm of nonnegative power value
        real(dp), intent(in) :: seed !! seed used to construct positive infinity
        real(dp) :: value
        if (log_value > log(huge(1.0_dp))) then
            value = positive_infinity(seed)
        else if (log_value < log(tiny(1.0_dp))) then
            value = 0.0_dp
        else
            value = exp(log_value)
        end if
    end function safe_power_from_log

    pure elemental function valid_shape(b) result(ok)
        real(dp), intent(in) :: b !! candidate shape parameter
        logical :: ok
        ok = ieee_is_finite(b) .and. b > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! requested location
        real(dp), intent(in), optional :: scale !! requested scale
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure elemental function affine(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! standardized quantile
        real(dp), intent(in) :: seed !! seed used to construct positive infinity
        real(dp) :: x
        if (.not. ieee_is_finite(z)) then
            x = positive_infinity(seed)
        else if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine
end module scifort_exponpow
