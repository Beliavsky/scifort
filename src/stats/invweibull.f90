! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Inverse Weibull distribution matching scipy.stats.invweibull.
! In standardized coordinates z=(x-loc)/scale > 0,
! f(z)=c*z**(-c-1)*exp(-z**(-c)), c>0.
module scifort_invweibull
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private
    public :: invweibull_pdf, invweibull_logpdf, invweibull_cdf, invweibull_sf
    public :: invweibull_logcdf, invweibull_logsf, invweibull_ppf, invweibull_isf
contains
    pure elemental function invweibull_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lz, mu, sigma, s, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            lz = log(z)
            s = -shape_times_log(c, lz, x)
            if (s > log(huge(1.0_dp))) then
                y = negative_infinity(x)
            else
                if (s < log(tiny(1.0_dp))) then
                    t = 0.0_dp
                else
                    t = exp(s)
                end if
                if (abs(lz) > huge(1.0_dp) / (c + 1.0_dp)) then
                    y = negative_infinity(x)
                else
                    y = log(c) - (c + 1.0_dp) * lz - t - log(sigma)
                end if
            end if
        end if
    end function invweibull_logpdf

    pure elemental function invweibull_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = invweibull_logpdf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function invweibull_pdf

    pure elemental function invweibull_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lz, mu, sigma, s, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            lz = log(z)
            s = -shape_times_log(c, lz, x)
            if (s > log(huge(1.0_dp))) then
                y = negative_infinity(x)
            else if (s < log(tiny(1.0_dp))) then
                y = 0.0_dp
            else
                y = -exp(s)
            end if
        end if
    end function invweibull_logcdf

    pure elemental function invweibull_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability
        real(dp), intent(in) :: c !! positive finite inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = invweibull_logcdf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function invweibull_cdf

    pure elemental function invweibull_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = invweibull_logcdf(x, c, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 0.0_dp
        else if (lc == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(-expm1_safe(lc))
        end if
    end function invweibull_logsf

    pure elemental function invweibull_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability
        real(dp), intent(in) :: c !! positive finite inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = invweibull_logcdf(x, c, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 1.0_dp
        else
            y = -expm1_safe(lc)
        end if
    end function invweibull_sf

    pure elemental function invweibull_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive finite inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: e, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            e = -log(-log(p)) / c
            z = exp_quantile(e, p)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function invweibull_ppf

    pure elemental function invweibull_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive finite inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: e, mu, q, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            q = -log1p_safe(-p)
            e = -log(q) / c
            z = exp_quantile(e, p)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function invweibull_isf

    pure elemental function exp_quantile(e, seed) result(z)
        real(dp), intent(in) :: e !! logarithm of standardized quantile
        real(dp), intent(in) :: seed !! seed for infinity construction
        real(dp) :: z
        if (e > log(huge(1.0_dp))) then
            z = positive_infinity(seed)
        else if (e < log(tiny(1.0_dp))) then
            z = 0.0_dp
        else
            z = exp(e)
        end if
    end function exp_quantile

    pure elemental function shape_times_log(c, lz, seed) result(t)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in) :: lz !! logarithm of positive standardized variate
        real(dp), intent(in) :: seed !! value used to construct infinity when needed
        real(dp) :: t
        if (lz == 0.0_dp) then
            t = 0.0_dp
        else if (abs(lz) > huge(1.0_dp) / c) then
            if (lz > 0.0_dp) then
                t = positive_infinity(seed)
            else
                t = negative_infinity(seed)
            end if
        else
            t = c * lz
        end if
    end function shape_times_log

    pure elemental function valid_shape(c) result(ok)
        real(dp), intent(in) :: c !! candidate shape parameter
        logical :: ok
        ok = ieee_is_finite(c) .and. c > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to form infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z) .or. z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive
end module scifort_invweibull
