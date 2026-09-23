! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Burr Type III distribution matching scipy.stats.burr.
! In standardized coordinates z=(x-loc)/scale >= 0,
! f(z)=c*d*z**(-c-1)*(1+z**(-c))**(-d-1), c,d>0.
module scifort_burr
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, log1pexp, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: burr_cdf, burr_isf, burr_logcdf, burr_logpdf
    public :: burr_logsf, burr_pdf, burr_ppf, burr_sf

contains

    pure elemental function burr_logpdf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: cd, l, lz, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c) .and. valid_shape(d))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z == 0.0_dp) then
                cd = c * d
                if (cd < 1.0_dp) then
                    y = positive_infinity(x)
                else if (cd == 1.0_dp) then
                    y = log(c) + log(d) - log(sigma)
                else
                    y = negative_infinity(x)
                end if
            else
                lz = log(z)
                t = shape_times_log(-c, lz, x)
                l = log1pexp(t)
                y = log(c) + log(d) - (c + 1.0_dp) * lz - &
                    (d + 1.0_dp) * l - log(sigma)
            end if
        end if
    end function burr_logpdf

    pure elemental function burr_pdf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = burr_logpdf(x, c, d, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function burr_pdf

    pure elemental function burr_logcdf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: l, lz, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c) .and. valid_shape(d))) then
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
            t = shape_times_log(-c, lz, x)
            l = log1pexp(t)
            if (d > 1.0_dp .and. l > huge(1.0_dp) / d) then
                y = negative_infinity(x)
            else
                y = -d * l
            end if
        end if
    end function burr_logcdf

    pure elemental function burr_cdf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = burr_logcdf(x, c, d, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function burr_cdf

    pure elemental function burr_logsf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lcdf

        lcdf = burr_logcdf(x, c, d, loc, scale)
        if (ieee_is_nan(lcdf)) then
            y = lcdf
        else if (lcdf == negative_infinity(lcdf)) then
            y = 0.0_dp
        else if (lcdf == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(-expm1_safe(lcdf))
        end if
    end function burr_logsf

    pure elemental function burr_sf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lcdf

        lcdf = burr_logcdf(x, c, d, loc, scale)
        if (ieee_is_nan(lcdf)) then
            y = lcdf
        else if (lcdf == negative_infinity(lcdf)) then
            y = 1.0_dp
        else if (lcdf == 0.0_dp) then
            y = 0.0_dp
        else
            y = -expm1_safe(lcdf)
        end if
    end function burr_sf

    pure elemental function burr_ppf(p, c, d, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, q, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c) .and. valid_shape(d))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            q = -log(p) / d
            z = burr_quantile(q, c, p)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function burr_ppf

    pure elemental function burr_isf(p, c, d, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, q, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c) .and. valid_shape(d))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            q = -log1p_safe(-p) / d
            z = burr_quantile(q, c, p)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function burr_isf

    pure elemental function burr_quantile(q, c, seed) result(z)
        real(dp), intent(in) :: q !! positive logarithmic probability transform
        real(dp), intent(in) :: c !! positive finite Burr shape parameter
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: z
        real(dp) :: lbase, lz

        if (.not. ieee_is_finite(q)) then
            z = 0.0_dp
            return
        end if
        if (q > 40.0_dp) then
            lbase = q + log1p_safe(-exp(-q))
        else
            lbase = log(expm1_safe(q))
        end if
        lz = -lbase / c
        if (lz > log(huge(1.0_dp))) then
            z = positive_infinity(seed)
        else if (lz < log(tiny(1.0_dp))) then
            z = 0.0_dp
        else
            z = exp(lz)
        end if
    end function burr_quantile

    pure elemental function shape_times_log(c, lz, seed) result(t)
        real(dp), intent(in) :: c !! finite multiplier applied to a logarithm
        real(dp), intent(in) :: lz !! logarithm of a positive standardized variate
        real(dp), intent(in) :: seed !! value used to construct infinity when needed
        real(dp) :: t

        if (lz == 0.0_dp .or. c == 0.0_dp) then
            t = 0.0_dp
        else if (abs(c) > 1.0_dp .and. abs(lz) > huge(1.0_dp) / abs(c)) then
            if ((c > 0.0_dp) .eqv. (lz > 0.0_dp)) then
                t = positive_infinity(seed)
            else
                t = negative_infinity(seed)
            end if
        else
            t = c * lz
        end if
    end function shape_times_log

    pure elemental logical function valid_shape(a) result(ok)
        real(dp), intent(in) :: a !! candidate positive shape parameter
        ok = ieee_is_finite(a) .and. a > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! resolved location, 0 when absent
        real(dp), intent(out) :: sigma !! resolved scale, 1 when absent
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to construct positive infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z) .or. z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive

end module scifort_burr
