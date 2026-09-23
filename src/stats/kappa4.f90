! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Four-parameter kappa distribution matching scipy.stats.kappa4.
module scifort_kappa4
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, log1pexp, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: kappa4_cdf, kappa4_isf, kappa4_logcdf, kappa4_logpdf
    public :: kappa4_logsf, kappa4_pdf, kappa4_ppf, kappa4_sf

contains

    pure elemental function kappa4_logpdf(x, h, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lower, upper, logbase, logt, logv, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(h, k) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            call support_bounds(h, k, lower, upper)
            if (z < lower .or. z > upper) then
                y = negative_infinity(x)
            else if (z == lower .and. ieee_is_finite(lower)) then
                y = lower_endpoint_logpdf(h, k, x) - log(sigma)
            else if (z == upper .and. ieee_is_finite(upper)) then
                y = upper_endpoint_logpdf(k, x) - log(sigma)
            else
                call transform_terms(z, h, k, logbase, logt, logv, t)
                if (h == 0.0_dp) then
                    if (.not. ieee_is_finite(t)) then
                        y = negative_infinity(x)
                    else
                        y = logt - logbase - t - log(sigma)
                    end if
                else
                    y = logt - logbase + (1.0_dp / h - 1.0_dp) * logv - log(sigma)
                end if
            end if
        end if
    end function kappa4_logpdf

    pure elemental function kappa4_pdf(x, h, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = kappa4_logpdf(x, h, k, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly) .or. ly >= log(huge(1.0_dp))) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function kappa4_pdf

    pure elemental function kappa4_logcdf(x, h, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lower, upper, logbase, logt, logv, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(h, k) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = negative_infinity(x)
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            call support_bounds(h, k, lower, upper)
            if (z <= lower) then
                y = negative_infinity(x)
            else if (z >= upper) then
                y = 0.0_dp
            else
                call transform_terms(z, h, k, logbase, logt, logv, t)
                if (h == 0.0_dp) then
                    if (.not. ieee_is_finite(t)) then
                        y = negative_infinity(x)
                    else
                        y = -t
                    end if
                else
                    y = logv / h
                end if
            end if
        end if
    end function kappa4_logcdf

    pure elemental function kappa4_cdf(x, h, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = kappa4_logcdf(x, h, k, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function kappa4_cdf

    pure elemental function kappa4_sf(x, h, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = kappa4_logcdf(x, h, k, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 1.0_dp
        else if (lc == 0.0_dp) then
            y = 0.0_dp
        else
            y = -expm1_safe(lc)
        end if
    end function kappa4_sf

    pure elemental function kappa4_logsf(x, h, k, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = kappa4_logcdf(x, h, k, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 0.0_dp
        else if (lc == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(-expm1_safe(lc))
        end if
    end function kappa4_logsf

    pure elemental function kappa4_ppf(p, h, k, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, lower, upper, logp, logt, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(h, k) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            call support_bounds(h, k, lower, upper)
            if (p == 0.0_dp) then
                x = affine_value(mu, sigma, lower, p)
            else if (p == 1.0_dp) then
                x = affine_value(mu, sigma, upper, p)
            else
                logp = log(p)
                logt = inverse_outer_logt(logp, h)
                z = inverse_inner(logt, k, p)
                x = affine_value(mu, sigma, z, p)
            end if
        end if
    end function kappa4_ppf

    pure elemental function kappa4_isf(p, h, k, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, lower, upper, logp, logt, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(h, k) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            call support_bounds(h, k, lower, upper)
            if (p == 1.0_dp) then
                x = affine_value(mu, sigma, lower, p)
            else if (p == 0.0_dp) then
                x = affine_value(mu, sigma, upper, p)
            else
                logp = log1p_safe(-p)
                logt = inverse_outer_logt(logp, h)
                z = inverse_inner(logt, k, p)
                x = affine_value(mu, sigma, z, p)
            end if
        end if
    end function kappa4_isf

    pure elemental subroutine transform_terms(z, h, k, logbase, logt, logv, t)
        real(dp), intent(in) :: z !! standardized point strictly inside support
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(out) :: logbase !! log(1-k*z), or zero for k=0
        real(dp), intent(out) :: logt !! logarithm of the inner transform
        real(dp), intent(out) :: logv !! log(1-h*T), or zero for h=0
        real(dp), intent(out) :: t !! inner transform T, or infinity when it overflows
        if (k == 0.0_dp) then
            logbase = 0.0_dp
            logt = -z
        else
            logbase = log1p_safe(-k * z)
            logt = logbase / k
        end if
        if (logt > log(huge(1.0_dp))) then
            t = positive_infinity(z)
        else if (logt < log(tiny(1.0_dp))) then
            t = 0.0_dp
        else
            t = exp(logt)
        end if
        if (h == 0.0_dp) then
            logv = 0.0_dp
        else if (h > 0.0_dp) then
            logv = log1p_safe(-h * t)
        else
            logv = log1pexp(log(-h) + logt)
        end if
    end subroutine transform_terms

    pure elemental function inverse_outer_logt(logp, h) result(logt)
        real(dp), intent(in) :: logp !! log of lower-tail probability in (-inf,0)
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp) :: logt, q
        if (h == 0.0_dp) then
            logt = log(-logp)
        else
            q = h * logp
            logt = log_abs_expm1(q) - log(abs(h))
        end if
    end function inverse_outer_logt

    pure elemental function inverse_inner(logt, k, seed) result(z)
        real(dp), intent(in) :: logt !! logarithm of positive inner transform
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in) :: seed !! probability used to construct infinity
        real(dp) :: z, q
        if (k == 0.0_dp) then
            z = -logt
        else
            q = k * logt
            if (q > log(huge(1.0_dp))) then
                if (k > 0.0_dp) then
                    z = negative_infinity(seed)
                else
                    z = positive_infinity(seed)
                end if
            else
                z = -expm1_safe(q) / k
            end if
        end if
    end function inverse_inner

    pure elemental function log_abs_expm1(x) result(y)
        real(dp), intent(in) :: x !! nonzero exponent in log(abs(exp(x)-1))
        real(dp) :: y
        if (x > 0.5_dp) then
            y = x + log1p_safe(-exp(-x))
        else if (x > 0.0_dp) then
            y = log(expm1_safe(x))
        else
            y = log(-expm1_safe(x))
        end if
    end function log_abs_expm1

    pure elemental subroutine support_bounds(h, k, lower, upper)
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(out) :: lower !! standardized lower support endpoint
        real(dp), intent(out) :: upper !! standardized upper support endpoint
        if (h > 0.0_dp) then
            lower = inverse_inner(-log(h), k, h)
        else if (k < 0.0_dp) then
            lower = 1.0_dp / k
        else
            lower = negative_infinity(h)
        end if
        if (k > 0.0_dp) then
            upper = 1.0_dp / k
        else
            upper = positive_infinity(k)
        end if
    end subroutine support_bounds

    pure elemental function lower_endpoint_logpdf(h, k, seed) result(y)
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in) :: seed !! value used to construct infinity
        real(dp) :: y, hk
        if (h > 0.0_dp) then
            if (h < 1.0_dp) then
                y = negative_infinity(seed)
            else if (h == 1.0_dp) then
                y = 0.0_dp
            else
                y = positive_infinity(seed)
            end if
        else if (h < 0.0_dp .and. k < 0.0_dp) then
            hk = h * k
            if (hk < 1.0_dp) then
                y = negative_infinity(seed)
            else if (hk == 1.0_dp) then
                y = (1.0_dp / h - 1.0_dp) * log(-h)
            else
                y = positive_infinity(seed)
            end if
        else
            y = negative_infinity(seed)
        end if
    end function lower_endpoint_logpdf

    pure elemental function upper_endpoint_logpdf(k, seed) result(y)
        real(dp), intent(in) :: k !! positive finite second kappa shape parameter
        real(dp), intent(in) :: seed !! value used to construct infinity
        real(dp) :: y
        if (k < 1.0_dp) then
            y = negative_infinity(seed)
        else if (k == 1.0_dp) then
            y = 0.0_dp
        else
            y = positive_infinity(seed)
        end if
    end function upper_endpoint_logpdf

    pure elemental logical function valid_shapes(h, k) result(ok)
        real(dp), intent(in) :: h !! candidate first kappa shape parameter
        real(dp), intent(in) :: k !! candidate second kappa shape parameter
        ok = ieee_is_finite(h) .and. ieee_is_finite(k)
    end function valid_shapes

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

    pure elemental function affine_value(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! standardized value
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
    end function affine_value

end module scifort_kappa4
