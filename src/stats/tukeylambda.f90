! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Tukey lambda distribution matching scipy.stats.tukeylambda.
module scifort_tukeylambda
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: tukeylambda_cdf, tukeylambda_isf, tukeylambda_logcdf
    public :: tukeylambda_logpdf, tukeylambda_logsf, tukeylambda_pdf
    public :: tukeylambda_ppf, tukeylambda_sf

contains

    pure elemental function tukeylambda_ppf(p, lam, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(lam) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            if (lam > 0.0_dp) then
                x = mu - sigma / lam
            else
                x = negative_infinity(p)
            end if
        else if (p == 1.0_dp) then
            if (lam > 0.0_dp) then
                x = mu + sigma / lam
            else
                x = positive_infinity(p)
            end if
        else
            z = standard_ppf(p, lam)
            x = mu + sigma * z
        end if
    end function tukeylambda_ppf

    pure elemental function tukeylambda_isf(p, lam, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(lam) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            if (lam > 0.0_dp) then
                x = mu - sigma / lam
            else
                x = negative_infinity(p)
            end if
        else if (p == 0.0_dp) then
            if (lam > 0.0_dp) then
                x = mu + sigma / lam
            else
                x = positive_infinity(p)
            end if
        else
            z = -standard_ppf(p, lam)
            x = mu + sigma * z
        end if
    end function tukeylambda_isf

    pure elemental function tukeylambda_cdf(x, lam, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z, lo, hi, mid, qmid, bound
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(lam) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (lam > 0.0_dp) then
                bound = 1.0_dp / lam
                if (z <= -bound) then
                    y = 0.0_dp
                    return
                else if (z >= bound) then
                    y = 1.0_dp
                    return
                end if
            else if (z == negative_infinity(z)) then
                y = 0.0_dp
                return
            else if (z == positive_infinity(z)) then
                y = 1.0_dp
                return
            end if
            lo = 0.0_dp
            hi = 1.0_dp
            do iter = 1, 110
                mid = 0.5_dp * (lo + hi)
                qmid = standard_ppf(mid, lam)
                if (qmid < z) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            y = 0.5_dp * (lo + hi)
        end if
    end function tukeylambda_cdf

    pure elemental function tukeylambda_sf(x, lam, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(lam) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            y = tukeylambda_cdf(-z, lam)
        end if
    end function tukeylambda_sf

    pure elemental function tukeylambda_logcdf(x, lam, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = tukeylambda_cdf(x, lam, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function tukeylambda_logcdf

    pure elemental function tukeylambda_logsf(x, lam, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = tukeylambda_sf(x, lam, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function tukeylambda_logsf

    pure elemental function tukeylambda_logpdf(x, lam, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, p, lp, lq, a, b, m

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. ieee_is_finite(lam) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        p = tukeylambda_cdf(x, lam, mu, sigma)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp .or. p >= 1.0_dp) then
            y = negative_infinity(x)
        else
            lp = log(p)
            lq = log1p_safe(-p)
            a = (lam - 1.0_dp) * lp
            b = (lam - 1.0_dp) * lq
            m = max(a, b)
            y = -(m + log(exp(a - m) + exp(b - m))) - log(sigma)
        end if
    end function tukeylambda_logpdf

    pure elemental function tukeylambda_pdf(x, lam, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = tukeylambda_logpdf(x, lam, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function tukeylambda_pdf

    pure elemental function standard_ppf(p, lam) result(z)
        real(dp), intent(in) :: p !! probability strictly between zero and one
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp) :: z, lp, lq
        lp = log(p)
        lq = log1p_safe(-p)
        if (abs(lam) < 1.0e-7_dp) then
            z = (lp - lq) + 0.5_dp * lam * (lp * lp - lq * lq) + &
                lam * lam * (lp**3 - lq**3) / 6.0_dp
        else
            z = (expm1_safe(lam * lp) - expm1_safe(lam * lq)) / lam
        end if
    end function standard_ppf

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

end module scifort_tukeylambda
