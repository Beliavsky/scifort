! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Double Pareto-lognormal distribution matching scipy.stats.dpareto_lognorm.
module scifort_dpareto_lognorm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_logcdf, normal_logpdf, normal_logsf
    implicit none
    private

    public :: dpareto_lognorm_cdf, dpareto_lognorm_isf, dpareto_lognorm_logcdf
    public :: dpareto_lognorm_logpdf, dpareto_lognorm_logsf, dpareto_lognorm_pdf
    public :: dpareto_lognorm_ppf, dpareto_lognorm_sf

contains

    pure elemental function dpareto_lognorm_logpdf(x, u, s, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, zstd, logx, z, x1, x2, lr1, lr2

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(u, s, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            zstd = (x - mu) / sigma
            if (zstd <= 0.0_dp .or. .not. ieee_is_finite(zstd)) then
                y = negative_infinity(x)
            else
                logx = log(zstd)
                z = (logx - u) / s
                x1 = a * s - z
                x2 = b * s + z
                lr1 = mills_log(x1)
                lr2 = mills_log(x2)
                y = log(a) + log(b) - log(a + b) - logx + normal_logpdf(z) + &
                    logaddexp(lr1, lr2) - log(sigma)
            end if
        end if
    end function dpareto_lognorm_logpdf

    pure elemental function dpareto_lognorm_pdf(x, u, s, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: y, ly

        ly = dpareto_lognorm_logpdf(x, u, s, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function dpareto_lognorm_pdf

    pure elemental function dpareto_lognorm_cdf(x, u, s, a, b, loc, scale) result(p)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: p
        real(dp) :: mu, sigma, zstd

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(u, s, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            p = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            p = quiet_nan(x)
        else
            zstd = (x - mu) / sigma
            if (zstd <= 0.0_dp) then
                p = 0.0_dp
            else if (.not. ieee_is_finite(zstd)) then
                p = 1.0_dp
            else
                p = standard_tail_logx(log(zstd), u, s, a, b, .false.)
            end if
        end if
    end function dpareto_lognorm_cdf

    pure elemental function dpareto_lognorm_sf(x, u, s, a, b, loc, scale) result(p)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: p
        real(dp) :: mu, sigma, zstd

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(u, s, a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            p = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            p = quiet_nan(x)
        else
            zstd = (x - mu) / sigma
            if (zstd <= 0.0_dp) then
                p = 1.0_dp
            else if (.not. ieee_is_finite(zstd)) then
                p = 0.0_dp
            else
                p = standard_tail_logx(log(zstd), u, s, a, b, .true.)
            end if
        end if
    end function dpareto_lognorm_sf

    pure elemental function dpareto_lognorm_logcdf(x, u, s, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: y, p

        p = dpareto_lognorm_cdf(x, u, s, a, b, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else if (p > 0.5_dp) then
            y = log1p_safe(-dpareto_lognorm_sf(x, u, s, a, b, loc, scale))
        else
            y = log(p)
        end if
    end function dpareto_lognorm_logcdf

    pure elemental function dpareto_lognorm_logsf(x, u, s, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: y, p

        p = dpareto_lognorm_sf(x, u, s, a, b, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else if (p > 0.5_dp) then
            y = log1p_safe(-dpareto_lognorm_cdf(x, u, s, a, b, loc, scale))
        else
            y = log(p)
        end if
    end function dpareto_lognorm_logsf

    pure elemental function dpareto_lognorm_ppf(p, u, s, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, t

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(u, s, a, b) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            t = inverse_logx(p, u, s, a, b, .false.)
            x = shifted_exp(t, mu, sigma)
        end if
    end function dpareto_lognorm_ppf

    pure elemental function dpareto_lognorm_isf(p, u, s, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, t

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(u, s, a, b) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            t = inverse_logx(p, u, s, a, b, .true.)
            x = shifted_exp(t, mu, sigma)
        end if
    end function dpareto_lognorm_isf

    pure elemental function standard_tail_logx(t, u, s, a, b, upper) result(p)
        real(dp), intent(in) :: t !! log of the positive standardized variate
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        logical, intent(in) :: upper !! true for survival probability
        real(dp) :: p
        real(dp) :: z, x1, x2, loga, logb, logc, logd

        z = (t - u) / s
        x1 = a * s - z
        x2 = b * s + z
        loga = normal_logcdf(z)
        logd = normal_logsf(z)
        logb = log(a) - log(a + b) + b * s * z + 0.5_dp * (b * s)**2 + normal_logsf(x2)
        logc = log(b) - log(a + b) - a * s * z + 0.5_dp * (a * s)**2 + normal_logsf(x1)
        if (upper) then
            if (z >= 0.0_dp) then
                p = signed_probability(logd, logc, logb)
            else
                p = 1.0_dp - signed_probability(loga, logb, logc)
            end if
        else
            if (z <= 0.0_dp) then
                p = signed_probability(loga, logb, logc)
            else
                p = 1.0_dp - signed_probability(logd, logc, logb)
            end if
        end if
        p = min(1.0_dp, max(0.0_dp, p))
    end function standard_tail_logx

    pure elemental function inverse_logx(p, u, s, a, b, upper) result(t)
        real(dp), intent(in) :: p !! requested lower or upper probability in (0,1)
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        logical, intent(in) :: upper !! true to invert the survival probability
        real(dp) :: t
        real(dp) :: lo, hi, mid, step
        integer :: iter

        lo = u
        hi = u
        step = max(1.0_dp, s)
        if (upper) then
            do while (standard_tail_logx(lo, u, s, a, b, .true.) < p)
                lo = lo - step
                step = 2.0_dp * step
            end do
            step = max(1.0_dp, s)
            do while (standard_tail_logx(hi, u, s, a, b, .true.) > p)
                hi = hi + step
                step = 2.0_dp * step
            end do
            do iter = 1, 100
                mid = 0.5_dp * (lo + hi)
                if (standard_tail_logx(mid, u, s, a, b, .true.) > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        else
            do while (standard_tail_logx(lo, u, s, a, b, .false.) > p)
                lo = lo - step
                step = 2.0_dp * step
            end do
            step = max(1.0_dp, s)
            do while (standard_tail_logx(hi, u, s, a, b, .false.) < p)
                hi = hi + step
                step = 2.0_dp * step
            end do
            do iter = 1, 100
                mid = 0.5_dp * (lo + hi)
                if (standard_tail_logx(mid, u, s, a, b, .false.) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        end if
        t = 0.5_dp * (lo + hi)
    end function inverse_logx

    pure elemental function shifted_exp(t, loc, scale) result(x)
        real(dp), intent(in) :: t !! logarithm of standardized positive variate
        real(dp), intent(in) :: loc !! outer location parameter
        real(dp), intent(in) :: scale !! positive outer scale
        real(dp) :: x

        if (t > log(huge(1.0_dp) / scale)) then
            x = positive_infinity(t)
        else if (t < log(tiny(1.0_dp))) then
            x = loc
        else
            x = loc + scale * exp(t)
        end if
    end function shifted_exp

    pure elemental function signed_probability(logx, logy, logz) result(p)
        real(dp), intent(in) :: logx !! log first positive term
        real(dp), intent(in) :: logy !! log second positive term
        real(dp), intent(in) :: logz !! log term to subtract
        real(dp) :: p
        real(dp) :: m, v

        m = max(logx, max(logy, logz))
        if (m == negative_infinity(m)) then
            p = 0.0_dp
        else
            v = exp(logx - m) + exp(logy - m) - exp(logz - m)
            p = exp(m) * max(0.0_dp, v)
        end if
    end function signed_probability

    pure elemental function mills_log(x) result(y)
        real(dp), intent(in) :: x !! argument of the normal Mills ratio
        real(dp) :: y

        y = normal_logsf(x) - normal_logpdf(x)
    end function mills_log

    pure elemental function logaddexp(x, y) result(z)
        real(dp), intent(in) :: x !! first logarithm
        real(dp), intent(in) :: y !! second logarithm
        real(dp) :: z
        real(dp) :: m

        m = max(x, y)
        if (m == negative_infinity(m)) then
            z = m
        else
            z = m + log(exp(x - m) + exp(y - m))
        end if
    end function logaddexp

    pure elemental logical function valid_shapes(u, s, a, b) result(ok)
        real(dp), intent(in) :: u !! candidate finite lognormal location shape
        real(dp), intent(in) :: s !! candidate positive lognormal scale shape
        real(dp), intent(in) :: a !! candidate positive upper-Pareto shape
        real(dp), intent(in) :: b !! candidate positive lower-Pareto shape
        ok = ieee_is_finite(u) .and. ieee_is_finite(s) .and. s > 0.0_dp .and. &
            ieee_is_finite(a) .and. a > 0.0_dp .and. ieee_is_finite(b) .and. b > 0.0_dp
    end function valid_shapes

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

end module scifort_dpareto_lognorm
