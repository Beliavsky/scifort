! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized exponential distribution matching scipy.stats.genexpon.
module scifort_genexpon
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: genexpon_cdf, genexpon_isf, genexpon_logcdf, genexpon_logpdf
    public :: genexpon_logsf, genexpon_pdf, genexpon_ppf, genexpon_sf

contains

    pure elemental function genexpon_logpdf(x, a, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive exponential-decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z, u, h

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, b, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                u = -expm1_safe(-c * z)
                h = a + b * u
                y = log(h) - (a + b) * z + b * u / c - log(sigma)
            end if
        end if
    end function genexpon_logpdf

    pure elemental function genexpon_pdf(x, a, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive exponential-decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = genexpon_logpdf(x, a, b, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function genexpon_pdf

    pure elemental function genexpon_logsf(x, a, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive exponential-decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z, u

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, b, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (.not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                u = -expm1_safe(-c * z)
                y = -(a + b) * z + b * u / c
            end if
        end if
    end function genexpon_logsf

    pure elemental function genexpon_sf(x, a, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive exponential-decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = genexpon_logsf(x, a, b, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function genexpon_sf

    pure elemental function genexpon_logcdf(x, a, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive exponential-decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ls

        ls = genexpon_logsf(x, a, b, c, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == 0.0_dp) then
            y = negative_infinity(x)
        else if (ls == negative_infinity(ls)) then
            y = 0.0_dp
        else
            y = log(-expm1_safe(ls))
        end if
    end function genexpon_logcdf

    pure elemental function genexpon_cdf(x, a, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive exponential-decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ls

        ls = genexpon_logsf(x, a, b, c, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == 0.0_dp) then
            y = 0.0_dp
        else if (ls == negative_infinity(ls)) then
            y = 1.0_dp
        else
            y = -expm1_safe(ls)
        end if
    end function genexpon_cdf

    pure elemental function genexpon_ppf(p, a, b, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive exponential-decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, s, t, arg, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, b, c) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            s = a + b
            t = (b - c * log1p_safe(-p)) / s
            arg = -(b / s) * exp(-t)
            z = (t + lambert_w0_negative(arg)) / c
            x = mu + sigma * z
        end if
    end function genexpon_ppf

    pure elemental function genexpon_isf(p, a, b, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive exponential-decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, s, t, arg, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, b, c) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            s = a + b
            t = (b - c * log(p)) / s
            arg = -(b / s) * exp(-t)
            z = (t + lambert_w0_negative(arg)) / c
            x = mu + sigma * z
        end if
    end function genexpon_isf

    pure elemental function lambert_w0_negative(x) result(w)
        real(dp), intent(in) :: x !! argument on the principal real branch, -1/e <= x <= 0
        real(dp) :: w
        real(dp) :: ew, f, denom, delta, q
        integer :: iter

        if (x == 0.0_dp) then
            w = 0.0_dp
            return
        end if
        if (x <= -exp(-1.0_dp)) then
            w = -1.0_dp
            return
        end if
        if (x > -0.05_dp) then
            w = x * (1.0_dp - x + 1.5_dp * x * x)
        else
            q = sqrt(max(0.0_dp, 2.0_dp * (1.0_dp + exp(1.0_dp) * x)))
            w = -1.0_dp + q - q * q / 3.0_dp + 11.0_dp * q**3 / 72.0_dp
        end if
        do iter = 1, 20
            ew = exp(w)
            f = w * ew - x
            denom = ew * (w + 1.0_dp) - (w + 2.0_dp) * f / (2.0_dp * (w + 1.0_dp))
            delta = f / denom
            w = w - delta
            if (abs(delta) <= 8.0_dp * epsilon(w) * max(1.0_dp, abs(w))) exit
        end do
    end function lambert_w0_negative

    pure elemental logical function valid_shape(a, b, c) result(ok)
        real(dp), intent(in) :: a !! first shape parameter
        real(dp), intent(in) :: b !! second shape parameter
        real(dp), intent(in) :: c !! third shape parameter
        ok = ieee_is_finite(a) .and. ieee_is_finite(b) .and. ieee_is_finite(c) .and. &
            a > 0.0_dp .and. b > 0.0_dp .and. c > 0.0_dp
    end function valid_shape

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

end module scifort_genexpon
