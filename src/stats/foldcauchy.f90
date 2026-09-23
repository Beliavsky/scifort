! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Folded-Cauchy distribution with nonnegative shape c, matching scipy.stats.foldcauchy.

module scifort_foldcauchy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: foldcauchy_cdf, foldcauchy_isf, foldcauchy_logcdf, foldcauchy_logpdf
    public :: foldcauchy_logsf, foldcauchy_pdf, foldcauchy_ppf, foldcauchy_sf

contains

    pure elemental function foldcauchy_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, logy
        logy = foldcauchy_logpdf(x, c, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function foldcauchy_pdf

    pure elemental function foldcauchy_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, a, b, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                a = inv_one_plus_square(z - c)
                b = inv_one_plus_square(z + c)
                y = log(a + b) - log(scifort_pi) - log(sigma)
            end if
        end if
    end function foldcauchy_logpdf

    pure elemental function foldcauchy_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            y = (atan(z - c) + atan(z + c)) / scifort_pi
            y = max(0.0_dp, min(1.0_dp, y))
        end if
    end function foldcauchy_cdf

    pure elemental function foldcauchy_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            y = (atan(1.0_dp, z - c) + atan(1.0_dp, z + c)) / scifort_pi
            y = max(0.0_dp, min(1.0_dp, y))
        end if
    end function foldcauchy_sf

    pure elemental function foldcauchy_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p, q
        p = foldcauchy_cdf(x, c, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(x)
        else if (p < 0.5_dp) then
            y = log(p)
        else
            q = foldcauchy_sf(x, c, loc, scale)
            y = log1p_safe(-q)
        end if
    end function foldcauchy_logcdf

    pure elemental function foldcauchy_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p, q
        q = foldcauchy_sf(x, c, loc, scale)
        if (ieee_is_nan(q)) then
            y = q
        else if (q == 0.0_dp) then
            y = negative_infinity(x)
        else if (q < 0.5_dp) then
            y = log(q)
        else
            p = foldcauchy_cdf(x, c, loc, scale)
            y = log1p_safe(-p)
        end if
    end function foldcauchy_logsf

    pure elemental function foldcauchy_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            z = standard_quantile(p, c, .true.)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function foldcauchy_ppf

    pure elemental function foldcauchy_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            z = standard_quantile(p, c, .false.)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function foldcauchy_isf

    pure elemental function standard_quantile(p, c, lower_tail) result(z)
        real(dp), intent(in) :: p !! requested tail probability strictly between zero and one
        real(dp), intent(in) :: c !! nonnegative folding offset
        logical, intent(in) :: lower_tail !! true for CDF inversion, false for SF inversion
        real(dp) :: z, lo, hi, mid, value
        integer :: iter
        lo = 0.0_dp
        hi = max(1.0_dp, c + 1.0_dp)
        do
            if (lower_tail) then
                value = (atan(hi - c) + atan(hi + c)) / scifort_pi
                if (value >= p) exit
            else
                value = (atan(1.0_dp, hi - c) + atan(1.0_dp, hi + c)) / scifort_pi
                if (value <= p) exit
            end if
            if (hi > 0.25_dp * huge(hi)) exit
            hi = 2.0_dp * hi
        end do
        do iter = 1, 100
            mid = 0.5_dp * (lo + hi)
            if (lower_tail) then
                value = (atan(mid - c) + atan(mid + c)) / scifort_pi
                if (value < p) then
                    lo = mid
                else
                    hi = mid
                end if
            else
                value = (atan(1.0_dp, mid - c) + atan(1.0_dp, mid + c)) / scifort_pi
                if (value > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end if
        end do
        z = 0.5_dp * (lo + hi)
    end function standard_quantile

    pure elemental function inv_one_plus_square(x) result(y)
        real(dp), intent(in) :: x !! real argument
        real(dp) :: y, a, r
        a = abs(x)
        if (a <= 1.0_dp) then
            y = 1.0_dp / (1.0_dp + x * x)
        else
            r = 1.0_dp / a
            y = r * r / (1.0_dp + r * r)
        end if
    end function inv_one_plus_square

    pure elemental logical function valid_shape(c) result(ok)
        real(dp), intent(in) :: c !! candidate folded-Cauchy shape parameter
        ok = ieee_is_finite(c) .and. c >= 0.0_dp
    end function valid_shape

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

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when absent
        real(dp), intent(out) :: sigma !! scale, 1 when absent
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_foldcauchy
