! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! ARGUS distribution matching scipy.stats.argus.
module scifort_argus
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi, scifort_log_two
    use scifort_incomplete_gamma, only : log_gammainc
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: argus_cdf, argus_isf, argus_logcdf, argus_logpdf
    public :: argus_logsf, argus_pdf, argus_ppf, argus_sf

contains

    pure elemental function argus_logpdf(x, chi, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: y, mu, sigma, z, one_minus_z2
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(chi) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp .or. z >= 1.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                one_minus_z2 = (1.0_dp - z) * (1.0_dp + z)
                y = 3.0_dp * log(chi) - scifort_log_sqrt_two_pi - argus_logphi(chi) + &
                    log(z) + 0.5_dp * log(one_minus_z2) - &
                    0.5_dp * chi * chi * one_minus_z2 - log(sigma)
            end if
        end if
    end function argus_logpdf

    pure elemental function argus_pdf(x, chi, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: y, ly
        ly = argus_logpdf(x, chi, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function argus_pdf

    pure elemental function argus_logsf(x, chi, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: y, mu, sigma, z, t
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(chi) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (z >= 1.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                t = chi * sqrt((1.0_dp - z) * (1.0_dp + z))
                y = argus_logphi(t) - argus_logphi(chi)
            end if
        end if
    end function argus_logsf

    pure elemental function argus_sf(x, chi, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: y, ly
        ly = argus_logsf(x, chi, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function argus_sf

    pure elemental function argus_logcdf(x, chi, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: y, ls
        ls = argus_logsf(x, chi, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == 0.0_dp) then
            y = negative_infinity(x)
        else if (ls == negative_infinity(ls)) then
            y = 0.0_dp
        else
            y = log(-expm1_safe(ls))
        end if
    end function argus_logcdf

    pure elemental function argus_cdf(x, chi, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: y, ls
        ls = argus_logsf(x, chi, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == 0.0_dp) then
            y = 0.0_dp
        else if (ls == negative_infinity(ls)) then
            y = 1.0_dp
        else
            y = -expm1_safe(ls)
        end if
    end function argus_cdf

    pure elemental function argus_ppf(p, chi, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(chi) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = mu + sigma
        else
            z = standard_ppf(p, chi)
            x = mu + sigma * z
        end if
    end function argus_ppf

    pure elemental function argus_isf(p, chi, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(chi) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = mu + sigma
        else
            z = standard_isf(p, chi)
            x = mu + sigma * z
        end if
    end function argus_isf

    pure elemental function standard_ppf(p, chi) result(z)
        real(dp), intent(in) :: p !! lower-tail probability strictly between zero and one
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp) :: z, lo, hi, mid, q
        integer :: iter
        lo = 0.0_dp
        hi = 1.0_dp
        if (p <= 0.5_dp) then
            do iter = 1, 80
                mid = 0.5_dp * (lo + hi)
                if (standard_cdf(mid, chi) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        else
            q = 1.0_dp - p
            do iter = 1, 80
                mid = 0.5_dp * (lo + hi)
                if (standard_sf(mid, chi) > q) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        end if
        z = 0.5_dp * (lo + hi)
    end function standard_ppf

    pure elemental function standard_isf(p, chi) result(z)
        real(dp), intent(in) :: p !! upper-tail probability strictly between zero and one
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp) :: z, lo, hi, mid, q
        integer :: iter
        lo = 0.0_dp
        hi = 1.0_dp
        if (p <= 0.5_dp) then
            do iter = 1, 80
                mid = 0.5_dp * (lo + hi)
                if (standard_sf(mid, chi) > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        else
            q = 1.0_dp - p
            do iter = 1, 80
                mid = 0.5_dp * (lo + hi)
                if (standard_cdf(mid, chi) < q) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        end if
        z = 0.5_dp * (lo + hi)
    end function standard_isf

    pure elemental function standard_sf(z, chi) result(y)
        real(dp), intent(in) :: z !! standardized point in (0,1)
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp) :: y, t
        t = chi * sqrt((1.0_dp - z) * (1.0_dp + z))
        y = exp(argus_logphi(t) - argus_logphi(chi))
    end function standard_sf

    pure elemental function standard_cdf(z, chi) result(y)
        real(dp), intent(in) :: z !! standardized point in (0,1)
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp) :: y, t, ls
        t = chi * sqrt((1.0_dp - z) * (1.0_dp + z))
        ls = argus_logphi(t) - argus_logphi(chi)
        y = -expm1_safe(ls)
    end function standard_cdf

    pure elemental function argus_logphi(chi) result(y)
        real(dp), intent(in) :: chi !! nonnegative ARGUS normalization argument
        real(dp) :: y
        if (chi == 0.0_dp) then
            y = negative_infinity(chi)
        else
            y = log_gammainc(1.5_dp, 0.5_dp * chi * chi) - scifort_log_two
        end if
    end function argus_logphi

    pure elemental logical function valid_shape(chi) result(ok)
        real(dp), intent(in) :: chi !! candidate ARGUS shape parameter
        ok = ieee_is_finite(chi) .and. chi > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! resolved lower support endpoint
        real(dp), intent(out) :: sigma !! resolved support width
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_argus
