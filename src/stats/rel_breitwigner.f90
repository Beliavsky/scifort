! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
!
! The complex CDF antiderivative is adapted from SciPy 1.17.0,
! scipy/stats/_continuous_distns.py, rel_breitwigner_gen._cdf.
! SciPy copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Licensed under BSD-3-Clause; see THIRD_PARTY_LICENSES.md and CODE_PROVENANCE.md.

! Relativistic Breit-Wigner distribution matching scipy.stats.rel_breitwigner.
module scifort_rel_breitwigner
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: rel_breitwigner_cdf, rel_breitwigner_isf, rel_breitwigner_logcdf
    public :: rel_breitwigner_logpdf, rel_breitwigner_logsf
    public :: rel_breitwigner_pdf, rel_breitwigner_ppf, rel_breitwigner_sf

contains

    pure elemental function rel_breitwigner_logpdf(x, rho, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: c, mu, q, sigma, u, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(rho) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            u = sqrt(1.0_dp + 1.0_dp / (rho * rho))
            c = 2.0_dp / scifort_pi * sqrt(2.0_dp * u * u / (1.0_dp + u))
            q = ((z - rho) * (z + rho)) / rho
            if (.not. ieee_is_finite(q)) then
                y = negative_infinity(x)
            else
                y = log(c) - log1p_safe(q * q) - log(sigma)
            end if
        end if
    end function rel_breitwigner_logpdf

    pure elemental function rel_breitwigner_pdf(x, rho, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = rel_breitwigner_logpdf(x, rho, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function rel_breitwigner_pdf

    pure elemental function rel_breitwigner_cdf(x, rho, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: c, mu, sigma, z
        complex(dp) :: a, d, w

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(rho) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (x == positive_infinity(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            c = sqrt(2.0_dp / (1.0_dp + sqrt(1.0_dp + 1.0_dp / (rho * rho)))) / &
                scifort_pi
            a = sqrt(cmplx(-1.0_dp, 1.0_dp / rho, kind=dp))
            d = sqrt(cmplx(-rho * rho, -rho, kind=dp))
            w = a * atan(cmplx(z, 0.0_dp, kind=dp) / d)
            y = 2.0_dp * c * aimag(w)
            y = min(1.0_dp, max(0.0_dp, y))
        end if
    end function rel_breitwigner_cdf

    pure elemental function rel_breitwigner_sf(x, rho, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: c, mu, sigma, z
        complex(dp) :: a, d, w

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(rho) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 1.0_dp
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            c = sqrt(2.0_dp / (1.0_dp + sqrt(1.0_dp + 1.0_dp / (rho * rho)))) / &
                scifort_pi
            a = sqrt(cmplx(-1.0_dp, 1.0_dp / rho, kind=dp))
            d = sqrt(cmplx(-rho * rho, -rho, kind=dp))
            w = a * atan(d / cmplx(z, 0.0_dp, kind=dp))
            y = 2.0_dp * c * aimag(w)
            y = min(1.0_dp, max(0.0_dp, y))
        end if
    end function rel_breitwigner_sf

    pure elemental function rel_breitwigner_logcdf(x, rho, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = rel_breitwigner_cdf(x, rho, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function rel_breitwigner_logcdf

    pure elemental function rel_breitwigner_logsf(x, rho, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = rel_breitwigner_sf(x, rho, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function rel_breitwigner_logsf

    pure elemental function rel_breitwigner_ppf(p, rho, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(rho) .or. .not. valid_loc_scale(mu, sigma) .or. &
                .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else if (p <= 0.5_dp) then
            z = lower_quantile(p, rho)
            x = mu + sigma * z
        else
            z = upper_quantile(1.0_dp - p, rho)
            x = mu + sigma * z
        end if
    end function rel_breitwigner_ppf

    pure elemental function rel_breitwigner_isf(p, rho, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(rho) .or. .not. valid_loc_scale(mu, sigma) .or. &
                .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p <= 0.5_dp) then
            z = upper_quantile(p, rho)
            x = mu + sigma * z
        else
            z = lower_quantile(1.0_dp - p, rho)
            x = mu + sigma * z
        end if
    end function rel_breitwigner_isf

    pure function lower_quantile(p, rho) result(z)
        real(dp), intent(in) :: p !! lower-tail probability strictly between 0 and 1
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp) :: z
        integer :: iter
        real(dp) :: hi, lo, mid

        lo = 0.0_dp
        hi = max(1.0_dp, rho)
        do iter = 1, 200
            if (rel_breitwigner_cdf(hi, rho) >= p) exit
            if (hi >= 0.25_dp * huge(hi)) then
                z = positive_infinity(p)
                return
            end if
            hi = 2.0_dp * hi
        end do
        do iter = 1, 140
            mid = 0.5_dp * (lo + hi)
            if (rel_breitwigner_cdf(mid, rho) < p) then
                lo = mid
            else
                hi = mid
            end if
            if (hi - lo <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        z = 0.5_dp * (lo + hi)
    end function lower_quantile

    pure function upper_quantile(p, rho) result(z)
        real(dp), intent(in) :: p !! upper-tail probability strictly between 0 and 1
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp) :: z
        integer :: iter
        real(dp) :: hi, lo, mid

        lo = 0.0_dp
        hi = max(1.0_dp, rho)
        do iter = 1, 200
            if (rel_breitwigner_sf(hi, rho) <= p) exit
            if (hi >= 0.25_dp * huge(hi)) then
                z = positive_infinity(p)
                return
            end if
            hi = 2.0_dp * hi
        end do
        do iter = 1, 140
            mid = 0.5_dp * (lo + hi)
            if (rel_breitwigner_sf(mid, rho) > p) then
                lo = mid
            else
                hi = mid
            end if
            if (hi - lo <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        z = 0.5_dp * (lo + hi)
    end function upper_quantile

    pure elemental logical function valid_shape(rho) result(valid)
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        valid = rho > 0.0_dp .and. ieee_is_finite(rho)
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location parameter
        real(dp), intent(in), optional :: scale !! optional positive scale parameter
        real(dp), intent(out) :: mu !! resolved location, default 0
        real(dp), intent(out) :: sigma !! resolved scale, default 1
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_rel_breitwigner
