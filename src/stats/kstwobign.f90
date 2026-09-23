! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Asymptotic two-sided Kolmogorov distribution matching scipy.stats.kstwobign.
module scifort_kstwobign
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: kstwobign_cdf, kstwobign_isf, kstwobign_logcdf, kstwobign_logpdf
    public :: kstwobign_logsf, kstwobign_pdf, kstwobign_ppf, kstwobign_sf
    public :: kstwobign_logpdf_derivative

contains

    pure elemental function kstwobign_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else
                y = standard_pdf(z) / sigma
            end if
        end if
    end function kstwobign_pdf

    pure elemental function kstwobign_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p

        p = kstwobign_pdf(x, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function kstwobign_logpdf

    pure elemental function kstwobign_cdf(x, loc, scale) result(p)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: p
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            p = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            p = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                p = 0.0_dp
            else
                p = standard_cdf(z)
            end if
        end if
    end function kstwobign_cdf

    pure elemental function kstwobign_sf(x, loc, scale) result(p)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: p
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            p = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            p = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                p = 1.0_dp
            else
                p = standard_sf(z)
            end if
        end if
    end function kstwobign_sf

    pure elemental function kstwobign_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p

        p = kstwobign_cdf(x, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else if (p > 0.5_dp) then
            y = log1p_safe(-kstwobign_sf(x, loc, scale))
        else
            y = log(p)
        end if
    end function kstwobign_logcdf

    pure elemental function kstwobign_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p

        p = kstwobign_sf(x, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else if (p > 0.5_dp) then
            y = log1p_safe(-kstwobign_cdf(x, loc, scale))
        else
            y = log(p)
        end if
    end function kstwobign_logsf

    pure elemental function kstwobign_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma) .or. .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            lo = 0.0_dp
            hi = 1.0_dp
            do while (standard_cdf(hi) < p)
                hi = 2.0_dp * hi
            end do
            do iter = 1, 80
                mid = 0.5_dp * (lo + hi)
                if (standard_cdf(mid) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function kstwobign_ppf

    pure elemental function kstwobign_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma) .or. .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            lo = 0.0_dp
            hi = 1.0_dp
            do while (standard_sf(hi) > p)
                hi = 2.0_dp * hi
            end do
            do iter = 1, 80
                mid = 0.5_dp * (lo + hi)
                if (standard_sf(mid) > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function kstwobign_isf

    pure elemental function kstwobign_logpdf_derivative(z) result(g)
        real(dp), intent(in) :: z !! positive standardized variate
        real(dp) :: g
        real(dp) :: f, fp

        call standard_pdf_and_derivative(z, f, fp)
        if (f > 0.0_dp) then
            g = fp / f
        else
            g = quiet_nan(z)
        end if
    end function kstwobign_logpdf_derivative

    pure elemental function standard_cdf(z) result(p)
        real(dp), intent(in) :: z !! positive standardized variate
        real(dp) :: p
        real(dp) :: a, term, total
        integer :: k

        if (z < 0.82_dp) then
            total = 0.0_dp
            do k = 1, 100000
                a = real(2 * k - 1, dp)**2 * scifort_pi**2 / (8.0_dp * z * z)
                term = exp(-a)
                total = total + term
                if (term <= epsilon(1.0_dp) * max(total, tiny(1.0_dp))) exit
            end do
            p = sqrt(2.0_dp * scifort_pi) * total / z
        else
            p = 1.0_dp - standard_sf(z)
        end if
        p = min(1.0_dp, max(0.0_dp, p))
    end function standard_cdf

    pure elemental function standard_sf(z) result(p)
        real(dp), intent(in) :: z !! positive standardized variate
        real(dp) :: p
        real(dp) :: term, total
        integer :: k

        if (z < 0.82_dp) then
            p = 1.0_dp - standard_cdf_small(z)
        else
            total = 0.0_dp
            do k = 1, 100000
                term = 2.0_dp * (-1.0_dp)**(k - 1) * exp(-2.0_dp * real(k * k, dp) * z * z)
                total = total + term
                if (abs(term) <= epsilon(1.0_dp) * max(abs(total), tiny(1.0_dp))) exit
            end do
            p = total
        end if
        p = min(1.0_dp, max(0.0_dp, p))
    end function standard_sf

    pure elemental function standard_cdf_small(z) result(p)
        real(dp), intent(in) :: z !! positive standardized variate below the series switch
        real(dp) :: p
        real(dp) :: a, term, total
        integer :: k

        total = 0.0_dp
        do k = 1, 100000
            a = real(2 * k - 1, dp)**2 * scifort_pi**2 / (8.0_dp * z * z)
            term = exp(-a)
            total = total + term
            if (term <= epsilon(1.0_dp) * max(total, tiny(1.0_dp))) exit
        end do
        p = sqrt(2.0_dp * scifort_pi) * total / z
        p = min(1.0_dp, max(0.0_dp, p))
    end function standard_cdf_small

    pure elemental function standard_pdf(z) result(f)
        real(dp), intent(in) :: z !! positive standardized variate
        real(dp) :: f, fp

        call standard_pdf_and_derivative(z, f, fp)
    end function standard_pdf

    pure elemental subroutine standard_pdf_and_derivative(z, f, fp)
        real(dp), intent(in) :: z !! positive standardized variate
        real(dp), intent(out) :: f !! standardized density
        real(dp), intent(out) :: fp !! derivative of standardized density
        real(dp) :: a, e, e0, e1, e2, s2, s4, signv
        integer :: k

        if (z < 0.82_dp) then
            e0 = 0.0_dp
            e1 = 0.0_dp
            e2 = 0.0_dp
            do k = 1, 100000
                a = real(2 * k - 1, dp)**2 * scifort_pi**2 / 8.0_dp
                e = exp(-a / (z * z))
                e0 = e0 + e
                e1 = e1 + a * e
                e2 = e2 + a * a * e
                if (e <= epsilon(1.0_dp) * max(e0, tiny(1.0_dp))) exit
            end do
            f = sqrt(2.0_dp * scifort_pi) * &
                (-e0 / (z * z) + 2.0_dp * e1 / (z**4))
            fp = sqrt(2.0_dp * scifort_pi) * &
                (2.0_dp * e0 / (z**3) - 10.0_dp * e1 / (z**5) + 4.0_dp * e2 / (z**7))
        else
            s2 = 0.0_dp
            s4 = 0.0_dp
            do k = 1, 100000
                if (mod(k, 2) == 1) then
                    signv = 1.0_dp
                else
                    signv = -1.0_dp
                end if
                e = exp(-2.0_dp * real(k * k, dp) * z * z)
                s2 = s2 + signv * real(k * k, dp) * e
                s4 = s4 + signv * real(k * k, dp)**2 * e
                if (e * real(k * k, dp) <= epsilon(1.0_dp) * max(abs(s2), tiny(1.0_dp))) exit
            end do
            f = 8.0_dp * z * s2
            fp = 8.0_dp * s2 - 32.0_dp * z * z * s4
        end if
        f = max(0.0_dp, f)
    end subroutine standard_pdf_and_derivative

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

end module scifort_kstwobign
