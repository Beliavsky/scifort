! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Rice distribution matching scipy.stats.rice.
module scifort_rice
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_incomplete_gamma, only : gammainc, gammaincc
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_gamma_one_plus
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    use scifort_special_elementary, only : i0e, i1e
    implicit none
    private

    public :: rice_cdf, rice_isf, rice_logcdf, rice_logpdf
    public :: rice_logsf, rice_pdf, rice_ppf, rice_sf
    public :: rice_i1_i0_ratio

contains

    pure elemental function rice_logpdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z, t, i0s

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                t = z * b
                i0s = i0e(t)
                y = log(z) - 0.5_dp * (z - b)**2 + log(i0s) - log(sigma)
            end if
        end if
    end function rice_logpdf

    pure elemental function rice_pdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = rice_logpdf(x, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function rice_pdf

    pure elemental function rice_cdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (.not. ieee_is_finite(z)) then
                y = 1.0_dp
            else
                y = rice_tail(z, b, .false.)
            end if
        end if
    end function rice_cdf

    pure elemental function rice_sf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 1.0_dp
            else if (.not. ieee_is_finite(z)) then
                y = 0.0_dp
            else
                y = rice_tail(z, b, .true.)
            end if
        end if
    end function rice_sf

    pure elemental function rice_logcdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = rice_cdf(x, b, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function rice_logcdf

    pure elemental function rice_logsf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = rice_sf(x, b, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function rice_logsf

    pure elemental function rice_ppf(p, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(b) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            lo = 0.0_dp
            hi = max(1.0_dp, b + 8.0_dp)
            do while (rice_tail(hi, b, .false.) < p)
                hi = 2.0_dp * hi
            end do
            do iter = 1, 110
                mid = 0.5_dp * (lo + hi)
                if (rice_tail(mid, b, .false.) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function rice_ppf

    pure elemental function rice_isf(p, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(b) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            lo = 0.0_dp
            hi = max(1.0_dp, b + 8.0_dp)
            do while (rice_tail(hi, b, .true.) > p)
                hi = 2.0_dp * hi
            end do
            do iter = 1, 110
                mid = 0.5_dp * (lo + hi)
                if (rice_tail(mid, b, .true.) > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function rice_isf

    pure elemental function rice_tail(z, b, upper) result(prob)
        real(dp), intent(in) :: z !! positive standardized variate
        real(dp), intent(in) :: b !! nonnegative Rice shape
        logical, intent(in) :: upper !! true for upper tail, false for lower tail
        real(dp) :: prob
        real(dp) :: m, t, weight, term, total, tol
        integer :: k0, k

        m = 0.5_dp * b * b
        t = 0.5_dp * z * z
        if (m == 0.0_dp) then
            if (upper) then
                prob = gammaincc(1.0_dp, t)
            else
                prob = gammainc(1.0_dp, t)
            end if
            return
        end if
        k0 = int(floor(m))
        weight = exp(-m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp)))
        if (upper) then
            term = weight * gammaincc(real(k0 + 1, dp), t)
        else
            term = weight * gammainc(real(k0 + 1, dp), t)
        end if
        total = term
        tol = 8.0_dp * epsilon(1.0_dp)

        weight = exp(-m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp)))
        do k = k0 - 1, 0, -1
            weight = weight * real(k + 1, dp) / m
            if (upper) then
                term = weight * gammaincc(real(k + 1, dp), t)
            else
                term = weight * gammainc(real(k + 1, dp), t)
            end if
            total = total + term
            if (weight < tol * max(total, tiny(1.0_dp)) .and. k < k0 / 2) exit
        end do

        weight = exp(-m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp)))
        do k = k0 + 1, k0 + 100000
            weight = weight * m / real(k, dp)
            if (upper) then
                term = weight * gammaincc(real(k + 1, dp), t)
            else
                term = weight * gammainc(real(k + 1, dp), t)
            end if
            total = total + term
            if (weight < tol * max(total, tiny(1.0_dp)) .and. k > k0 + 16) exit
        end do
        prob = min(1.0_dp, max(0.0_dp, total))
    end function rice_tail

    pure elemental function rice_i1_i0_ratio(x) result(r)
        real(dp), intent(in) :: x !! nonnegative Bessel argument
        real(dp) :: r
        if (x == 0.0_dp) then
            r = 0.0_dp
        else
            r = i1e(x) / i0e(x)
        end if
    end function rice_i1_i0_ratio

    pure elemental logical function valid_shape(b) result(ok)
        real(dp), intent(in) :: b !! Rice shape parameter
        ok = ieee_is_finite(b) .and. b >= 0.0_dp
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

end module scifort_rice
