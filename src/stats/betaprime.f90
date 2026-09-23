! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Beta-prime distribution matching scipy.stats.betaprime.
! In standardized coordinates z=(x-loc)/scale > 0,
! f(z)=z**(a-1)*(1+z)**(-a-b)/B(a,b), a,b>0.
module scifort_betaprime
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_beta, only : beta_inverse_xy, incomplete_beta_xy
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private
    public :: betaprime_pdf, betaprime_logpdf, betaprime_cdf, betaprime_sf
    public :: betaprime_logcdf, betaprime_logsf, betaprime_ppf, betaprime_isf
contains
    pure elemental function betaprime_logpdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive finite first shape parameter
        real(dp), intent(in) :: b !! positive finite second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z == 0.0_dp) then
                if (a < 1.0_dp) then
                    y = positive_infinity(x)
                else if (a == 1.0_dp) then
                    y = log(b) - log(sigma)
                else
                    y = negative_infinity(x)
                end if
            else
                y = (a - 1.0_dp) * log(z) - (a + b) * log1p_safe(z) - &
                    log_beta(a, b) - log(sigma)
            end if
        end if
    end function betaprime_logpdf

    pure elemental function betaprime_pdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive finite first shape parameter
        real(dp), intent(in) :: b !! positive finite second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = betaprime_logpdf(x, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function betaprime_pdf

    pure elemental subroutine beta_prime_tails(z, a, b, p, q, logp, logq)
        real(dp), intent(in) :: z !! positive standardized variate
        real(dp), intent(in) :: a !! positive first shape
        real(dp), intent(in) :: b !! positive second shape
        real(dp), intent(out) :: p !! lower-tail probability
        real(dp), intent(out) :: q !! upper-tail probability
        real(dp), intent(out) :: logp !! logarithm of lower tail
        real(dp), intent(out) :: logq !! logarithm of upper tail
        real(dp) :: u, v, pu, qu, lpu, lqu
        if (z <= 1.0_dp) then
            u = z / (1.0_dp + z)
            v = 1.0_dp / (1.0_dp + z)
            call incomplete_beta_xy(a, b, u, v, p, q, logp, logq)
        else
            v = 1.0_dp / (1.0_dp + z)
            u = 1.0_dp - v
            call incomplete_beta_xy(b, a, v, u, pu, qu, lpu, lqu)
            p = qu
            q = pu
            logp = lqu
            logq = lpu
        end if
    end subroutine beta_prime_tails

    pure elemental function betaprime_cdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability
        real(dp), intent(in) :: a !! positive finite first shape parameter
        real(dp), intent(in) :: b !! positive finite second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, q, lp, lq
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            call beta_prime_tails(z, a, b, y, q, lp, lq)
        end if
    end function betaprime_cdf

    pure elemental function betaprime_sf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability
        real(dp), intent(in) :: a !! positive finite first shape parameter
        real(dp), intent(in) :: b !! positive finite second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, p, lp, lq
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            call beta_prime_tails(z, a, b, p, y, lp, lq)
        end if
    end function betaprime_sf

    pure elemental function betaprime_logcdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive finite first shape parameter
        real(dp), intent(in) :: b !! positive finite second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, p, q, lq
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            call beta_prime_tails(z, a, b, p, q, y, lq)
        end if
    end function betaprime_logcdf

    pure elemental function betaprime_logsf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive finite first shape parameter
        real(dp), intent(in) :: b !! positive finite second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, p, q, lp
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            call beta_prime_tails(z, a, b, p, q, lp, y)
        end if
    end function betaprime_logsf

    pure elemental function betaprime_ppf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive finite first shape parameter
        real(dp), intent(in) :: b !! positive finite second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, u, v, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            call beta_inverse_xy(a, b, p, 1.0_dp - p, u, v)
            if (v == 0.0_dp) then
                z = positive_infinity(p)
            else
                z = u / v
            end if
            x = affine_positive(mu, sigma, z, p)
        end if
    end function betaprime_ppf

    pure elemental function betaprime_isf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive finite first shape parameter
        real(dp), intent(in) :: b !! positive finite second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, u, v, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            call beta_inverse_xy(b, a, p, 1.0_dp - p, v, u)
            if (v == 0.0_dp) then
                z = positive_infinity(p)
            else
                z = u / v
            end if
            x = affine_positive(mu, sigma, z, p)
        end if
    end function betaprime_isf

    pure elemental function valid_shape(a) result(ok)
        real(dp), intent(in) :: a !! candidate shape parameter
        logical :: ok
        ok = ieee_is_finite(a) .and. a > 0.0_dp
    end function valid_shape

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

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to form infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z) .or. z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive
end module scifort_betaprime
