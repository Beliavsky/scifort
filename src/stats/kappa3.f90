! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Three-parameter kappa distribution matching scipy.stats.kappa3.
module scifort_kappa3
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, log1pexp, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: kappa3_cdf, kappa3_isf, kappa3_logcdf, kappa3_logpdf
    public :: kappa3_logsf, kappa3_pdf, kappa3_ppf, kappa3_sf

contains

    pure elemental function kappa3_logpdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, logden, logz, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z == 0.0_dp) then
                y = -log(a) / a - log(sigma)
            else
                logz = log(z)
                logden = log(a) + log1pexp(a * logz - log(a))
                y = log(a) - (1.0_dp + 1.0_dp / a) * logden - log(sigma)
            end if
        end if
    end function kappa3_logpdf

    pure elemental function kappa3_pdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = kappa3_logpdf(x, a, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function kappa3_pdf

    pure elemental function kappa3_logcdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, logz, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            logz = log(z)
            y = -log1pexp(log(a) - a * logz) / a
        end if
    end function kappa3_logcdf

    pure elemental function kappa3_cdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = kappa3_logcdf(x, a, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function kappa3_cdf

    pure elemental function kappa3_sf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = kappa3_logcdf(x, a, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 1.0_dp
        else if (lc == 0.0_dp) then
            y = 0.0_dp
        else
            y = -expm1_safe(lc)
        end if
    end function kappa3_sf

    pure elemental function kappa3_logsf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = kappa3_logcdf(x, a, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 0.0_dp
        else if (lc == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(-expm1_safe(lc))
        end if
    end function kappa3_logsf

    pure elemental function kappa3_ppf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            z = quantile_from_logprob(log(p), a, p)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function kappa3_ppf

    pure elemental function kappa3_isf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            z = quantile_from_logprob(log1p_safe(-p), a, p)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function kappa3_isf

    pure elemental function quantile_from_logprob(logp, a, seed) result(z)
        real(dp), intent(in) :: logp !! logarithm of a probability strictly between zero and one
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: z, t, logden, logz
        t = -a * logp
        if (t > log(huge(1.0_dp))) then
            z = 0.0_dp
        else
            logden = log(expm1_safe(t))
            logz = (log(a) - logden) / a
            if (logz > log(huge(1.0_dp))) then
                z = positive_infinity(seed)
            else if (logz < log(tiny(1.0_dp))) then
                z = 0.0_dp
            else
                z = exp(logz)
            end if
        end if
    end function quantile_from_logprob

    pure elemental logical function valid_shape(a) result(ok)
        real(dp), intent(in) :: a !! candidate kappa shape parameter
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
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z) .or. z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive

end module scifort_kappa3
