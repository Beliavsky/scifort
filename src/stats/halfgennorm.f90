! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Upper half of the generalized normal distribution, matching scipy.stats.halfgennorm.
module scifort_halfgennorm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv, log_gammainc, log_gammaincc
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, &
        valid_loc_scale
    implicit none
    private

    public :: halfgennorm_cdf, halfgennorm_isf, halfgennorm_logcdf, halfgennorm_logpdf
    public :: halfgennorm_logsf, halfgennorm_pdf, halfgennorm_ppf, halfgennorm_sf

contains

    pure elemental function halfgennorm_logpdf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            t = positive_power(z, beta)
            if (.not. ieee_is_finite(t)) then
                y = negative_infinity(x)
            else
                y = log(beta) - log_gamma(1.0_dp / beta) - t - log(sigma)
            end if
        end if
    end function halfgennorm_logpdf

    pure elemental function halfgennorm_pdf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = halfgennorm_logpdf(x, beta, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function halfgennorm_pdf

    pure elemental function halfgennorm_cdf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            t = positive_power(z, beta)
            y = gammainc(1.0_dp / beta, t)
        end if
    end function halfgennorm_cdf

    pure elemental function halfgennorm_sf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            t = positive_power(z, beta)
            y = gammaincc(1.0_dp / beta, t)
        end if
    end function halfgennorm_sf

    pure elemental function halfgennorm_logcdf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            t = positive_power(z, beta)
            y = log_gammainc(1.0_dp / beta, t)
        end if
    end function halfgennorm_logcdf

    pure elemental function halfgennorm_logsf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            t = positive_power(z, beta)
            y = log_gammaincc(1.0_dp / beta, t)
        end if
    end function halfgennorm_logsf

    pure elemental function halfgennorm_ppf(p, beta, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            t = gammaincinv(1.0_dp / beta, p)
            z = positive_power(t, 1.0_dp / beta)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function halfgennorm_ppf

    pure elemental function halfgennorm_isf(p, beta, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            t = gammainccinv(1.0_dp / beta, p)
            z = positive_power(t, 1.0_dp / beta)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function halfgennorm_isf

    pure elemental function positive_power(x, p) result(y)
        real(dp), intent(in) :: x !! nonnegative base
        real(dp), intent(in) :: p !! positive finite exponent
        real(dp) :: y, logy
        if (x == 0.0_dp) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = positive_infinity(x)
        else
            logy = p * log(x)
            if (logy > log(huge(1.0_dp))) then
                y = positive_infinity(x)
            else if (logy < log(tiny(1.0_dp))) then
                y = 0.0_dp
            else
                y = exp(logy)
            end if
        end if
    end function positive_power

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive scale parameter
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z) .or. z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive

    pure elemental logical function valid_shape(beta) result(ok)
        real(dp), intent(in) :: beta !! candidate generalized-normal shape parameter
        ok = ieee_is_finite(beta) .and. beta > 0.0_dp
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

end module scifort_halfgennorm
