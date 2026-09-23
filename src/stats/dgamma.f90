! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Double-gamma (reflected gamma) distribution. In standardized coordinates
! z=(x-loc)/scale, f(z)=|z|**(a-1)*exp(-|z|)/(2*Gamma(a)), a>0.
! This matches scipy.stats.dgamma.

module scifort_dgamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv, log_gammainc, log_gammaincc
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_special_elementary, only : xlogy
    implicit none
    private

    public :: dgamma_cdf
    public :: dgamma_isf
    public :: dgamma_logcdf
    public :: dgamma_logpdf
    public :: dgamma_logsf
    public :: dgamma_pdf
    public :: dgamma_ppf
    public :: dgamma_sf

contains

    pure elemental function dgamma_pdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        y = exp(dgamma_logpdf(x, a, loc, scale))
    end function dgamma_pdf

    pure elemental function dgamma_logpdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: r
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
            return
        end if
        z = (x - mu) / sigma
        r = abs(z)
        y = xlogy(a - 1.0_dp, r) - r - log(2.0_dp) - log_gamma(a) - log(sigma)
    end function dgamma_logpdf

    pure elemental function dgamma_cdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = 0.0_dp
            else
                y = 1.0_dp
            end if
        else
            z = (x - mu) / sigma
            if (z > 0.0_dp) then
                y = 0.5_dp + 0.5_dp * gammainc(a, z)
            else
                y = 0.5_dp * gammaincc(a, -z)
            end if
        end if
    end function dgamma_cdf

    pure elemental function dgamma_sf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = 1.0_dp
            else
                y = 0.0_dp
            end if
        else
            z = (x - mu) / sigma
            if (z > 0.0_dp) then
                y = 0.5_dp * gammaincc(a, z)
            else
                y = 0.5_dp + 0.5_dp * gammainc(a, -z)
            end if
        end if
    end function dgamma_sf

    pure elemental function dgamma_logcdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: loghalf
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = negative_infinity(x)
            else
                y = 0.0_dp
            end if
        else
            loghalf = -log(2.0_dp)
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = loghalf + log_gammaincc(a, -z)
            else
                y = log1p_safe(-0.5_dp * gammaincc(a, z))
            end if
        end if
    end function dgamma_logcdf

    pure elemental function dgamma_logsf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: loghalf
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = 0.0_dp
            else
                y = negative_infinity(x)
            end if
        else
            loghalf = -log(2.0_dp)
            z = (x - mu) / sigma
            if (z >= 0.0_dp) then
                y = loghalf + log_gammaincc(a, z)
            else
                y = log1p_safe(-0.5_dp * gammaincc(a, -z))
            end if
        end if
    end function dgamma_logsf

    pure elemental function dgamma_ppf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = negative_infinity(p)
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
        else if (p > 0.5_dp) then
            z = gammaincinv(a, 2.0_dp * p - 1.0_dp)
            x = mu + sigma * z
        else if (p < 0.5_dp) then
            z = gammainccinv(a, 2.0_dp * p)
            x = mu - sigma * z
        else
            x = mu
        end if
    end function dgamma_ppf

    pure elemental function dgamma_isf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = positive_infinity(p)
        else if (p >= 1.0_dp) then
            x = negative_infinity(p)
        else if (p > 0.5_dp) then
            z = gammaincinv(a, 2.0_dp * p - 1.0_dp)
            x = mu - sigma * z
        else if (p < 0.5_dp) then
            z = gammainccinv(a, 2.0_dp * p)
            x = mu + sigma * z
        else
            x = mu
        end if
    end function dgamma_isf

    pure elemental logical function valid_shape(a) result(valid)
        real(dp), intent(in) :: a !! shape parameter to validate
        valid = ieee_is_finite(a) .and. a > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location argument
        real(dp), intent(in), optional :: scale !! optional scale argument
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_dgamma
