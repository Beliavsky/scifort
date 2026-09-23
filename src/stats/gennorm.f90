! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized normal distribution with shape beta > 0.
! In standardized coordinates z = (x-loc)/scale,
! f(z) = beta * exp(-abs(z)**beta) / (2*Gamma(1/beta)).
! This matches scipy.stats.gennorm.

module scifort_gennorm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_incomplete_gamma, only : gammaincc, gammainccinv, log_gammaincc
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: gennorm_cdf
    public :: gennorm_isf
    public :: gennorm_logcdf
    public :: gennorm_logpdf
    public :: gennorm_logsf
    public :: gennorm_pdf
    public :: gennorm_ppf
    public :: gennorm_sf

contains

    pure elemental function gennorm_logpdf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(beta))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = abs((x - mu) / sigma)
            t = abs_power(z, beta)
            if (.not. ieee_is_finite(t)) then
                y = negative_infinity(x)
            else
                y = log(beta) - scifort_log_two - log_gamma(1.0_dp / beta) - t - log(sigma)
            end if
        end if
    end function gennorm_logpdf

    pure elemental function gennorm_pdf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy
        logy = gennorm_logpdf(x, beta, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function gennorm_pdf

    pure elemental function gennorm_logcdf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: a
        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(beta))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = negative_infinity(x)
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z == 0.0_dp) then
                y = -scifort_log_two
            else
                a = 1.0_dp / beta
                t = abs_power(abs(z), beta)
                if (z < 0.0_dp) then
                    y = -scifort_log_two + log_gammaincc(a, t)
                else
                    q = gammaincc(a, t)
                    y = log1p_safe(-0.5_dp * q)
                end if
            end if
        end if
    end function gennorm_logcdf

    pure elemental function gennorm_logsf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: a
        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(beta))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 0.0_dp
        else if (x == positive_infinity(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z == 0.0_dp) then
                y = -scifort_log_two
            else
                a = 1.0_dp / beta
                t = abs_power(abs(z), beta)
                if (z > 0.0_dp) then
                    y = -scifort_log_two + log_gammaincc(a, t)
                else
                    q = gammaincc(a, t)
                    y = log1p_safe(-0.5_dp * q)
                end if
            end if
        end if
    end function gennorm_logsf

    pure elemental function gennorm_cdf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy
        logy = gennorm_logcdf(x, beta, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function gennorm_cdf

    pure elemental function gennorm_sf(x, beta, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy
        logy = gennorm_logsf(x, beta, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function gennorm_sf

    pure elemental function gennorm_ppf(p, beta, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(beta))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else if (p == 0.5_dp) then
            x = mu
        else
            z = standardized_quantile(p, beta)
            if (.not. ieee_is_finite(z) .or. abs(z) > (huge(1.0_dp) - abs(mu)) / sigma) then
                if (z < 0.0_dp) then
                    x = negative_infinity(p)
                else
                    x = positive_infinity(p)
                end if
            else
                x = mu + sigma * z
            end if
        end if
    end function gennorm_ppf

    pure elemental function gennorm_isf(p, beta, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: beta !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(beta))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 0.5_dp) then
            x = mu
        else
            z = -standardized_quantile(p, beta)
            if (.not. ieee_is_finite(z) .or. abs(z) > (huge(1.0_dp) - abs(mu)) / sigma) then
                if (z < 0.0_dp) then
                    x = negative_infinity(p)
                else
                    x = positive_infinity(p)
                end if
            else
                x = mu + sigma * z
            end if
        end if
    end function gennorm_isf

    pure elemental function standardized_quantile(p, beta) result(z)
        real(dp), intent(in) :: p !! lower-tail probability strictly between zero and one
        real(dp), intent(in) :: beta !! positive shape parameter
        real(dp) :: z
        real(dp) :: q
        real(dp) :: t
        if (p < 0.5_dp) then
            q = 2.0_dp * p
            t = gammainccinv(1.0_dp / beta, q)
            z = -positive_power(t, 1.0_dp / beta)
        else
            q = 2.0_dp * exp(log1p_safe(-p))
            t = gammainccinv(1.0_dp / beta, q)
            z = positive_power(t, 1.0_dp / beta)
        end if
    end function standardized_quantile

    pure elemental function abs_power(z, beta) result(t)
        real(dp), intent(in) :: z !! nonnegative magnitude
        real(dp), intent(in) :: beta !! positive exponent
        real(dp) :: t
        real(dp) :: lt
        if (z == 0.0_dp) then
            t = 0.0_dp
        else
            lt = beta * log(z)
            if (lt > log(huge(1.0_dp))) then
                t = positive_infinity(z)
            else
                t = exp(lt)
            end if
        end if
    end function abs_power

    pure elemental function positive_power(t, exponent) result(z)
        real(dp), intent(in) :: t !! nonnegative base
        real(dp), intent(in) :: exponent !! positive exponent
        real(dp) :: z
        real(dp) :: lz
        if (t == 0.0_dp) then
            z = 0.0_dp
        else if (.not. ieee_is_finite(t)) then
            z = positive_infinity(t)
        else
            lz = exponent * log(t)
            if (lz > log(huge(1.0_dp))) then
                z = positive_infinity(t)
            else
                z = exp(lz)
            end if
        end if
    end function positive_power

    pure elemental logical function valid_shape(beta) result(valid)
        real(dp), intent(in) :: beta !! shape parameter to check
        valid = ieee_is_finite(beta) .and. beta > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location
        real(dp), intent(in), optional :: scale !! optional scale
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_gennorm
