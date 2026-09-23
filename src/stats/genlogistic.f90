! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized logistic distribution with shape c > 0.
! In standardized coordinates z = (x-loc)/scale,
! F(z) = (1 + exp(-z))**(-c). This matches scipy.stats.genlogistic.

module scifort_genlogistic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, log1pexp, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: genlogistic_cdf
    public :: genlogistic_isf
    public :: genlogistic_logcdf
    public :: genlogistic_logpdf
    public :: genlogistic_logsf
    public :: genlogistic_pdf
    public :: genlogistic_ppf
    public :: genlogistic_sf

contains

    pure elemental function genlogistic_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            y = log(c) - z - (c + 1.0_dp) * log1pexp(-z) - log(sigma)
        end if
    end function genlogistic_logpdf

    pure elemental function genlogistic_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy
        logy = genlogistic_logpdf(x, c, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function genlogistic_pdf

    pure elemental function genlogistic_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = negative_infinity(x)
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            y = -c * log1pexp(-z)
        end if
    end function genlogistic_logcdf

    pure elemental function genlogistic_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lcdf
        lcdf = genlogistic_logcdf(x, c, loc, scale)
        if (ieee_is_nan(lcdf)) then
            y = lcdf
        else if (lcdf == 0.0_dp) then
            y = negative_infinity(x)
        else if (lcdf == negative_infinity(lcdf)) then
            y = 0.0_dp
        else
            y = log(-expm1_safe(lcdf))
        end if
    end function genlogistic_logsf

    pure elemental function genlogistic_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy
        logy = genlogistic_logcdf(x, c, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function genlogistic_cdf

    pure elemental function genlogistic_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lcdf
        lcdf = genlogistic_logcdf(x, c, loc, scale)
        if (ieee_is_nan(lcdf)) then
            y = lcdf
        else if (lcdf == 0.0_dp) then
            y = 0.0_dp
        else if (lcdf == negative_infinity(lcdf)) then
            y = 1.0_dp
        else
            y = -expm1_safe(lcdf)
        end if
    end function genlogistic_sf

    pure elemental function genlogistic_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            t = -log(p) / c
            z = neg_log_expm1_positive(t)
            x = mu + sigma * z
        end if
    end function genlogistic_ppf

    pure elemental function genlogistic_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            t = -log1p_safe(-p) / c
            z = neg_log_expm1_positive(t)
            x = mu + sigma * z
        end if
    end function genlogistic_isf

    pure elemental function neg_log_expm1_positive(t) result(z)
        real(dp), intent(in) :: t !! strictly positive transformed probability
        real(dp) :: z
        if (t > 0.5_dp) then
            z = -t - log1p_safe(-exp(-t))
        else
            z = -log(expm1_safe(t))
        end if
    end function neg_log_expm1_positive

    pure elemental logical function valid_shape(c) result(valid)
        real(dp), intent(in) :: c !! shape parameter to check
        valid = ieee_is_finite(c) .and. c > 0.0_dp
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

end module scifort_genlogistic
