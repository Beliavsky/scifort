! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Log-uniform (reciprocal) distribution with 0<a<b. In standardized
! coordinates z=(x-loc)/scale, f(z)=1/(z*log(b/a)) on [a,b].
! This matches scipy.stats.loguniform / scipy.stats.reciprocal.

module scifort_loguniform
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: loguniform_cdf
    public :: loguniform_isf
    public :: loguniform_logcdf
    public :: loguniform_logpdf
    public :: loguniform_logsf
    public :: loguniform_pdf
    public :: loguniform_ppf
    public :: loguniform_sf

contains

    pure elemental function loguniform_pdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = exp(loguniform_logpdf(x, a, b, loc, scale))
    end function loguniform_pdf

    pure elemental function loguniform_logpdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: width
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (z < a .or. z > b .or. .not. ieee_is_finite(z)) then
            y = negative_infinity(x)
        else
            width = log(b) - log(a)
            y = -log(z) - log(width) - log(sigma)
        end if
    end function loguniform_logpdf

    pure elemental function loguniform_cdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: width
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= a) then
                y = 0.0_dp
            else if (z >= b) then
                y = 1.0_dp
            else
                width = log(b) - log(a)
                y = (log(z) - log(a)) / width
            end if
        end if
    end function loguniform_cdf

    pure elemental function loguniform_sf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: width
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= a) then
                y = 1.0_dp
            else if (z >= b) then
                y = 0.0_dp
            else
                width = log(b) - log(a)
                y = (log(b) - log(z)) / width
            end if
        end if
    end function loguniform_sf

    pure elemental function loguniform_logcdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: cdf_value

        cdf_value = loguniform_cdf(x, a, b, loc, scale)
        if (ieee_is_nan(cdf_value)) then
            y = cdf_value
        else if (cdf_value <= 0.0_dp) then
            y = negative_infinity(cdf_value)
        else
            y = log(cdf_value)
        end if
    end function loguniform_logcdf

    pure elemental function loguniform_logsf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: sf_value

        sf_value = loguniform_sf(x, a, b, loc, scale)
        if (ieee_is_nan(sf_value)) then
            y = sf_value
        else if (sf_value <= 0.0_dp) then
            y = negative_infinity(sf_value)
        else
            y = log(sf_value)
        end if
    end function loguniform_logsf

    pure elemental function loguniform_ppf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            z = exp(log(a) + p * (log(b) - log(a)))
            x = mu + sigma * z
        end if
    end function loguniform_ppf

    pure elemental function loguniform_isf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            z = exp(log(b) - p * (log(b) - log(a)))
            x = mu + sigma * z
        end if
    end function loguniform_isf

    pure elemental logical function valid_shapes(a, b) result(valid)
        real(dp), intent(in) :: a !! lower shape parameter to validate
        real(dp), intent(in) :: b !! upper shape parameter to validate
        valid = ieee_is_finite(a) .and. ieee_is_finite(b) .and. a > 0.0_dp .and. b > a
    end function valid_shapes

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

end module scifort_loguniform
