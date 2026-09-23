! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Truncated exponential distribution with shape b > 0 on
! [loc, loc + scale*b]. In standardized coordinates z = (x-loc)/scale,
! f(z) = exp(-z) / (1 - exp(-b)). This matches scipy.stats.truncexpon.

module scifort_truncexpon
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: truncexpon_cdf
    public :: truncexpon_isf
    public :: truncexpon_logcdf
    public :: truncexpon_logpdf
    public :: truncexpon_logsf
    public :: truncexpon_pdf
    public :: truncexpon_ppf
    public :: truncexpon_sf

contains

    pure elemental function truncexpon_pdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: b !! positive finite standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp .or. z > b) then
                y = 0.0_dp
            else
                y = exp(-z) / (-expm1_safe(-b) * sigma)
            end if
        end if
    end function truncexpon_pdf

    pure elemental function truncexpon_logpdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: b !! positive finite standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp .or. z > b) then
                y = negative_infinity(x)
            else
                y = -z - log(-expm1_safe(-b)) - log(sigma)
            end if
        end if
    end function truncexpon_logpdf

    pure elemental function truncexpon_cdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: b !! positive finite standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            if (z >= b) then
                y = 1.0_dp
            else
                y = expm1_safe(-z) / expm1_safe(-b)
            end if
        end if
    end function truncexpon_cdf

    pure elemental function truncexpon_sf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: b !! positive finite standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z >= b) then
                y = 0.0_dp
            else
                y = exp(-z) * (-expm1_safe(-(b - z))) / (-expm1_safe(-b))
            end if
        end if
    end function truncexpon_sf

    pure elemental function truncexpon_logcdf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: b !! positive finite standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z >= b) then
                y = 0.0_dp
            else
                y = log(-expm1_safe(-z)) - log(-expm1_safe(-b))
            end if
        end if
    end function truncexpon_logcdf

    pure elemental function truncexpon_logsf(x, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: b !! positive finite standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z >= b) then
                y = negative_infinity(x)
            else
                y = -z + log(-expm1_safe(-(b - z))) - log(-expm1_safe(-b))
            end if
        end if
    end function truncexpon_logsf

    pure elemental function truncexpon_ppf(p, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: b !! positive finite standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: denom
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = affine_nonnegative(mu, sigma, b, p)
        else
            denom = -expm1_safe(-b)
            z = -log1p_safe(-p * denom)
            x = affine_nonnegative(mu, sigma, z, p)
        end if
    end function truncexpon_ppf

    pure elemental function truncexpon_isf(p, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: b !! positive finite standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: denom
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: tail
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = affine_nonnegative(mu, sigma, b, p)
        else
            denom = -expm1_safe(-b)
            tail = exp(-b) + p * denom
            if (tail <= 0.0_dp) then
                z = b
            else
                z = -log(tail)
            end if
            x = affine_nonnegative(mu, sigma, min(z, b), p)
        end if
    end function truncexpon_isf

    pure elemental function affine_nonnegative(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to construct positive infinity if needed
        real(dp) :: x

        if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_nonnegative

    pure elemental logical function valid_shape(b) result(valid)
        real(dp), intent(in) :: b !! shape parameter to validate

        valid = ieee_is_finite(b) .and. b > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when loc is absent
        real(dp), intent(out) :: sigma !! scale, 1 when scale is absent

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_truncexpon
