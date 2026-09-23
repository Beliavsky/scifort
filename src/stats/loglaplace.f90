! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Log-Laplace distribution with shape c > 0, location loc, and scale > 0.
! For z = (x - loc) / scale > 0,
! f(x) = c * exp(-c*abs(log(z))) / (2*z*scale).
! This matches scipy.stats.loglaplace.

module scifort_loglaplace
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: loglaplace_cdf
    public :: loglaplace_isf
    public :: loglaplace_logcdf
    public :: loglaplace_logpdf
    public :: loglaplace_logsf
    public :: loglaplace_pdf
    public :: loglaplace_ppf
    public :: loglaplace_sf

contains

    pure elemental function loglaplace_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
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
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp .or. .not. ieee_is_finite(z)) then
                y = 0.0_dp
            else if (z == 0.0_dp) then
                if (c < 1.0_dp) then
                    y = positive_infinity(x)
                else if (c == 1.0_dp) then
                    y = 0.5_dp * c / sigma
                else
                    y = 0.0_dp
                end if
            else
                y = exp(log(c) - scifort_log_two - log(z) - &
                    c * abs(log(z)) - log(sigma))
            end if
        end if
    end function loglaplace_pdf

    pure elemental function loglaplace_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
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
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else if (z == 0.0_dp) then
                if (c < 1.0_dp) then
                    y = positive_infinity(x)
                else if (c == 1.0_dp) then
                    y = log(c) - scifort_log_two - log(sigma)
                else
                    y = negative_infinity(x)
                end if
            else
                y = log(c) - scifort_log_two - log(z) - &
                    c * abs(log(z)) - log(sigma)
            end if
        end if
    end function loglaplace_logpdf

    pure elemental function loglaplace_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logtail
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            if (z < 1.0_dp) then
                y = 0.5_dp * exp(c * log(z))
            else
                logtail = -scifort_log_two - c * log(z)
                y = -expm1_safe(logtail)
            end if
        end if
    end function loglaplace_cdf

    pure elemental function loglaplace_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: loghead
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z < 1.0_dp) then
                loghead = -scifort_log_two + c * log(z)
                y = -expm1_safe(loghead)
            else
                y = 0.5_dp * exp(-c * log(z))
            end if
        end if
    end function loglaplace_sf

    pure elemental function loglaplace_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logtail
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z < 1.0_dp) then
                y = -scifort_log_two + c * log(z)
            else
                logtail = -scifort_log_two - c * log(z)
                y = log1p_safe(-exp(logtail))
            end if
        end if
    end function loglaplace_logcdf

    pure elemental function loglaplace_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: loghead
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z < 1.0_dp) then
                loghead = -scifort_log_two + c * log(z)
                y = log1p_safe(-exp(loghead))
            else
                y = -scifort_log_two - c * log(z)
            end if
        end if
    end function loglaplace_logsf

    pure elemental function loglaplace_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: exponent
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else if (p <= 0.5_dp) then
            exponent = (scifort_log_two + log(p)) / c
            z = exp(exponent)
            x = mu + sigma * z
        else
            exponent = -(scifort_log_two + log1p_safe(-p)) / c
            if (exponent > log(huge(1.0_dp))) then
                x = positive_infinity(p)
            else
                z = exp(exponent)
                x = affine_positive(mu, sigma, z, p)
            end if
        end if
    end function loglaplace_ppf

    pure elemental function loglaplace_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: exponent
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p <= 0.5_dp) then
            exponent = -(scifort_log_two + log(p)) / c
            if (exponent > log(huge(1.0_dp))) then
                x = positive_infinity(p)
            else
                z = exp(exponent)
                x = affine_positive(mu, sigma, z, p)
            end if
        else
            exponent = (scifort_log_two + log1p_safe(-p)) / c
            z = exp(exponent)
            x = mu + sigma * z
        end if
    end function loglaplace_isf

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to construct positive infinity if needed
        real(dp) :: x

        if (.not. ieee_is_finite(z) .or. z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive

    pure elemental logical function valid_shape(c) result(valid)
        real(dp), intent(in) :: c !! shape parameter to check

        valid = ieee_is_finite(c) .and. c > 0.0_dp
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

end module scifort_loglaplace
