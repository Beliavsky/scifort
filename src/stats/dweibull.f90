! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Double-Weibull distribution with shape c > 0, location loc, and scale > 0.
! In standardized coordinates z = (x-loc)/scale,
! f(z) = c*abs(z)**(c-1)*exp(-abs(z)**c)/2. This matches scipy.stats.dweibull.

module scifort_dweibull
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: dweibull_cdf
    public :: dweibull_isf
    public :: dweibull_logcdf
    public :: dweibull_logpdf
    public :: dweibull_logsf
    public :: dweibull_pdf
    public :: dweibull_ppf
    public :: dweibull_sf

contains

    pure elemental function dweibull_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logy

        logy = dweibull_logpdf(x, c, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else if (logy == positive_infinity(logy)) then
            y = positive_infinity(x)
        else
            y = exp(logy)
        end if
    end function dweibull_pdf

    pure elemental function dweibull_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: log_power
        real(dp) :: logr
        real(dp) :: mu
        real(dp) :: power
        real(dp) :: r
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
            r = abs(z)
            if (r == 0.0_dp) then
                if (c < 1.0_dp) then
                    y = positive_infinity(x)
                else if (c == 1.0_dp) then
                    y = -scifort_log_two - log(sigma)
                else
                    y = negative_infinity(x)
                end if
            else
                logr = log(r)
                log_power = shape_times_log(c, logr, x)
                if (log_power > log(huge(1.0_dp))) then
                    y = negative_infinity(x)
                else
                    power = exp(log_power)
                    y = log(c) - scifort_log_two - log(sigma) + &
                        log_power - logr - power
                end if
            end if
        end if
    end function dweibull_logpdf

    pure elemental function dweibull_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: power
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
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
            power = abs_power(z, c, x)
            if (z < 0.0_dp) then
                if (power == positive_infinity(power)) then
                    y = 0.0_dp
                else
                    y = 0.5_dp * exp(-power)
                end if
            else
                if (power == positive_infinity(power)) then
                    y = 1.0_dp
                else
                    y = 1.0_dp - 0.5_dp * exp(-power)
                end if
            end if
        end if
    end function dweibull_cdf

    pure elemental function dweibull_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: power
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
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
            power = abs_power(z, c, x)
            if (z < 0.0_dp) then
                if (power == positive_infinity(power)) then
                    y = 1.0_dp
                else
                    y = 1.0_dp - 0.5_dp * exp(-power)
                end if
            else
                if (power == positive_infinity(power)) then
                    y = 0.0_dp
                else
                    y = 0.5_dp * exp(-power)
                end if
            end if
        end if
    end function dweibull_sf

    pure elemental function dweibull_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: power
        real(dp) :: sigma
        real(dp) :: tail
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
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
            z = (x - mu) / sigma
            power = abs_power(z, c, x)
            if (z < 0.0_dp) then
                if (power == positive_infinity(power)) then
                    y = negative_infinity(x)
                else
                    y = -scifort_log_two - power
                end if
            else
                if (power == positive_infinity(power)) then
                    y = 0.0_dp
                else
                    tail = 0.5_dp * exp(-power)
                    y = log1p_safe(-tail)
                end if
            end if
        end if
    end function dweibull_logcdf

    pure elemental function dweibull_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: power
        real(dp) :: sigma
        real(dp) :: tail
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
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
            z = (x - mu) / sigma
            power = abs_power(z, c, x)
            if (z < 0.0_dp) then
                if (power == positive_infinity(power)) then
                    y = 0.0_dp
                else
                    tail = 0.5_dp * exp(-power)
                    y = log1p_safe(-tail)
                end if
            else
                if (power == positive_infinity(power)) then
                    y = negative_infinity(x)
                else
                    y = -scifort_log_two - power
                end if
            end if
        end if
    end function dweibull_logsf

    pure elemental function dweibull_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma

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
            q = standardized_ppf(p, c)
            x = affine_signed(mu, sigma, q, p)
        end if
    end function dweibull_ppf

    pure elemental function dweibull_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma

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
            q = -standardized_ppf(p, c)
            x = affine_signed(mu, sigma, q, p)
        end if
    end function dweibull_isf

    pure elemental function standardized_ppf(p, c) result(q)
        real(dp), intent(in) :: p !! probability strictly between zero and one
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp) :: q

        real(dp) :: a
        real(dp) :: exponent
        real(dp) :: r

        if (p == 0.5_dp) then
            q = 0.0_dp
        else if (p < 0.5_dp) then
            a = -(scifort_log_two + log(p))
            exponent = log(a) / c
            if (exponent > log(huge(1.0_dp))) then
                q = negative_infinity(p)
            else
                r = exp(exponent)
                q = -r
            end if
        else
            a = -(scifort_log_two + log1p_safe(-p))
            exponent = log(a) / c
            if (exponent > log(huge(1.0_dp))) then
                q = positive_infinity(p)
            else
                r = exp(exponent)
                q = r
            end if
        end if
    end function standardized_ppf

    pure elemental function abs_power(z, c, seed) result(value)
        real(dp), intent(in) :: z !! finite standardized variate
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in) :: seed !! value used to construct infinity when needed
        real(dp) :: value

        real(dp) :: log_power
        real(dp) :: r

        r = abs(z)
        if (r == 0.0_dp) then
            value = 0.0_dp
        else
            log_power = shape_times_log(c, log(r), seed)
            if (log_power > log(huge(1.0_dp))) then
                value = positive_infinity(seed)
            else if (log_power < log(tiny(1.0_dp))) then
                value = 0.0_dp
            else
                value = exp(log_power)
            end if
        end if
    end function abs_power

    pure elemental function shape_times_log(c, logr, seed) result(value)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in) :: logr !! logarithm of positive absolute standardized variate
        real(dp), intent(in) :: seed !! value used to construct infinity when needed
        real(dp) :: value

        if (logr == 0.0_dp) then
            value = 0.0_dp
        else if (abs(logr) > huge(1.0_dp) / c) then
            if (logr > 0.0_dp) then
                value = positive_infinity(seed)
            else
                value = negative_infinity(seed)
            end if
        else
            value = c * logr
        end if
    end function shape_times_log

    pure elemental function affine_signed(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! standardized quantile, possibly infinite
        real(dp), intent(in) :: seed !! value used to construct infinity when needed
        real(dp) :: x

        if (.not. ieee_is_finite(z)) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else if (abs(z) > (huge(1.0_dp) - abs(mu)) / sigma) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else
            x = mu + sigma * z
        end if
    end function affine_signed

    pure elemental logical function valid_shape(c) result(valid)
        real(dp), intent(in) :: c !! shape parameter to validate

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

end module scifort_dweibull
