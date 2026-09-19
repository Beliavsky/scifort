! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Rayleigh distribution with location loc and scale > 0:
! with z = (x - loc) / scale >= 0, f(x) = (z / scale) * exp(-z**2 / 2).
! This matches the parameterization of scipy.stats.rayleigh.

module scifort_rayleigh
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: rayleigh_cdf
    public :: rayleigh_isf
    public :: rayleigh_logcdf
    public :: rayleigh_logpdf
    public :: rayleigh_logsf
    public :: rayleigh_pdf
    public :: rayleigh_ppf
    public :: rayleigh_sf

contains

    pure elemental function rayleigh_pdf(x, loc, scale) result(y)
        ! Evaluates the probability density function of the Rayleigh distribution.
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else
            t = 0.5_dp * z * z
            if (t > 709.0_dp) then
                y = 0.0_dp
            else
                y = (z / sigma) * exp(-t)
            end if
        end if
    end function rayleigh_pdf

    pure elemental function rayleigh_logpdf(x, loc, scale) result(y)
        ! Evaluates the log probability density function of the Rayleigh distribution.
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = negative_infinity(x)
        else
            t = 0.5_dp * z * z
            y = log(z) - log(sigma) - t
        end if
    end function rayleigh_logpdf

    pure elemental function rayleigh_cdf(x, loc, scale) result(y)
        ! Evaluates the cumulative distribution function of the Rayleigh distribution.
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else
            t = 0.5_dp * z * z
            if (t > 709.0_dp) then
                y = 1.0_dp
            else
                y = -expm1_safe(-t)
            end if
        end if
    end function rayleigh_cdf

    pure elemental function rayleigh_sf(x, loc, scale) result(y)
        ! Evaluates the survival function P(X > x) of the Rayleigh distribution.
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 1.0_dp
        else
            t = 0.5_dp * z * z
            if (t > 709.0_dp) then
                y = 0.0_dp
            else
                y = exp(-t)
            end if
        end if
    end function rayleigh_sf

    pure elemental function rayleigh_logcdf(x, loc, scale) result(y)
        ! Evaluates the logarithm of the cumulative distribution function of the Rayleigh distribution.
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = negative_infinity(x)
        else
            t = 0.5_dp * z * z
            if (t < scifort_log_two) then
                y = log(-expm1_safe(-t))
            else if (t > 709.0_dp) then
                y = 0.0_dp
            else
                y = log1p_safe(-exp(-t))
            end if
        end if
    end function rayleigh_logcdf

    pure elemental function rayleigh_logsf(x, loc, scale) result(y)
        ! Evaluates the logarithm of the survival function of the Rayleigh distribution.
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else
            y = -0.5_dp * z * z
        end if
    end function rayleigh_logsf

    pure elemental function rayleigh_ppf(p, loc, scale) result(x)
        ! Evaluates the percent point function (quantile) of the Rayleigh distribution.
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: v
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
            return
        end if

        if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            v = -log1p_safe(-p)
            z = sqrt(2.0_dp * v)
            x = mu + sigma * z
        end if
    end function rayleigh_ppf

    pure elemental function rayleigh_isf(p, loc, scale) result(x)
        ! Evaluates the inverse survival function of the Rayleigh distribution.
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: v
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
            return
        end if

        if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = mu
        else
            v = -log(p)
            z = sqrt(2.0_dp * v)
            x = mu + sigma * z
        end if
    end function rayleigh_isf

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        ! Extracts location and scale arguments, substituting default values if absent.
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when loc is absent
        real(dp), intent(out) :: sigma !! scale, 1 when scale is absent

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_rayleigh
