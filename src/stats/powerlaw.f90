! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Power-function distribution with shape a > 0. In standardized coordinates
! z = (x - loc) / scale on [0, 1], f(z) = a*z**(a - 1).

module scifort_powerlaw
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: powerlaw_cdf
    public :: powerlaw_isf
    public :: powerlaw_logcdf
    public :: powerlaw_logpdf
    public :: powerlaw_logsf
    public :: powerlaw_pdf
    public :: powerlaw_ppf
    public :: powerlaw_sf

contains

    pure elemental function powerlaw_pdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: logp
        real(dp) :: mu
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

        z = (x - mu) / sigma
        if (z < 0.0_dp .or. z > 1.0_dp) then
            y = 0.0_dp
        else if (z == 0.0_dp) then
            if (a == 1.0_dp) then
                y = 1.0_dp / sigma
            else
                y = 0.0_dp
            end if
        else
            logp = log(a) - log(sigma) + (a - 1.0_dp) * log(z)
            if (logp > log(huge(1.0_dp))) then
                y = positive_infinity(logp)
            else
                y = exp(logp)
            end if
        end if
    end function powerlaw_pdf

    pure elemental function powerlaw_logpdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
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

        z = (x - mu) / sigma
        if (z < 0.0_dp .or. z > 1.0_dp) then
            y = negative_infinity(x)
        else if (z == 0.0_dp) then
            if (a == 1.0_dp) then
                y = -log(sigma)
            else
                y = negative_infinity(x)
            end if
        else
            y = log(a) - log(sigma) + (a - 1.0_dp) * log(z)
        end if
    end function powerlaw_logpdf

    pure elemental function powerlaw_cdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
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

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else if (z >= 1.0_dp) then
            y = 1.0_dp
        else
            y = exp(a * log(z))
        end if
    end function powerlaw_cdf

    pure elemental function powerlaw_sf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
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

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 1.0_dp
        else if (z >= 1.0_dp) then
            y = 0.0_dp
        else
            t = a * log(z)
            y = -expm1_safe(t)
        end if
    end function powerlaw_sf

    pure elemental function powerlaw_logcdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
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

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = negative_infinity(x)
        else if (z >= 1.0_dp) then
            y = 0.0_dp
        else
            y = a * log(z)
        end if
    end function powerlaw_logcdf

    pure elemental function powerlaw_logsf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
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

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else if (z >= 1.0_dp) then
            y = negative_infinity(x)
        else
            t = a * log(z)
            y = log(-expm1_safe(t))
        end if
    end function powerlaw_logsf

    pure elemental function powerlaw_ppf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = mu + sigma
        else
            x = mu + sigma * exp(log(p) / a)
        end if
    end function powerlaw_ppf

    pure elemental function powerlaw_isf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu + sigma
        else if (p == 1.0_dp) then
            x = mu
        else
            x = mu + sigma * exp(log1p_safe(-p) / a)
        end if
    end function powerlaw_isf

    pure elemental logical function valid_shape(a) result(valid)
        real(dp), intent(in) :: a !! candidate shape parameter

        valid = ieee_is_finite(a) .and. a > 0.0_dp
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

end module scifort_powerlaw
