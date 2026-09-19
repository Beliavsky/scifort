! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_logistic
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, log1pexp, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: logistic_cdf
    public :: logistic_isf
    public :: logistic_logcdf
    public :: logistic_logpdf
    public :: logistic_logsf
    public :: logistic_pdf
    public :: logistic_ppf
    public :: logistic_sf

contains

    pure elemental function logistic_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        y = exp(logistic_logpdf(x, loc, scale))
    end function logistic_pdf

    pure elemental function logistic_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: a
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        a = abs(z)
        y = -log(sigma) - a - 2.0_dp * log1p_safe(exp(-a))
    end function logistic_logpdf

    pure elemental function logistic_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: e
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z >= 0.0_dp) then
            y = 1.0_dp / (1.0_dp + exp(-z))
        else
            e = exp(z)
            y = e / (1.0_dp + e)
        end if
    end function logistic_cdf

    pure elemental function logistic_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: e
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z >= 0.0_dp) then
            e = exp(-z)
            y = e / (1.0_dp + e)
        else
            y = 1.0_dp / (1.0_dp + exp(z))
        end if
    end function logistic_sf

    pure elemental function logistic_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        y = -log1pexp(-z)
    end function logistic_logcdf

    pure elemental function logistic_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        y = -log1pexp(z)
    end function logistic_logsf

    pure elemental function logistic_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = negative_infinity(p)
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
        else
            x = mu + sigma * (log(p) - log1p_safe(-p))
        end if
    end function logistic_ppf

    pure elemental function logistic_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = positive_infinity(p)
        else if (p >= 1.0_dp) then
            x = negative_infinity(p)
        else
            x = mu + sigma * (log1p_safe(-p) - log(p))
        end if
    end function logistic_isf

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

end module scifort_logistic
