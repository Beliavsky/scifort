! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Triangular distribution on [loc, loc + scale], with standardized mode c in [0, 1].

module scifort_triang
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: triang_cdf
    public :: triang_isf
    public :: triang_logcdf
    public :: triang_logpdf
    public :: triang_logsf
    public :: triang_pdf
    public :: triang_ppf
    public :: triang_sf

contains

    pure elemental function triang_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
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
        else if (c == 0.0_dp) then
            y = 2.0_dp * (1.0_dp - z) / sigma
        else if (c == 1.0_dp) then
            y = 2.0_dp * z / sigma
        else if (z < c) then
            y = 2.0_dp * z / (c * sigma)
        else
            y = 2.0_dp * (1.0_dp - z) / ((1.0_dp - c) * sigma)
        end if
    end function triang_pdf

    pure elemental function triang_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: density

        density = triang_pdf(x, c, loc, scale)
        if (ieee_is_nan(density)) then
            y = density
        else if (density <= 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(density)
        end if
    end function triang_logpdf

    pure elemental function triang_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
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
        else if (c == 0.0_dp) then
            y = z * (2.0_dp - z)
        else if (c == 1.0_dp .or. z <= c) then
            y = z * z / c
        else
            y = 1.0_dp - (1.0_dp - z) ** 2 / (1.0_dp - c)
        end if
    end function triang_cdf

    pure elemental function triang_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
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
        else if (c == 1.0_dp) then
            y = (1.0_dp - z) * (1.0_dp + z)
        else if (c == 0.0_dp .or. z >= c) then
            y = (1.0_dp - z) ** 2 / (1.0_dp - c)
        else
            y = 1.0_dp - z * z / c
        end if
    end function triang_sf

    pure elemental function triang_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
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
        else if (c == 0.0_dp) then
            y = log(z) + log(2.0_dp - z)
        else if (c == 1.0_dp .or. z <= c) then
            y = 2.0_dp * log(z) - log(c)
        else
            y = log1p_safe(-(1.0_dp - z) ** 2 / (1.0_dp - c))
        end if
    end function triang_logcdf

    pure elemental function triang_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
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
        else if (c == 1.0_dp) then
            y = log1p_safe(-z) + log1p_safe(z)
        else if (c == 0.0_dp .or. z >= c) then
            y = 2.0_dp * log1p_safe(-z) - log(1.0_dp - c)
        else
            y = log1p_safe(-z * z / c)
        end if
    end function triang_logsf

    pure elemental function triang_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: root
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = mu + sigma
        else if (p < c) then
            z = sqrt(p * c)
            x = mu + sigma * z
        else
            root = sqrt((1.0_dp - p) * (1.0_dp - c))
            z = (p + c - p * c) / (1.0_dp + root)
            x = mu + sigma * z
        end if
    end function triang_ppf

    pure elemental function triang_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu + sigma
        else if (p == 1.0_dp) then
            x = mu
        else if (p > 1.0_dp - c) then
            z = sqrt(c * (1.0_dp - p))
            x = mu + sigma * z
        else
            z = 1.0_dp - sqrt(p * (1.0_dp - c))
            x = mu + sigma * z
        end if
    end function triang_isf

    pure elemental logical function valid_shape(c) result(valid)
        real(dp), intent(in) :: c !! candidate standardized mode

        valid = ieee_is_finite(c) .and. c >= 0.0_dp .and. c <= 1.0_dp
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

end module scifort_triang
