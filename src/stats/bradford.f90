! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Bradford distribution with shape c > 0 on [loc, loc + scale].
! In standardized coordinates z = (x - loc) / scale in [0, 1],
! f(z) = c / ((1 + c*z) * log(1 + c)). This matches scipy.stats.bradford.

module scifort_bradford
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: bradford_cdf
    public :: bradford_isf
    public :: bradford_logcdf
    public :: bradford_logpdf
    public :: bradford_logsf
    public :: bradford_pdf
    public :: bradford_ppf
    public :: bradford_sf

contains

    pure elemental function bradford_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
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
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp .or. z > 1.0_dp) then
                y = 0.0_dp
            else
                y = c / ((1.0_dp + c * z) * log1p_safe(c) * sigma)
            end if
        end if
    end function bradford_pdf

    pure elemental function bradford_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
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
            if (z < 0.0_dp .or. z > 1.0_dp) then
                y = negative_infinity(x)
            else
                y = log(c) - log1p_safe(c * z) - log(log1p_safe(c)) - log(sigma)
            end if
        end if
    end function bradford_logpdf

    pure elemental function bradford_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y

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
            if (z >= 1.0_dp) then
                y = 1.0_dp
            else
                y = log1p_safe(c * z) / log1p_safe(c)
            end if
        end if
    end function bradford_cdf

    pure elemental function bradford_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: ratio
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
            if (z >= 1.0_dp) then
                y = 0.0_dp
            else
                ratio = c * (1.0_dp - z) / (1.0_dp + c * z)
                y = log1p_safe(ratio) / log1p_safe(c)
            end if
        end if
    end function bradford_sf

    pure elemental function bradford_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y

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
            if (z >= 1.0_dp) then
                y = 0.0_dp
            else
                y = log(log1p_safe(c * z)) - log(log1p_safe(c))
            end if
        end if
    end function bradford_logcdf

    pure elemental function bradford_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: ratio
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
            if (z >= 1.0_dp) then
                y = negative_infinity(x)
            else
                ratio = c * (1.0_dp - z) / (1.0_dp + c * z)
                y = log(log1p_safe(ratio)) - log(log1p_safe(c))
            end if
        end if
    end function bradford_logsf

    pure elemental function bradford_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
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
            x = mu
        else if (p == 1.0_dp) then
            x = mu + sigma
        else
            z = expm1_safe(p * log1p_safe(c)) / c
            x = mu + sigma * z
        end if
    end function bradford_ppf

    pure elemental function bradford_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x

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
            x = mu + sigma
        else
            z = 1.0_dp + (1.0_dp + c) * expm1_safe(-p * log1p_safe(c)) / c
            x = mu + sigma * z
        end if
    end function bradford_isf

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

end module scifort_bradford
