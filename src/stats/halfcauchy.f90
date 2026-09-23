! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Half-Cauchy distribution with lower support endpoint loc and scale > 0.
! With z = (x - loc) / scale >= 0, f(x) = 2 / (pi*scale*(1 + z**2)).
! This matches scipy.stats.halfcauchy.

module scifort_halfcauchy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_log_two, scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: halfcauchy_cdf
    public :: halfcauchy_isf
    public :: halfcauchy_logcdf
    public :: halfcauchy_logpdf
    public :: halfcauchy_logsf
    public :: halfcauchy_pdf
    public :: halfcauchy_ppf
    public :: halfcauchy_sf

contains

    pure elemental function halfcauchy_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z > sqrt(huge(1.0_dp))) then
                y = 0.0_dp
            else
                y = 2.0_dp / (scifort_pi * sigma * (1.0_dp + z * z))
            end if
        end if
    end function halfcauchy_pdf

    pure elemental function halfcauchy_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: a
        real(dp) :: log_one_plus_z2
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
        if (x < mu) then
            y = negative_infinity(x)
            return
        end if

        z = (x - mu) / sigma
        a = abs(z)
        if (a <= sqrt(huge(1.0_dp))) then
            log_one_plus_z2 = log(1.0_dp + a * a)
        else
            log_one_plus_z2 = 2.0_dp * log(a)
        end if
        y = scifort_log_two - log(scifort_pi) - log(sigma) - log_one_plus_z2
    end function halfcauchy_logpdf

    pure elemental function halfcauchy_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z <= 1.0_dp) then
                y = 2.0_dp * atan(z) / scifort_pi
            else
                y = 1.0_dp - 2.0_dp * atan(1.0_dp / z) / scifort_pi
            end if
        end if
    end function halfcauchy_cdf

    pure elemental function halfcauchy_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            y = 2.0_dp * atan(1.0_dp / z) / scifort_pi
        end if
    end function halfcauchy_sf

    pure elemental function halfcauchy_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: cdf_value
        real(dp) :: mu
        real(dp) :: sf_value
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z <= 1.0_dp) then
                cdf_value = 2.0_dp * atan(z) / scifort_pi
                y = log(cdf_value)
            else
                sf_value = 2.0_dp * atan(1.0_dp / z) / scifort_pi
                y = log1p_safe(-sf_value)
            end if
        end if
    end function halfcauchy_logcdf

    pure elemental function halfcauchy_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: sf_value
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            sf_value = 2.0_dp * atan(1.0_dp / z) / scifort_pi
            if (sf_value > 0.0_dp) then
                y = log(sf_value)
            else
                y = negative_infinity(x)
            end if
        end if
    end function halfcauchy_logsf

    pure elemental function halfcauchy_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else if (p <= 0.5_dp) then
            z = tan(0.5_dp * scifort_pi * p)
            x = mu + sigma * z
        else
            z = tan(0.5_dp * scifort_pi * (1.0_dp - p))
            if (z == 0.0_dp .or. sigma > huge(1.0_dp) * z) then
                x = positive_infinity(p)
            else
                x = mu + sigma / z
            end if
        end if
    end function halfcauchy_ppf

    pure elemental function halfcauchy_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p < 0.5_dp) then
            z = tan(0.5_dp * scifort_pi * p)
            if (z == 0.0_dp .or. sigma > huge(1.0_dp) * z) then
                x = positive_infinity(p)
            else
                x = mu + sigma / z
            end if
        else
            z = tan(0.5_dp * scifort_pi * (1.0_dp - p))
            x = mu + sigma * z
        end if
    end function halfcauchy_isf

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

end module scifort_halfcauchy
