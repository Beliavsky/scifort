! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_uniform
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: uniform_cdf
    public :: uniform_isf
    public :: uniform_logcdf
    public :: uniform_logpdf
    public :: uniform_logsf
    public :: uniform_pdf
    public :: uniform_ppf
    public :: uniform_sf

contains

    pure elemental function uniform_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
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
        if (ieee_is_nan(z)) then
            y = quiet_nan(x)
        else if (z < 0.0_dp .or. z > 1.0_dp) then
            y = 0.0_dp
        else
            y = 1.0_dp / sigma
        end if
    end function uniform_pdf

    pure elemental function uniform_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
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
        if (ieee_is_nan(z)) then
            y = quiet_nan(x)
        else if (z < 0.0_dp .or. z > 1.0_dp) then
            y = negative_infinity(x)
        else
            y = -log(sigma)
        end if
    end function uniform_logpdf

    pure elemental function uniform_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
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
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else if (z >= 1.0_dp) then
            y = 1.0_dp
        else
            y = z
        end if
    end function uniform_cdf

    pure elemental function uniform_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
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
        if (z <= 0.0_dp) then
            y = 1.0_dp
        else if (z >= 1.0_dp) then
            y = 0.0_dp
        else
            y = 1.0_dp - z
        end if
    end function uniform_sf

    pure elemental function uniform_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: cdf
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        cdf = uniform_cdf(x, mu, sigma)
        if (cdf <= 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(cdf)
        end if
    end function uniform_logcdf

    pure elemental function uniform_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sf
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        sf = uniform_sf(x, mu, sigma)
        if (sf <= 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(sf)
        end if
    end function uniform_logsf

    pure elemental function uniform_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            x = mu + sigma * p
        end if
    end function uniform_ppf

    pure elemental function uniform_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            x = mu + sigma * (1.0_dp - p)
        end if
    end function uniform_isf

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

end module scifort_uniform
