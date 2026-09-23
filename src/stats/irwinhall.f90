! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Irwin-Hall distribution matching scipy.stats.irwinhall.
module scifort_irwinhall
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: irwinhall_cdf, irwinhall_isf, irwinhall_logcdf, irwinhall_logpdf
    public :: irwinhall_logsf, irwinhall_pdf, irwinhall_ppf, irwinhall_sf
    public :: irwinhall_logpdf_derivative

contains

    pure elemental function irwinhall_pdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z
        integer :: ni

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        ni = int(n)
        z = (x - mu) / sigma
        if (z < 0.0_dp .or. z > n) then
            y = 0.0_dp
        else if (ni == 1) then
            y = 1.0_dp / sigma
        else
            y = cardinal_bspline(ni, z) / sigma
        end if
    end function irwinhall_pdf

    pure elemental function irwinhall_logpdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: y, p
        p = irwinhall_pdf(x, n, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function irwinhall_logpdf

    pure elemental function irwinhall_cdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z
        integer :: ni

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        ni = int(n)
        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else if (z >= n) then
            y = 1.0_dp
        else if (z <= 0.5_dp * n) then
            y = standard_cdf(z, ni)
        else
            y = 1.0_dp - standard_cdf(n - z, ni)
        end if
        y = max(0.0_dp, min(1.0_dp, y))
    end function irwinhall_cdf

    pure elemental function irwinhall_sf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z
        integer :: ni

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        ni = int(n)
        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 1.0_dp
        else if (z >= n) then
            y = 0.0_dp
        else if (z >= 0.5_dp * n) then
            y = standard_cdf(n - z, ni)
        else
            y = 1.0_dp - standard_cdf(z, ni)
        end if
        y = max(0.0_dp, min(1.0_dp, y))
    end function irwinhall_sf

    pure elemental function irwinhall_logcdf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: y, p
        p = irwinhall_cdf(x, n, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function irwinhall_logcdf

    pure elemental function irwinhall_logsf(x, n, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: y, p
        p = irwinhall_sf(x, n, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function irwinhall_logsf

    pure elemental function irwinhall_ppf(p, n, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = mu
        else if (p >= 1.0_dp) then
            x = mu + sigma * n
        else if (p == 0.5_dp) then
            x = mu + 0.5_dp * sigma * n
        else if (p > 0.5_dp) then
            x = irwinhall_isf(1.0_dp - p, n, mu, sigma)
        else
            lo = 0.0_dp
            hi = n
            do iter = 1, 90
                mid = 0.5_dp * (lo + hi)
                if (irwinhall_cdf(mu + sigma * mid, n, mu, sigma) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function irwinhall_ppf

    pure elemental function irwinhall_isf(p, n, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(n) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p >= 1.0_dp) then
            x = mu
        else if (p <= 0.0_dp) then
            x = mu + sigma * n
        else if (p == 0.5_dp) then
            x = mu + 0.5_dp * sigma * n
        else
            lo = 0.0_dp
            hi = n
            do iter = 1, 90
                mid = 0.5_dp * (lo + hi)
                if (irwinhall_sf(mu + sigma * mid, n, mu, sigma) > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function irwinhall_isf

    pure elemental function irwinhall_logpdf_derivative(x, n) result(y)
        real(dp), intent(in) :: x !! standardized interior point
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp) :: y
        real(dp) :: f, fp
        integer :: ni

        if (.not. valid_shape(n)) then
            y = quiet_nan(x)
            return
        end if
        ni = int(n)
        if (x <= 0.0_dp .or. x >= n) then
            y = quiet_nan(x)
        else if (ni == 1) then
            y = 0.0_dp
        else
            f = cardinal_bspline(ni, x)
            fp = cardinal_bspline(ni - 1, x) - cardinal_bspline(ni - 1, x - 1.0_dp)
            if (f <= 0.0_dp) then
                y = quiet_nan(x)
            else
                y = fp / f
            end if
        end if
    end function irwinhall_logpdf_derivative

    pure function standard_cdf(x, n) result(p)
        real(dp), intent(in) :: x !! standardized point in (0,n/2]
        integer, intent(in) :: n !! positive integer order
        real(dp) :: p
        integer :: j, last

        p = 0.0_dp
        last = min(int(floor(x)), n)
        do j = 0, last
            p = p + cardinal_bspline(n + 1, x - real(j, dp))
        end do
    end function standard_cdf

    pure function cardinal_bspline(order, x) result(y)
        integer, intent(in) :: order !! cardinal B-spline order, >= 1
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: y
        real(dp) :: old(0:order + 1), new(0:order + 1), t
        integer :: k, r

        if (order < 1 .or. x < 0.0_dp .or. x > real(order, dp)) then
            y = 0.0_dp
            return
        end if
        if (order == 1) then
            if (x >= 0.0_dp .and. x <= 1.0_dp) then
                y = 1.0_dp
            else
                y = 0.0_dp
            end if
            return
        end if

        old = 0.0_dp
        do k = 0, order + 1
            t = x - real(k, dp)
            if (t >= 0.0_dp .and. t < 1.0_dp) old(k) = 1.0_dp
        end do
        do r = 2, order
            new = 0.0_dp
            do k = 0, order + 1 - r
                t = x - real(k, dp)
                new(k) = t * old(k) / real(r - 1, dp) + &
                    (real(r, dp) - t) * old(k + 1) / real(r - 1, dp)
            end do
            old = new
        end do
        y = max(0.0_dp, old(0))
    end function cardinal_bspline

    pure elemental function valid_shape(n) result(ok)
        real(dp), intent(in) :: n !! candidate integer shape
        logical :: ok
        ok = ieee_is_finite(n) .and. n >= 1.0_dp .and. n == floor(n)
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! requested location, absent for zero
        real(dp), intent(in), optional :: scale !! requested scale, absent for one
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_irwinhall
