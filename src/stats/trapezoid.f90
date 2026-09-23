! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Trapezoidal distribution, matching scipy.stats.trapezoid.

module scifort_trapezoid
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: trapezoid_cdf, trapezoid_isf, trapezoid_logcdf, trapezoid_logpdf
    public :: trapezoid_logsf, trapezoid_pdf, trapezoid_ppf, trapezoid_sf

contains

    pure elemental function trapezoid_pdf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! standardized left edge of the plateau in [0,1]
        real(dp), intent(in) :: d !! standardized right edge of the plateau in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, mu, sigma, z, h
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, d) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            h = 2.0_dp / (1.0_dp + d - c) / sigma
            if (z < 0.0_dp .or. z > 1.0_dp .or. .not. ieee_is_finite(z)) then
                y = 0.0_dp
            else if (z < c) then
                if (c == 0.0_dp) then
                    y = h
                else
                    y = h * z / c
                end if
            else if (z <= d) then
                y = h
            else
                if (d == 1.0_dp) then
                    y = h
                else
                    y = h * (1.0_dp - z) / (1.0_dp - d)
                end if
            end if
        end if
    end function trapezoid_pdf

    pure elemental function trapezoid_logpdf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! standardized left edge of the plateau in [0,1]
        real(dp), intent(in) :: d !! standardized right edge of the plateau in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, p
        p = trapezoid_pdf(x, c, d, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(p)
        end if
    end function trapezoid_logpdf

    pure elemental function trapezoid_cdf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! standardized left edge of the plateau in [0,1]
        real(dp), intent(in) :: d !! standardized right edge of the plateau in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, mu, sigma, z, den
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, d) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            den = 1.0_dp + d - c
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (z >= 1.0_dp) then
                y = 1.0_dp
            else if (z < c) then
                y = z * z / (c * den)
            else if (z <= d) then
                y = (c + 2.0_dp * (z - c)) / den
            else
                y = 1.0_dp - (1.0_dp - z) * (1.0_dp - z) / ((1.0_dp - d) * den)
            end if
        end if
    end function trapezoid_cdf

    pure elemental function trapezoid_sf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! standardized left edge of the plateau in [0,1]
        real(dp), intent(in) :: d !! standardized right edge of the plateau in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, mu, sigma, z, den
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, d) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            den = 1.0_dp + d - c
            if (z <= 0.0_dp) then
                y = 1.0_dp
            else if (z >= 1.0_dp) then
                y = 0.0_dp
            else if (z < c) then
                y = 1.0_dp - z * z / (c * den)
            else if (z <= d) then
                y = 1.0_dp - (c + 2.0_dp * (z - c)) / den
            else
                y = (1.0_dp - z) * (1.0_dp - z) / ((1.0_dp - d) * den)
            end if
        end if
    end function trapezoid_sf

    pure elemental function trapezoid_logcdf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! standardized left edge of the plateau in [0,1]
        real(dp), intent(in) :: d !! standardized right edge of the plateau in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, p, q
        p = trapezoid_cdf(x, c, d, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(x)
        else if (p < 0.5_dp) then
            y = log(p)
        else
            q = trapezoid_sf(x, c, d, loc, scale)
            y = log1p_safe(-q)
        end if
    end function trapezoid_logcdf

    pure elemental function trapezoid_logsf(x, c, d, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! standardized left edge of the plateau in [0,1]
        real(dp), intent(in) :: d !! standardized right edge of the plateau in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, p, q
        q = trapezoid_sf(x, c, d, loc, scale)
        if (ieee_is_nan(q)) then
            y = q
        else if (q == 0.0_dp) then
            y = negative_infinity(x)
        else if (q < 0.5_dp) then
            y = log(q)
        else
            p = trapezoid_cdf(x, c, d, loc, scale)
            y = log1p_safe(-p)
        end if
    end function trapezoid_logsf

    pure elemental function trapezoid_ppf(p, c, d, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! standardized left edge of the plateau in [0,1]
        real(dp), intent(in) :: d !! standardized right edge of the plateau in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x, mu, sigma, den, qc, qd, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, d) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            den = 1.0_dp + d - c
            qc = c / den
            qd = (2.0_dp * d - c) / den
            if (p == 0.0_dp) then
                z = 0.0_dp
            else if (p < qc) then
                z = sqrt(p * c * den)
            else if (p <= qd) then
                z = 0.5_dp * (p * den + c)
            else if (p < 1.0_dp) then
                z = 1.0_dp - sqrt((1.0_dp - p) * den * (1.0_dp - d))
            else
                z = 1.0_dp
            end if
            x = affine_bounded(mu, sigma, z)
        end if
    end function trapezoid_ppf

    pure elemental function trapezoid_isf(p, c, d, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! standardized left edge of the plateau in [0,1]
        real(dp), intent(in) :: d !! standardized right edge of the plateau in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x, mu, sigma, den, sc, sd, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(c, d) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            den = 1.0_dp + d - c
            sc = 1.0_dp - c / den
            sd = (1.0_dp - d) / den
            if (p == 1.0_dp) then
                z = 0.0_dp
            else if (p > sc) then
                z = sqrt((1.0_dp - p) * c * den)
            else if (p >= sd) then
                z = 0.5_dp * ((1.0_dp - p) * den + c)
            else if (p > 0.0_dp) then
                z = 1.0_dp - sqrt(p * den * (1.0_dp - d))
            else
                z = 1.0_dp
            end if
            x = affine_bounded(mu, sigma, z)
        end if
    end function trapezoid_isf

    pure elemental logical function valid_shapes(c, d) result(ok)
        real(dp), intent(in) :: c !! candidate left plateau edge
        real(dp), intent(in) :: d !! candidate right plateau edge
        ok = ieee_is_finite(c) .and. ieee_is_finite(d) .and. &
            c >= 0.0_dp .and. c <= d .and. d <= 1.0_dp
    end function valid_shapes

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

    pure elemental function affine_bounded(mu, sigma, z) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! standardized quantile in [0,1]
        real(dp) :: x
        if (mu > 0.0_dp .and. sigma * z > huge(1.0_dp) - mu) then
            x = huge(1.0_dp)
        else
            x = mu + sigma * z
        end if
    end function affine_bounded

end module scifort_trapezoid
