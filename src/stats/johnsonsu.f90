! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Johnson SU distribution, matching scipy.stats.johnsonsu.

module scifort_johnsonsu
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, normal_logsf, &
        normal_ppf, normal_sf
    implicit none
    private

    public :: johnsonsu_cdf, johnsonsu_isf, johnsonsu_logcdf, johnsonsu_logpdf
    public :: johnsonsu_logsf, johnsonsu_pdf, johnsonsu_ppf, johnsonsu_sf

contains

    pure elemental function johnsonsu_pdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = johnsonsu_logpdf(x, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function johnsonsu_pdf

    pure elemental function johnsonsu_logpdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, h, t, q
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            h = hypot(z, 1.0_dp)
            t = asinh(z)
            q = a + b * t
            y = log(b) - log(h) - 0.5_dp * q * q - &
                scifort_log_sqrt_two_pi - log(sigma)
        end if
    end function johnsonsu_logpdf

    pure elemental function johnsonsu_cdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
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
            y = normal_cdf(a + b * asinh(z))
        end if
    end function johnsonsu_cdf

    pure elemental function johnsonsu_sf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
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
            y = normal_sf(a + b * asinh(z))
        end if
    end function johnsonsu_sf

    pure elemental function johnsonsu_logcdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
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
            y = normal_logcdf(a + b * asinh(z))
        end if
    end function johnsonsu_logcdf

    pure elemental function johnsonsu_logsf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
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
            y = normal_logsf(a + b * asinh(z))
        end if
    end function johnsonsu_logsf

    pure elemental function johnsonsu_ppf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            t = (normal_ppf(p) - a) / b
            z = safe_sinh(t, p)
            x = affine(mu, sigma, z, p)
        end if
    end function johnsonsu_ppf

    pure elemental function johnsonsu_isf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            t = (normal_isf(p) - a) / b
            z = safe_sinh(t, p)
            x = affine(mu, sigma, z, p)
        end if
    end function johnsonsu_isf

    pure elemental logical function valid_shapes(a, b) result(ok)
        real(dp), intent(in) :: a !! candidate first Johnson SU shape parameter
        real(dp), intent(in) :: b !! candidate second Johnson SU shape parameter
        ok = ieee_is_finite(a) .and. ieee_is_finite(b) .and. b > 0.0_dp
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

    pure elemental function safe_sinh(t, seed) result(z)
        real(dp), intent(in) :: t !! real hyperbolic-sine argument
        real(dp), intent(in) :: seed !! value used to construct signed infinity
        real(dp) :: z, limit
        limit = log(huge(1.0_dp)) + log(2.0_dp)
        if (t > limit) then
            z = positive_infinity(seed)
        else if (t < -limit) then
            z = negative_infinity(seed)
        else
            z = sinh(t)
        end if
    end function safe_sinh

    pure elemental function affine(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! standardized quantile
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
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
    end function affine

end module scifort_johnsonsu
