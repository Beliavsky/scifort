! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Johnson SB distribution, matching scipy.stats.johnsonsb.

module scifort_johnsonsb
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, normal_logsf, &
        normal_ppf, normal_sf
    use scifort_special_elementary, only : expit, logit
    implicit none
    private

    public :: johnsonsb_cdf, johnsonsb_isf, johnsonsb_logcdf, johnsonsb_logpdf
    public :: johnsonsb_logsf, johnsonsb_pdf, johnsonsb_ppf, johnsonsb_sf

contains

    pure elemental function johnsonsb_pdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, ly
        ly = johnsonsb_logpdf(x, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function johnsonsb_pdf

    pure elemental function johnsonsb_logpdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, mu, sigma, z, t, q
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp .or. z >= 1.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                t = logit(z)
                q = a + b * t
                y = log(b) - log(z) - log(1.0_dp - z) - &
                    0.5_dp * q * q - scifort_log_sqrt_two_pi - log(sigma)
            end if
        end if
    end function johnsonsb_logpdf

    pure elemental function johnsonsb_cdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (z >= 1.0_dp) then
                y = 1.0_dp
            else
                y = normal_cdf(a + b * logit(z))
            end if
        end if
    end function johnsonsb_cdf

    pure elemental function johnsonsb_sf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 1.0_dp
            else if (z >= 1.0_dp) then
                y = 0.0_dp
            else
                y = normal_sf(a + b * logit(z))
            end if
        end if
    end function johnsonsb_sf

    pure elemental function johnsonsb_logcdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = negative_infinity(x)
            else if (z >= 1.0_dp) then
                y = 0.0_dp
            else
                y = normal_logcdf(a + b * logit(z))
            end if
        end if
    end function johnsonsb_logcdf

    pure elemental function johnsonsb_logsf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (z >= 1.0_dp) then
                y = negative_infinity(x)
            else
                y = normal_logsf(a + b * logit(z))
            end if
        end if
    end function johnsonsb_logsf

    pure elemental function johnsonsb_ppf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = affine_bounded(mu, sigma, 1.0_dp, p)
        else
            z = expit((normal_ppf(p) - a) / b)
            x = affine_bounded(mu, sigma, z, p)
        end if
    end function johnsonsb_ppf

    pure elemental function johnsonsb_isf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive finite second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = affine_bounded(mu, sigma, 1.0_dp, p)
        else
            z = expit((normal_isf(p) - a) / b)
            x = affine_bounded(mu, sigma, z, p)
        end if
    end function johnsonsb_isf

    pure elemental logical function valid_shapes(a, b) result(ok)
        real(dp), intent(in) :: a !! candidate first Johnson SB shape parameter
        real(dp), intent(in) :: b !! candidate second Johnson SB shape parameter
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

    pure elemental function affine_bounded(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! standardized quantile in [0,1]
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: x
        if (mu > 0.0_dp .and. sigma * z > huge(1.0_dp) - mu) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_bounded

end module scifort_johnsonsb
