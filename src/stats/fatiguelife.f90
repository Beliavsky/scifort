! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Birnbaum-Saunders (fatigue-life) distribution with shape c > 0.
! In standardized coordinates z = (x-loc)/scale > 0,
! F(z) = Phi((sqrt(z)-1/sqrt(z))/c). This matches scipy.stats.fatiguelife.

module scifort_fatiguelife
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, normal_logsf, &
        normal_ppf, normal_sf
    implicit none
    private

    public :: fatiguelife_cdf
    public :: fatiguelife_isf
    public :: fatiguelife_logcdf
    public :: fatiguelife_logpdf
    public :: fatiguelife_logsf
    public :: fatiguelife_pdf
    public :: fatiguelife_ppf
    public :: fatiguelife_sf

contains

    pure elemental function fatiguelife_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy
        logy = fatiguelife_logpdf(x, c, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function fatiguelife_pdf

    pure elemental function fatiguelife_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: h
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            h = z - 2.0_dp + 1.0_dp / z
            y = log(z + 1.0_dp) - log(2.0_dp * c) - scifort_log_sqrt_two_pi - &
                1.5_dp * log(z) - 0.5_dp * h / (c * c) - log(sigma)
        end if
    end function fatiguelife_logpdf

    pure elemental function standardized_w(x, c, mu, sigma) result(w)
        real(dp), intent(in) :: x !! finite point above the lower support endpoint
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive scale parameter
        real(dp) :: w
        real(dp) :: root_z
        real(dp) :: z
        z = (x - mu) / sigma
        root_z = sqrt(z)
        w = (root_z - 1.0_dp / root_z) / c
    end function standardized_w

    pure elemental function fatiguelife_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu
        real(dp) :: sigma
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
            y = normal_cdf(standardized_w(x, c, mu, sigma))
        end if
    end function fatiguelife_cdf

    pure elemental function fatiguelife_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu
        real(dp) :: sigma
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
            y = normal_sf(standardized_w(x, c, mu, sigma))
        end if
    end function fatiguelife_sf

    pure elemental function fatiguelife_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu
        real(dp) :: sigma
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
            y = normal_logcdf(standardized_w(x, c, mu, sigma))
        end if
    end function fatiguelife_logcdf

    pure elemental function fatiguelife_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu
        real(dp) :: sigma
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
            y = normal_logsf(standardized_w(x, c, mu, sigma))
        end if
    end function fatiguelife_logsf

    pure elemental function fatiguelife_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: y
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            y = c * normal_ppf(p)
            x = affine_positive(mu, sigma, positive_root_square(y), p)
        end if
    end function fatiguelife_ppf

    pure elemental function fatiguelife_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: y
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            y = c * normal_isf(p)
            x = affine_positive(mu, sigma, positive_root_square(y), p)
        end if
    end function fatiguelife_isf

    pure elemental function positive_root_square(y) result(z)
        real(dp), intent(in) :: y !! transformed normal quantile
        real(dp) :: z
        real(dp) :: h
        real(dp) :: root
        h = hypot(y, 2.0_dp)
        if (y >= 0.0_dp) then
            root = 0.5_dp * (y + h)
        else
            root = 2.0_dp / (h - y)
        end if
        z = root * root
    end function positive_root_square

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to construct positive infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z) .or. z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive

    pure elemental logical function valid_shape(c) result(valid)
        real(dp), intent(in) :: c !! shape parameter to check
        valid = ieee_is_finite(c) .and. c > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location
        real(dp), intent(in), optional :: scale !! optional scale
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_fatiguelife
