! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Reciprocal inverse-Gaussian distribution, matching scipy.stats.recipinvgauss.

module scifort_recipinvgauss
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_invgauss, only : invgauss_cdf, invgauss_isf, invgauss_logcdf, &
        invgauss_logsf, invgauss_ppf, invgauss_sf
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: recipinvgauss_cdf, recipinvgauss_isf, recipinvgauss_logcdf
    public :: recipinvgauss_logpdf, recipinvgauss_logsf, recipinvgauss_pdf
    public :: recipinvgauss_ppf, recipinvgauss_sf

contains

    pure elemental function recipinvgauss_pdf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, logy
        logy = recipinvgauss_logpdf(x, mu_shape, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function recipinvgauss_pdf

    pure elemental function recipinvgauss_logpdf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, delta
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(mu_shape) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                delta = 1.0_dp - mu_shape * z
                y = -0.5_dp * delta * delta / (z * mu_shape * mu_shape) - &
                    0.5_dp * log(z) - scifort_log_sqrt_two_pi - log(sigma)
            end if
        end if
    end function recipinvgauss_logpdf

    pure elemental function recipinvgauss_cdf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(mu_shape) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            y = invgauss_sf(1.0_dp / z, mu_shape)
        end if
    end function recipinvgauss_cdf

    pure elemental function recipinvgauss_sf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(mu_shape) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            y = invgauss_cdf(1.0_dp / z, mu_shape)
        end if
    end function recipinvgauss_sf

    pure elemental function recipinvgauss_logcdf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(mu_shape) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            y = invgauss_logsf(1.0_dp / z, mu_shape)
        end if
    end function recipinvgauss_logcdf

    pure elemental function recipinvgauss_logsf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(mu_shape) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            y = invgauss_logcdf(1.0_dp / z, mu_shape)
        end if
    end function recipinvgauss_logsf

    pure elemental function recipinvgauss_ppf(p, mu_shape, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, q, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(mu_shape) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            q = invgauss_isf(p, mu_shape)
            z = 1.0_dp / q
            x = affine_positive(mu, sigma, z, p)
        end if
    end function recipinvgauss_ppf

    pure elemental function recipinvgauss_isf(p, mu_shape, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, q, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(mu_shape) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            q = invgauss_ppf(p, mu_shape)
            z = 1.0_dp / q
            x = affine_positive(mu, sigma, z, p)
        end if
    end function recipinvgauss_isf

    pure elemental logical function valid_shape(mu_shape) result(ok)
        real(dp), intent(in) :: mu_shape !! candidate reciprocal inverse-Gaussian shape
        ok = ieee_is_finite(mu_shape) .and. mu_shape > 0.0_dp
    end function valid_shape

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

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when absent
        real(dp), intent(out) :: sigma !! scale, 1 when absent
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_recipinvgauss
