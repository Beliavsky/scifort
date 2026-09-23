! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Asymmetric Laplace distribution with shape kappa > 0. In standardized
! coordinates z=(x-loc)/scale, density is exp(-kappa*z)/(kappa+1/kappa)
! for z>=0 and exp(z/kappa)/(kappa+1/kappa) for z<0.
! This matches scipy.stats.laplace_asymmetric.

module scifort_laplace_asymmetric
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: laplace_asymmetric_cdf
    public :: laplace_asymmetric_isf
    public :: laplace_asymmetric_logcdf
    public :: laplace_asymmetric_logpdf
    public :: laplace_asymmetric_logsf
    public :: laplace_asymmetric_pdf
    public :: laplace_asymmetric_ppf
    public :: laplace_asymmetric_sf

contains

    pure elemental function laplace_asymmetric_pdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = exp(laplace_asymmetric_logpdf(x, kappa, loc, scale))
    end function laplace_asymmetric_pdf

    pure elemental function laplace_asymmetric_logpdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: invk
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(kappa))) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        invk = 1.0_dp / kappa
        z = (x - mu) / sigma
        if (z >= 0.0_dp) then
            y = -kappa * z - log(kappa + invk) - log(sigma)
        else
            y = invk * z - log(kappa + invk) - log(sigma)
        end if
    end function laplace_asymmetric_logpdf

    pure elemental function laplace_asymmetric_cdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: invk
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(kappa))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            invk = 1.0_dp / kappa
            z = (x - mu) / sigma
            if (z >= 0.0_dp) then
                y = 1.0_dp - exp(-kappa * z) * invk / (kappa + invk)
            else
                y = exp(invk * z) * kappa / (kappa + invk)
            end if
        end if
    end function laplace_asymmetric_cdf

    pure elemental function laplace_asymmetric_sf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: invk
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(kappa))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            invk = 1.0_dp / kappa
            z = (x - mu) / sigma
            if (z >= 0.0_dp) then
                y = exp(-kappa * z) * invk / (kappa + invk)
            else
                y = 1.0_dp - exp(invk * z) * kappa / (kappa + invk)
            end if
        end if
    end function laplace_asymmetric_sf

    pure elemental function laplace_asymmetric_logcdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: invk
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: tail
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(kappa))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            invk = 1.0_dp / kappa
            z = (x - mu) / sigma
            if (z >= 0.0_dp) then
                tail = exp(-kappa * z) * invk / (kappa + invk)
                y = log1p_safe(-tail)
            else
                y = invk * z + log(kappa) - log(kappa + invk)
            end if
        end if
    end function laplace_asymmetric_logcdf

    pure elemental function laplace_asymmetric_logsf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: invk
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: lower
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(kappa))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            invk = 1.0_dp / kappa
            z = (x - mu) / sigma
            if (z >= 0.0_dp) then
                y = -kappa * z + log(invk) - log(kappa + invk)
            else
                lower = exp(invk * z) * kappa / (kappa + invk)
                y = log1p_safe(-lower)
            end if
        end if
    end function laplace_asymmetric_logsf

    pure elemental function laplace_asymmetric_ppf(p, kappa, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: cutoff
        real(dp) :: invk
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(kappa))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = negative_infinity(p)
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
        else
            invk = 1.0_dp / kappa
            cutoff = kappa / (kappa + invk)
            if (p >= cutoff) then
                z = -log((1.0_dp - p) * (kappa + invk) * kappa) * invk
            else
                z = log(p * (kappa + invk) / kappa) * kappa
            end if
            x = mu + sigma * z
        end if
    end function laplace_asymmetric_ppf

    pure elemental function laplace_asymmetric_isf(p, kappa, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: cutoff
        real(dp) :: invk
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(kappa))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = positive_infinity(p)
        else if (p >= 1.0_dp) then
            x = negative_infinity(p)
        else
            invk = 1.0_dp / kappa
            cutoff = invk / (kappa + invk)
            if (p <= cutoff) then
                z = -log(p * (kappa + invk) * kappa) * invk
            else
                z = log((1.0_dp - p) * (kappa + invk) / kappa) * kappa
            end if
            x = mu + sigma * z
        end if
    end function laplace_asymmetric_isf

    pure elemental logical function valid_shape(kappa) result(valid)
        real(dp), intent(in) :: kappa !! asymmetry parameter to validate
        valid = ieee_is_finite(kappa) .and. kappa > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location argument
        real(dp), intent(in), optional :: scale !! optional scale argument
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_laplace_asymmetric
