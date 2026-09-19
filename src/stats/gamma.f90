! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Gamma distribution with shape a > 0, location loc, and scale > 0:
! f(x) = z**(a - 1) exp(-z) / (Gamma(a) scale), z = (x - loc) / scale >= 0.
! Probabilities use the regularized incomplete gamma functions in
! scifort_incomplete_gamma. This matches the parameterization of
! scipy.stats.gamma.

module scifort_gamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv, log_gamma_kernel, log_gammainc, log_gammaincc
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, &
        valid_loc_scale
    implicit none
    private

    public :: gamma_cdf
    public :: gamma_isf
    public :: gamma_logcdf
    public :: gamma_logpdf
    public :: gamma_logsf
    public :: gamma_pdf
    public :: gamma_ppf
    public :: gamma_sf

contains

    pure elemental function gamma_pdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z < 0.0_dp .or. z > huge(z)) then
            y = 0.0_dp
        else if (z == 0.0_dp) then
            if (a < 1.0_dp) then
                y = positive_infinity(z)
            else if (a == 1.0_dp) then
                y = 1.0_dp / sigma
            else
                y = 0.0_dp
            end if
        else
            y = exp(standard_logpdf(a, z)) / sigma
        end if
    end function gamma_pdf

    pure elemental function gamma_logpdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z < 0.0_dp .or. z > huge(z)) then
            y = negative_infinity(z)
        else if (z == 0.0_dp) then
            if (a < 1.0_dp) then
                y = positive_infinity(z)
            else if (a == 1.0_dp) then
                y = -log(sigma)
            else
                y = negative_infinity(z)
            end if
        else
            y = standard_logpdf(a, z) - log(sigma)
        end if
    end function gamma_logpdf

    pure elemental function gamma_cdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = 0.0_dp
        else
            y = gammainc(a, z)
        end if
    end function gamma_cdf

    pure elemental function gamma_sf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = 1.0_dp
        else
            y = gammaincc(a, z)
        end if
    end function gamma_sf

    pure elemental function gamma_logcdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = negative_infinity(z)
        else
            y = log_gammainc(a, z)
        end if
    end function gamma_logcdf

    pure elemental function gamma_logsf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = 0.0_dp
        else
            y = log_gammaincc(a, z)
        end if
    end function gamma_logsf

    pure elemental function gamma_ppf(p, a, loc, scale) result(y)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: a
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(p)
            return
        end if

        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else
            y = mu + sigma * gammaincinv(a, p)
        end if
    end function gamma_ppf

    pure elemental function gamma_isf(p, a, loc, scale) result(y)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: a
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(p)
            return
        end if

        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else
            y = mu + sigma * gammainccinv(a, p)
        end if
    end function gamma_isf

    ! log of the standard density for finite z > 0. For large a the density is
    ! formed as a * kernel(a, z) / z, where kernel = z**a exp(-z) / Gamma(a + 1)
    ! is evaluated without the cancellation of (a - 1) log(z) - z - log(Gamma(a)).
    pure elemental function standard_logpdf(a, z) result(y)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: z
        real(dp) :: y

        if (a < 10.0_dp) then
            y = (a - 1.0_dp) * log(z) - z - log_gamma(a)
        else
            y = log(a) + log_gamma_kernel(a, z) - log(z)
        end if
    end function standard_logpdf

    pure elemental logical function valid_shape(a) result(valid)
        real(dp), intent(in) :: a

        valid = ieee_is_finite(a) .and. a > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp), intent(out) :: mu
        real(dp), intent(out) :: sigma

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_gamma
