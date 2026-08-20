! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_cauchy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_log_two, scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: cauchy_cdf
    public :: cauchy_isf
    public :: cauchy_logcdf
    public :: cauchy_logpdf
    public :: cauchy_logsf
    public :: cauchy_pdf
    public :: cauchy_ppf
    public :: cauchy_sf

contains

    pure elemental function cauchy_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        y = exp(cauchy_logpdf(x, loc, scale))
    end function cauchy_pdf

    pure elemental function cauchy_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
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

        z = (x - mu) / sigma
        a = abs(z)
        if (a <= sqrt(huge(1.0_dp))) then
            log_one_plus_z2 = log(1.0_dp + a * a)
        else
            log_one_plus_z2 = 2.0_dp * log(a)
        end if
        y = -log(scifort_pi) - log(sigma) - log_one_plus_z2
    end function cauchy_logpdf

    pure elemental function cauchy_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
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
        else if (z < 0.0_dp) then
            y = atan(-1.0_dp / z) / scifort_pi
        else if (z > 0.0_dp) then
            y = 1.0_dp - atan(1.0_dp / z) / scifort_pi
        else
            y = 0.5_dp
        end if
    end function cauchy_cdf

    pure elemental function cauchy_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
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
        else if (z > 0.0_dp) then
            y = atan(1.0_dp / z) / scifort_pi
        else if (z < 0.0_dp) then
            y = 1.0_dp - atan(-1.0_dp / z) / scifort_pi
        else
            y = 0.5_dp
        end if
    end function cauchy_sf

    pure elemental function cauchy_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: tail
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = quiet_nan(x)
        else if (z < 0.0_dp) then
            tail = atan(-1.0_dp / z) / scifort_pi
            if (tail <= 0.0_dp) then
                y = negative_infinity(x)
            else
                y = log(tail)
            end if
        else if (z > 0.0_dp) then
            tail = atan(1.0_dp / z) / scifort_pi
            y = log1p_safe(-tail)
        else
            y = -scifort_log_two
        end if
    end function cauchy_logcdf

    pure elemental function cauchy_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: tail
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = quiet_nan(x)
        else if (z > 0.0_dp) then
            tail = atan(1.0_dp / z) / scifort_pi
            if (tail <= 0.0_dp) then
                y = negative_infinity(x)
            else
                y = log(tail)
            end if
        else if (z < 0.0_dp) then
            tail = atan(-1.0_dp / z) / scifort_pi
            y = log1p_safe(-tail)
        else
            y = -scifort_log_two
        end if
    end function cauchy_logsf

    pure elemental function cauchy_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
            return
        end if

        z = standard_cauchy_ppf(p)
        x = mu + sigma * z
    end function cauchy_ppf

    pure elemental function cauchy_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
            return
        end if

        z = standard_cauchy_ppf(p)
        x = mu - sigma * z
    end function cauchy_isf

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

    pure elemental function standard_cauchy_ppf(p) result(z)
        real(dp), intent(in) :: p
        real(dp) :: z

        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            z = quiet_nan(p)
        else if (p <= 0.0_dp) then
            z = negative_infinity(p)
        else if (p >= 1.0_dp) then
            z = positive_infinity(p)
        else if (p < 0.5_dp) then
            z = -1.0_dp / tan(scifort_pi * p)
        else if (p > 0.5_dp) then
            z = 1.0_dp / tan(scifort_pi * (1.0_dp - p))
        else
            z = 0.0_dp
        end if
    end function standard_cauchy_ppf

end module scifort_cauchy
