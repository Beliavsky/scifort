! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_exponential
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: exponential_cdf
    public :: exponential_isf
    public :: exponential_logcdf
    public :: exponential_logpdf
    public :: exponential_logsf
    public :: exponential_pdf
    public :: exponential_ppf
    public :: exponential_sf

contains

    pure elemental function exponential_pdf(x, loc, scale) result(y)
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
        if (z < 0.0_dp) then
            y = 0.0_dp
        else
            y = exp(-z) / sigma
        end if
    end function exponential_pdf

    pure elemental function exponential_logpdf(x, loc, scale) result(y)
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
        if (z < 0.0_dp) then
            y = negative_infinity(x)
        else
            y = -z - log(sigma)
        end if
    end function exponential_logpdf

    pure elemental function exponential_cdf(x, loc, scale) result(y)
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
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else
            y = -expm1_safe(-z)
        end if
    end function exponential_cdf

    pure elemental function exponential_sf(x, loc, scale) result(y)
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
        if (z < 0.0_dp) then
            y = 1.0_dp
        else
            y = exp(-z)
        end if
    end function exponential_sf

    pure elemental function exponential_logcdf(x, loc, scale) result(y)
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
        if (z <= 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(-expm1_safe(-z))
        end if
    end function exponential_logcdf

    pure elemental function exponential_logsf(x, loc, scale) result(y)
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
        if (z < 0.0_dp) then
            y = 0.0_dp
        else
            y = -z
        end if
    end function exponential_logsf

    pure elemental function exponential_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
        else
            x = mu - sigma * log1p_safe(-p)
        end if
    end function exponential_ppf

    pure elemental function exponential_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = positive_infinity(p)
        else
            x = mu - sigma * log(p)
        end if
    end function exponential_isf

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

end module scifort_exponential
