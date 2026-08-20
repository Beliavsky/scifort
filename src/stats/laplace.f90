! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_laplace
    use scifort_constants, only : scifort_log_two
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: laplace_cdf
    public :: laplace_isf
    public :: laplace_logcdf
    public :: laplace_logpdf
    public :: laplace_logsf
    public :: laplace_pdf
    public :: laplace_ppf
    public :: laplace_sf

contains

    pure elemental function laplace_pdf(x, loc, scale) result(y)
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
        y = 0.5_dp * exp(-abs(z)) / sigma
    end function laplace_pdf

    pure elemental function laplace_logpdf(x, loc, scale) result(y)
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
        y = -scifort_log_two - log(sigma) - abs(z)
    end function laplace_logpdf

    pure elemental function laplace_cdf(x, loc, scale) result(y)
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
            y = 0.5_dp * exp(z)
        else
            y = 1.0_dp - 0.5_dp * exp(-z)
        end if
    end function laplace_cdf

    pure elemental function laplace_sf(x, loc, scale) result(y)
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
            y = 1.0_dp - 0.5_dp * exp(z)
        else
            y = 0.5_dp * exp(-z)
        end if
    end function laplace_sf

    pure elemental function laplace_logcdf(x, loc, scale) result(y)
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
            y = -scifort_log_two + z
        else
            y = log1p_safe(-0.5_dp * exp(-z))
        end if
    end function laplace_logcdf

    pure elemental function laplace_logsf(x, loc, scale) result(y)
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
            y = log1p_safe(-0.5_dp * exp(z))
        else
            y = -scifort_log_two - z
        end if
    end function laplace_logsf

    pure elemental function laplace_ppf(p, loc, scale) result(x)
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
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
            return
        else if (p <= 0.0_dp) then
            x = negative_infinity(p)
            return
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
            return
        end if

        if (p < 0.5_dp) then
            z = scifort_log_two + log(p)
        else if (p > 0.5_dp) then
            z = -scifort_log_two - log1p_safe(-p)
        else
            z = 0.0_dp
        end if
        x = mu + sigma * z
    end function laplace_ppf

    pure elemental function laplace_isf(p, loc, scale) result(x)
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
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
            return
        else if (p <= 0.0_dp) then
            x = positive_infinity(p)
            return
        else if (p >= 1.0_dp) then
            x = negative_infinity(p)
            return
        end if

        if (p < 0.5_dp) then
            z = -scifort_log_two - log(p)
        else if (p > 0.5_dp) then
            z = scifort_log_two + log1p_safe(-p)
        else
            z = 0.0_dp
        end if
        x = mu + sigma * z
    end function laplace_isf

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

end module scifort_laplace
