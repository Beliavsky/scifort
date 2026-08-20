! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_normal
    use scifort_constants, only : scifort_inv_sqrt_two_pi, &
        scifort_log_sqrt_two_pi, scifort_sqrt_two
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: normal_cdf
    public :: normal_isf
    public :: normal_logcdf
    public :: normal_logpdf
    public :: normal_logsf
    public :: normal_pdf
    public :: normal_ppf
    public :: normal_sf

contains

    pure elemental function normal_pdf(x, loc, scale) result(y)
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
        y = scifort_inv_sqrt_two_pi * exp(-0.5_dp * z * z) / sigma
    end function normal_pdf

    pure elemental function normal_logpdf(x, loc, scale) result(y)
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
        y = -0.5_dp * z * z - log(sigma) - scifort_log_sqrt_two_pi
    end function normal_logpdf

    pure elemental function normal_cdf(x, loc, scale) result(y)
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
        y = standard_normal_cdf(z)
    end function normal_cdf

    pure elemental function normal_sf(x, loc, scale) result(y)
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
        y = standard_normal_sf(z)
    end function normal_sf

    pure elemental function normal_logcdf(x, loc, scale) result(y)
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
        y = standard_normal_logcdf(z)
    end function normal_logcdf

    pure elemental function normal_logsf(x, loc, scale) result(y)
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
        y = standard_normal_logcdf(-z)
    end function normal_logsf

    pure elemental function normal_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
            return
        end if

        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            x = mu + sigma * standard_normal_ppf(p)
        end if
    end function normal_ppf

    pure elemental function normal_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
            return
        end if

        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            x = mu - sigma * standard_normal_ppf(p)
        end if
    end function normal_isf

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

    pure elemental function standard_normal_cdf(z) result(p)
        real(dp), intent(in) :: z
        real(dp) :: p

        p = 0.5_dp * erfc(-z / scifort_sqrt_two)
    end function standard_normal_cdf

    pure elemental function standard_normal_sf(z) result(p)
        real(dp), intent(in) :: z
        real(dp) :: p

        p = 0.5_dp * erfc(z / scifort_sqrt_two)
    end function standard_normal_sf

    pure elemental function standard_normal_logcdf(z) result(y)
        real(dp), intent(in) :: z
        real(dp) :: y

        real(dp) :: inv_z2
        real(dp) :: series
        real(dp) :: tail

        if (z < -20.0_dp) then
            inv_z2 = 1.0_dp / (z * z)
            series = 1.0_dp - inv_z2 * (1.0_dp - inv_z2 * &
                (3.0_dp - inv_z2 * (15.0_dp - inv_z2 * &
                (105.0_dp - inv_z2 * (945.0_dp - 10395.0_dp * inv_z2)))))
            y = -0.5_dp * z * z - log(-z) - scifort_log_sqrt_two_pi + log(series)
        else if (z <= 0.0_dp) then
            y = log(standard_normal_cdf(z))
        else
            tail = standard_normal_sf(z)
            y = log1p_safe(-tail)
        end if
    end function standard_normal_logcdf

    pure elemental function standard_normal_ppf(p) result(x)
        real(dp), intent(in) :: p
        real(dp) :: x

        integer :: iteration
        real(dp) :: lower
        real(dp) :: midpoint
        real(dp) :: q
        real(dp) :: target
        real(dp) :: upper

        if (p <= 0.0_dp) then
            x = negative_infinity(p)
            return
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
            return
        end if

        if (p < 0.5_dp) then
            q = p
        else if (p > 0.5_dp) then
            q = 1.0_dp - p
        else
            x = 0.0_dp
            return
        end if
        lower = -40.0_dp
        upper = 0.0_dp

        if (q < 1.0e-100_dp) then
            target = log(q)
            do iteration = 1, 100
                midpoint = 0.5_dp * (lower + upper)
                if (midpoint <= lower .or. midpoint >= upper) exit
                if (standard_normal_logcdf(midpoint) < target) then
                    lower = midpoint
                else
                    upper = midpoint
                end if
            end do
        else
            target = q
            do iteration = 1, 100
                midpoint = 0.5_dp * (lower + upper)
                if (midpoint <= lower .or. midpoint >= upper) exit
                if (standard_normal_cdf(midpoint) < target) then
                    lower = midpoint
                else
                    upper = midpoint
                end if
            end do
        end if

        x = 0.5_dp * (lower + upper)
        if (p > 0.5_dp) x = -x
    end function standard_normal_ppf

end module scifort_normal
