! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Left-skewed Gumbel distribution. With z = (x - loc) / scale,
! f(x) = exp(z - exp(z)) / scale.

module scifort_gumbel_l
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: gumbel_l_cdf
    public :: gumbel_l_isf
    public :: gumbel_l_logcdf
    public :: gumbel_l_logpdf
    public :: gumbel_l_logsf
    public :: gumbel_l_pdf
    public :: gumbel_l_ppf
    public :: gumbel_l_sf

contains

    pure elemental function gumbel_l_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z > log(huge(1.0_dp))) then
            y = negative_infinity(x)
        else
            t = exp(z)
            y = z - t - log(sigma)
        end if
    end function gumbel_l_logpdf

    pure elemental function gumbel_l_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logp

        logp = gumbel_l_logpdf(x, loc, scale)
        if (ieee_is_nan(logp)) then
            y = logp
        else if (logp == negative_infinity(logp)) then
            y = 0.0_dp
        else
            y = exp(logp)
        end if
    end function gumbel_l_pdf

    pure elemental function gumbel_l_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z > log(huge(1.0_dp))) then
            y = negative_infinity(x)
        else
            y = -exp(z)
        end if
    end function gumbel_l_logsf

    pure elemental function gumbel_l_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logsf

        logsf = gumbel_l_logsf(x, loc, scale)
        if (ieee_is_nan(logsf)) then
            y = logsf
        else if (logsf == negative_infinity(logsf)) then
            y = 0.0_dp
        else
            y = exp(logsf)
        end if
    end function gumbel_l_sf

    pure elemental function gumbel_l_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z > log(huge(1.0_dp))) then
            y = 1.0_dp
        else
            t = exp(z)
            y = -expm1_safe(-t)
        end if
    end function gumbel_l_cdf

    pure elemental function gumbel_l_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z > log(huge(1.0_dp))) then
            y = 0.0_dp
        else
            t = exp(z)
            if (t == 0.0_dp) then
                y = z
            else
                y = log(-expm1_safe(-t))
            end if
        end if
    end function gumbel_l_logcdf

    pure elemental function gumbel_l_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            x = mu + sigma * log(-log1p_safe(-p))
        end if
    end function gumbel_l_ppf

    pure elemental function gumbel_l_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else
            x = mu + sigma * log(-log(p))
        end if
    end function gumbel_l_isf

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when loc is absent
        real(dp), intent(out) :: sigma !! scale, 1 when scale is absent

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_gumbel_l
