! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Real-line von Mises distribution matching scipy.stats.vonmises_line.
module scifort_vonmises_line
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan, valid_loc_scale
    use scifort_vonmises, only : vonmises_line_cdf_standard, vonmises_line_ppf_standard, &
        vonmises_standard_logpdf
    implicit none
    private

    public :: vonmises_line_cdf, vonmises_line_isf, vonmises_line_logcdf
    public :: vonmises_line_logpdf, vonmises_line_logsf, vonmises_line_pdf
    public :: vonmises_line_ppf, vonmises_line_sf

contains

    pure elemental function vonmises_line_logpdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! center of the finite support (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(kappa) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (.not. ieee_is_finite(z) .or. z < -scifort_pi .or. z > scifort_pi) then
                y = negative_infinity(x)
            else
                y = vonmises_standard_logpdf(z, kappa) - log(sigma)
            end if
        end if
    end function vonmises_line_logpdf

    pure elemental function vonmises_line_pdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! center of the finite support (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y, ly

        ly = vonmises_line_logpdf(x, kappa, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function vonmises_line_pdf

    pure elemental function vonmises_line_cdf(x, kappa, loc, scale) result(p)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! center of the finite support (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: p
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(kappa) .or. .not. valid_loc_scale(mu, sigma)) then
            p = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            p = quiet_nan(x)
        else
            z = (x - mu) / sigma
            p = vonmises_line_cdf_standard(z, kappa)
        end if
    end function vonmises_line_cdf

    pure elemental function vonmises_line_sf(x, kappa, loc, scale) result(p)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! center of the finite support (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: p
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(kappa) .or. .not. valid_loc_scale(mu, sigma)) then
            p = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            p = quiet_nan(x)
        else
            z = (x - mu) / sigma
            p = vonmises_line_cdf_standard(-z, kappa)
        end if
    end function vonmises_line_sf

    pure elemental function vonmises_line_logcdf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! center of the finite support (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y, p

        p = vonmises_line_cdf(x, kappa, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else if (p > 0.5_dp) then
            y = log1p_safe(-vonmises_line_sf(x, kappa, loc, scale))
        else
            y = log(p)
        end if
    end function vonmises_line_logcdf

    pure elemental function vonmises_line_logsf(x, kappa, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! center of the finite support (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: y, p

        p = vonmises_line_sf(x, kappa, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p == 0.0_dp) then
            y = negative_infinity(p)
        else if (p > 0.5_dp) then
            y = log1p_safe(-vonmises_line_cdf(x, kappa, loc, scale))
        else
            y = log(p)
        end if
    end function vonmises_line_logsf

    pure elemental function vonmises_line_ppf(p, kappa, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! center of the finite support (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(kappa) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            x = mu + sigma * vonmises_line_ppf_standard(p, kappa)
        end if
    end function vonmises_line_ppf

    pure elemental function vonmises_line_isf(p, kappa, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! center of the finite support (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(kappa) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            x = mu - sigma * vonmises_line_ppf_standard(p, kappa)
        end if
    end function vonmises_line_isf

    pure elemental logical function valid_shape(kappa) result(ok)
        real(dp), intent(in) :: kappa !! candidate concentration parameter
        ok = ieee_is_finite(kappa) .and. kappa >= 0.0_dp
    end function valid_shape

    pure elemental subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location parameter
        real(dp), intent(in), optional :: scale !! optional positive scale parameter
        real(dp), intent(out) :: mu !! resolved location parameter
        real(dp), intent(out) :: sigma !! resolved scale parameter

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_vonmises_line
