! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! R distribution with shape c > 0, standardized support [-1, 1]. It is the
! affine image 2*Y-1 for Y ~ Beta(c/2, c/2), matching scipy.stats.rdist.

module scifort_rdist
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_beta, only : beta_cdf, beta_isf, beta_logcdf, beta_logpdf, &
        beta_logsf, beta_pdf, beta_ppf, beta_sf
    use scifort_constants, only : scifort_log_two
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: rdist_cdf
    public :: rdist_isf
    public :: rdist_logcdf
    public :: rdist_logpdf
    public :: rdist_logsf
    public :: rdist_pdf
    public :: rdist_ppf
    public :: rdist_sf

contains

    pure elemental function rdist_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: y
        real(dp) :: logy

        logy = rdist_logpdf(x, c, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function rdist_pdf

    pure elemental function rdist_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: y
        real(dp) :: a, mu, sigma, u, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < -1.0_dp .or. z > 1.0_dp) then
                y = negative_infinity(x)
            else
                a = 0.5_dp * c
                u = 0.5_dp * (z + 1.0_dp)
                y = beta_logpdf(u, a, a) - scifort_log_two - log(sigma)
            end if
        end if
    end function rdist_logpdf

    pure elemental function rdist_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: y
        real(dp) :: a, mu, sigma, u, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= -1.0_dp) then
                y = 0.0_dp
            else if (z >= 1.0_dp) then
                y = 1.0_dp
            else
                a = 0.5_dp * c
                u = 0.5_dp * (z + 1.0_dp)
                y = beta_cdf(u, a, a)
            end if
        end if
    end function rdist_cdf

    pure elemental function rdist_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: y
        real(dp) :: a, mu, sigma, u, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= -1.0_dp) then
                y = 1.0_dp
            else if (z >= 1.0_dp) then
                y = 0.0_dp
            else
                a = 0.5_dp * c
                u = 0.5_dp * (z + 1.0_dp)
                y = beta_sf(u, a, a)
            end if
        end if
    end function rdist_sf

    pure elemental function rdist_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: y
        real(dp) :: a, mu, sigma, u, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= -1.0_dp) then
                y = negative_infinity(x)
            else if (z >= 1.0_dp) then
                y = 0.0_dp
            else
                a = 0.5_dp * c
                u = 0.5_dp * (z + 1.0_dp)
                y = beta_logcdf(u, a, a)
            end if
        end if
    end function rdist_logcdf

    pure elemental function rdist_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: y
        real(dp) :: a, mu, sigma, u, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= -1.0_dp) then
                y = 0.0_dp
            else if (z >= 1.0_dp) then
                y = negative_infinity(x)
            else
                a = 0.5_dp * c
                u = 0.5_dp * (z + 1.0_dp)
                y = beta_logsf(u, a, a)
            end if
        end if
    end function rdist_logsf

    pure elemental function rdist_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: x
        real(dp) :: a, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            a = 0.5_dp * c
            z = 2.0_dp * beta_ppf(p, a, a) - 1.0_dp
            x = mu + sigma * z
        end if
    end function rdist_ppf

    pure elemental function rdist_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: x
        real(dp) :: a, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            a = 0.5_dp * c
            z = 2.0_dp * beta_isf(p, a, a) - 1.0_dp
            x = mu + sigma * z
        end if
    end function rdist_isf

    pure elemental logical function valid_shape(c) result(ok)
        real(dp), intent(in) :: c !! candidate R-distribution shape parameter
        ok = ieee_is_finite(c) .and. c > 0.0_dp
    end function valid_shape

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

end module scifort_rdist
