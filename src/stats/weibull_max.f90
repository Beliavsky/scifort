! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Weibull maximum distribution with shape c > 0. This is the reflection of
! scipy.stats.weibull_min and matches scipy.stats.weibull_max.

module scifort_weibull_max
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, &
        valid_loc_scale
    use scifort_weibull, only : weibull_cdf, weibull_isf, weibull_logcdf, &
        weibull_logpdf, weibull_logsf, weibull_ppf, weibull_sf
    implicit none
    private

    public :: weibull_max_cdf
    public :: weibull_max_isf
    public :: weibull_max_logcdf
    public :: weibull_max_logpdf
    public :: weibull_max_logsf
    public :: weibull_max_pdf
    public :: weibull_max_ppf
    public :: weibull_max_sf

contains

    pure elemental function weibull_max_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy

        logy = weibull_max_logpdf(x, c, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function weibull_max_pdf

    pure elemental function weibull_max_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            y = weibull_logpdf(-z, c) - log(sigma)
        end if
    end function weibull_max_logpdf

    pure elemental function weibull_max_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x >= mu) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            y = weibull_sf(-z, c)
        end if
    end function weibull_max_cdf

    pure elemental function weibull_max_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x >= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            y = weibull_cdf(-z, c)
        end if
    end function weibull_max_sf

    pure elemental function weibull_max_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x >= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            y = weibull_logsf(-z, c)
        end if
    end function weibull_max_logcdf

    pure elemental function weibull_max_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x >= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            y = weibull_logcdf(-z, c)
        end if
    end function weibull_max_logsf

    pure elemental function weibull_max_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = mu
        else
            z = -weibull_isf(p, c)
            x = affine(mu, sigma, z, p)
        end if
    end function weibull_max_ppf

    pure elemental function weibull_max_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(c) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else if (p == 0.0_dp) then
            x = mu
        else
            z = -weibull_ppf(p, c)
            x = affine(mu, sigma, z, p)
        end if
    end function weibull_max_isf

    pure elemental logical function valid_shape(c) result(ok)
        real(dp), intent(in) :: c !! candidate Weibull shape parameter
        ok = ieee_is_finite(c) .and. c > 0.0_dp
    end function valid_shape

    pure elemental function affine(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! standardized quantile
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: x

        if (.not. ieee_is_finite(z)) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else if (abs(z) > (huge(1.0_dp) - abs(mu)) / sigma) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else
            x = mu + sigma * z
        end if
    end function affine

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

end module scifort_weibull_max
