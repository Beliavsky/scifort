! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Log-gamma distribution matching scipy.stats.loggamma.
! Standard density: exp(c*z-exp(z))/Gamma(c), c>0, z real.
module scifort_loggamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_gamma, only : gammaincinv, gammainccinv, log_gammainc, log_gammaincc
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private
    public :: loggamma_pdf, loggamma_logpdf, loggamma_cdf, loggamma_sf
    public :: loggamma_logcdf, loggamma_logsf, loggamma_ppf, loggamma_isf
contains
    pure elemental function loggamma_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, ez
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z > log(huge(1.0_dp))) then
                y = negative_infinity(x)
            else
                ez = exp(z)
                y = c * z - ez - log_gamma(c) - log(sigma)
            end if
        end if
    end function loggamma_logpdf

    pure elemental function loggamma_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = loggamma_logpdf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function loggamma_pdf

    pure elemental function loggamma_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, ez
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = negative_infinity(x)
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z < log(tiny(1.0_dp))) then
                y = c * z - log_gamma(c + 1.0_dp)
            else if (z > log(huge(1.0_dp))) then
                y = 0.0_dp
            else
                ez = exp(z)
                y = log_gammainc(c, ez)
            end if
        end if
    end function loggamma_logcdf

    pure elemental function loggamma_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, ez, lc
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 0.0_dp
        else if (x == positive_infinity(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z < log(tiny(1.0_dp))) then
                lc = c * z - log_gamma(c + 1.0_dp)
                y = log1p_safe(-exp(lc))
            else if (z > log(huge(1.0_dp))) then
                y = negative_infinity(x)
            else
                ez = exp(z)
                y = log_gammaincc(c, ez)
            end if
        end if
    end function loggamma_logsf

    pure elemental function loggamma_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of P(X <= x)
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = loggamma_logcdf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function loggamma_cdf

    pure elemental function loggamma_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of P(X > x)
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = loggamma_logsf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function loggamma_sf

    pure elemental function loggamma_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, g, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            g = gammaincinv(c, p)
            if (g > 0.0_dp) then
                z = log(g)
            else
                z = (log(p) + log_gamma(c + 1.0_dp)) / c
            end if
            x = affine(mu, sigma, z, p)
        end if
    end function loggamma_ppf

    pure elemental function loggamma_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, g, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            g = gammainccinv(c, p)
            if (g > 0.0_dp) then
                z = log(g)
            else
                z = (log1p_safe(-p) + log_gamma(c + 1.0_dp)) / c
            end if
            x = affine(mu, sigma, z, p)
        end if
    end function loggamma_isf

    pure elemental function valid_shape(c) result(ok)
        real(dp), intent(in) :: c !! candidate shape parameter
        logical :: ok
        ok = ieee_is_finite(c) .and. c > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        if (present(loc)) mu = loc
        sigma = 1.0_dp
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure elemental function affine(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! standardized quantile
        real(dp), intent(in) :: seed !! value used to form infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z) .or. abs(z) > (huge(1.0_dp) - abs(mu)) / sigma) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else
            x = mu + sigma * z
        end if
    end function affine
end module scifort_loggamma
