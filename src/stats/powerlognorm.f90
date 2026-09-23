! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Power log-normal distribution matching scipy.stats.powerlognorm.
! Standard density is c/(z*s)*phi(log(z)/s)*Phi(-log(z)/s)**(c-1).
module scifort_powerlognorm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_logcdf, normal_logpdf, normal_ppf
    implicit none
    private
    public :: powerlognorm_pdf, powerlognorm_logpdf, powerlognorm_cdf, powerlognorm_sf
    public :: powerlognorm_logcdf, powerlognorm_logsf, powerlognorm_ppf, powerlognorm_isf
contains
    pure elemental function powerlognorm_logpdf(x, c, s, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite power shape parameter
        real(dp), intent(in) :: s !! positive finite lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lz, mu, sigma, w, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c) .and. valid_shape(s))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            lz = log(z)
            w = lz / s
            y = log(c) - lz - log(s) + normal_logpdf(w) + &
                (c - 1.0_dp) * normal_logcdf(-w) - log(sigma)
        end if
    end function powerlognorm_logpdf

    pure elemental function powerlognorm_pdf(x, c, s, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite power shape parameter
        real(dp), intent(in) :: s !! positive finite lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = powerlognorm_logpdf(x, c, s, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function powerlognorm_pdf

    pure elemental function powerlognorm_logsf(x, c, s, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite power shape parameter
        real(dp), intent(in) :: s !! positive finite lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, w, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c) .and. valid_shape(s))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            w = log(z) / s
            y = c * normal_logcdf(-w)
        end if
    end function powerlognorm_logsf

    pure elemental function powerlognorm_sf(x, c, s, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability
        real(dp), intent(in) :: c !! positive finite power shape parameter
        real(dp), intent(in) :: s !! positive finite lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = powerlognorm_logsf(x, c, s, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function powerlognorm_sf

    pure elemental function powerlognorm_cdf(x, c, s, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability
        real(dp), intent(in) :: c !! positive finite power shape parameter
        real(dp), intent(in) :: s !! positive finite lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ls
        ls = powerlognorm_logsf(x, c, s, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == negative_infinity(ls)) then
            y = 1.0_dp
        else
            y = -expm1_safe(ls)
        end if
    end function powerlognorm_cdf

    pure elemental function powerlognorm_logcdf(x, c, s, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite power shape parameter
        real(dp), intent(in) :: s !! positive finite lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ls
        ls = powerlognorm_logsf(x, c, s, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == 0.0_dp) then
            y = negative_infinity(x)
        else if (ls == negative_infinity(ls)) then
            y = 0.0_dp
        else
            y = log(-expm1_safe(ls))
        end if
    end function powerlognorm_logcdf

    pure elemental function powerlognorm_ppf(p, c, s, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive finite power shape parameter
        real(dp), intent(in) :: s !! positive finite lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, q, sigma, w, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c) .and. valid_shape(s)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            q = exp(log1p_safe(-p) / c)
            w = normal_ppf(q)
            z = exp_quantile(-s * w, p)
            x = affine(mu, sigma, z, p)
        end if
    end function powerlognorm_ppf

    pure elemental function powerlognorm_isf(p, c, s, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive finite power shape parameter
        real(dp), intent(in) :: s !! positive finite lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, q, sigma, w, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c) .and. valid_shape(s)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            q = exp(log(p) / c)
            w = normal_ppf(q)
            z = exp_quantile(-s * w, p)
            x = affine(mu, sigma, z, p)
        end if
    end function powerlognorm_isf

    pure elemental function exp_quantile(logz, seed) result(z)
        real(dp), intent(in) :: logz !! logarithm of positive standardized quantile
        real(dp), intent(in) :: seed !! seed used to construct positive infinity
        real(dp) :: z
        if (logz > log(huge(1.0_dp))) then
            z = positive_infinity(seed)
        else if (logz < log(tiny(1.0_dp))) then
            z = 0.0_dp
        else
            z = exp(logz)
        end if
    end function exp_quantile

    pure elemental function valid_shape(c) result(ok)
        real(dp), intent(in) :: c !! candidate shape parameter
        logical :: ok
        ok = ieee_is_finite(c) .and. c > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! requested location
        real(dp), intent(in), optional :: scale !! requested scale
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure elemental function affine(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! standardized quantile
        real(dp), intent(in) :: seed !! seed used to construct positive infinity
        real(dp) :: x
        if (.not. ieee_is_finite(z)) then
            x = positive_infinity(seed)
        else if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine
end module scifort_powerlognorm
