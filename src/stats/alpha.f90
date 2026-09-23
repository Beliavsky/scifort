! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Alpha distribution with shape a > 0, location loc, and scale > 0.
! In standardized coordinates z = (x-loc)/scale > 0,
! F(z) = Phi(a - 1/z) / Phi(a). This matches scipy.stats.alpha.

module scifort_alpha
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_logcdf, normal_logpdf, normal_ppf
    implicit none
    private

    public :: alpha_cdf
    public :: alpha_isf
    public :: alpha_logcdf
    public :: alpha_logpdf
    public :: alpha_logsf
    public :: alpha_pdf
    public :: alpha_ppf
    public :: alpha_sf

contains

    pure elemental function alpha_pdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy

        logy = alpha_logpdf(x, a, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function alpha_pdf

    pure elemental function alpha_logpdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            w = a - 1.0_dp / z
            y = normal_logpdf(w) - normal_logcdf(a) - 2.0_dp * log(z) - log(sigma)
        end if
    end function alpha_logpdf

    pure elemental function alpha_logcdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            y = min(0.0_dp, normal_logcdf(a - 1.0_dp / z) - normal_logcdf(a))
        end if
    end function alpha_logcdf

    pure elemental function alpha_logsf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: la
        real(dp) :: lw
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            la = normal_logcdf(a)
            lw = normal_logcdf(a - 1.0_dp / z)
            y = logdiffexp_pair(la, lw) - la
        end if
    end function alpha_logsf

    pure elemental function alpha_cdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy

        logy = alpha_logcdf(x, a, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function alpha_cdf

    pure elemental function alpha_sf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logy

        logy = alpha_logsf(x, a, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function alpha_sf

    pure elemental function alpha_ppf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: denom
        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            q = exp(log(p) + normal_logcdf(a))
            denom = a - normal_ppf(q)
            x = affine_positive(mu, sigma, 1.0_dp / denom, p)
        end if
    end function alpha_ppf

    pure elemental function alpha_isf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: denom
        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            q = exp(log1p_safe(-p) + normal_logcdf(a))
            denom = a - normal_ppf(q)
            x = affine_positive(mu, sigma, 1.0_dp / denom, p)
        end if
    end function alpha_isf

    pure elemental function logdiffexp_pair(a, b) result(y)
        real(dp), intent(in) :: a !! log of the positive leading term
        real(dp), intent(in) :: b !! log of the nonnegative subtracted term
        real(dp) :: y

        if (b == negative_infinity(b)) then
            y = a
        else if (b >= a) then
            y = negative_infinity(a)
        else
            y = a + log1p_safe(-exp(b - a))
        end if
    end function logdiffexp_pair

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to construct positive infinity if needed
        real(dp) :: x

        if (.not. ieee_is_finite(z) .or. z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive

    pure elemental logical function valid_shape(a) result(valid)
        real(dp), intent(in) :: a !! shape parameter to check
        valid = ieee_is_finite(a) .and. a > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location
        real(dp), intent(in), optional :: scale !! optional scale
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_alpha
