! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Inverse-Gaussian distribution in the scipy.stats.invgauss parameterization.
! With shape mu_shape > 0 and z = (x - loc) / scale > 0,
! f(x) = exp(-(z-mu_shape)**2/(2*mu_shape**2*z)) /
!        (sqrt(2*pi*z**3) * scale).

module scifort_invgauss
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_logcdf
    implicit none
    private

    public :: invgauss_cdf
    public :: invgauss_isf
    public :: invgauss_logcdf
    public :: invgauss_logpdf
    public :: invgauss_logsf
    public :: invgauss_pdf
    public :: invgauss_ppf
    public :: invgauss_sf

contains

    pure elemental function invgauss_pdf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logy

        logy = invgauss_logpdf(x, mu_shape, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function invgauss_pdf

    pure elemental function invgauss_logpdf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(mu_shape))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                y = -scifort_log_sqrt_two_pi - 1.5_dp * log(z) - log(sigma) - &
                    0.5_dp * (z - mu_shape) * (z - mu_shape) / &
                    (mu_shape * mu_shape * z)
            end if
        end if
    end function invgauss_logpdf

    pure elemental function invgauss_cdf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logp

        logp = invgauss_logcdf(x, mu_shape, loc, scale)
        if (ieee_is_nan(logp)) then
            y = logp
        else if (logp == negative_infinity(logp)) then
            y = 0.0_dp
        else if (logp >= 0.0_dp) then
            y = 1.0_dp
        else
            y = exp(logp)
        end if
    end function invgauss_cdf

    pure elemental function invgauss_sf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logp

        logp = invgauss_logsf(x, mu_shape, loc, scale)
        if (ieee_is_nan(logp)) then
            y = logp
        else if (logp == negative_infinity(logp)) then
            y = 0.0_dp
        else if (logp >= 0.0_dp) then
            y = 1.0_dp
        else
            y = exp(logp)
        end if
    end function invgauss_sf

    pure elemental function invgauss_logcdf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(mu_shape))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            y = standard_logcdf(z, mu_shape)
        end if
    end function invgauss_logcdf

    pure elemental function invgauss_logsf(x, mu_shape, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(mu_shape))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            y = standard_logsf(z, mu_shape)
        end if
    end function invgauss_logsf

    pure elemental function invgauss_ppf(p, mu_shape, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(mu_shape))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            z = standard_quantile(p, mu_shape, .true.)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function invgauss_ppf

    pure elemental function invgauss_isf(p, mu_shape, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(mu_shape))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            z = standard_quantile(p, mu_shape, .false.)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function invgauss_isf

    pure elemental function standard_logcdf(z, mu_shape) result(y)
        real(dp), intent(in) :: z !! positive standardized value
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp) :: y

        real(dp) :: a
        real(dp) :: b
        real(dp) :: root_z
        real(dp) :: term1
        real(dp) :: term2

        root_z = sqrt(z)
        a = root_z / mu_shape - 1.0_dp / root_z
        b = -root_z / mu_shape - 1.0_dp / root_z
        term1 = normal_logcdf(a)
        term2 = 2.0_dp / mu_shape + normal_logcdf(b)
        y = min(0.0_dp, logaddexp_pair(term1, term2))
    end function standard_logcdf

    pure elemental function standard_logsf(z, mu_shape) result(y)
        real(dp), intent(in) :: z !! positive standardized value
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        real(dp) :: y

        real(dp) :: a
        real(dp) :: b
        real(dp) :: root_z
        real(dp) :: term1
        real(dp) :: term2

        root_z = sqrt(z)
        a = root_z / mu_shape - 1.0_dp / root_z
        b = -root_z / mu_shape - 1.0_dp / root_z
        term1 = normal_logcdf(-a)
        term2 = 2.0_dp / mu_shape + normal_logcdf(b)
        y = min(0.0_dp, logdiffexp_pair(term1, term2))
    end function standard_logsf

    pure function standard_quantile(p, mu_shape, lower_tail) result(z)
        real(dp), intent(in) :: p !! requested tail probability in (0, 1)
        real(dp), intent(in) :: mu_shape !! positive SciPy shape parameter
        logical, intent(in) :: lower_tail !! true for CDF inversion, false for SF inversion
        real(dp) :: z

        integer :: iteration
        real(dp) :: lower
        real(dp) :: logp
        real(dp) :: logvalue
        real(dp) :: midpoint
        real(dp) :: upper

        logp = log(p)
        lower = 0.0_dp
        upper = max(1.0_dp, mu_shape)

        do
            if (lower_tail) then
                logvalue = standard_logcdf(upper, mu_shape)
                if (logvalue >= logp) exit
            else
                logvalue = standard_logsf(upper, mu_shape)
                if (logvalue <= logp) exit
            end if
            lower = upper
            if (upper > 0.5_dp * huge(1.0_dp)) then
                z = positive_infinity(p)
                return
            end if
            upper = 2.0_dp * upper
        end do

        do iteration = 1, 220
            midpoint = lower + 0.5_dp * (upper - lower)
            if (midpoint == lower .or. midpoint == upper) exit
            if (lower_tail) then
                logvalue = standard_logcdf(midpoint, mu_shape)
                if (logvalue < logp) then
                    lower = midpoint
                else
                    upper = midpoint
                end if
            else
                logvalue = standard_logsf(midpoint, mu_shape)
                if (logvalue > logp) then
                    lower = midpoint
                else
                    upper = midpoint
                end if
            end if
        end do
        z = lower + 0.5_dp * (upper - lower)
    end function standard_quantile

    pure elemental function logaddexp_pair(a, b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y

        real(dp) :: hi
        real(dp) :: lo

        hi = max(a, b)
        lo = min(a, b)
        if (hi == negative_infinity(hi)) then
            y = hi
        else
            y = hi + log1p_safe(exp(lo - hi))
        end if
    end function logaddexp_pair

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

    pure elemental logical function valid_shape(mu_shape) result(valid)
        real(dp), intent(in) :: mu_shape !! shape parameter to check

        valid = ieee_is_finite(mu_shape) .and. mu_shape > 0.0_dp
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

end module scifort_invgauss
