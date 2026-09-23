! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Moyal distribution. In standardized coordinates z = (x - loc) / scale,
! f(z) = exp(-(z + exp(-z))/2) / sqrt(2*pi).
! This matches scipy.stats.moyal.

module scifort_moyal
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi, scifort_log_two, &
        scifort_pi, scifort_sqrt_two
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_ppf
    implicit none
    private

    public :: moyal_cdf
    public :: moyal_isf
    public :: moyal_logcdf
    public :: moyal_logpdf
    public :: moyal_logsf
    public :: moyal_pdf
    public :: moyal_ppf
    public :: moyal_sf

contains

    pure elemental function moyal_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: e
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
        if (-z > log(huge(1.0_dp))) then
            y = negative_infinity(x)
        else
            e = exp(-z)
            y = -0.5_dp * (z + e) - scifort_log_sqrt_two_pi - log(sigma)
        end if
    end function moyal_logpdf

    pure elemental function moyal_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logp

        logp = moyal_logpdf(x, loc, scale)
        if (ieee_is_nan(logp)) then
            y = logp
        else if (logp == negative_infinity(logp)) then
            y = 0.0_dp
        else
            y = exp(logp)
        end if
    end function moyal_pdf

    pure elemental function moyal_cdf(x, loc, scale) result(y)
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
        if (-z > log(huge(1.0_dp))) then
            y = 0.0_dp
        else
            t = exp(-0.5_dp * z)
            y = erfc(t / scifort_sqrt_two)
        end if
    end function moyal_cdf

    pure elemental function moyal_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
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
        if (-z > log(huge(1.0_dp))) then
            y = 1.0_dp
        else
            t = exp(-0.5_dp * z)
            y = erf(t / scifort_sqrt_two)
        end if
    end function moyal_sf

    pure elemental function moyal_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: inv_t2
        real(dp) :: mu
        real(dp) :: series
        real(dp) :: sf_value
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: t2
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
        if (z > 0.0_dp) then
            sf_value = moyal_sf(x, mu, sigma)
            y = log1p_safe(-sf_value)
        else if (-z > log(huge(1.0_dp))) then
            y = negative_infinity(x)
        else
            t = exp(-0.5_dp * z)
            if (t <= 20.0_dp) then
                y = log(erfc(t / scifort_sqrt_two))
            else
                t2 = exp(-z)
                inv_t2 = exp(z)
                series = 1.0_dp - inv_t2 * (1.0_dp - inv_t2 * &
                    (3.0_dp - inv_t2 * (15.0_dp - inv_t2 * &
                    (105.0_dp - inv_t2 * (945.0_dp - 10395.0_dp * inv_t2)))))
                y = scifort_log_two - 0.5_dp * t2 + 0.5_dp * z - &
                    scifort_log_sqrt_two_pi + log(series)
            end if
        end if
    end function moyal_logcdf

    pure elemental function moyal_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sf_value
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
        if (z > 40.0_dp) then
            y = 0.5_dp * (scifort_log_two - log(scifort_pi)) - 0.5_dp * z
        else
            sf_value = moyal_sf(x, mu, sigma)
            if (sf_value >= 1.0_dp) then
                y = 0.0_dp
            else if (sf_value <= 0.0_dp) then
                y = negative_infinity(x)
            else
                y = log(sf_value)
            end if
        end if
    end function moyal_logsf

    pure elemental function moyal_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = negative_infinity(p)
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
        else
            q = 0.5_dp * p
            if (q == 0.0_dp) then
                z = lower_quantile_from_log(log(p))
            else
                z = -2.0_dp * log(-normal_ppf(q))
            end if
            x = mu + sigma * z
        end if
    end function moyal_ppf

    pure elemental function moyal_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: log_c
        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = positive_infinity(p)
        else if (p >= 1.0_dp) then
            x = negative_infinity(p)
        else if (p < 1.0e-8_dp) then
            log_c = 0.5_dp * (scifort_log_two - log(scifort_pi))
            z = -2.0_dp * (log(p) - log_c)
            x = mu + sigma * z
        else
            q = 0.5_dp * (1.0_dp - p)
            t = -normal_ppf(q)
            z = -2.0_dp * log(t)
            x = mu + sigma * z
        end if
    end function moyal_isf

    pure elemental function lower_quantile_from_log(logp) result(z)
        real(dp), intent(in) :: logp !! natural logarithm of a tiny lower-tail probability
        real(dp) :: z

        integer :: i
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid

        lo = -50.0_dp
        hi = 0.0_dp
        do i = 1, 100
            mid = 0.5_dp * (lo + hi)
            if (standard_logcdf(mid) < logp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        z = 0.5_dp * (lo + hi)
    end function lower_quantile_from_log

    pure elemental function standard_logcdf(z) result(y)
        real(dp), intent(in) :: z !! standardized Moyal variate
        real(dp) :: y

        real(dp) :: inv_t2
        real(dp) :: series
        real(dp) :: t
        real(dp) :: t2

        if (-z > log(huge(1.0_dp))) then
            y = negative_infinity(z)
            return
        end if
        t = exp(-0.5_dp * z)
        if (t <= 20.0_dp) then
            y = log(erfc(t / scifort_sqrt_two))
        else
            t2 = exp(-z)
            inv_t2 = exp(z)
            series = 1.0_dp - inv_t2 * (1.0_dp - inv_t2 * &
                (3.0_dp - inv_t2 * (15.0_dp - inv_t2 * &
                (105.0_dp - inv_t2 * (945.0_dp - 10395.0_dp * inv_t2)))))
            y = scifort_log_two - 0.5_dp * t2 + 0.5_dp * z - &
                scifort_log_sqrt_two_pi + log(series)
        end if
    end function standard_logcdf

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

end module scifort_moyal
