! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Real modified Bessel K for positive arguments.  The implementation uses the
! defining integral K_nu(x) = integral_0^inf exp(-x cosh(t)) cosh(nu*t) dt,
! evaluated after logarithmic scaling with composite 16-point Gauss-Legendre
! quadrature.  Large x with modest order uses the standard asymptotic series.
module scifort_bessel_k
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: log_besselk
    public :: log_besselk_scaled
    public :: log_besselk_order_derivative
    public :: besselk_log_derivative_x

    real(dp), parameter :: gl_x(8) = [ &
        0.0950125098376374401853193354250_dp, &
        0.281603550779258913230460501460_dp, &
        0.458016777657227386342419442984_dp, &
        0.617876244402643748446671764049_dp, &
        0.755404408355003033895101194847_dp, &
        0.865631202387831743880467897712_dp, &
        0.944575023073232576077988415535_dp, &
        0.989400934991649932596154173450_dp]
    real(dp), parameter :: gl_w(8) = [ &
        0.189450610455068496285396723208_dp, &
        0.182603415044923588866763667969_dp, &
        0.169156519395002538189312079030_dp, &
        0.149595988816576732081501730547_dp, &
        0.124628971255533872052476282192_dp, &
        0.0951585116824927848099251076022_dp, &
        0.0622535239386478928628438369944_dp, &
        0.0271524594117540948517805724560_dp]

contains

    pure elemental function log_besselk(nu, x) result(y)
        real(dp), intent(in) :: nu !! real order of K_nu
        real(dp), intent(in) :: x !! positive real argument
        real(dp) :: y
        real(dp) :: anu, tpeak, tmax, hmax, integral

        if (ieee_is_nan(nu) .or. ieee_is_nan(x) .or. .not. ieee_is_finite(nu) .or. x < 0.0_dp) then
            y = quiet_nan(x + nu)
            return
        end if
        if (x == 0.0_dp) then
            y = positive_infinity(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
            return
        end if

        anu = abs(nu)
        if (x >= max(12.0_dp, 0.5_dp * anu * anu + 8.0_dp)) then
            y = log_besselk_asymptotic(anu, x)
            return
        end if

        if (anu > 0.0_dp) then
            tpeak = asinh(anu / x)
        else
            tpeak = 0.0_dp
        end if
        hmax = max(log_integrand(0.0_dp, anu, x), log_integrand(tpeak, anu, x))
        tmax = max(12.0_dp, tpeak + 10.0_dp, log(2.0_dp / min(x, 1.0_dp)) + 10.0_dp)
        do while (log_integrand(tmax, anu, x) - hmax > -45.0_dp .and. tmax < 100.0_dp)
            tmax = min(100.0_dp, tmax + 5.0_dp)
        end do
        integral = integrate_scaled(anu, x, hmax, tmax)
        if (integral <= 0.0_dp) then
            y = quiet_nan(x)
        else
            y = hmax + log(integral)
        end if
    end function log_besselk

    pure elemental function log_besselk_scaled(nu, x) result(y)
        real(dp), intent(in) :: nu !! real order of exp(x)*K_nu(x)
        real(dp), intent(in) :: x !! positive real argument
        real(dp) :: y
        real(dp) :: anu, tpeak, tmax, hmax, integral

        if (ieee_is_nan(nu) .or. ieee_is_nan(x) .or. .not. ieee_is_finite(nu) .or. x < 0.0_dp) then
            y = quiet_nan(x + nu)
            return
        end if
        if (x == 0.0_dp) then
            y = positive_infinity(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
            return
        end if

        anu = abs(nu)
        if (x >= max(12.0_dp, 0.5_dp * anu * anu + 8.0_dp)) then
            y = log_besselk_asymptotic_scaled(anu, x)
            return
        end if

        if (anu > 0.0_dp) then
            tpeak = asinh(anu / x)
        else
            tpeak = 0.0_dp
        end if
        hmax = max(log_integrand_scaled(0.0_dp, anu, x), log_integrand_scaled(tpeak, anu, x))
        tmax = max(12.0_dp, tpeak + 10.0_dp, log(2.0_dp / min(x, 1.0_dp)) + 10.0_dp)
        do while (log_integrand_scaled(tmax, anu, x) - hmax > -45.0_dp .and. tmax < 100.0_dp)
            tmax = min(100.0_dp, tmax + 5.0_dp)
        end do
        integral = integrate_scaled_besselk(anu, x, hmax, tmax)
        if (integral <= 0.0_dp) then
            y = quiet_nan(x)
        else
            y = hmax + log(integral)
        end if
    end function log_besselk_scaled

    pure elemental function log_besselk_order_derivative(nu, x) result(y)
        real(dp), intent(in) :: nu !! real order at which d log(K_nu(x))/d nu is evaluated
        real(dp), intent(in) :: x !! positive real argument
        real(dp) :: y, anu, tpeak, tmax, hmax, denominator, numerator

        if (ieee_is_nan(nu) .or. ieee_is_nan(x) .or. .not. ieee_is_finite(nu) .or. &
            x <= 0.0_dp .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x + nu)
            return
        end if
        if (nu == 0.0_dp) then
            y = 0.0_dp
            return
        end if

        anu = abs(nu)
        tpeak = asinh(anu / x)
        hmax = max(log_integrand(0.0_dp, anu, x), log_integrand(tpeak, anu, x))
        tmax = max(12.0_dp, tpeak + 10.0_dp, log(2.0_dp / min(x, 1.0_dp)) + 10.0_dp)
        do while (log_integrand(tmax, anu, x) - hmax > -45.0_dp .and. tmax < 100.0_dp)
            tmax = min(100.0_dp, tmax + 5.0_dp)
        end do
        denominator = integrate_scaled(anu, x, hmax, tmax)
        numerator = integrate_scaled_order_moment(anu, x, hmax, tmax)
        y = sign(1.0_dp, nu) * numerator / denominator
    end function log_besselk_order_derivative

    pure elemental function besselk_log_derivative_x(nu, x) result(y)
        real(dp), intent(in) :: nu !! real order at which d log(K_nu(x))/dx is evaluated
        real(dp), intent(in) :: x !! positive real argument
        real(dp) :: y, lp, lm, l0
        l0 = log_besselk(nu, x)
        lp = log_besselk(nu + 1.0_dp, x)
        lm = log_besselk(nu - 1.0_dp, x)
        y = -0.5_dp * (exp(lp - l0) + exp(lm - l0))
    end function besselk_log_derivative_x

    pure elemental function log_besselk_asymptotic(nu, x) result(y)
        real(dp), intent(in) :: nu !! nonnegative real order
        real(dp), intent(in) :: x !! positive large argument
        real(dp) :: y
        real(dp) :: mu, term, sumv
        integer :: k

        mu = 4.0_dp * nu * nu
        term = 1.0_dp
        sumv = 1.0_dp
        do k = 1, 20
            term = term * (mu - real((2*k - 1)**2, dp)) / (real(k, dp) * 8.0_dp * x)
            if (abs(term) > abs(sumv) .and. k > 6) exit
            sumv = sumv + term
            if (abs(term) <= 2.0_dp * epsilon(sumv) * abs(sumv)) exit
        end do
        if (sumv <= 0.0_dp) then
            y = log_besselk_integral(nu, x)
        else
            y = -x + 0.5_dp * (log(scifort_pi) - log(2.0_dp * x)) + log(sumv)
        end if
    end function log_besselk_asymptotic

    pure elemental function log_besselk_asymptotic_scaled(nu, x) result(y)
        real(dp), intent(in) :: nu !! nonnegative real order
        real(dp), intent(in) :: x !! positive large argument
        real(dp) :: y
        real(dp) :: mu, term, sumv
        integer :: k

        mu = 4.0_dp * nu * nu
        term = 1.0_dp
        sumv = 1.0_dp
        do k = 1, 20
            term = term * (mu - real((2*k - 1)**2, dp)) / (real(k, dp) * 8.0_dp * x)
            if (abs(term) > abs(sumv) .and. k > 6) exit
            sumv = sumv + term
            if (abs(term) <= 2.0_dp * epsilon(sumv) * abs(sumv)) exit
        end do
        if (sumv <= 0.0_dp) then
            y = log_besselk(nu, x) + x
        else
            y = 0.5_dp * (log(scifort_pi) - log(2.0_dp * x)) + log(sumv)
        end if
    end function log_besselk_asymptotic_scaled

    pure elemental function log_besselk_integral(nu, x) result(y)
        real(dp), intent(in) :: nu !! nonnegative real order
        real(dp), intent(in) :: x !! positive real argument
        real(dp) :: y, tpeak, tmax, hmax, integral
        if (nu > 0.0_dp) then
            tpeak = asinh(nu / x)
        else
            tpeak = 0.0_dp
        end if
        hmax = max(log_integrand(0.0_dp, nu, x), log_integrand(tpeak, nu, x))
        tmax = max(12.0_dp, tpeak + 10.0_dp)
        integral = integrate_scaled(nu, x, hmax, tmax)
        y = hmax + log(integral)
    end function log_besselk_integral

    pure function integrate_scaled(nu, x, hmax, tmax) result(value)
        real(dp), intent(in) :: nu !! nonnegative Bessel order
        real(dp), intent(in) :: x !! positive Bessel argument
        real(dp), intent(in) :: hmax !! logarithmic scale for the integrand
        real(dp), intent(in) :: tmax !! finite truncation of the integral
        real(dp) :: value
        real(dp) :: a, b, mid, half, t1, t2
        real(dp) :: segment_width
        integer :: nseg, j, i

        segment_width = min(0.5_dp, 4.0_dp / sqrt(max(1.0_dp, hypot(x, nu))))
        nseg = max(1, ceiling(tmax / segment_width))
        value = 0.0_dp
        do j = 1, nseg
            a = tmax * real(j - 1, dp) / real(nseg, dp)
            b = tmax * real(j, dp) / real(nseg, dp)
            mid = 0.5_dp * (a + b)
            half = 0.5_dp * (b - a)
            do i = 1, 8
                t1 = mid - half * gl_x(i)
                t2 = mid + half * gl_x(i)
                value = value + half * gl_w(i) * ( &
                    exp(log_integrand(t1, nu, x) - hmax) + &
                    exp(log_integrand(t2, nu, x) - hmax))
            end do
        end do
    end function integrate_scaled

    pure function integrate_scaled_besselk(nu, x, hmax, tmax) result(value)
        real(dp), intent(in) :: nu !! nonnegative Bessel order
        real(dp), intent(in) :: x !! positive Bessel argument
        real(dp), intent(in) :: hmax !! logarithmic scale for the scaled integrand
        real(dp), intent(in) :: tmax !! finite truncation of the integral
        real(dp) :: value
        real(dp) :: a, b, mid, half, t1, t2
        real(dp) :: segment_width
        integer :: nseg, j, i

        segment_width = min(0.5_dp, 4.0_dp / sqrt(max(1.0_dp, hypot(x, nu))))
        nseg = max(1, ceiling(tmax / segment_width))
        value = 0.0_dp
        do j = 1, nseg
            a = tmax * real(j - 1, dp) / real(nseg, dp)
            b = tmax * real(j, dp) / real(nseg, dp)
            mid = 0.5_dp * (a + b)
            half = 0.5_dp * (b - a)
            do i = 1, 8
                t1 = mid - half * gl_x(i)
                t2 = mid + half * gl_x(i)
                value = value + half * gl_w(i) * ( &
                    exp(log_integrand_scaled(t1, nu, x) - hmax) + &
                    exp(log_integrand_scaled(t2, nu, x) - hmax))
            end do
        end do
    end function integrate_scaled_besselk

    pure function integrate_scaled_order_moment(nu, x, hmax, tmax) result(value)
        real(dp), intent(in) :: nu !! nonnegative Bessel order
        real(dp), intent(in) :: x !! positive Bessel argument
        real(dp), intent(in) :: hmax !! logarithmic scale for the integrand
        real(dp), intent(in) :: tmax !! finite truncation of the integral
        real(dp) :: value
        real(dp) :: a, b, mid, half, t1, t2, f1, f2
        real(dp) :: segment_width
        integer :: nseg, j, i

        segment_width = min(0.5_dp, 4.0_dp / sqrt(max(1.0_dp, hypot(x, nu))))
        nseg = max(1, ceiling(tmax / segment_width))
        value = 0.0_dp
        do j = 1, nseg
            a = tmax * real(j - 1, dp) / real(nseg, dp)
            b = tmax * real(j, dp) / real(nseg, dp)
            mid = 0.5_dp * (a + b)
            half = 0.5_dp * (b - a)
            do i = 1, 8
                t1 = mid - half * gl_x(i)
                t2 = mid + half * gl_x(i)
                f1 = exp(log_integrand(t1, nu, x) - hmax) * t1 * tanh(nu * t1)
                f2 = exp(log_integrand(t2, nu, x) - hmax) * t2 * tanh(nu * t2)
                value = value + half * gl_w(i) * (f1 + f2)
            end do
        end do
    end function integrate_scaled_order_moment

    pure elemental function log_integrand(t, nu, x) result(y)
        real(dp), intent(in) :: t !! nonnegative integration coordinate
        real(dp), intent(in) :: nu !! nonnegative Bessel order
        real(dp), intent(in) :: x !! positive Bessel argument
        real(dp) :: y
        y = -x * cosh(t) + log_cosh(nu * t)
    end function log_integrand

    pure elemental function log_integrand_scaled(t, nu, x) result(y)
        real(dp), intent(in) :: t !! nonnegative integration coordinate
        real(dp), intent(in) :: nu !! nonnegative integration order
        real(dp), intent(in) :: x !! positive Bessel argument
        real(dp) :: y
        y = -x * (cosh(t) - 1.0_dp) + log_cosh(nu * t)
    end function log_integrand_scaled

    pure elemental function log_cosh(x) result(y)
        real(dp), intent(in) :: x !! real argument
        real(dp) :: y, ax
        ax = abs(x)
        y = ax + log(0.5_dp * (1.0_dp + exp(-2.0_dp * ax)))
    end function log_cosh

end module scifort_bessel_k
