! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Elementary real-valued special functions with SciPy-style semantics where
! practical. The implementation is independent and built from Fortran
! intrinsics plus existing SciFort kernels.

module scifort_special_elementary
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi, scifort_sqrt_two
    use scifort_digamma, only : digamma_positive
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan
    use scifort_normal, only : normal_cdf, normal_logcdf, normal_ppf
    implicit none
    private

    public :: betaln
    public :: boxcox
    public :: boxcox1p
    public :: cosm1
    public :: digamma
    public :: entr
    public :: expit
    public :: exprel
    public :: gammaln
    public :: i0e
    public :: i1e
    public :: inv_boxcox
    public :: inv_boxcox1p
    public :: kl_div
    public :: log_expit
    public :: log_ndtr
    public :: logit
    public :: ndtr
    public :: ndtri
    public :: psi
    public :: rel_entr
    public :: xlog1py
    public :: xlogy
    public :: erf
    public :: erfc
    public :: erfcinv
    public :: erfinv

    interface erf
        module procedure erf_real
    end interface erf

    interface erfc
        module procedure erfc_real
    end interface erfc

contains

    pure elemental function erf_real(x) result(y)
        real(dp), intent(in) :: x !! real argument of the error function
        real(dp) :: y

        intrinsic :: erf

        y = erf(x)
    end function erf_real

    pure elemental function erfc_real(x) result(y)
        real(dp), intent(in) :: x !! real argument of the complementary error function
        real(dp) :: y

        intrinsic :: erfc

        y = erfc(x)
    end function erfc_real

    pure elemental function erfinv(y_arg) result(x)
        real(dp), intent(in) :: y_arg !! value in [-1, 1] whose inverse error function is required
        real(dp) :: x

        real(dp) :: p

        if (ieee_is_nan(y_arg) .or. y_arg < -1.0_dp .or. y_arg > 1.0_dp) then
            x = quiet_nan(y_arg)
        else if (y_arg == -1.0_dp) then
            x = negative_infinity(y_arg)
        else if (y_arg == 1.0_dp) then
            x = positive_infinity(y_arg)
        else if (y_arg < 0.0_dp) then
            p = 0.5_dp * (1.0_dp + y_arg)
            x = normal_ppf(p) / scifort_sqrt_two
        else if (y_arg > 0.0_dp) then
            p = 0.5_dp * (1.0_dp - y_arg)
            x = -normal_ppf(p) / scifort_sqrt_two
        else
            x = 0.0_dp
        end if
    end function erfinv

    pure elemental function erfcinv(y_arg) result(x)
        real(dp), intent(in) :: y_arg !! value in [0, 2] whose inverse erfc is required
        real(dp) :: x

        real(dp) :: p

        if (ieee_is_nan(y_arg) .or. y_arg < 0.0_dp .or. y_arg > 2.0_dp) then
            x = quiet_nan(y_arg)
        else if (y_arg == 0.0_dp) then
            x = positive_infinity(y_arg)
        else if (y_arg == 2.0_dp) then
            x = negative_infinity(y_arg)
        else if (y_arg <= 1.0_dp) then
            p = 0.5_dp * y_arg
            x = -normal_ppf(p) / scifort_sqrt_two
        else
            p = 0.5_dp * (2.0_dp - y_arg)
            x = normal_ppf(p) / scifort_sqrt_two
        end if
    end function erfcinv

    pure elemental function gammaln(x) result(y)
        real(dp), intent(in) :: x !! real argument; returns log(abs(Gamma(x)))
        real(dp) :: y

        real(dp) :: s

        if (ieee_is_nan(x)) then
            y = x
        else if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = positive_infinity(x)
            else
                y = negative_infinity(x)
            end if
        else if (x > 0.0_dp) then
            y = log_gamma(x)
        else if (x == 0.0_dp .or. x == anint(x)) then
            y = positive_infinity(x)
        else
            s = sinpi_reduced(x)
            if (s == 0.0_dp) then
                y = positive_infinity(x)
            else
                y = log(scifort_pi) - log(abs(s)) - log_gamma(1.0_dp - x)
            end if
        end if
    end function gammaln

    pure elemental function betaln(a, b) result(y)
        real(dp), intent(in) :: a !! first beta argument; positive finite values are supported
        real(dp), intent(in) :: b !! second beta argument; positive finite values are supported
        real(dp) :: y

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            y = quiet_nan(a)
        else if (a > 0.0_dp .and. b > 0.0_dp .and. &
                ieee_is_finite(a) .and. ieee_is_finite(b)) then
            y = log_beta(a, b)
        else if (a == 0.0_dp .or. b == 0.0_dp) then
            y = positive_infinity(a + b)
        else
            ! SciFort's public betaln currently promises the positive-shape
            ! statistical domain only. Do not manufacture values outside it.
            y = quiet_nan(a)
        end if
    end function betaln

    pure elemental function digamma(x) result(y)
        real(dp), intent(in) :: x !! real argument of the digamma function
        real(dp) :: y

        real(dp) :: c
        real(dp) :: s

        if (ieee_is_nan(x)) then
            y = x
        else if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = positive_infinity(x)
            else
                y = quiet_nan(x)
            end if
        else if (x > 0.0_dp) then
            y = digamma_positive(x)
        else if (x == 0.0_dp) then
            y = negative_infinity(x)
        else if (x == anint(x)) then
            y = quiet_nan(x)
        else
            s = sinpi_reduced(x)
            c = cospi_reduced(x)
            y = digamma_positive(1.0_dp - x) - scifort_pi * c / s
        end if
    end function digamma

    pure elemental function i0e(x) result(y)
        real(dp), intent(in) :: x !! real argument; returns exp(-abs(x))*I0(x)
        real(dp) :: y
        real(dp) :: ax, term, sumv, u, product
        integer :: k

        if (ieee_is_nan(x)) then
            y = x
            return
        end if
        ax = abs(x)
        if (.not. ieee_is_finite(ax)) then
            y = 0.0_dp
        else if (ax <= 100.0_dp) then
            u = 0.25_dp * ax * ax
            term = 1.0_dp
            sumv = 1.0_dp
            do k = 1, 10000
                term = term * u / real(k * k, dp)
                sumv = sumv + term
                if (abs(term) <= epsilon(sumv) * abs(sumv)) exit
            end do
            y = exp(-ax) * sumv
        else
            sumv = 1.0_dp
            product = 1.0_dp
            do k = 1, 12
                product = product * (-real((2 * k - 1)**2, dp)) / &
                    (real(k, dp) * 8.0_dp * ax)
                sumv = sumv + (-1.0_dp)**k * product
            end do
            y = sumv / sqrt(2.0_dp * scifort_pi * ax)
        end if
    end function i0e

    pure elemental function i1e(x) result(y)
        real(dp), intent(in) :: x !! real argument; returns exp(-abs(x))*I1(x)
        real(dp) :: y
        real(dp) :: ax, term, sumv, u, product
        integer :: k

        if (ieee_is_nan(x)) then
            y = x
            return
        end if
        ax = abs(x)
        if (.not. ieee_is_finite(ax)) then
            y = 0.0_dp
        else if (ax == 0.0_dp) then
            y = 0.0_dp
        else if (ax <= 100.0_dp) then
            u = 0.25_dp * ax * ax
            term = 0.5_dp * ax
            sumv = term
            do k = 1, 10000
                term = term * u / (real(k, dp) * real(k + 1, dp))
                sumv = sumv + term
                if (abs(term) <= epsilon(sumv) * abs(sumv)) exit
            end do
            y = exp(-ax) * sumv
        else
            sumv = 1.0_dp
            product = 1.0_dp
            do k = 1, 12
                product = product * (4.0_dp - real((2 * k - 1)**2, dp)) / &
                    (real(k, dp) * 8.0_dp * ax)
                sumv = sumv + (-1.0_dp)**k * product
            end do
            y = sumv / sqrt(2.0_dp * scifort_pi * ax)
        end if
        if (x < 0.0_dp) y = -y
    end function i1e

    pure elemental function psi(x) result(y)
        real(dp), intent(in) :: x !! real argument; alias of digamma
        real(dp) :: y

        y = digamma(x)
    end function psi

    pure elemental function ndtr(x) result(y)
        real(dp), intent(in) :: x !! standard-normal variate
        real(dp) :: y

        y = normal_cdf(x)
    end function ndtr

    pure elemental function log_ndtr(x) result(y)
        real(dp), intent(in) :: x !! standard-normal variate
        real(dp) :: y

        y = normal_logcdf(x)
    end function log_ndtr

    pure elemental function ndtri(p) result(x)
        real(dp), intent(in) :: p !! probability in [0, 1]
        real(dp) :: x

        x = normal_ppf(p)
    end function ndtri

    pure elemental function expit(x) result(y)
        real(dp), intent(in) :: x !! logistic-transform argument
        real(dp) :: y

        real(dp) :: e

        if (ieee_is_nan(x)) then
            y = x
        else if (x >= 0.0_dp) then
            y = 1.0_dp / (1.0_dp + exp(-x))
        else
            e = exp(x)
            y = e / (1.0_dp + e)
        end if
    end function expit

    pure elemental function logit(p) result(x)
        real(dp), intent(in) :: p !! probability in [0, 1]
        real(dp) :: x

        if (ieee_is_nan(p)) then
            x = p
        else if (p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            x = log(p) - log1p_safe(-p)
        end if
    end function logit

    pure elemental function log_expit(x) result(y)
        real(dp), intent(in) :: x !! logistic-transform argument
        real(dp) :: y

        if (ieee_is_nan(x)) then
            y = x
        else if (x >= 0.0_dp) then
            y = -log1p_safe(exp(-x))
        else
            y = x - log1p_safe(exp(x))
        end if
    end function log_expit

    pure elemental function xlogy(x, y_arg) result(value)
        real(dp), intent(in) :: x !! multiplier
        real(dp), intent(in) :: y_arg !! logarithm argument
        real(dp) :: value

        if (ieee_is_nan(x) .or. ieee_is_nan(y_arg)) then
            value = quiet_nan(x)
        else if (x == 0.0_dp) then
            value = 0.0_dp
        else if (y_arg < 0.0_dp) then
            value = quiet_nan(y_arg)
        else if (y_arg == 0.0_dp) then
            if (x > 0.0_dp) then
                value = negative_infinity(x)
            else
                value = positive_infinity(x)
            end if
        else
            value = x * log(y_arg)
        end if
    end function xlogy

    pure elemental function xlog1py(x, y_arg) result(value)
        real(dp), intent(in) :: x !! multiplier
        real(dp), intent(in) :: y_arg !! argument in log(1 + y_arg), >= -1
        real(dp) :: value

        if (ieee_is_nan(x) .or. ieee_is_nan(y_arg)) then
            value = quiet_nan(x)
        else if (x == 0.0_dp) then
            value = 0.0_dp
        else if (y_arg < -1.0_dp) then
            value = quiet_nan(y_arg)
        else if (y_arg == -1.0_dp) then
            if (x > 0.0_dp) then
                value = negative_infinity(x)
            else
                value = positive_infinity(x)
            end if
        else
            value = x * log1p_safe(y_arg)
        end if
    end function xlog1py

    pure elemental function entr(x) result(y)
        real(dp), intent(in) :: x !! nonnegative argument of -x*log(x)
        real(dp) :: y

        if (ieee_is_nan(x)) then
            y = x
        else if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else if (x == 0.0_dp) then
            y = 0.0_dp
        else
            y = -x * log(x)
        end if
    end function entr

    pure elemental function rel_entr(x, y_arg) result(value)
        real(dp), intent(in) :: x !! first relative-entropy argument
        real(dp), intent(in) :: y_arg !! second relative-entropy argument
        real(dp) :: value

        if (ieee_is_nan(x) .or. ieee_is_nan(y_arg)) then
            value = quiet_nan(x)
        else if (x == 0.0_dp .and. y_arg >= 0.0_dp) then
            value = 0.0_dp
        else if (x <= 0.0_dp .or. y_arg <= 0.0_dp) then
            value = positive_infinity(x)
        else if (.not. ieee_is_finite(x) .and. .not. ieee_is_finite(y_arg)) then
            value = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            value = positive_infinity(x)
        else if (.not. ieee_is_finite(y_arg)) then
            value = negative_infinity(y_arg)
        else
            value = x * (log(x) - log(y_arg))
        end if
    end function rel_entr

    pure elemental function kl_div(x, y_arg) result(value)
        real(dp), intent(in) :: x !! first Kullback-Leibler argument
        real(dp), intent(in) :: y_arg !! second Kullback-Leibler argument
        real(dp) :: value

        if (ieee_is_nan(x) .or. ieee_is_nan(y_arg)) then
            value = quiet_nan(x)
        else if (x == 0.0_dp .and. y_arg >= 0.0_dp) then
            value = y_arg
        else if (x <= 0.0_dp .or. y_arg <= 0.0_dp) then
            value = positive_infinity(x)
        else if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y_arg)) then
            value = quiet_nan(x)
        else
            value = x * (log(x) - log(y_arg)) - x + y_arg
        end if
    end function kl_div

    pure elemental function exprel(x) result(y)
        real(dp), intent(in) :: x !! argument of (exp(x) - 1) / x
        real(dp) :: y

        real(dp), parameter :: log_huge = log(huge(1.0_dp))

        if (ieee_is_nan(x)) then
            y = x
        else if (x == 0.0_dp) then
            y = 1.0_dp
        else if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = positive_infinity(x)
            else
                y = 0.0_dp
            end if
        else if (x > log_huge) then
            y = positive_infinity(x)
        else
            y = expm1_safe(x) / x
        end if
    end function exprel

    pure elemental function cosm1(x) result(y)
        real(dp), intent(in) :: x !! argument of cos(x) - 1
        real(dp) :: y

        if (.not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else
            y = -2.0_dp * sin(0.5_dp * x) ** 2
        end if
    end function cosm1

    pure elemental function boxcox(x, lambda) result(y)
        real(dp), intent(in) :: x !! nonnegative Box-Cox input
        real(dp), intent(in) :: lambda !! Box-Cox power parameter
        real(dp) :: y

        if (ieee_is_nan(x) .or. ieee_is_nan(lambda) .or. x < 0.0_dp) then
            y = quiet_nan(x)
        else if (lambda == 0.0_dp) then
            if (x == 0.0_dp) then
                y = negative_infinity(x)
            else
                y = log(x)
            end if
        else if (x == 0.0_dp) then
            if (lambda > 0.0_dp) then
                y = -1.0_dp / lambda
            else
                y = negative_infinity(x)
            end if
        else
            y = expm1_safe(lambda * log(x)) / lambda
        end if
    end function boxcox

    pure elemental function boxcox1p(x, lambda) result(y)
        real(dp), intent(in) :: x !! Box-Cox input shifted by one, x >= -1
        real(dp), intent(in) :: lambda !! Box-Cox power parameter
        real(dp) :: y

        if (ieee_is_nan(x) .or. ieee_is_nan(lambda) .or. x < -1.0_dp) then
            y = quiet_nan(x)
        else if (lambda == 0.0_dp) then
            y = log1p_safe(x)
        else if (x == -1.0_dp) then
            if (lambda > 0.0_dp) then
                y = -1.0_dp / lambda
            else
                y = negative_infinity(x)
            end if
        else
            y = expm1_safe(lambda * log1p_safe(x)) / lambda
        end if
    end function boxcox1p

    pure elemental function inv_boxcox(y_arg, lambda) result(x)
        real(dp), intent(in) :: y_arg !! transformed Box-Cox value
        real(dp), intent(in) :: lambda !! Box-Cox power parameter
        real(dp) :: x

        real(dp) :: exponent_value
        real(dp) :: t

        if (ieee_is_nan(y_arg) .or. ieee_is_nan(lambda)) then
            x = quiet_nan(y_arg)
            return
        end if
        if (lambda == 0.0_dp) then
            x = safe_exp(y_arg)
            return
        end if

        t = lambda * y_arg
        if (t < -1.0_dp) then
            x = quiet_nan(t)
        else if (t == -1.0_dp) then
            if (lambda > 0.0_dp) then
                x = 0.0_dp
            else
                x = positive_infinity(t)
            end if
        else
            exponent_value = log1p_safe(t) / lambda
            x = safe_exp(exponent_value)
        end if
    end function inv_boxcox

    pure elemental function inv_boxcox1p(y_arg, lambda) result(x)
        real(dp), intent(in) :: y_arg !! transformed Box-Cox value
        real(dp), intent(in) :: lambda !! Box-Cox power parameter
        real(dp) :: x

        real(dp) :: exponent_value
        real(dp) :: t

        if (ieee_is_nan(y_arg) .or. ieee_is_nan(lambda)) then
            x = quiet_nan(y_arg)
            return
        end if
        if (lambda == 0.0_dp) then
            x = safe_expm1(y_arg)
            return
        end if

        t = lambda * y_arg
        if (t < -1.0_dp) then
            x = quiet_nan(t)
        else if (t == -1.0_dp) then
            if (lambda > 0.0_dp) then
                x = -1.0_dp
            else
                x = positive_infinity(t)
            end if
        else
            exponent_value = log1p_safe(t) / lambda
            x = safe_expm1(exponent_value)
        end if
    end function inv_boxcox1p

    pure elemental function sinpi_reduced(x) result(y)
        real(dp), intent(in) :: x !! finite real argument of sin(pi*x)
        real(dp) :: y

        real(dp) :: r

        r = modulo(x, 2.0_dp)
        if (r > 1.0_dp) r = r - 2.0_dp
        y = sin(scifort_pi * r)
    end function sinpi_reduced

    pure elemental function cospi_reduced(x) result(y)
        real(dp), intent(in) :: x !! finite real argument of cos(pi*x)
        real(dp) :: y

        real(dp) :: r

        r = modulo(x, 2.0_dp)
        if (r > 1.0_dp) r = r - 2.0_dp
        y = cos(scifort_pi * r)
    end function cospi_reduced

    pure elemental function safe_exp(x) result(y)
        real(dp), intent(in) :: x !! real exponent, possibly infinite
        real(dp) :: y

        real(dp), parameter :: log_huge = log(huge(1.0_dp))

        if (ieee_is_nan(x)) then
            y = x
        else if (x > log_huge) then
            y = positive_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = positive_infinity(x)
            else
                y = 0.0_dp
            end if
        else
            y = exp(x)
        end if
    end function safe_exp

    pure elemental function safe_expm1(x) result(y)
        real(dp), intent(in) :: x !! real exponent in exp(x) - 1
        real(dp) :: y

        real(dp), parameter :: log_huge = log(huge(1.0_dp))

        if (ieee_is_nan(x)) then
            y = x
        else if (x > log_huge) then
            y = positive_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = positive_infinity(x)
            else
                y = -1.0_dp
            end if
        else
            y = expm1_safe(x)
        end if
    end function safe_expm1

end module scifort_special_elementary
