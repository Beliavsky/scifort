! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Elementary real-valued special functions with SciPy-style semantics where
! practical. The implementation is independent and built from Fortran
! intrinsics plus existing SciFort kernels.

module scifort_special_elementary
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi, scifort_pi, scifort_sqrt_two
    use scifort_digamma, only : digamma_positive
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta
    use scifort_math, only : expm1_safe, log1p_safe, log1pexp, negative_infinity, &
        positive_infinity, quiet_nan
    use scifort_normal, only : normal_cdf, normal_logcdf, normal_ppf
    implicit none
    private

    public :: agm
    public :: betaln
    public :: boxcox
    public :: boxcox1p
    public :: cbrt
    public :: cosdg
    public :: cotdg
    public :: diric
    public :: exp10
    public :: exp2
    public :: expm1
    public :: cosm1
    public :: digamma
    public :: entr
    public :: expit
    public :: exprel
    public :: gammaln
    public :: huber
    public :: i0e
    public :: i1e
    public :: inv_boxcox
    public :: inv_boxcox1p
    public :: kl_div
    public :: log_expit
    public :: log1p
    public :: log_ndtr
    public :: logit
    public :: ndtr
    public :: ndtri_exp
    public :: ndtri
    public :: psi
    public :: pseudo_huber
    public :: powm1
    public :: radian
    public :: round
    public :: sinc
    public :: sindg
    public :: softplus
    public :: tandg
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
        real(dp), intent(in) :: a !! first real beta-function argument
        real(dp), intent(in) :: b !! second real beta-function argument
        real(dp) :: y

        real(dp) :: m
        real(dp) :: p

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            y = quiet_nan(a)
            return
        end if

        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            if (.not. ieee_is_finite(a) .and. .not. ieee_is_finite(b)) then
                if (a > 0.0_dp .and. b > 0.0_dp) then
                    y = quiet_nan(a)
                else
                    y = positive_infinity(a)
                end if
            else if (.not. ieee_is_finite(a)) then
                if (a < 0.0_dp .or. b <= 0.0_dp) then
                    y = positive_infinity(a)
                else
                    y = negative_infinity(a)
                end if
            else
                if (b < 0.0_dp .or. a <= 0.0_dp) then
                    y = positive_infinity(b)
                else
                    y = negative_infinity(b)
                end if
            end if
            return
        end if

        ! Preserve the well-conditioned positive-shape implementation used by
        ! SciFort's probability distributions.
        if (a > 0.0_dp .and. b > 0.0_dp) then
            y = log_beta(a, b)
            return
        end if

        ! At a nonpositive-integer pole, a positive integer in the other
        ! argument can cancel the gamma pole exactly.  Use the finite-product
        ! identity rather than an infinity-minus-infinity gamma-log expression.
        if (is_nonpositive_integer_beta(a)) then
            if (is_positive_integer_beta(b) .and. b <= -a) then
                p = -a
                m = b
                y = gammaln(m) - gammaln(p + 1.0_dp) + gammaln(p - m + 1.0_dp)
            else
                y = positive_infinity(a)
            end if
            return
        end if
        if (is_nonpositive_integer_beta(b)) then
            if (is_positive_integer_beta(a) .and. a <= -b) then
                p = -b
                m = a
                y = gammaln(m) - gammaln(p + 1.0_dp) + gammaln(p - m + 1.0_dp)
            else
                y = positive_infinity(b)
            end if
            return
        end if

        if (is_nonpositive_integer_beta(a + b)) then
            y = negative_infinity(a + b)
        else
            y = gammaln(a) + gammaln(b) - gammaln(a + b)
        end if
    end function betaln

    pure elemental logical function is_nonpositive_integer_beta(x) result(is_integer)
        real(dp), intent(in) :: x !! finite value to classify for beta poles

        is_integer = x <= 0.0_dp .and. x == anint(x)
    end function is_nonpositive_integer_beta

    pure elemental logical function is_positive_integer_beta(x) result(is_integer)
        real(dp), intent(in) :: x !! finite value to classify for beta cancellation

        is_integer = x > 0.0_dp .and. x == anint(x)
    end function is_positive_integer_beta

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
            if (sign(1.0_dp, x) < 0.0_dp) then
                y = positive_infinity(x)
            else
                y = negative_infinity(x)
            end if
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


    pure elemental function ndtri_exp(y) result(x)
        real(dp), intent(in) :: y !! logarithm of a probability in (-infinity, 0]
        real(dp) :: x

        real(dp), parameter :: log_half = -0.69314718055994530942_dp
        integer :: iteration
        real(dp) :: derivative
        real(dp) :: f
        real(dp) :: inv_x2
        real(dp) :: logcdf
        real(dp) :: logpdf
        real(dp) :: q
        real(dp) :: series
        real(dp) :: step

        if (ieee_is_nan(y) .or. y > 0.0_dp) then
            x = quiet_nan(y)
            return
        else if (y == 0.0_dp) then
            x = positive_infinity(y)
            return
        else if (.not. ieee_is_finite(y)) then
            x = negative_infinity(y)
            return
        end if

        if (y > log_half) then
            q = -expm1_safe(y)
            x = -normal_ppf(q)
            return
        else if (y > -700.0_dp) then
            x = normal_ppf(exp(y))
            return
        end if

        x = -scifort_sqrt_two * sqrt(-y)
        do iteration = 1, 12
            if (abs(x) > sqrt(huge(1.0_dp))) then
                inv_x2 = 0.0_dp
            else
                inv_x2 = 1.0_dp / (x * x)
            end if
            series = 1.0_dp - inv_x2 * (1.0_dp - inv_x2 * &
                (3.0_dp - inv_x2 * (15.0_dp - inv_x2 * &
                (105.0_dp - inv_x2 * (945.0_dp - 10395.0_dp * inv_x2)))))
            logpdf = (-0.5_dp * x) * x - scifort_log_sqrt_two_pi
            logcdf = logpdf - log(-x) + log(series)
            f = logcdf - y
            derivative = exp(logpdf - logcdf)
            step = f / derivative
            x = x - step
            if (abs(step) <= 4.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x))) exit
        end do
    end function ndtri_exp

    pure elemental function agm(a, b) result(value)
        real(dp), intent(in) :: a !! first real arithmetic-geometric mean argument
        real(dp), intent(in) :: b !! second real arithmetic-geometric mean argument
        real(dp) :: value

        real(dp) :: an
        real(dp) :: bn
        real(dp) :: next_a
        real(dp) :: next_b
        real(dp) :: sign_value
        integer :: iteration

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            value = quiet_nan(a + b)
            return
        end if
        if ((a < 0.0_dp .and. b > 0.0_dp) .or. &
                (a > 0.0_dp .and. b < 0.0_dp)) then
            value = quiet_nan(a + b)
            return
        end if
        if (a == 0.0_dp .or. b == 0.0_dp) then
            if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
                value = quiet_nan(a + b)
            else
                value = 0.0_dp
            end if
            return
        end if

        sign_value = 1.0_dp
        an = a
        bn = b
        if (a < 0.0_dp .and. b < 0.0_dp) then
            sign_value = -1.0_dp
            an = -a
            bn = -b
        end if
        if (.not. ieee_is_finite(an) .or. .not. ieee_is_finite(bn)) then
            value = sign_value * positive_infinity(an + bn)
            return
        end if

        do iteration = 1, 64
            next_a = 0.5_dp * an + 0.5_dp * bn
            next_b = sqrt(an) * sqrt(bn)
            if (abs(next_a - next_b) <= 2.0_dp * epsilon(1.0_dp) * next_a) then
                value = sign_value * next_a
                return
            end if
            an = next_a
            bn = next_b
        end do
        value = sign_value * 0.5_dp * (an + bn)
    end function agm

    pure elemental function huber(delta, r) result(value)
        real(dp), intent(in) :: delta !! nonnegative transition between quadratic and linear loss
        real(dp), intent(in) :: r !! residual
        real(dp) :: value

        real(dp) :: ar
        real(dp) :: term

        if (delta < 0.0_dp) then
            value = positive_infinity(delta)
            return
        end if
        if (ieee_is_nan(delta) .or. ieee_is_nan(r)) then
            value = quiet_nan(delta + r)
            return
        end if

        ar = abs(r)
        if (ar <= delta) then
            if (.not. ieee_is_finite(ar)) then
                value = positive_infinity(ar)
            else if (ar > sqrt(2.0_dp) * sqrt(huge(1.0_dp))) then
                value = positive_infinity(ar)
            else
                value = (0.5_dp * ar) * ar
            end if
        else
            if (delta == 0.0_dp .and. .not. ieee_is_finite(ar)) then
                value = quiet_nan(ar)
                return
            end if
            term = ar - 0.5_dp * delta
            if (delta == 0.0_dp) then
                value = delta * term
            else if (term > huge(1.0_dp) / abs(delta)) then
                value = positive_infinity(term)
            else
                value = delta * term
            end if
        end if
    end function huber

    pure elemental function pseudo_huber(delta, r) result(value)
        real(dp), intent(in) :: delta !! nonnegative soft transition scale
        real(dp), intent(in) :: r !! residual
        real(dp) :: value

        real(dp) :: ar
        real(dp) :: factor
        real(dp) :: q
        real(dp) :: term

        if (delta < 0.0_dp) then
            value = positive_infinity(delta)
            return
        end if
        if (delta == 0.0_dp .or. r == 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (ieee_is_nan(delta) .or. ieee_is_nan(r)) then
            value = quiet_nan(delta + r)
            return
        end if
        if (.not. ieee_is_finite(delta)) then
            value = quiet_nan(delta)
            return
        end if
        if (.not. ieee_is_finite(r)) then
            value = positive_infinity(r)
            return
        end if

        ar = abs(r)
        if (ar <= delta) then
            q = ar / delta
            factor = sqrt(1.0_dp + q * q) + 1.0_dp
            term = ar / factor
            if (term > 1.0_dp) then
                if (ar > huge(1.0_dp) / term) then
                    value = positive_infinity(ar)
                else
                    value = term * ar
                end if
            else
                value = term * ar
            end if
        else
            q = delta / ar
            factor = sqrt(1.0_dp + q * q) - q
            term = delta * factor
            if (term > 1.0_dp) then
                if (ar > huge(1.0_dp) / term) then
                    value = positive_infinity(ar)
                else
                    value = term * ar
                end if
            else
                value = term * ar
            end if
        end if
    end function pseudo_huber

    pure elemental function cbrt(x) result(y)
        real(dp), intent(in) :: x !! real value whose real cube root is required
        real(dp) :: y

        real(dp) :: ax
        integer :: iteration

        if (ieee_is_nan(x) .or. x == 0.0_dp .or. .not. ieee_is_finite(x)) then
            y = x
            return
        end if

        ax = abs(x)
        y = exp(log(ax) / 3.0_dp)
        do iteration = 1, 2
            y = (2.0_dp * y + ax / (y * y)) / 3.0_dp
        end do
        y = sign(y, x)
    end function cbrt

    pure elemental function exp10(x) result(y)
        real(dp), intent(in) :: x !! real exponent in 10**x
        real(dp) :: y

        real(dp), parameter :: log_ten = log(10.0_dp)
        real(dp), parameter :: overflow_x = log(huge(1.0_dp)) / log_ten

        if (ieee_is_nan(x)) then
            y = x
        else if (x > overflow_x) then
            y = positive_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = positive_infinity(x)
            else
                y = 0.0_dp
            end if
        else
            y = 10.0_dp ** x
        end if
    end function exp10

    pure elemental function exp2(x) result(y)
        real(dp), intent(in) :: x !! real exponent in 2**x
        real(dp) :: y

        real(dp), parameter :: log_two = log(2.0_dp)
        real(dp), parameter :: overflow_x = log(huge(1.0_dp)) / log_two

        if (ieee_is_nan(x)) then
            y = x
        else if (x > overflow_x) then
            y = positive_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = positive_infinity(x)
            else
                y = 0.0_dp
            end if
        else
            y = 2.0_dp ** x
        end if
    end function exp2

    pure elemental function expm1(x) result(y)
        real(dp), intent(in) :: x !! real exponent in exp(x)-1
        real(dp) :: y

        y = safe_expm1(x)
    end function expm1

    pure elemental function log1p(x) result(y)
        real(dp), intent(in) :: x !! real argument in log(1+x), with x >= -1
        real(dp) :: y

        y = log1p_safe(x)
    end function log1p

    pure elemental function powm1(x, exponent_arg) result(y)
        real(dp), intent(in) :: x !! real base in x**exponent_arg - 1
        real(dp), intent(in) :: exponent_arg !! real exponent
        real(dp) :: y

        logical :: odd_integer
        real(dp) :: log_abs_x
        real(dp) :: t

        if (exponent_arg == 0.0_dp .or. x == 1.0_dp) then
            y = 0.0_dp
            return
        end if
        if (ieee_is_nan(x) .or. ieee_is_nan(exponent_arg)) then
            y = quiet_nan(x + exponent_arg)
            return
        end if

        if (x == 0.0_dp) then
            if (exponent_arg > 0.0_dp) then
                y = -1.0_dp
            else
                y = positive_infinity(exponent_arg)
            end if
            return
        end if

        if (.not. ieee_is_finite(exponent_arg)) then
            if (x < 0.0_dp) then
                y = quiet_nan(x)
            else if (x > 1.0_dp) then
                if (exponent_arg > 0.0_dp) then
                    y = positive_infinity(exponent_arg)
                else
                    y = -1.0_dp
                end if
            else
                if (exponent_arg > 0.0_dp) then
                    y = -1.0_dp
                else
                    y = positive_infinity(exponent_arg)
                end if
            end if
            return
        end if

        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                if (exponent_arg > 0.0_dp) then
                    y = positive_infinity(x)
                else
                    y = -1.0_dp
                end if
            else if (exponent_arg /= anint(exponent_arg)) then
                y = quiet_nan(x)
            else if (exponent_arg < 0.0_dp) then
                y = -1.0_dp
            else if (modulo(abs(exponent_arg), 2.0_dp) == 1.0_dp) then
                y = negative_infinity(x)
            else
                y = positive_infinity(x)
            end if
            return
        end if

        if (x > 0.0_dp) then
            log_abs_x = log(x)
            if (abs(log_abs_x) > 1.0_dp) then
                if (abs(exponent_arg) > huge(1.0_dp) / abs(log_abs_x)) then
                    if (sign(1.0_dp, exponent_arg) == sign(1.0_dp, log_abs_x)) then
                        y = positive_infinity(x)
                    else
                        y = -1.0_dp
                    end if
                    return
                end if
            end if
            t = exponent_arg * log_abs_x
            if (abs(t) <= 0.5_dp) then
                y = safe_expm1(t)
            else if (t > log(huge(1.0_dp))) then
                y = positive_infinity(x)
            else
                y = x ** exponent_arg - 1.0_dp
            end if
            return
        end if

        if (exponent_arg /= anint(exponent_arg)) then
            y = quiet_nan(x)
            return
        end if

        odd_integer = modulo(abs(exponent_arg), 2.0_dp) == 1.0_dp
        log_abs_x = log(abs(x))
        if (abs(log_abs_x) > 1.0_dp) then
            if (abs(exponent_arg) > huge(1.0_dp) / abs(log_abs_x)) then
                if (sign(1.0_dp, exponent_arg) == sign(1.0_dp, log_abs_x)) then
                    if (odd_integer) then
                        y = negative_infinity(x)
                    else
                        y = positive_infinity(x)
                    end if
                else
                    y = -1.0_dp
                end if
                return
            end if
        end if

        t = exponent_arg * log_abs_x
        if (.not. odd_integer) then
            if (abs(t) <= 0.5_dp) then
                y = safe_expm1(t)
            else if (t > log(huge(1.0_dp))) then
                y = positive_infinity(x)
            else
                y = abs(x) ** exponent_arg - 1.0_dp
            end if
        else if (t > log(huge(1.0_dp))) then
            y = negative_infinity(x)
        else
            y = -(abs(x) ** exponent_arg) - 1.0_dp
        end if
    end function powm1

    pure elemental function radian(d, m, s) result(value)
        real(dp), intent(in) :: d !! degrees component
        real(dp), intent(in) :: m !! arc-minutes component
        real(dp), intent(in) :: s !! arc-seconds component
        real(dp) :: value

        real(dp), parameter :: degree_scale = scifort_pi / 180.0_dp
        real(dp), parameter :: minute_scale = scifort_pi / 10800.0_dp
        real(dp), parameter :: second_scale = scifort_pi / 648000.0_dp
        logical :: has_negative_infinity
        logical :: has_positive_infinity

        if (ieee_is_nan(d) .or. ieee_is_nan(m) .or. ieee_is_nan(s)) then
            value = quiet_nan(d + m + s)
            return
        end if

        has_positive_infinity = (.not. ieee_is_finite(d) .and. d > 0.0_dp) .or. &
            (.not. ieee_is_finite(m) .and. m > 0.0_dp) .or. &
            (.not. ieee_is_finite(s) .and. s > 0.0_dp)
        has_negative_infinity = (.not. ieee_is_finite(d) .and. d < 0.0_dp) .or. &
            (.not. ieee_is_finite(m) .and. m < 0.0_dp) .or. &
            (.not. ieee_is_finite(s) .and. s < 0.0_dp)
        if (has_positive_infinity .and. has_negative_infinity) then
            value = quiet_nan(d)
        else if (has_positive_infinity) then
            value = positive_infinity(d)
        else if (has_negative_infinity) then
            value = negative_infinity(d)
        else
            value = d * degree_scale + m * minute_scale + s * second_scale
        end if
    end function radian

    pure elemental function round(x) result(y)
        real(dp), intent(in) :: x !! real value rounded to the nearest even integer-valued real
        real(dp) :: y

        real(dp), parameter :: all_integral_threshold = 4503599627370496.0_dp
        real(dp) :: ax
        real(dp) :: fraction_part
        real(dp) :: lower_integer
        real(dp) :: nearest_integer

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x) .or. x == 0.0_dp) then
            y = x
            return
        end if
        ax = abs(x)
        if (ax >= all_integral_threshold) then
            y = x
            return
        end if

        lower_integer = aint(ax)
        fraction_part = ax - lower_integer
        if (fraction_part < 0.5_dp) then
            nearest_integer = lower_integer
        else if (fraction_part > 0.5_dp) then
            nearest_integer = lower_integer + 1.0_dp
        else if (modulo(lower_integer, 2.0_dp) == 0.0_dp) then
            nearest_integer = lower_integer
        else
            nearest_integer = lower_integer + 1.0_dp
        end if

        if (nearest_integer == 0.0_dp) then
            y = 0.0_dp
        else
            y = sign(nearest_integer, x)
        end if
    end function round

    pure elemental function sinc(x) result(y)
        real(dp), intent(in) :: x !! argument of normalized sin(pi*x)/(pi*x)
        real(dp) :: y

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else if (x == 0.0_dp) then
            y = 1.0_dp
        else if (x == anint(x)) then
            y = 0.0_dp
        else
            y = (sinpi_reduced(x) / x) / scifort_pi
        end if
    end function sinc

    pure elemental function sindg(x) result(y)
        real(dp), intent(in) :: x !! angle in degrees
        real(dp) :: y

        real(dp) :: r

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
            return
        end if
        r = degree_remainder(x)
        if (r == 0.0_dp .or. abs(r) == 180.0_dp) then
            y = 0.0_dp
        else if (r == 90.0_dp) then
            y = 1.0_dp
        else if (r == -90.0_dp) then
            y = -1.0_dp
        else
            y = sin(scifort_pi * r / 180.0_dp)
        end if
    end function sindg

    pure elemental function cosdg(x) result(y)
        real(dp), intent(in) :: x !! angle in degrees
        real(dp) :: y

        real(dp) :: r

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
            return
        end if
        r = degree_remainder(x)
        if (r == 0.0_dp) then
            y = 1.0_dp
        else if (abs(r) == 180.0_dp) then
            y = -1.0_dp
        else if (abs(r) == 90.0_dp) then
            y = 0.0_dp
        else
            y = cos(scifort_pi * r / 180.0_dp)
        end if
    end function cosdg

    pure elemental function tandg(x) result(y)
        real(dp), intent(in) :: x !! angle in degrees
        real(dp) :: y

        real(dp) :: r

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
            return
        end if
        r = degree_remainder(x)
        if (r == 0.0_dp .or. abs(r) == 180.0_dp) then
            y = 0.0_dp
        else if (abs(r) == 90.0_dp) then
            y = positive_infinity(x)
        else
            y = tan(scifort_pi * r / 180.0_dp)
        end if
    end function tandg

    pure elemental function cotdg(x) result(y)
        real(dp), intent(in) :: x !! angle in degrees
        real(dp) :: y

        real(dp) :: r

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
            return
        end if
        r = degree_remainder(x)
        if (r == 0.0_dp .or. abs(r) == 180.0_dp) then
            y = positive_infinity(x)
        else if (abs(r) == 90.0_dp) then
            y = 0.0_dp
        else
            y = tan(scifort_pi * r / 180.0_dp)
            if (abs(y) < 1.0_dp / huge(1.0_dp)) then
                y = positive_infinity(x)
            else
                y = 1.0_dp / y
            end if
        end if
    end function cotdg

    pure elemental function softplus(x) result(y)
        real(dp), intent(in) :: x !! real argument of log(1+exp(x))
        real(dp) :: y

        y = log1pexp(x)
    end function softplus

    pure elemental function diric(x, n) result(y)
        real(dp), intent(in) :: x !! real phase argument in radians
        integer, intent(in) :: n !! positive integer Dirichlet-kernel order
        real(dp) :: y

        real(dp) :: half_r
        real(dp) :: r
        real(dp) :: sin_half
        real(dp), parameter :: two_pi = 2.0_dp * scifort_pi
        real(dp), parameter :: four_pi = 4.0_dp * scifort_pi

        if (n <= 0 .or. ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
            return
        end if

        r = modulo(x, four_pi)
        if (r > two_pi) r = r - four_pi
        if (r == 0.0_dp) then
            y = 1.0_dp
            return
        else if (abs(r) == two_pi) then
            if (modulo(n, 2) == 0) then
                y = -1.0_dp
            else
                y = 1.0_dp
            end if
            return
        end if

        half_r = 0.5_dp * r
        sin_half = sin(half_r)
        y = sin(real(n, dp) * half_r) / (real(n, dp) * sin_half)
    end function diric

    pure elemental function degree_remainder(x) result(r)
        real(dp), intent(in) :: x !! finite angle in degrees to reduce modulo 360
        real(dp) :: r

        r = modulo(x, 360.0_dp)
        if (r > 180.0_dp) r = r - 360.0_dp
    end function degree_remainder

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
