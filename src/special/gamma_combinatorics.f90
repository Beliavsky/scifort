! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Real gamma-derived special functions: beta, generalized binomial,
! Pochhammer/rising factorial, and integer-order polygamma.
!
! The implementation is independent of SciPy source code. It uses standard
! gamma-ratio identities, the polygamma/Hurwitz-zeta identity, recurrence, and
! reflection formulas from the NIST DLMF.

module scifort_gamma_combinatorics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_gamma_family, only : gammasgn
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_special_elementary, only : betaln, digamma, gammaln
    use scifort_zeta, only : hurwitz_zeta
    implicit none
    private

    public :: beta
    public :: binom
    public :: poch
    public :: polygamma

contains

    pure elemental function beta(a, b) result(y)
        real(dp), intent(in) :: a !! first real beta-function argument
        real(dp), intent(in) :: b !! second real beta-function argument
        real(dp) :: y

        real(dp) :: log_abs
        real(dp) :: sign_value

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            y = quiet_nan(a)
            return
        end if

        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            y = beta_infinite_argument(a, b)
            return
        end if

        if (is_nonpositive_integer(a)) then
            if (is_positive_integer(b) .and. b <= -a) then
                call beta_integer_cancellation(a, b, log_abs, sign_value)
                y = signed_exp(log_abs, sign_value, a)
            else
                y = positive_infinity(a)
            end if
            return
        end if
        if (is_nonpositive_integer(b)) then
            if (is_positive_integer(a) .and. a <= -b) then
                call beta_integer_cancellation(b, a, log_abs, sign_value)
                y = signed_exp(log_abs, sign_value, b)
            else
                y = positive_infinity(b)
            end if
            return
        end if

        if (is_nonpositive_integer(a + b)) then
            sign_value = gammasgn(a) * gammasgn(b)
            y = sign(0.0_dp, sign_value)
            return
        end if

        log_abs = betaln(a, b)
        sign_value = gammasgn(a) * gammasgn(b) * gammasgn(a + b)
        y = signed_exp(log_abs, sign_value, a)
    end function beta

    pure elemental function binom(n, k) result(y)
        real(dp), intent(in) :: n !! generalized upper binomial argument
        real(dp), intent(in) :: k !! generalized lower binomial argument
        real(dp) :: y

        real(dp) :: log_abs
        real(dp) :: sign_value

        if (ieee_is_nan(n) .or. ieee_is_nan(k)) then
            y = quiet_nan(n)
            return
        end if

        if (.not. ieee_is_finite(n) .or. .not. ieee_is_finite(k)) then
            y = binom_infinite_argument(n, k)
            return
        end if

        ! SciPy defines generalized binom as NaN for a negative-integer upper
        ! argument, even where a combinatorial continuation could be assigned.
        if (n < 0.0_dp .and. n == anint(n)) then
            y = quiet_nan(n)
            return
        end if

        if (k < 0.0_dp .and. k == anint(k)) then
            y = 0.0_dp
            return
        end if
        if (n >= 0.0_dp .and. n == anint(n) .and. k == anint(k)) then
            if (k > n) then
                y = 0.0_dp
                return
            end if
            if (n <= 1.0e6_dp .and. min(k, n - k) <= 64.0_dp) then
                y = binom_small_integer(n, k)
                return
            end if
        end if

        if (is_nonpositive_integer(k + 1.0_dp) .or. &
                is_nonpositive_integer(n - k + 1.0_dp)) then
            y = 0.0_dp
            return
        end if

        log_abs = gammaln(n + 1.0_dp) - gammaln(k + 1.0_dp) - &
            gammaln(n - k + 1.0_dp)
        sign_value = gammasgn(n + 1.0_dp) * gammasgn(k + 1.0_dp) * &
            gammasgn(n - k + 1.0_dp)
        y = signed_exp(log_abs, sign_value, n)
    end function binom

    pure elemental function binom_small_integer(n, k) result(y)
        real(dp), intent(in) :: n !! nonnegative integer upper argument, at most 1e6
        real(dp), intent(in) :: k !! integer lower argument in [0,n]
        real(dp) :: y

        integer :: j
        integer :: r
        real(dp) :: r_real

        r_real = min(k, n - k)
        r = int(r_real)
        y = 1.0_dp
        do j = 1, r
            y = y * (n - r_real + real(j, dp)) / real(j, dp)
        end do
    end function binom_small_integer

    pure elemental function poch(a, m) result(y)
        real(dp), intent(in) :: a !! initial value of the rising factorial
        real(dp), intent(in) :: m !! real increment/order of the rising factorial
        real(dp) :: y

        real(dp) :: log_abs
        real(dp) :: sign_value

        if (m == 0.0_dp) then
            y = 1.0_dp
            return
        end if
        if (ieee_is_nan(a) .or. ieee_is_nan(m)) then
            y = quiet_nan(a)
            return
        end if

        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(m)) then
            y = poch_infinite_argument(a, m)
            return
        end if

        if (m == anint(m) .and. a == anint(a)) then
            y = poch_integer_pair(a, m)
            return
        end if

        if (is_nonpositive_integer(a)) then
            y = 0.0_dp
            return
        end if
        if (is_nonpositive_integer(a + m)) then
            y = positive_infinity(a + m)
            return
        end if

        log_abs = gammaln(a + m) - gammaln(a)
        sign_value = gammasgn(a + m) * gammasgn(a)
        y = signed_exp(log_abs, sign_value, a)
    end function poch

    pure elemental function polygamma(n, x) result(y)
        integer, intent(in) :: n !! derivative order, nonnegative
        real(dp), intent(in) :: x !! real argument of the polygamma function
        real(dp) :: y

        integer :: shift_count
        real(dp) :: correction
        real(dp) :: q

        if (n < 0) then
            y = quiet_nan(x)
            return
        end if
        if (n == 0) then
            y = digamma(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = parity_signed_zero(n)
            else
                y = quiet_nan(x)
            end if
            return
        end if
        if (x <= 0.0_dp .and. x == anint(x)) then
            y = parity_signed_infinity(n, x)
            return
        end if
        if (x > 0.0_dp) then
            y = polygamma_positive(n, x)
            return
        end if

        ! For nearby negative arguments, recurrence avoids the cancellation and
        ! coefficient growth of a high-order cotangent reflection polynomial.
        if (x > -16.0_dp) then
            shift_count = int(floor(-x)) + 1
            q = x + real(shift_count, dp)
            correction = polygamma_recurrence_correction(n, x, shift_count)
            if (.not. ieee_is_finite(correction)) then
                y = correction
            else
                y = polygamma_positive(n, q) + correction
            end if
        else
            y = polygamma_reflected(n, x)
        end if
    end function polygamma

    pure function polygamma_positive(n, x) result(y)
        integer, intent(in) :: n !! strictly positive derivative order
        real(dp), intent(in) :: x !! strictly positive finite argument
        real(dp) :: y

        real(dp) :: log_abs
        real(dp) :: log_factorial
        real(dp) :: leading_log
        real(dp) :: s
        real(dp) :: zeta_value

        s = real(n + 1, dp)
        log_factorial = gammaln(real(n + 1, dp))
        leading_log = log_factorial - s * log(x)
        if (leading_log > log(huge(1.0_dp))) then
            y = parity_signed_infinity(n, x)
            return
        end if

        zeta_value = hurwitz_zeta(s, x)
        if (zeta_value == 0.0_dp) then
            y = parity_signed_zero(n)
            return
        end if
        if (.not. ieee_is_finite(zeta_value)) then
            y = parity_signed_infinity(n, x)
            return
        end if

        log_abs = log_factorial + log(zeta_value)
        y = signed_exp(log_abs, parity_sign(n), x)
    end function polygamma_positive

    pure function polygamma_recurrence_correction(n, x, shift_count) result(y)
        integer, intent(in) :: n !! strictly positive derivative order
        real(dp), intent(in) :: x !! negative noninteger starting argument
        integer, intent(in) :: shift_count !! number of unit recurrences to apply
        real(dp) :: y

        integer :: j
        real(dp) :: compensation
        real(dp) :: log_abs
        real(dp) :: log_factorial
        real(dp) :: term
        real(dp) :: term_sign
        real(dp) :: sum_temp
        real(dp) :: t
        real(dp) :: xj

        y = 0.0_dp
        compensation = 0.0_dp
        log_factorial = gammaln(real(n + 1, dp))
        do j = 0, shift_count - 1
            xj = x + real(j, dp)
            log_abs = log_factorial - real(n + 1, dp) * log(abs(xj))
            if (log_abs > log(huge(1.0_dp))) then
                term_sign = parity_sign(n)
                if (modulo(real(n + 1, dp), 2.0_dp) /= 0.0_dp .and. xj < 0.0_dp) then
                    term_sign = -term_sign
                end if
                y = sign(positive_infinity(x), term_sign)
                return
            end if
            term_sign = parity_sign(n)
            if (modulo(real(n + 1, dp), 2.0_dp) /= 0.0_dp .and. xj < 0.0_dp) then
                term_sign = -term_sign
            end if
            term = sign(exp(log_abs), term_sign)
            t = term - compensation
            sum_temp = y + t
            compensation = (sum_temp - y) - t
            y = sum_temp
        end do
    end function polygamma_recurrence_correction

    pure function polygamma_reflected(n, x) result(y)
        integer, intent(in) :: n !! strictly positive derivative order
        real(dp), intent(in) :: x !! finite negative noninteger argument
        real(dp) :: y

        real(dp) :: cot_value
        real(dp) :: derivative_value
        real(dp) :: positive_value
        real(dp) :: reflection

        cot_value = cospi_reduced(x) / sinpi_reduced(x)
        derivative_value = cot_derivative_polynomial(n, cot_value)
        reflection = (scifort_pi ** (n + 1)) * derivative_value
        positive_value = polygamma_positive(n, 1.0_dp - x)
        if (modulo(real(n, dp), 2.0_dp) == 0.0_dp) then
            y = positive_value - reflection
        else
            y = -positive_value - reflection
        end if
    end function polygamma_reflected

    pure function cot_derivative_polynomial(n, c) result(value)
        integer, intent(in) :: n !! derivative order with respect to the angle
        real(dp), intent(in) :: c !! cotangent of the reduced angle
        real(dp) :: value

        integer :: degree
        integer :: i
        integer :: j
        real(dp) :: coeff(0:n + 1)
        real(dp) :: next_coeff(0:n + 1)

        coeff = 0.0_dp
        coeff(1) = 1.0_dp
        degree = 1
        do i = 1, n
            next_coeff = 0.0_dp
            do j = 1, degree
                next_coeff(j - 1) = next_coeff(j - 1) - real(j, dp) * coeff(j)
                next_coeff(j + 1) = next_coeff(j + 1) - real(j, dp) * coeff(j)
            end do
            coeff = next_coeff
            degree = degree + 1
        end do

        value = coeff(degree)
        do j = degree - 1, 0, -1
            value = value * c + coeff(j)
        end do
    end function cot_derivative_polynomial

    pure elemental function beta_infinite_argument(a, b) result(y)
        real(dp), intent(in) :: a !! first possibly infinite beta argument
        real(dp), intent(in) :: b !! second possibly infinite beta argument
        real(dp) :: y

        real(dp) :: finite_arg

        if (a < 0.0_dp .or. b < 0.0_dp) then
            if ((.not. ieee_is_finite(a) .and. a < 0.0_dp) .or. &
                    (.not. ieee_is_finite(b) .and. b < 0.0_dp)) then
                y = positive_infinity(a)
                return
            end if
        end if
        if (.not. ieee_is_finite(a) .and. .not. ieee_is_finite(b)) then
            y = quiet_nan(a)
            return
        end if

        if (.not. ieee_is_finite(a)) then
            finite_arg = b
        else
            finite_arg = a
        end if

        if (finite_arg > 0.0_dp) then
            y = 0.0_dp
        else if (finite_arg == 0.0_dp .or. is_nonpositive_integer(finite_arg)) then
            y = positive_infinity(finite_arg)
        else
            y = sign(positive_infinity(finite_arg), gammasgn(finite_arg))
        end if
    end function beta_infinite_argument

    pure elemental function binom_infinite_argument(n, k) result(y)
        real(dp), intent(in) :: n !! possibly infinite generalized upper argument
        real(dp), intent(in) :: k !! possibly infinite generalized lower argument
        real(dp) :: y

        if (.not. ieee_is_finite(n)) then
            if (n < 0.0_dp) then
                y = quiet_nan(n)
            else if (.not. ieee_is_finite(k)) then
                if (k < 0.0_dp) then
                    y = 0.0_dp
                else
                    y = quiet_nan(k)
                end if
            else if (k == 0.0_dp) then
                y = 1.0_dp
            else if (k < 0.0_dp .and. k == anint(k)) then
                y = 0.0_dp
            else if (k > 0.0_dp .and. k == anint(k)) then
                y = positive_infinity(n)
            else
                y = quiet_nan(k)
            end if
        else
            y = quiet_nan(k)
        end if
    end function binom_infinite_argument

    pure elemental function poch_infinite_argument(a, m) result(y)
        real(dp), intent(in) :: a !! possibly infinite initial value
        real(dp), intent(in) :: m !! possibly infinite increment
        real(dp) :: y

        if (.not. ieee_is_finite(a)) then
            if (a > 0.0_dp) then
                if (m > 0.0_dp) then
                    y = positive_infinity(a)
                else if (m < 0.0_dp) then
                    y = 0.0_dp
                else
                    y = 1.0_dp
                end if
            else if (m == 1.0_dp) then
                y = negative_infinity(a)
            else if (m == -1.0_dp) then
                y = -0.0_dp
            else
                y = quiet_nan(a)
            end if
            return
        end if

        if (m > 0.0_dp) then
            if (is_nonpositive_integer(a)) then
                y = 0.0_dp
            else
                y = sign(positive_infinity(a), gammasgn(a))
            end if
        else
            if (is_nonpositive_integer(a)) then
                y = 0.0_dp
            else
                y = quiet_nan(m)
            end if
        end if
    end function poch_infinite_argument

    pure elemental function poch_integer_pair(a, m) result(y)
        real(dp), intent(in) :: a !! finite integer initial value
        real(dp), intent(in) :: m !! finite nonzero integer increment
        real(dp) :: y

        real(dp) :: log_abs
        real(dp) :: p
        real(dp) :: r
        real(dp) :: sign_value

        if (m > 0.0_dp) then
            if (a <= 0.0_dp) then
                p = -a
                if (m > p) then
                    y = 0.0_dp
                    return
                end if
                log_abs = gammaln(p + 1.0_dp) - gammaln(p - m + 1.0_dp)
                sign_value = merge(1.0_dp, -1.0_dp, modulo(m, 2.0_dp) == 0.0_dp)
            else
                log_abs = gammaln(a + m) - gammaln(a)
                sign_value = 1.0_dp
            end if
        else
            r = -m
            if (a >= 1.0_dp .and. a <= r) then
                y = positive_infinity(a)
                return
            end if
            if (a <= 0.0_dp) then
                log_abs = gammaln(-a + 1.0_dp) - gammaln(-a + r + 1.0_dp)
                sign_value = merge(1.0_dp, -1.0_dp, modulo(r, 2.0_dp) == 0.0_dp)
            else
                log_abs = gammaln(a - r) - gammaln(a)
                sign_value = 1.0_dp
            end if
        end if

        y = signed_exp(log_abs, sign_value, a)
    end function poch_integer_pair

    pure subroutine beta_integer_cancellation(pole_arg, positive_integer_arg, &
            log_abs, sign_value)
        real(dp), intent(in) :: pole_arg !! nonpositive integer beta argument
        real(dp), intent(in) :: positive_integer_arg !! positive integer cancellation argument
        real(dp), intent(out) :: log_abs !! logarithm of finite beta magnitude
        real(dp), intent(out) :: sign_value !! sign of the finite beta value

        real(dp) :: p
        real(dp) :: m

        p = -pole_arg
        m = positive_integer_arg
        log_abs = gammaln(m) - gammaln(p + 1.0_dp) + gammaln(p - m + 1.0_dp)
        sign_value = merge(1.0_dp, -1.0_dp, modulo(m, 2.0_dp) == 0.0_dp)
    end subroutine beta_integer_cancellation

    pure elemental logical function is_nonpositive_integer(x) result(is_integer)
        real(dp), intent(in) :: x !! finite value to classify

        is_integer = x <= 0.0_dp .and. x == anint(x)
    end function is_nonpositive_integer

    pure elemental logical function is_positive_integer(x) result(is_integer)
        real(dp), intent(in) :: x !! finite value to classify

        is_integer = x > 0.0_dp .and. x == anint(x)
    end function is_positive_integer

    pure elemental function parity_sign(n) result(sign_value)
        integer, intent(in) :: n !! positive derivative order
        real(dp) :: sign_value

        if (modulo(n, 2) == 0) then
            sign_value = -1.0_dp
        else
            sign_value = 1.0_dp
        end if
    end function parity_sign

    pure elemental function parity_signed_zero(n) result(y)
        integer, intent(in) :: n !! positive derivative order
        real(dp) :: y

        y = sign(0.0_dp, parity_sign(n))
    end function parity_signed_zero

    pure elemental function parity_signed_infinity(n, kind_arg) result(y)
        integer, intent(in) :: n !! positive derivative order
        real(dp), intent(in) :: kind_arg !! real kind selector for IEEE infinity
        real(dp) :: y

        y = sign(positive_infinity(kind_arg), parity_sign(n))
    end function parity_signed_infinity

    pure elemental function signed_exp(log_abs, sign_value, kind_arg) result(y)
        real(dp), intent(in) :: log_abs !! logarithm of the requested magnitude
        real(dp), intent(in) :: sign_value !! requested sign, normally -1 or +1
        real(dp), intent(in) :: kind_arg !! real kind selector for IEEE infinity
        real(dp) :: y

        if (log_abs > log(huge(1.0_dp))) then
            y = sign(positive_infinity(kind_arg), sign_value)
        else if (log_abs == negative_infinity(log_abs)) then
            y = sign(0.0_dp, sign_value)
        else
            y = sign(exp(log_abs), sign_value)
        end if
    end function signed_exp

    pure elemental function sinpi_reduced(x) result(y)
        real(dp), intent(in) :: x !! finite real argument of sin(pi*x)
        real(dp) :: y

        real(dp) :: r

        r = modulo(x, 2.0_dp)
        if (r > 1.0_dp) r = r - 2.0_dp
        if (r == 0.0_dp .or. abs(r) == 1.0_dp) then
            y = 0.0_dp
        else if (abs(r) == 0.5_dp) then
            y = sign(1.0_dp, r)
        else
            y = sin(scifort_pi * r)
        end if
    end function sinpi_reduced

    pure elemental function cospi_reduced(x) result(y)
        real(dp), intent(in) :: x !! finite real argument of cos(pi*x)
        real(dp) :: y

        real(dp) :: r

        r = modulo(x, 2.0_dp)
        if (r == 0.0_dp) then
            y = 1.0_dp
        else if (r == 1.0_dp) then
            y = -1.0_dp
        else if (r == 0.5_dp .or. r == 1.5_dp) then
            y = 0.0_dp
        else
            if (r > 1.0_dp) r = r - 2.0_dp
            if (abs(r) == 0.5_dp) then
                y = 0.0_dp
            else
                y = cos(scifort_pi * r)
            end if
        end if
    end function cospi_reduced

end module scifort_gamma_combinatorics
