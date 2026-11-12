! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized hypergeometric special functions.
!
! hyp0f1 is implemented independently from the defining generalized-
! hypergeometric series and differential equation in NIST DLMF 16.2 and 16.8.
! Real large arguments additionally use the DLMF 10.16.9 Bessel connection.
!
! hyp1f1 is implemented independently from the defining Kummer series, the
! differential equation in DLMF 13.2, Kummer's transformation, and the
! large-argument expansions in DLMF 13.7.
!
! hyp2f1 is implemented independently from the defining Gauss series and the
! hypergeometric differential equation in NIST DLMF 15.2 and 15.10. Analytic
! continuation uses local Taylor series along paths that avoid z=1.
!
! hyperu is implemented independently from the Kummer differential equation,
! connection formula, and large-positive-argument asymptotic expansion in
! NIST DLMF Chapter 13.

module scifort_hypergeometric
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_is_negative
    use scifort_bessel_jy, only : jv
    use scifort_bessel_order, only : iv
    use scifort_gamma_family, only : gammasgn
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_special_elementary, only : gammaln
    implicit none
    private

    public :: hyp0f1
    public :: hyp1f1
    public :: hyp2f1
    public :: hyperu

    interface hyp0f1
        module procedure hyp0f1_real
        module procedure hyp0f1_complex
    end interface hyp0f1

    interface hyp1f1
        module procedure hyp1f1_real
        module procedure hyp1f1_complex
    end interface hyp1f1

    interface hyp2f1
        module procedure hyp2f1_real
        module procedure hyp2f1_complex
    end interface hyp2f1

    real(dp), parameter :: origin_series_radius = 32.0_dp
    real(dp), parameter :: continuation_step_ratio = 0.25_dp
    integer, parameter :: max_origin_terms = 10000
    integer, parameter :: max_taylor_terms = 160

    real(dp), parameter :: hyp1f1_origin_radius = 6.0_dp
    real(dp), parameter :: hyp1f1_asymptotic_radius = 80.0_dp
    real(dp), parameter :: hyp1f1_step_ratio = 0.15_dp
    integer, parameter :: hyp1f1_max_taylor_terms = 240
    integer, parameter :: hyp1f1_max_asymptotic_terms = 160

    real(dp), parameter :: hyperu_asymptotic_start = 24.0_dp
    real(dp), parameter :: hyperu_asymptotic_target = 2.0e-14_dp
    real(dp), parameter :: hyperu_integer_offset = 1.0e-3_dp
    real(dp), parameter :: hyperu_near_integer_width = 1.0e-2_dp
    real(dp), parameter :: hyperu_small_x = 0.125_dp
    real(dp), parameter :: hyperu_step_ratio = 0.01_dp
    integer, parameter :: hyperu_max_anchor_expansions = 24
    integer, parameter :: hyperu_max_asymptotic_terms = 320
    integer, parameter :: hyperu_max_continuation_steps = 10000

    real(dp), parameter :: hyp2f1_origin_radius = 0.30_dp
    real(dp), parameter :: hyp2f1_step_ratio = 0.25_dp
    real(dp), parameter :: hyp2f1_detour_real_threshold = 0.75_dp
    integer, parameter :: hyp2f1_max_taylor_terms = 240
    integer, parameter :: hyp2f1_max_continuation_steps = 10000

contains

    pure elemental function hyperu(a, b, x) result(value)
        real(dp), intent(in) :: a !! first real parameter
        real(dp), intent(in) :: b !! second real parameter
        real(dp), intent(in) :: x !! real argument; x >= 0 is on the principal real branch
        real(dp) :: value

        if (ieee_is_nan(a) .or. ieee_is_nan(b) .or. ieee_is_nan(x)) then
            value = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            value = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp) then
            value = quiet_nan(x)
            return
        end if

        if (x == 0.0_dp) then
            value = hyperu_at_zero(a, b)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (a > 0.0_dp) then
                value = 0.0_dp
            else if (a == 0.0_dp) then
                value = 1.0_dp
            else
                value = positive_infinity(x)
            end if
            return
        end if

        if (a == 0.0_dp) then
            value = 1.0_dp
            return
        end if
        if (is_negative_integer(a)) then
            value = hyperu_negative_integer(a, b, x)
            return
        end if

        if (x < hyperu_small_x) then
            value = hyperu_small_positive_x(a, b, x)
        else
            value = hyperu_ode_continuation(a, b, x)
        end if
    end function hyperu

    pure function hyperu_at_zero(a, b) result(value)
        real(dp), intent(in) :: a !! finite first parameter
        real(dp), intent(in) :: b !! finite second parameter
        real(dp) :: value

        real(dp) :: log_abs
        real(dp) :: sign_value

        if (a == 0.0_dp) then
            value = 1.0_dp
            return
        end if
        if (is_negative_integer(a)) then
            value = hyperu_negative_integer(a, b, 0.0_dp)
            return
        end if

        if (b < 1.0_dp) then
            if (is_parameter_pole(1.0_dp + a - b)) then
                value = 0.0_dp
                return
            end if
            log_abs = gammaln(1.0_dp - b) - gammaln(1.0_dp + a - b)
            sign_value = gammasgn(1.0_dp + a - b)
            value = signed_exp_from_log(log_abs, sign_value)
        else
            sign_value = gammasgn(a)
            value = sign(positive_infinity(0.0_dp), sign_value)
        end if
    end function hyperu_at_zero

    pure function hyperu_negative_integer(a, b, x) result(value)
        real(dp), intent(in) :: a !! negative-integer first parameter
        real(dp), intent(in) :: b !! finite second parameter
        real(dp), intent(in) :: x !! nonnegative finite real argument
        real(dp) :: value

        integer :: k
        integer :: n
        real(dp) :: previous
        real(dp) :: current
        real(dp) :: next_value

        n = int(-anint(a))
        if (n == 0) then
            value = 1.0_dp
            return
        end if

        previous = 1.0_dp
        current = x - b
        if (n == 1) then
            value = current
            return
        end if

        do k = 1, n - 1
            next_value = (x - b - 2.0_dp * real(k, dp)) * current - &
                real(k, dp) * (real(k, dp) + b - 1.0_dp) * previous
            previous = current
            current = next_value
        end do
        value = current
    end function hyperu_negative_integer

    pure function hyperu_small_positive_x(a, b, x) result(value)
        real(dp), intent(in) :: a !! finite noninteger first parameter
        real(dp), intent(in) :: b !! finite second parameter
        real(dp), intent(in) :: x !! small positive argument
        real(dp) :: value

        real(dp) :: b_nearest

        b_nearest = anint(b)
        if (abs(b - b_nearest) < hyperu_near_integer_width) then
            value = hyperu_near_integer_connection(a, b, x)
        else
            value = hyperu_connection_formula(a, b, x)
        end if

        if (ieee_is_nan(value)) value = hyperu_ode_continuation(a, b, x)
    end function hyperu_small_positive_x

    pure function hyperu_near_integer_connection(a, b, x) result(value)
        real(dp), intent(in) :: a !! finite noninteger first parameter
        real(dp), intent(in) :: b !! finite second parameter near an integer
        real(dp), intent(in) :: x !! small positive argument
        real(dp) :: value

        real(dp) :: coarse
        real(dp) :: fine
        real(dp) :: h
        real(dp) :: raw_value
        real(dp) :: v1m
        real(dp) :: v1p
        real(dp) :: v2m
        real(dp) :: v2p

        h = hyperu_integer_offset
        if (abs((b + h) - anint(b + h)) <= 64.0_dp * epsilon(1.0_dp) .or. &
                abs((b - h) - anint(b - h)) <= 64.0_dp * epsilon(1.0_dp)) then
            h = 1.3_dp * hyperu_integer_offset
        end if

        v1p = hyperu_connection_formula(a, b + h, x)
        v1m = hyperu_connection_formula(a, b - h, x)
        v2p = hyperu_connection_formula(a, b + 0.5_dp * h, x)
        v2m = hyperu_connection_formula(a, b - 0.5_dp * h, x)

        if (ieee_is_finite(v1p) .and. ieee_is_finite(v1m) .and. &
                ieee_is_finite(v2p) .and. ieee_is_finite(v2m)) then
            coarse = 0.5_dp * (v1p + v1m)
            fine = 0.5_dp * (v2p + v2m)
            value = (4.0_dp * fine - coarse) / 3.0_dp
            return
        end if

        raw_value = hyperu_connection_formula(a, b, x)
        value = raw_value
    end function hyperu_near_integer_connection

    pure function hyperu_connection_formula(a, b, x) result(value)
        real(dp), intent(in) :: a !! finite noninteger first parameter
        real(dp), intent(in) :: b !! finite noninteger second parameter
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: value

        real(dp) :: first_m
        real(dp) :: first_term
        real(dp) :: second_m
        real(dp) :: second_term

        if (b == anint(b)) then
            value = quiet_nan(x)
            return
        end if

        first_m = hyp1f1_real(a, b, x)
        second_m = hyp1f1_real(1.0_dp + a - b, 2.0_dp - b, x)
        first_term = hyperu_gamma_power_term(1.0_dp - b, 1.0_dp + a - b, &
            0.0_dp, x, first_m)
        second_term = hyperu_gamma_power_term(b - 1.0_dp, a, 1.0_dp - b, &
            x, second_m)

        if ((.not. ieee_is_finite(first_term)) .and. (.not. ieee_is_finite(second_term)) .and. &
                sign(1.0_dp, first_term) /= sign(1.0_dp, second_term)) then
            value = quiet_nan(x)
        else
            value = first_term + second_term
        end if
    end function hyperu_connection_formula

    pure function hyperu_gamma_power_term(gamma_numerator, gamma_denominator, &
            power, x, multiplier) result(value)
        real(dp), intent(in) :: gamma_numerator !! numerator Gamma argument
        real(dp), intent(in) :: gamma_denominator !! denominator Gamma argument
        real(dp), intent(in) :: power !! exponent applied to positive x
        real(dp), intent(in) :: x !! positive finite base
        real(dp), intent(in) :: multiplier !! finite confluent-hypergeometric factor
        real(dp) :: value

        real(dp) :: log_abs
        real(dp) :: sign_value

        if (multiplier == 0.0_dp .or. is_parameter_pole(gamma_denominator)) then
            value = 0.0_dp
            return
        end if
        if (is_parameter_pole(gamma_numerator) .or. ieee_is_nan(multiplier)) then
            value = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(multiplier)) then
            value = sign(positive_infinity(x), multiplier)
            return
        end if

        log_abs = gammaln(gamma_numerator) - gammaln(gamma_denominator) + &
            power * log(x) + log(abs(multiplier))
        sign_value = gammasgn(gamma_numerator) * gammasgn(gamma_denominator) * &
            sign(1.0_dp, multiplier)
        value = signed_exp_from_log(log_abs, sign_value)
    end function hyperu_gamma_power_term

    pure function hyperu_ode_continuation(a, b, x) result(value)
        real(dp), intent(in) :: a !! finite noninteger first parameter
        real(dp), intent(in) :: b !! finite second parameter
        real(dp), intent(in) :: x !! positive finite target argument
        real(dp) :: value

        complex(dp) :: derivative_complex
        complex(dp) :: h
        complex(dp) :: value_complex
        real(dp) :: anchor
        real(dp) :: derivative
        real(dp) :: relative_term
        real(dp) :: remaining
        real(dp) :: step_size
        logical :: terminated
        integer :: expansion
        integer :: step

        anchor = max(x, hyperu_asymptotic_start)
        do expansion = 1, hyperu_max_anchor_expansions
            call hyperu_asymptotic_anchor(a, b, anchor, value, derivative, &
                relative_term, terminated)
            if (terminated .or. relative_term <= hyperu_asymptotic_target) exit
            anchor = 1.5_dp * anchor
        end do

        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(derivative)) return
        if (anchor == x) return

        value_complex = cmplx(value, 0.0_dp, kind=dp)
        derivative_complex = cmplx(derivative, 0.0_dp, kind=dp)

        do step = 1, hyperu_max_continuation_steps
            remaining = anchor - x
            if (remaining <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, x)) exit
            step_size = min(remaining, hyperu_step_ratio * anchor)
            h = cmplx(-step_size, 0.0_dp, kind=dp)
            call hyp1f1_taylor_step(a, b, cmplx(anchor, 0.0_dp, kind=dp), h, &
                value_complex, derivative_complex)
            anchor = anchor - step_size
            if (.not. complex_is_finite(value_complex) .or. &
                    .not. complex_is_finite(derivative_complex)) exit
        end do

        value = real(value_complex, dp)
    end function hyperu_ode_continuation

    pure subroutine hyperu_asymptotic_anchor(a, b, x, value, derivative, &
            relative_term, terminated)
        real(dp), intent(in) :: a !! finite noninteger first parameter
        real(dp), intent(in) :: b !! finite second parameter
        real(dp), intent(in) :: x !! positive asymptotic anchor
        real(dp), intent(out) :: value !! asymptotic U value at x
        real(dp), intent(out) :: derivative !! derivative with respect to x
        real(dp), intent(out) :: relative_term !! relative magnitude of best retained term
        logical, intent(out) :: terminated !! true when the asymptotic series terminates exactly

        real(dp) :: best_derivative_sum
        real(dp) :: best_sum
        real(dp) :: derivative_sum
        real(dp) :: factor
        real(dp) :: log_scale
        real(dp) :: magnitude
        real(dp) :: previous_magnitude
        real(dp) :: scale
        real(dp) :: sum_value
        real(dp) :: term
        integer :: k

        term = 1.0_dp
        sum_value = 1.0_dp
        derivative_sum = -a
        best_sum = sum_value
        best_derivative_sum = derivative_sum
        relative_term = huge(1.0_dp)
        previous_magnitude = huge(1.0_dp)
        terminated = .false.

        do k = 1, hyperu_max_asymptotic_terms
            factor = -((a + real(k - 1, dp)) * (1.0_dp + a - b + real(k - 1, dp))) / &
                (real(k, dp) * x)
            term = term * factor
            sum_value = sum_value + term
            derivative_sum = derivative_sum + (-a - real(k, dp)) * term
            magnitude = abs(term)

            if (magnitude / max(tiny(1.0_dp), abs(sum_value)) < relative_term) then
                relative_term = magnitude / max(tiny(1.0_dp), abs(sum_value))
                best_sum = sum_value
                best_derivative_sum = derivative_sum
            end if

            if (term == 0.0_dp) then
                relative_term = 0.0_dp
                best_sum = sum_value
                best_derivative_sum = derivative_sum
                terminated = .true.
                exit
            end if
            if (k >= 4 .and. magnitude > previous_magnitude) exit
            previous_magnitude = magnitude
            if (k >= 8 .and. magnitude <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(sum_value))) exit
        end do

        log_scale = -a * log(x)
        if (log_scale > log(huge(1.0_dp))) then
            value = sign(positive_infinity(x), best_sum)
            derivative = sign(positive_infinity(x), best_derivative_sum)
        else if (log_scale < log(tiny(1.0_dp))) then
            value = sign(0.0_dp, best_sum)
            derivative = sign(0.0_dp, best_derivative_sum)
        else
            scale = exp(log_scale)
            value = scale * best_sum
            derivative = (scale / x) * best_derivative_sum
        end if
    end subroutine hyperu_asymptotic_anchor

    pure elemental function hyp0f1_real(v, z) result(value)
        real(dp), intent(in) :: v !! denominator parameter; nonpositive integers are poles
        real(dp), intent(in) :: z !! real hypergeometric argument
        real(dp) :: value

        real(dp) :: bessel_value

        if (z == 0.0_dp) then
            if (ieee_is_nan(v)) then
                value = 1.0_dp
            else if (is_parameter_pole(v) .or. (.not. ieee_is_finite(v) .and. v < 0.0_dp)) then
                value = quiet_nan(z)
            else
                value = 1.0_dp
            end if
            return
        end if

        if (ieee_is_nan(z) .or. ieee_is_nan(v)) then
            value = quiet_nan(z)
            return
        end if
        if (.not. ieee_is_finite(z)) then
            value = quiet_nan(z)
            return
        end if
        if (.not. ieee_is_finite(v)) then
            if (v > 0.0_dp) then
                value = 1.0_dp
            else
                value = quiet_nan(z)
            end if
            return
        end if
        if (is_parameter_pole(v)) then
            value = quiet_nan(z)
            return
        end if

        if (abs(z) <= origin_series_radius) then
            value = hyp0f1_real_series(v, z)
            return
        end if

        if (z < 0.0_dp) then
            bessel_value = jv(v - 1.0_dp, 2.0_dp * sqrt(-z))
        else
            bessel_value = iv(v - 1.0_dp, 2.0_dp * sqrt(z))
        end if

        if (bessel_value == 0.0_dp .or. .not. ieee_is_finite(bessel_value)) then
            value = hyp0f1_real_series(v, z)
        else
            value = hyp0f1_from_real_bessel(v, abs(z), bessel_value)
        end if
    end function hyp0f1_real

    pure elemental function hyp0f1_complex(v, z) result(value)
        real(dp), intent(in) :: v !! denominator parameter; nonpositive integers are poles
        complex(dp), intent(in) :: z !! complex hypergeometric argument
        complex(dp) :: value

        complex(dp) :: derivative
        complex(dp) :: h
        complex(dp) :: z0
        real(dp) :: radius
        real(dp) :: remaining

        if (real(z, dp) == 0.0_dp .and. aimag(z) == 0.0_dp) then
            if (ieee_is_nan(v)) then
                value = cmplx(1.0_dp, 0.0_dp, kind=dp)
            else if (is_parameter_pole(v) .or. (.not. ieee_is_finite(v) .and. v < 0.0_dp)) then
                value = complex_nan(.false.)
            else
                value = cmplx(1.0_dp, 0.0_dp, kind=dp)
            end if
            return
        end if

        if (complex_has_nan(z) .or. ieee_is_nan(v)) then
            value = complex_nan(.true.)
            return
        end if
        if (.not. complex_is_finite(z)) then
            value = complex_nan(.true.)
            return
        end if
        if (.not. ieee_is_finite(v)) then
            if (v > 0.0_dp) then
                value = cmplx(1.0_dp, 0.0_dp, kind=dp)
            else
                value = complex_nan(.false.)
            end if
            return
        end if
        if (is_parameter_pole(v)) then
            value = complex_nan(.false.)
            return
        end if

        radius = abs(z)
        if (radius <= origin_series_radius) then
            call hyp0f1_origin_series(v, z, value, derivative)
            return
        end if

        z0 = z * (origin_series_radius / radius)
        call hyp0f1_origin_series(v, z0, value, derivative)

        do
            remaining = abs(z - z0)
            if (remaining <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, radius)) exit

            if (remaining <= continuation_step_ratio * abs(z0)) then
                h = z - z0
            else
                h = (z - z0) * (continuation_step_ratio * abs(z0) / remaining)
            end if

            call hyp0f1_taylor_step(v, z0, h, value, derivative)
            z0 = z0 + h
        end do
    end function hyp0f1_complex

    pure elemental function hyp1f1_real(a, b, z) result(value)
        real(dp), intent(in) :: a !! real numerator parameter
        real(dp), intent(in) :: b !! real denominator parameter
        real(dp), intent(in) :: z !! real hypergeometric argument
        real(dp) :: value

        complex(dp) :: continued

        if (ieee_is_nan(a) .or. ieee_is_nan(b) .or. ieee_is_nan(z)) then
            value = quiet_nan(z)
            return
        end if

        if (is_parameter_pole(b)) then
            if (real_pole_has_finite_polynomial(a, b)) then
                value = hyp1f1_polynomial_real(a, b, z)
            else
                value = positive_infinity(z)
            end if
            return
        end if

        if (.not. ieee_is_finite(a)) then
            value = quiet_nan(z)
            return
        end if
        if (.not. ieee_is_finite(b)) then
            if (b > 0.0_dp) then
                value = 1.0_dp
            else
                value = positive_infinity(z)
            end if
            return
        end if

        if (a == 0.0_dp) then
            value = 1.0_dp
            return
        end if
        if (a == b) then
            value = real_exp_safely(z)
            return
        end if
        if (is_negative_integer(a)) then
            value = hyp1f1_polynomial_real(a, b, z)
            return
        end if

        if (.not. ieee_is_finite(z)) then
            value = hyp1f1_real_infinite_argument(a, b, z)
            return
        end if
        if (z == 0.0_dp) then
            value = 1.0_dp
            return
        end if

        if (z >= hyp1f1_asymptotic_radius) then
            value = hyp1f1_real_asymptotic_positive(a, b, z)
        else if (z <= -hyp1f1_asymptotic_radius) then
            value = hyp1f1_real_asymptotic_negative(a, b, z)
        else if (abs(z) <= hyp1f1_origin_radius) then
            value = hyp1f1_real_series(a, b, z)
        else if (z < 0.0_dp) then
            continued = hyp1f1_continued(b - a, b, cmplx(-z, 0.0_dp, kind=dp))
            value = exp(z) * real(continued, dp)
        else
            continued = hyp1f1_continued(a, b, cmplx(z, 0.0_dp, kind=dp))
            value = real(continued, dp)
        end if
    end function hyp1f1_real

    pure elemental function hyp1f1_complex(a, b, z) result(value)
        real(dp), intent(in) :: a !! real numerator parameter
        real(dp), intent(in) :: b !! real denominator parameter
        complex(dp), intent(in) :: z !! complex hypergeometric argument
        complex(dp) :: value

        complex(dp) :: transformed
        real(dp) :: radius

        if (ieee_is_nan(a) .or. ieee_is_nan(b) .or. complex_has_nan(z)) then
            value = complex_nan(.true.)
            return
        end if

        if (is_parameter_pole(b)) then
            value = cmplx(positive_infinity(0.0_dp), 0.0_dp, kind=dp)
            return
        end if

        if (.not. ieee_is_finite(a)) then
            value = complex_nan(.true.)
            return
        end if
        if (.not. ieee_is_finite(b)) then
            value = cmplx(1.0_dp, 0.0_dp, kind=dp)
            return
        end if
        if (.not. complex_is_finite(z)) then
            value = complex_nan(.true.)
            return
        end if

        if (a == 0.0_dp) then
            value = cmplx(1.0_dp, 0.0_dp, kind=dp)
            return
        end if
        if (a == b) then
            value = exp(z)
            return
        end if
        if (is_negative_integer(a)) then
            value = hyp1f1_polynomial_complex(a, b, z)
            return
        end if
        if (aimag(z) == 0.0_dp) then
            value = cmplx(hyp1f1_real(a, b, real(z, dp)), 0.0_dp, kind=dp)
            return
        end if

        radius = abs(z)
        if (radius <= hyp1f1_origin_radius) then
            call hyp1f1_origin_series(a, b, z, value, transformed)
            return
        end if
        if (radius > hyp1f1_asymptotic_radius) then
            value = hyp1f1_complex_asymptotic(a, b, z)
            return
        end if

        if (real(z, dp) < 0.0_dp) then
            transformed = hyp1f1_continued(b - a, b, -z)
            value = exp(z) * transformed
        else
            value = hyp1f1_continued(a, b, z)
        end if
    end function hyp1f1_complex

    pure elemental function hyp2f1_real(a, b, c, z) result(value)
        real(dp), intent(in) :: a !! first real numerator parameter
        real(dp), intent(in) :: b !! second real numerator parameter
        real(dp), intent(in) :: c !! real denominator parameter; nonpositive integers are poles
        real(dp), intent(in) :: z !! real Gauss-hypergeometric argument
        real(dp) :: value

        complex(dp) :: continued
        real(dp) :: imaginary_tolerance

        if (z == 0.0_dp .or. a == 0.0_dp .or. b == 0.0_dp) then
            value = 1.0_dp
            return
        end if

        if (ieee_is_nan(a) .or. ieee_is_nan(b) .or. ieee_is_nan(c) .or. ieee_is_nan(z)) then
            value = quiet_nan(z)
            return
        end if

        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            value = quiet_nan(z)
            return
        end if
        if (.not. ieee_is_finite(c)) then
            if (c > 0.0_dp) then
                value = 1.0_dp
            else
                value = quiet_nan(z)
            end if
            return
        end if

        if (a == c) then
            value = hyp2f1_cancelled_real(b, z)
            return
        end if
        if (b == c) then
            value = hyp2f1_cancelled_real(a, z)
            return
        end if

        if (is_parameter_pole(c)) then
            if (hyp2f1_pole_has_finite_polynomial(a, b, c)) then
                value = hyp2f1_polynomial_real(a, b, c, z)
            else
                value = positive_infinity(z)
            end if
            return
        end if

        if (is_negative_integer(a) .or. is_negative_integer(b)) then
            value = hyp2f1_polynomial_real(a, b, c, z)
            return
        end if

        if (.not. ieee_is_finite(z)) then
            value = quiet_nan(z)
            return
        end if

        if (z == 1.0_dp) then
            value = hyp2f1_at_one(a, b, c)
            return
        end if

        if (abs(z) <= hyp2f1_origin_radius) then
            value = hyp2f1_real_series(a, b, c, z)
            return
        end if

        continued = hyp2f1_continued(a, b, c, cmplx(z, 0.0_dp, kind=dp))
        if (z < 1.0_dp) then
            value = real(continued, dp)
            return
        end if

        imaginary_tolerance = 128.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, abs(real(continued, dp)))
        if (abs(aimag(continued)) <= imaginary_tolerance) then
            value = real(continued, dp)
        else
            value = positive_infinity(z)
        end if
    end function hyp2f1_real

    pure elemental function hyp2f1_complex(a, b, c, z) result(value)
        real(dp), intent(in) :: a !! first real numerator parameter
        real(dp), intent(in) :: b !! second real numerator parameter
        real(dp), intent(in) :: c !! real denominator parameter; nonpositive integers are poles
        complex(dp), intent(in) :: z !! complex Gauss-hypergeometric argument
        complex(dp) :: value

        if ((real(z, dp) == 0.0_dp .and. aimag(z) == 0.0_dp) .or. &
                a == 0.0_dp .or. b == 0.0_dp) then
            value = cmplx(1.0_dp, 0.0_dp, kind=dp)
            return
        end if

        if (ieee_is_nan(a) .or. ieee_is_nan(b) .or. ieee_is_nan(c) .or. complex_has_nan(z)) then
            value = complex_nan(.true.)
            return
        end if

        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            value = complex_nan(.true.)
            return
        end if
        if (.not. ieee_is_finite(c)) then
            if (c > 0.0_dp) then
                value = cmplx(1.0_dp, 0.0_dp, kind=dp)
            else
                value = complex_nan(.true.)
            end if
            return
        end if

        if (a == c) then
            value = hyp2f1_cancelled_complex(b, z)
            return
        end if
        if (b == c) then
            value = hyp2f1_cancelled_complex(a, z)
            return
        end if

        if (is_parameter_pole(c)) then
            if (hyp2f1_pole_has_finite_polynomial(a, b, c)) then
                value = hyp2f1_polynomial_complex(a, b, c, z)
            else
                value = cmplx(positive_infinity(0.0_dp), 0.0_dp, kind=dp)
            end if
            return
        end if

        if (is_negative_integer(a) .or. is_negative_integer(b)) then
            value = hyp2f1_polynomial_complex(a, b, c, z)
            return
        end if

        if (.not. complex_is_finite(z)) then
            value = complex_nan(.true.)
            return
        end if

        if (real(z, dp) == 1.0_dp .and. aimag(z) == 0.0_dp) then
            value = cmplx(hyp2f1_at_one(a, b, c), 0.0_dp, kind=dp)
            return
        end if

        value = hyp2f1_continued(a, b, c, z)
    end function hyp2f1_complex

    pure function hyp2f1_real_series(a, b, c, z) result(value)
        real(dp), intent(in) :: a !! finite first numerator parameter
        real(dp), intent(in) :: b !! finite second numerator parameter
        real(dp), intent(in) :: c !! finite non-pole denominator parameter
        real(dp), intent(in) :: z !! finite real argument within the origin disk
        real(dp) :: value

        real(dp) :: compensation
        real(dp) :: corrected
        real(dp) :: term
        real(dp) :: updated
        integer :: n

        term = 1.0_dp
        value = 1.0_dp
        compensation = 0.0_dp

        do n = 1, max_origin_terms
            term = term * ((a + real(n - 1, dp)) * (b + real(n - 1, dp))) / &
                ((c + real(n - 1, dp)) * real(n, dp)) * z
            corrected = term - compensation
            updated = value + corrected
            compensation = (updated - value) - corrected
            value = updated

            if (n >= 12 .and. abs(term) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value))) exit
            if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(term)) exit
        end do
    end function hyp2f1_real_series

    pure subroutine hyp2f1_origin_series(a, b, c, z, value, derivative)
        real(dp), intent(in) :: a !! finite first numerator parameter
        real(dp), intent(in) :: b !! finite second numerator parameter
        real(dp), intent(in) :: c !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z !! finite nonzero argument within the origin disk
        complex(dp), intent(out) :: value !! 2F1 value at z
        complex(dp), intent(out) :: derivative !! derivative with respect to z

        complex(dp) :: compensation
        complex(dp) :: corrected
        complex(dp) :: derivative_term
        complex(dp) :: term
        complex(dp) :: updated
        integer :: n

        term = cmplx(1.0_dp, 0.0_dp, kind=dp)
        value = term
        derivative = cmplx(0.0_dp, 0.0_dp, kind=dp)
        compensation = cmplx(0.0_dp, 0.0_dp, kind=dp)

        do n = 1, max_origin_terms
            term = term * ((a + real(n - 1, dp)) * (b + real(n - 1, dp))) / &
                ((c + real(n - 1, dp)) * real(n, dp)) * z
            corrected = term - compensation
            updated = value + corrected
            compensation = (updated - value) - corrected
            value = updated

            derivative_term = real(n, dp) * term / z
            derivative = derivative + derivative_term

            if (n >= 12 .and. abs(term) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value)) .and. &
                    abs(derivative_term) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(derivative))) exit
        end do
    end subroutine hyp2f1_origin_series

    pure function hyp2f1_continued(a, b, c, z) result(value)
        real(dp), intent(in) :: a !! finite first numerator parameter
        real(dp), intent(in) :: b !! finite second numerator parameter
        real(dp), intent(in) :: c !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z !! finite target away from the singular point z=1
        complex(dp) :: value

        complex(dp) :: derivative
        complex(dp) :: first_target
        complex(dp) :: waypoint
        complex(dp) :: z0
        real(dp) :: first_radius
        real(dp) :: side

        if (abs(z) <= hyp2f1_origin_radius) then
            call hyp2f1_origin_series(a, b, c, z, value, derivative)
            return
        end if

        if (real(z, dp) >= 1.0_dp .or. &
                (real(z, dp) > hyp2f1_detour_real_threshold .and. aimag(z) /= 0.0_dp)) then
            if (ieee_is_negative(aimag(z))) then
                side = -1.0_dp
            else
                side = 1.0_dp
            end if
            waypoint = cmplx(0.45_dp, 0.55_dp * side, kind=dp)
            first_target = waypoint
        else
            waypoint = z
            first_target = z
        end if

        first_radius = abs(first_target)
        z0 = first_target * (hyp2f1_origin_radius / first_radius)
        call hyp2f1_origin_series(a, b, c, z0, value, derivative)

        call hyp2f1_continue_segment(a, b, c, z0, first_target, value, derivative)
        if (abs(first_target - z) > &
                8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(z))) then
            call hyp2f1_continue_segment(a, b, c, z0, z, value, derivative)
        end if
    end function hyp2f1_continued

    pure subroutine hyp2f1_continue_segment(a, b, c, z0, target, value, derivative)
        real(dp), intent(in) :: a !! finite first numerator parameter
        real(dp), intent(in) :: b !! finite second numerator parameter
        real(dp), intent(in) :: c !! finite non-pole denominator parameter
        complex(dp), intent(inout) :: z0 !! current continuation point, updated to target
        complex(dp), intent(in) :: target !! endpoint of a segment avoiding z=0 and z=1
        complex(dp), intent(inout) :: value !! input/output function value at z0
        complex(dp), intent(inout) :: derivative !! input/output derivative at z0

        complex(dp) :: h
        real(dp) :: convergence_radius
        real(dp) :: remaining
        real(dp) :: step_size
        integer :: step

        do step = 1, hyp2f1_max_continuation_steps
            remaining = abs(target - z0)
            if (remaining <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(target))) then
                z0 = target
                return
            end if

            convergence_radius = min(abs(z0), abs(1.0_dp - z0))
            if (convergence_radius <= 0.0_dp) then
                value = complex_nan(.true.)
                derivative = value
                return
            end if
            step_size = hyp2f1_step_ratio * convergence_radius
            if (remaining <= step_size) then
                h = target - z0
            else
                h = (target - z0) * (step_size / remaining)
            end if

            call hyp2f1_taylor_step(a, b, c, z0, h, value, derivative)
            z0 = z0 + h
        end do

        value = complex_nan(.true.)
        derivative = value
    end subroutine hyp2f1_continue_segment

    pure subroutine hyp2f1_taylor_step(a, b, c, z0, h, value, derivative)
        real(dp), intent(in) :: a !! finite first numerator parameter
        real(dp), intent(in) :: b !! finite second numerator parameter
        real(dp), intent(in) :: c !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z0 !! expansion center away from z=0 and z=1
        complex(dp), intent(in) :: h !! step inside the local convergence disk
        complex(dp), intent(inout) :: value !! input/output function value
        complex(dp), intent(inout) :: derivative !! input/output derivative with respect to z

        complex(dp) :: c_nm1
        complex(dp) :: c_nm2
        complex(dp) :: c_next
        complex(dp) :: derivative_new
        complex(dp) :: derivative_term
        complex(dp) :: h_power
        complex(dp) :: value_new
        complex(dp) :: value_term
        integer :: n
        real(dp) :: n_real

        c_nm2 = value
        c_nm1 = derivative
        value_new = value + derivative * h
        derivative_new = derivative
        h_power = h

        do n = 0, hyp2f1_max_taylor_terms - 3
            n_real = real(n, dp)
            c_next = ((n_real + a) * (n_real + b) * c_nm2 - &
                (n_real + 1.0_dp) * &
                (n_real + c - z0 * (2.0_dp * n_real + a + b + 1.0_dp)) * c_nm1) / &
                (z0 * (1.0_dp - z0) * (n_real + 2.0_dp) * (n_real + 1.0_dp))
            h_power = h_power * h
            value_term = c_next * h_power
            derivative_term = (n_real + 2.0_dp) * c_next * (h_power / h)
            value_new = value_new + value_term
            derivative_new = derivative_new + derivative_term

            if (n >= 10 .and. abs(value_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value_new)) .and. &
                    abs(derivative_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(derivative_new))) exit

            c_nm2 = c_nm1
            c_nm1 = c_next
        end do

        value = value_new
        derivative = derivative_new
    end subroutine hyp2f1_taylor_step

    pure function hyp2f1_polynomial_real(a, b, c, z) result(value)
        real(dp), intent(in) :: a !! first numerator parameter, with at least one terminating numerator
        real(dp), intent(in) :: b !! second numerator parameter, with at least one terminating numerator
        real(dp), intent(in) :: c !! denominator parameter with no pole before termination
        real(dp), intent(in) :: z !! real argument
        real(dp) :: value

        integer :: degree
        integer :: n
        real(dp) :: term

        degree = hyp2f1_polynomial_degree(a, b)
        term = 1.0_dp
        value = 1.0_dp
        do n = 1, degree
            term = term * ((a + real(n - 1, dp)) * (b + real(n - 1, dp))) / &
                ((c + real(n - 1, dp)) * real(n, dp)) * z
            value = value + term
        end do
    end function hyp2f1_polynomial_real

    pure function hyp2f1_polynomial_complex(a, b, c, z) result(value)
        real(dp), intent(in) :: a !! first numerator parameter, with at least one terminating numerator
        real(dp), intent(in) :: b !! second numerator parameter, with at least one terminating numerator
        real(dp), intent(in) :: c !! denominator parameter with no pole before termination
        complex(dp), intent(in) :: z !! complex argument
        complex(dp) :: value

        integer :: degree
        integer :: n
        complex(dp) :: term

        degree = hyp2f1_polynomial_degree(a, b)
        term = cmplx(1.0_dp, 0.0_dp, kind=dp)
        value = term
        do n = 1, degree
            term = term * ((a + real(n - 1, dp)) * (b + real(n - 1, dp))) / &
                ((c + real(n - 1, dp)) * real(n, dp)) * z
            value = value + term
        end do
    end function hyp2f1_polynomial_complex

    pure function hyp2f1_cancelled_real(exponent, z) result(value)
        real(dp), intent(in) :: exponent !! remaining numerator parameter in (1-z)^(-exponent)
        real(dp), intent(in) :: z !! real argument
        real(dp) :: value

        complex(dp) :: complex_value
        real(dp) :: imaginary_tolerance

        if (z < 1.0_dp) then
            value = exp(-exponent * log(1.0_dp - z))
            return
        end if
        if (z == 1.0_dp) then
            if (exponent > 0.0_dp) then
                value = positive_infinity(z)
            else if (exponent < 0.0_dp) then
                value = 0.0_dp
            else
                value = 1.0_dp
            end if
            return
        end if
        if (.not. ieee_is_finite(z)) then
            value = quiet_nan(z)
            return
        end if

        complex_value = hyp2f1_cancelled_complex(exponent, cmplx(z, 0.0_dp, kind=dp))
        imaginary_tolerance = 128.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, abs(real(complex_value, dp)))
        if (abs(aimag(complex_value)) <= imaginary_tolerance) then
            value = real(complex_value, dp)
        else
            value = positive_infinity(z)
        end if
    end function hyp2f1_cancelled_real

    pure function hyp2f1_cancelled_complex(exponent, z) result(value)
        real(dp), intent(in) :: exponent !! remaining numerator parameter in (1-z)^(-exponent)
        complex(dp), intent(in) :: z !! complex argument
        complex(dp) :: value

        complex(dp) :: one_minus_z

        if (complex_has_nan(z)) then
            value = complex_nan(.true.)
            return
        end if
        if (.not. complex_is_finite(z)) then
            value = complex_nan(.true.)
            return
        end if

        one_minus_z = cmplx(1.0_dp - real(z, dp), -aimag(z), kind=dp)
        if (real(one_minus_z, dp) == 0.0_dp .and. aimag(one_minus_z) == 0.0_dp) then
            if (exponent > 0.0_dp) then
                value = cmplx(positive_infinity(0.0_dp), 0.0_dp, kind=dp)
            else if (exponent < 0.0_dp) then
                value = cmplx(0.0_dp, 0.0_dp, kind=dp)
            else
                value = cmplx(1.0_dp, 0.0_dp, kind=dp)
            end if
            return
        end if

        value = exp(-exponent * log(one_minus_z))
    end function hyp2f1_cancelled_complex

    pure elemental function hyp2f1_polynomial_degree(a, b) result(degree)
        real(dp), intent(in) :: a !! first numerator parameter
        real(dp), intent(in) :: b !! second numerator parameter
        integer :: degree

        integer :: degree_a
        integer :: degree_b

        degree_a = huge(1)
        degree_b = huge(1)
        if (is_negative_integer(a)) degree_a = int(-anint(a))
        if (is_negative_integer(b)) degree_b = int(-anint(b))
        degree = min(degree_a, degree_b)
    end function hyp2f1_polynomial_degree

    pure elemental function hyp2f1_pole_has_finite_polynomial(a, b, c) result(is_finite_case)
        real(dp), intent(in) :: a !! first numerator parameter
        real(dp), intent(in) :: b !! second numerator parameter
        real(dp), intent(in) :: c !! nonpositive-integer denominator parameter
        logical :: is_finite_case

        integer :: degree

        if (.not. (is_negative_integer(a) .or. is_negative_integer(b))) then
            is_finite_case = .false.
            return
        end if

        degree = hyp2f1_polynomial_degree(a, b)
        is_finite_case = degree <= int(-anint(c))
    end function hyp2f1_pole_has_finite_polynomial

    pure function hyp2f1_at_one(a, b, c) result(value)
        real(dp), intent(in) :: a !! finite nonzero first numerator parameter
        real(dp), intent(in) :: b !! finite nonzero second numerator parameter
        real(dp), intent(in) :: c !! finite non-pole denominator parameter
        real(dp) :: value

        real(dp) :: delta
        real(dp) :: log_abs
        real(dp) :: sign_value

        delta = c - a - b
        if (delta <= 0.0_dp) then
            sign_value = gammasgn(c) * gammasgn(a) * gammasgn(b)
            value = sign(positive_infinity(0.0_dp), sign_value)
            return
        end if

        if (is_parameter_pole(c - a) .or. is_parameter_pole(c - b)) then
            value = 0.0_dp
            return
        end if

        log_abs = gammaln(c) + gammaln(delta) - gammaln(c - a) - gammaln(c - b)
        sign_value = gammasgn(c) * gammasgn(delta) * gammasgn(c - a) * gammasgn(c - b)
        value = signed_exp_from_log(log_abs, sign_value)
    end function hyp2f1_at_one

    pure function hyp1f1_real_series(a, b, z) result(value)
        real(dp), intent(in) :: a !! finite numerator parameter
        real(dp), intent(in) :: b !! finite non-pole denominator parameter
        real(dp), intent(in) :: z !! finite argument in the origin region
        real(dp) :: value

        real(dp) :: compensation
        real(dp) :: corrected
        real(dp) :: term
        real(dp) :: updated
        integer :: n

        term = 1.0_dp
        value = 1.0_dp
        compensation = 0.0_dp

        do n = 1, max_origin_terms
            term = term * ((a + real(n - 1, dp)) / (b + real(n - 1, dp))) * &
                (z / real(n, dp))
            corrected = term - compensation
            updated = value + corrected
            compensation = (updated - value) - corrected
            value = updated

            if (n >= 12 .and. abs(term) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value))) exit
            if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(term)) exit
        end do
    end function hyp1f1_real_series

    pure subroutine hyp1f1_origin_series(a, b, z, value, derivative)
        real(dp), intent(in) :: a !! finite numerator parameter
        real(dp), intent(in) :: b !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z !! finite nonzero argument in the origin disk
        complex(dp), intent(out) :: value !! 1F1 value at z
        complex(dp), intent(out) :: derivative !! derivative with respect to z

        complex(dp) :: compensation
        complex(dp) :: corrected
        complex(dp) :: derivative_term
        complex(dp) :: term
        complex(dp) :: updated
        integer :: n

        term = cmplx(1.0_dp, 0.0_dp, kind=dp)
        value = term
        derivative = cmplx(0.0_dp, 0.0_dp, kind=dp)
        compensation = cmplx(0.0_dp, 0.0_dp, kind=dp)

        do n = 1, max_origin_terms
            term = term * ((a + real(n - 1, dp)) / (b + real(n - 1, dp))) * &
                (z / real(n, dp))
            corrected = term - compensation
            updated = value + corrected
            compensation = (updated - value) - corrected
            value = updated

            derivative_term = real(n, dp) * term / z
            derivative = derivative + derivative_term

            if (n >= 12 .and. abs(term) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value)) .and. &
                    abs(derivative_term) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(derivative))) exit
        end do
    end subroutine hyp1f1_origin_series

    pure function hyp1f1_continued(a, b, z) result(value)
        real(dp), intent(in) :: a !! finite numerator parameter
        real(dp), intent(in) :: b !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z !! finite nonzero continuation target
        complex(dp) :: value

        complex(dp) :: derivative
        complex(dp) :: h
        complex(dp) :: z0
        real(dp) :: radius
        real(dp) :: remaining
        real(dp) :: step_size

        radius = abs(z)
        if (radius <= hyp1f1_origin_radius) then
            call hyp1f1_origin_series(a, b, z, value, derivative)
            return
        end if

        z0 = z * (hyp1f1_origin_radius / radius)
        call hyp1f1_origin_series(a, b, z0, value, derivative)

        do
            remaining = abs(z - z0)
            if (remaining <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, radius)) exit
            step_size = hyp1f1_step_ratio * max(abs(z0), 1.0_dp)
            if (remaining <= step_size) then
                h = z - z0
            else
                h = (z - z0) * (step_size / remaining)
            end if

            call hyp1f1_taylor_step(a, b, z0, h, value, derivative)
            z0 = z0 + h
        end do
    end function hyp1f1_continued

    pure subroutine hyp1f1_taylor_step(a, b, z0, h, value, derivative)
        real(dp), intent(in) :: a !! finite numerator parameter
        real(dp), intent(in) :: b !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z0 !! nonzero expansion center
        complex(dp), intent(in) :: h !! continuation step with abs(h) < abs(z0)
        complex(dp), intent(inout) :: value !! input/output function value
        complex(dp), intent(inout) :: derivative !! input/output derivative with respect to z

        complex(dp) :: c_nm1
        complex(dp) :: c_nm2
        complex(dp) :: c_next
        complex(dp) :: derivative_new
        complex(dp) :: derivative_term
        complex(dp) :: h_power
        complex(dp) :: value_new
        complex(dp) :: value_term
        integer :: n
        real(dp) :: n_real

        c_nm2 = value
        c_nm1 = derivative
        value_new = value + derivative * h
        derivative_new = derivative
        h_power = h

        do n = 0, hyp1f1_max_taylor_terms - 3
            n_real = real(n, dp)
            c_next = ((n_real + a) * c_nm2 - &
                (n_real + 1.0_dp) * (n_real + b - z0) * c_nm1) / &
                (z0 * (n_real + 2.0_dp) * (n_real + 1.0_dp))
            h_power = h_power * h
            value_term = c_next * h_power
            derivative_term = (n_real + 2.0_dp) * c_next * (h_power / h)
            value_new = value_new + value_term
            derivative_new = derivative_new + derivative_term

            if (n >= 10 .and. abs(value_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value_new)) .and. &
                    abs(derivative_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(derivative_new))) exit

            c_nm2 = c_nm1
            c_nm1 = c_next
        end do

        value = value_new
        derivative = derivative_new
    end subroutine hyp1f1_taylor_step

    pure function hyp1f1_real_asymptotic_positive(a, b, z) result(value)
        real(dp), intent(in) :: a !! finite noninteger-or-positive numerator parameter
        real(dp), intent(in) :: b !! finite non-pole denominator parameter
        real(dp), intent(in) :: z !! large positive argument
        real(dp) :: value

        complex(dp) :: series_complex
        real(dp) :: log_abs
        real(dp) :: series_value
        real(dp) :: sign_value

        call hyp1f1_asymptotic_series1(a, b, cmplx(z, 0.0_dp, kind=dp), series_complex)
        series_value = real(series_complex, dp)
        if (series_value == 0.0_dp) then
            value = real(hyp1f1_continued(a, b, cmplx(z, 0.0_dp, kind=dp)), dp)
            return
        end if

        log_abs = gammaln(b) - gammaln(a) + z + (a - b) * log(z) + log(abs(series_value))
        sign_value = gammasgn(b) * gammasgn(a) * sign(1.0_dp, series_value)
        value = signed_exp_from_log(log_abs, sign_value)
    end function hyp1f1_real_asymptotic_positive

    pure function hyp1f1_real_asymptotic_negative(a, b, z) result(value)
        real(dp), intent(in) :: a !! finite numerator parameter
        real(dp), intent(in) :: b !! finite non-pole denominator parameter
        real(dp), intent(in) :: z !! large negative argument
        real(dp) :: value

        complex(dp) :: series_complex
        real(dp) :: ap
        real(dp) :: log_abs
        real(dp) :: series_value
        real(dp) :: sign_value

        ap = b - a
        if (is_parameter_pole(ap)) then
            value = exp(z) * hyp1f1_polynomial_real(ap, b, -z)
            return
        end if

        call hyp1f1_asymptotic_series2(a, b, cmplx(z, 0.0_dp, kind=dp), series_complex)
        series_value = real(series_complex, dp)
        if (series_value == 0.0_dp) then
            value = real(hyp1f1_continued(a, b, cmplx(z, 0.0_dp, kind=dp)), dp)
            return
        end if

        log_abs = gammaln(b) - gammaln(b - a) - a * log(-z) + log(abs(series_value))
        sign_value = gammasgn(b) * gammasgn(b - a) * sign(1.0_dp, series_value)
        value = signed_exp_from_log(log_abs, sign_value)
    end function hyp1f1_real_asymptotic_negative

    pure function hyp1f1_complex_asymptotic(a, b, z) result(value)
        real(dp), intent(in) :: a !! finite numerator parameter
        real(dp), intent(in) :: b !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z !! large finite nonreal argument
        complex(dp) :: value

        complex(dp) :: exponent1
        complex(dp) :: exponent2
        complex(dp) :: series1
        complex(dp) :: series2
        complex(dp) :: term1
        complex(dp) :: term2
        real(dp) :: phase_sign
        real(dp), parameter :: pi = acos(-1.0_dp)

        call hyp1f1_asymptotic_series1(a, b, z, series1)
        exponent1 = z + (a - b) * log(z) + &
            cmplx(gammaln(b) - gammaln(a), 0.0_dp, kind=dp)
        term1 = (gammasgn(b) * gammasgn(a)) * exp(exponent1) * series1

        if (is_parameter_pole(b - a)) then
            term2 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        else
            call hyp1f1_asymptotic_series2(a, b, z, series2)
            if (aimag(z) >= 0.0_dp) then
                phase_sign = 1.0_dp
            else
                phase_sign = -1.0_dp
            end if
            exponent2 = -a * log(z) + &
                cmplx(gammaln(b) - gammaln(b - a), phase_sign * pi * a, kind=dp)
            term2 = (gammasgn(b) * gammasgn(b - a)) * exp(exponent2) * series2
        end if

        value = term1 + term2
    end function hyp1f1_complex_asymptotic

    pure subroutine hyp1f1_asymptotic_series1(a, b, z, series_value)
        real(dp), intent(in) :: a !! finite numerator parameter
        real(dp), intent(in) :: b !! finite denominator parameter
        complex(dp), intent(in) :: z !! large nonzero argument
        complex(dp), intent(out) :: series_value !! optimally truncated first asymptotic series

        complex(dp) :: term
        real(dp) :: previous_magnitude
        real(dp) :: term_magnitude
        integer :: n

        series_value = cmplx(1.0_dp, 0.0_dp, kind=dp)
        term = series_value
        previous_magnitude = abs(term)

        do n = 1, hyp1f1_max_asymptotic_terms
            term = term * ((b - a + real(n - 1, dp)) * &
                (1.0_dp - a + real(n - 1, dp))) / (real(n, dp) * z)
            term_magnitude = abs(term)
            if (n > 3 .and. term_magnitude > previous_magnitude) exit
            series_value = series_value + term
            if (term_magnitude <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(series_value))) exit
            previous_magnitude = term_magnitude
        end do
    end subroutine hyp1f1_asymptotic_series1

    pure subroutine hyp1f1_asymptotic_series2(a, b, z, series_value)
        real(dp), intent(in) :: a !! finite numerator parameter
        real(dp), intent(in) :: b !! finite denominator parameter
        complex(dp), intent(in) :: z !! large nonzero argument
        complex(dp), intent(out) :: series_value !! optimally truncated second asymptotic series

        complex(dp) :: term
        real(dp) :: previous_magnitude
        real(dp) :: term_magnitude
        integer :: n

        series_value = cmplx(1.0_dp, 0.0_dp, kind=dp)
        term = series_value
        previous_magnitude = abs(term)

        do n = 1, hyp1f1_max_asymptotic_terms
            term = term * ((a + real(n - 1, dp)) * &
                (1.0_dp + a - b + real(n - 1, dp))) / (real(n, dp) * (-z))
            term_magnitude = abs(term)
            if (n > 3 .and. term_magnitude > previous_magnitude) exit
            series_value = series_value + term
            if (term_magnitude <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(series_value))) exit
            previous_magnitude = term_magnitude
        end do
    end subroutine hyp1f1_asymptotic_series2

    pure function hyp1f1_polynomial_real(a, b, z) result(value)
        real(dp), intent(in) :: a !! nonpositive integer numerator parameter
        real(dp), intent(in) :: b !! denominator parameter with no encountered pole
        real(dp), intent(in) :: z !! real argument
        real(dp) :: value

        integer :: m
        integer :: n
        real(dp) :: term

        m = max(0, int(-anint(a)))
        term = 1.0_dp
        value = 1.0_dp
        do n = 1, m
            term = term * ((a + real(n - 1, dp)) / (b + real(n - 1, dp))) * &
                (z / real(n, dp))
            value = value + term
        end do
    end function hyp1f1_polynomial_real

    pure function hyp1f1_polynomial_complex(a, b, z) result(value)
        real(dp), intent(in) :: a !! negative integer numerator parameter
        real(dp), intent(in) :: b !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z !! complex argument
        complex(dp) :: value

        integer :: m
        integer :: n
        complex(dp) :: term

        m = max(0, int(-anint(a)))
        term = cmplx(1.0_dp, 0.0_dp, kind=dp)
        value = term
        do n = 1, m
            term = term * ((a + real(n - 1, dp)) / (b + real(n - 1, dp))) * &
                (z / real(n, dp))
            value = value + term
        end do
    end function hyp1f1_polynomial_complex

    pure elemental function real_pole_has_finite_polynomial(a, b) result(is_finite_case)
        real(dp), intent(in) :: a !! numerator parameter
        real(dp), intent(in) :: b !! nonpositive-integer denominator parameter
        logical :: is_finite_case

        is_finite_case = is_negative_integer(a) .and. a >= b
    end function real_pole_has_finite_polynomial

    pure function hyp1f1_real_infinite_argument(a, b, z) result(value)
        real(dp), intent(in) :: a !! finite numerator parameter
        real(dp), intent(in) :: b !! finite non-pole denominator parameter
        real(dp), intent(in) :: z !! positive or negative infinity
        real(dp) :: value

        real(dp) :: sign_value

        if (z > 0.0_dp) then
            sign_value = gammasgn(b) * gammasgn(a)
            value = sign(positive_infinity(z), sign_value)
        else if (is_parameter_pole(b - a)) then
            value = 0.0_dp
        else if (a > 0.0_dp) then
            value = 0.0_dp
        else
            sign_value = gammasgn(b) * gammasgn(b - a)
            value = sign(positive_infinity(z), sign_value)
        end if
    end function hyp1f1_real_infinite_argument

    pure function signed_exp_from_log(log_abs, sign_value) result(value)
        real(dp), intent(in) :: log_abs !! natural logarithm of absolute result
        real(dp), intent(in) :: sign_value !! sign carrier
        real(dp) :: value

        if (log_abs > log(huge(1.0_dp))) then
            value = sign(positive_infinity(log_abs), sign_value)
        else if (log_abs < log(tiny(1.0_dp))) then
            value = sign(0.0_dp, sign_value)
        else
            value = sign(1.0_dp, sign_value) * exp(log_abs)
        end if
    end function signed_exp_from_log

    pure function real_exp_safely(z) result(value)
        real(dp), intent(in) :: z !! real exponential argument
        real(dp) :: value

        if (z > log(huge(1.0_dp))) then
            value = positive_infinity(z)
        else
            value = exp(z)
        end if
    end function real_exp_safely

    pure function hyp0f1_real_series(v, z) result(value)
        real(dp), intent(in) :: v !! finite non-pole denominator parameter
        real(dp), intent(in) :: z !! finite real argument
        real(dp) :: value

        integer :: n
        real(dp) :: compensation
        real(dp) :: corrected
        real(dp) :: term
        real(dp) :: updated

        term = 1.0_dp
        value = 1.0_dp
        compensation = 0.0_dp

        do n = 1, max_origin_terms
            term = (term * (z / real(n, dp))) / (v + real(n - 1, dp))
            corrected = term - compensation
            updated = value + corrected
            compensation = (updated - value) - corrected
            value = updated

            if (n >= 12 .and. abs(term) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value))) exit
            if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(term)) exit
        end do
    end function hyp0f1_real_series

    pure subroutine hyp0f1_origin_series(v, z, value, derivative)
        real(dp), intent(in) :: v !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z !! finite nonzero argument inside the origin disk
        complex(dp), intent(out) :: value !! 0F1 value at z
        complex(dp), intent(out) :: derivative !! derivative with respect to z

        complex(dp) :: compensation
        complex(dp) :: corrected
        complex(dp) :: derivative_term
        complex(dp) :: term
        complex(dp) :: updated
        integer :: n

        term = cmplx(1.0_dp, 0.0_dp, kind=dp)
        value = term
        derivative = cmplx(0.0_dp, 0.0_dp, kind=dp)
        compensation = cmplx(0.0_dp, 0.0_dp, kind=dp)

        do n = 1, max_origin_terms
            term = (term * (z / real(n, dp))) / (v + real(n - 1, dp))
            corrected = term - compensation
            updated = value + corrected
            compensation = (updated - value) - corrected
            value = updated

            derivative_term = real(n, dp) * term / z
            derivative = derivative + derivative_term

            if (n >= 12 .and. abs(term) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value)) .and. &
                    abs(derivative_term) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(derivative))) exit
        end do
    end subroutine hyp0f1_origin_series

    pure subroutine hyp0f1_taylor_step(v, z0, h, value, derivative)
        real(dp), intent(in) :: v !! finite non-pole denominator parameter
        complex(dp), intent(in) :: z0 !! nonzero expansion center
        complex(dp), intent(in) :: h !! outward continuation step with abs(h) < abs(z0)
        complex(dp), intent(inout) :: value !! input/output function value
        complex(dp), intent(inout) :: derivative !! input/output derivative with respect to z

        complex(dp) :: a_nm1
        complex(dp) :: a_nm2
        complex(dp) :: a_next
        complex(dp) :: derivative_new
        complex(dp) :: derivative_term
        complex(dp) :: h_power
        complex(dp) :: value_new
        complex(dp) :: value_term
        integer :: n
        real(dp) :: n_real

        a_nm2 = value
        a_nm1 = derivative
        value_new = value + derivative * h
        derivative_new = derivative
        h_power = h

        do n = 0, max_taylor_terms - 3
            n_real = real(n, dp)
            a_next = (a_nm2 - (n_real + 1.0_dp) * (n_real + v) * a_nm1) / &
                (z0 * (n_real + 2.0_dp) * (n_real + 1.0_dp))
            h_power = h_power * h
            value_term = a_next * h_power
            derivative_term = (n_real + 2.0_dp) * a_next * (h_power / h)
            value_new = value_new + value_term
            derivative_new = derivative_new + derivative_term

            if (n >= 10 .and. abs(value_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value_new)) .and. &
                    abs(derivative_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(derivative_new))) exit

            a_nm2 = a_nm1
            a_nm1 = a_next
        end do

        value = value_new
        derivative = derivative_new
    end subroutine hyp0f1_taylor_step

    pure function hyp0f1_from_real_bessel(v, magnitude_z, bessel_value) result(value)
        real(dp), intent(in) :: v !! finite non-pole denominator parameter
        real(dp), intent(in) :: magnitude_z !! absolute real argument, > 0
        real(dp), intent(in) :: bessel_value !! finite nonzero Bessel J or I value
        real(dp) :: value

        real(dp) :: log_abs
        real(dp) :: sign_value

        log_abs = gammaln(v) + 0.5_dp * (1.0_dp - v) * log(magnitude_z) + &
            log(abs(bessel_value))
        sign_value = gammasgn(v) * sign(1.0_dp, bessel_value)

        if (log_abs > log(huge(1.0_dp))) then
            value = sign(positive_infinity(magnitude_z), sign_value)
        else
            value = sign_value * exp(log_abs)
        end if
    end function hyp0f1_from_real_bessel

    pure elemental function is_parameter_pole(v) result(is_pole)
        real(dp), intent(in) :: v !! denominator parameter tested for a nonpositive-integer pole
        logical :: is_pole

        is_pole = ieee_is_finite(v) .and. v <= 0.0_dp .and. v == anint(v)
    end function is_parameter_pole

    pure elemental function is_negative_integer(v) result(is_integer)
        real(dp), intent(in) :: v !! value tested for a negative integer
        logical :: is_integer

        is_integer = ieee_is_finite(v) .and. v < 0.0_dp .and. v == anint(v)
    end function is_negative_integer

    pure elemental function complex_has_nan(z) result(has_nan)
        complex(dp), intent(in) :: z !! complex value tested componentwise
        logical :: has_nan

        has_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
    end function complex_has_nan

    pure elemental function complex_is_finite(z) result(is_finite)
        complex(dp), intent(in) :: z !! complex value tested componentwise
        logical :: is_finite

        is_finite = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function complex_is_finite

    pure elemental function complex_nan(both_components) result(value)
        logical, intent(in) :: both_components !! true to return NaN in both components
        complex(dp) :: value

        real(dp) :: qnan

        qnan = quiet_nan(0.0_dp)
        if (both_components) then
            value = cmplx(qnan, qnan, kind=dp)
        else
            value = cmplx(qnan, 0.0_dp, kind=dp)
        end if
    end function complex_nan

end module scifort_hypergeometric
