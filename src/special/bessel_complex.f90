! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Complex cylindrical and modified Bessel functions for real order.
! J/Y are reconstructed from scaled Hankel pairs, with direct defining
! series used where large-order cancellation would dominate. I likewise
! uses its defining series in the high-order/turning regime and otherwise
! uses the DLMF 10.27 connection to J; K uses the DLMF 10.27 connection
! to scaled Hankel functions. Derivatives use the standard shifted-order
! finite sums from DLMF 10.6.7 and 10.29.5.

module scifort_bessel_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_hankel, only : hankel1e, hankel2e
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    complex(dp), parameter :: imaginary_unit = (0.0_dp, 1.0_dp)
    integer, parameter :: derivative_weight_limit = 1000
    integer, parameter :: family_i = 1
    integer, parameter :: family_j = 2
    integer, parameter :: family_k = 3
    integer, parameter :: family_y = 4

    public :: iv_complex
    public :: ive_complex
    public :: ivp_complex
    public :: jn_complex
    public :: jv_complex
    public :: jve_complex
    public :: jvp_complex
    public :: kv_complex
    public :: kve_complex
    public :: kvp_complex
    public :: yv_complex
    public :: yve_complex
    public :: yvp_complex

contains

    pure elemental function jn_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        value = jv_complex(v, z)
    end function jn_complex

    pure elemental function jv_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        if (invalid_inputs(v, z)) then
            value = complex_nan()
        else if (complex_is_zero(z)) then
            value = j_or_i_at_zero(v)
        else if (use_j_series(v, z)) then
            value = jv_series_order(v, z)
        else
            value = scale_complex_real(jve_complex(v, z), abs(aimag(z)))
        end if
    end function jv_complex

    pure elemental function yv_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument on the principal branch
        complex(dp) :: value

        if (invalid_inputs(v, z)) then
            value = complex_nan()
        else if (complex_is_zero(z)) then
            value = y_at_zero(v)
        else
            value = scale_complex_real(yve_complex(v, z), abs(aimag(z)))
        end if
    end function yv_complex

    pure elemental function jve_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real order of exp(-abs(aimag(z)))*J_v(z)
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        complex(dp) :: term1
        complex(dp) :: term2
        real(dp) :: x
        real(dp) :: y

        if (invalid_inputs(v, z)) then
            value = complex_nan()
            return
        end if
        if (complex_is_zero(z)) then
            value = j_or_i_at_zero(v)
            return
        end if
        if (use_j_series(v, z)) then
            value = exp(-abs(aimag(z))) * jv_series_order(v, z)
            return
        end if

        x = real(z, dp)
        y = aimag(z)
        if (y >= 0.0_dp) then
            term1 = hankel1e(v, z) * exp(-2.0_dp * y) * unit_phase(x)
            term2 = hankel2e(v, z) * unit_phase(-x)
        else
            term1 = hankel1e(v, z) * unit_phase(x)
            term2 = hankel2e(v, z) * exp(2.0_dp * y) * unit_phase(-x)
        end if

        if (complex_is_nan(term1) .or. complex_is_nan(term2)) then
            value = complex_nan()
        else
            value = 0.5_dp * (term1 + term2)
        end if
    end function jve_complex

    pure elemental function yve_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real order of exp(-abs(aimag(z)))*Y_v(z)
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        complex(dp) :: value

        complex(dp) :: term1
        complex(dp) :: term2
        real(dp) :: x
        real(dp) :: y

        if (invalid_inputs(v, z) .or. complex_is_zero(z)) then
            value = complex_nan()
            return
        end if

        x = real(z, dp)
        y = aimag(z)
        if (y >= 0.0_dp) then
            term1 = hankel1e(v, z) * exp(-2.0_dp * y) * unit_phase(x)
            term2 = hankel2e(v, z) * unit_phase(-x)
        else
            term1 = hankel1e(v, z) * unit_phase(x)
            term2 = hankel2e(v, z) * exp(2.0_dp * y) * unit_phase(-x)
        end if

        if (complex_is_nan(term1) .or. complex_is_nan(term2)) then
            value = complex_nan()
        else
            value = (term1 - term2) / (2.0_dp * imaginary_unit)
        end if
    end function yve_complex

    pure elemental function iv_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-I order
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        if (invalid_inputs(v, z)) then
            value = complex_nan()
        else if (complex_is_zero(z)) then
            value = j_or_i_at_zero(v)
        else if (use_i_series(v, z)) then
            value = iv_series_nonnegative(v, z)
        else
            value = scale_complex_real(ive_complex(v, z), abs(real(z, dp)))
        end if
    end function iv_complex

    pure elemental function ive_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real order of exp(-abs(real(z)))*I_v(z)
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        if (invalid_inputs(v, z)) then
            value = complex_nan()
        else if (complex_is_zero(z)) then
            value = j_or_i_at_zero(v)
        else if (use_i_series(v, z)) then
            value = exp(-abs(real(z, dp))) * iv_series_nonnegative(v, z)
        else if (aimag(z) >= 0.0_dp) then
            value = phase_half_pi(v, 1) * jve_complex(v, -imaginary_unit * z)
        else
            value = phase_half_pi(v, -1) * jve_complex(v, imaginary_unit * z)
        end if
    end function ive_complex

    pure elemental function kv_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-K order
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        complex(dp) :: value

        if (invalid_inputs(v, z) .or. complex_is_zero(z)) then
            value = complex_nan()
        else
            value = multiply_complex_exp(kve_complex(v, z), -z)
        end if
    end function kv_complex

    pure elemental function kve_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real order of exp(z)*K_v(z)
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        complex(dp) :: value

        if (invalid_inputs(v, z) .or. complex_is_zero(z)) then
            value = complex_nan()
        else if (aimag(z) >= 0.0_dp) then
            value = -0.5_dp * scifort_pi * imaginary_unit * phase_half_pi(v, -1) * &
                hankel2e(v, -imaginary_unit * z)
        else
            value = 0.5_dp * scifort_pi * imaginary_unit * phase_half_pi(v, 1) * &
                hankel1e(v, imaginary_unit * z)
        end if
    end function kve_complex

    pure elemental function jvp_complex(v, z, n) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order before differentiation
        complex(dp), intent(in) :: z !! finite complex argument
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value
        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n
        if (derivative_order < 0) then
            value = complex_nan()
        else if (derivative_order == 0) then
            value = jv_complex(v, z)
        else
            value = bessel_derivative_sum(v, z, derivative_order, family_j, -1)
        end if
    end function jvp_complex

    pure elemental function yvp_complex(v, z, n) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order before differentiation
        complex(dp), intent(in) :: z !! finite complex argument on the principal branch
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value
        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n
        if (derivative_order < 0) then
            value = complex_nan()
        else if (derivative_order == 0) then
            value = yv_complex(v, z)
        else
            value = bessel_derivative_sum(v, z, derivative_order, family_y, -1)
        end if
    end function yvp_complex

    pure elemental function ivp_complex(v, z, n) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-I order before differentiation
        complex(dp), intent(in) :: z !! finite complex argument
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value
        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n
        if (derivative_order < 0) then
            value = complex_nan()
        else if (derivative_order == 0) then
            value = iv_complex(v, z)
        else
            value = bessel_derivative_sum(v, z, derivative_order, family_i, 1)
        end if
    end function ivp_complex

    pure elemental function kvp_complex(v, z, n) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-K order before differentiation
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value
        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n
        if (derivative_order < 0) then
            value = complex_nan()
        else if (derivative_order == 0) then
            value = kv_complex(v, z)
        else
            value = bessel_derivative_sum(v, z, derivative_order, family_k, 1)
            if (mod(derivative_order, 2) /= 0) value = -value
        end if
    end function kvp_complex

    pure function bessel_derivative_sum(v, z, n, family, coefficient_phase) result(value)
        real(dp), intent(in) :: v !! real Bessel order before differentiation
        complex(dp), intent(in) :: z !! complex Bessel argument
        integer, intent(in) :: n !! positive derivative order
        integer, intent(in) :: family !! internal Bessel-family selector
        integer, intent(in) :: coefficient_phase !! shifted-order coefficient phase, -1 or +1
        complex(dp) :: value

        complex(dp) :: base_value
        complex(dp) :: scaled_sum
        complex(dp) :: term
        integer :: k
        real(dp) :: abs_term
        real(dp) :: log_weight
        real(dp) :: max_term
        real(dp) :: order
        real(dp) :: weight

        max_term = 0.0_dp
        scaled_sum = (0.0_dp, 0.0_dp)

        if (n <= derivative_weight_limit) then
            weight = scale(1.0_dp, -n)
            do k = 0, n
                if (k > 0) then
                    weight = real(coefficient_phase, dp) * weight * &
                        real(n - k + 1, dp) / real(k, dp)
                end if
                if (weight == 0.0_dp) cycle
                order = v - real(n, dp) + 2.0_dp * real(k, dp)
                base_value = bessel_family_value(family, order, z)
                if (.not. complex_is_finite(base_value)) then
                    value = complex_nan()
                    return
                end if
                term = weight * base_value
                if (term == (0.0_dp, 0.0_dp)) cycle
                abs_term = abs(term)
                if (max_term == 0.0_dp) then
                    max_term = abs_term
                    scaled_sum = term / abs_term
                else if (abs_term > max_term) then
                    scaled_sum = scaled_sum * (max_term / abs_term) + term / abs_term
                    max_term = abs_term
                else
                    scaled_sum = scaled_sum + term / max_term
                end if
            end do
        else
            do k = 0, n
                log_weight = log_gamma(real(n, dp) + 1.0_dp) - &
                    log_gamma(real(k, dp) + 1.0_dp) - &
                    log_gamma(real(n - k, dp) + 1.0_dp) - real(n, dp) * log(2.0_dp)
                weight = exp(log_weight)
                if (coefficient_phase < 0 .and. mod(k, 2) /= 0) weight = -weight
                if (weight == 0.0_dp) cycle
                order = v - real(n, dp) + 2.0_dp * real(k, dp)
                base_value = bessel_family_value(family, order, z)
                if (.not. complex_is_finite(base_value)) then
                    value = complex_nan()
                    return
                end if
                term = weight * base_value
                if (term == (0.0_dp, 0.0_dp)) cycle
                abs_term = abs(term)
                if (max_term == 0.0_dp) then
                    max_term = abs_term
                    scaled_sum = term / abs_term
                else if (abs_term > max_term) then
                    scaled_sum = scaled_sum * (max_term / abs_term) + term / abs_term
                    max_term = abs_term
                else
                    scaled_sum = scaled_sum + term / max_term
                end if
            end do
        end if

        if (max_term == 0.0_dp .or. scaled_sum == (0.0_dp, 0.0_dp)) then
            value = (0.0_dp, 0.0_dp)
        else
            value = max_term * scaled_sum
        end if
    end function bessel_derivative_sum

    pure function bessel_family_value(family, v, z) result(value)
        integer, intent(in) :: family !! internal Bessel-family selector
        real(dp), intent(in) :: v !! shifted real Bessel order
        complex(dp), intent(in) :: z !! complex Bessel argument
        complex(dp) :: value

        select case (family)
        case (family_i)
            value = iv_complex(v, z)
        case (family_j)
            value = jv_complex(v, z)
        case (family_k)
            value = kv_complex(v, z)
        case (family_y)
            value = yv_complex(v, z)
        case default
            value = complex_nan()
        end select
    end function bessel_family_value

    pure elemental logical function use_j_series(v, z) result(use_series)
        real(dp), intent(in) :: v !! finite real cylindrical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        real(dp) :: nu

        use_series = .false.
        if (v >= 0.0_dp) then
            nu = v
        else if (is_integer_value(v)) then
            nu = -v
        else
            return
        end if
        use_series = abs(z) <= 0.5_dp * (nu + 1.0_dp)
    end function use_j_series

    pure elemental logical function use_i_series(v, z) result(use_series)
        real(dp), intent(in) :: v !! finite real modified-Bessel-I order
        complex(dp), intent(in) :: z !! finite nonzero complex argument

        use_series = v >= 0.0_dp .and. abs(z) <= 2.25_dp * (v + 1.0_dp)
    end function use_i_series

    pure elemental function iv_series_nonnegative(v, z) result(value)
        real(dp), intent(in) :: v !! finite nonnegative modified-Bessel-I order
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        complex(dp) :: value

        integer, parameter :: max_terms = 100000
        complex(dp) :: factor
        complex(dp) :: series_sum
        complex(dp) :: term
        complex(dp) :: prefactor
        integer :: k
        real(dp) :: tolerance

        tolerance = 4.0_dp * epsilon(1.0_dp)
        factor = 0.25_dp * z * z
        term = (1.0_dp, 0.0_dp)
        series_sum = term
        do k = 1, max_terms
            term = term * factor / (real(k, dp) * (v + real(k, dp)))
            series_sum = series_sum + term
            if (abs(term) <= tolerance * max(1.0_dp, abs(series_sum))) exit
        end do

        if (v == 0.0_dp) then
            prefactor = (1.0_dp, 0.0_dp)
        else
            prefactor = exp(v * log(0.5_dp * z) - log_gamma(v + 1.0_dp))
        end if
        value = prefactor * series_sum
    end function iv_series_nonnegative

    pure elemental function jv_series_order(v, z) result(value)
        real(dp), intent(in) :: v !! finite real order eligible for direct series evaluation
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        complex(dp) :: value

        real(dp) :: nu
        real(dp) :: parity

        if (v < 0.0_dp) then
            nu = -v
            parity = 1.0_dp
            if (odd_integer_value(nu)) parity = -1.0_dp
            value = parity * jv_series_nonnegative(nu, z)
        else
            value = jv_series_nonnegative(v, z)
        end if
    end function jv_series_order

    pure elemental function jv_series_nonnegative(v, z) result(value)
        real(dp), intent(in) :: v !! finite nonnegative cylindrical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        complex(dp) :: value

        integer, parameter :: max_terms = 100000
        complex(dp) :: factor
        complex(dp) :: series_sum
        complex(dp) :: term
        complex(dp) :: prefactor
        integer :: k
        real(dp) :: tolerance

        tolerance = 4.0_dp * epsilon(1.0_dp)
        factor = -0.25_dp * z * z
        term = (1.0_dp, 0.0_dp)
        series_sum = term
        do k = 1, max_terms
            term = term * factor / (real(k, dp) * (v + real(k, dp)))
            series_sum = series_sum + term
            if (abs(term) <= tolerance * max(1.0_dp, abs(series_sum))) exit
        end do

        if (v == 0.0_dp) then
            prefactor = (1.0_dp, 0.0_dp)
        else
            prefactor = exp(v * log(0.5_dp * z) - log_gamma(v + 1.0_dp))
        end if
        value = prefactor * series_sum
    end function jv_series_nonnegative

    pure elemental function j_or_i_at_zero(v) result(value)
        real(dp), intent(in) :: v !! finite real J/I order at zero
        complex(dp) :: value

        if (v == 0.0_dp) then
            value = (1.0_dp, 0.0_dp)
        else if (v > 0.0_dp .or. is_integer_value(v)) then
            value = (0.0_dp, 0.0_dp)
        else
            value = complex_nan()
        end if
    end function j_or_i_at_zero

    pure elemental function y_at_zero(v) result(value)
        real(dp), intent(in) :: v !! finite real Y order at zero
        complex(dp) :: value
        real(dp) :: infinity_value

        if (v < 0.0_dp .and. .not. is_integer_value(v)) then
            value = complex_nan()
            return
        end if

        infinity_value = negative_infinity(0.0_dp)
        if (v < 0.0_dp .and. odd_integer_value(v)) then
            infinity_value = positive_infinity(0.0_dp)
        end if
        value = cmplx(infinity_value, 0.0_dp, kind=dp)
    end function y_at_zero

    pure elemental function phase_half_pi(v, sign_phase) result(value)
        real(dp), intent(in) :: v !! finite real order in exp(sign*i*pi*v/2)
        integer, intent(in) :: sign_phase !! phase sign, -1 or +1
        complex(dp) :: value
        real(dp) :: angle
        real(dp) :: reduced_v

        reduced_v = modulo(v, 4.0_dp)
        angle = 0.5_dp * real(sign_phase, dp) * scifort_pi * reduced_v
        value = unit_phase(angle)
    end function phase_half_pi

    pure elemental function unit_phase(angle) result(value)
        real(dp), intent(in) :: angle !! finite phase angle in radians
        complex(dp) :: value

        value = cmplx(cos(angle), sin(angle), kind=dp)
    end function unit_phase

    pure elemental function scale_complex_real(value, exponent) result(scaled)
        complex(dp), intent(in) :: value !! finite complex value to multiply by exp(exponent)
        real(dp), intent(in) :: exponent !! nonnegative real logarithmic scale
        complex(dp) :: scaled

        real(dp) :: factor
        real(dp) :: infinity_value
        real(dp) :: log_huge
        real(dp) :: xr
        real(dp) :: yi

        if (complex_is_nan(value) .or. ieee_is_nan(exponent)) then
            scaled = complex_nan()
            return
        end if
        if (.not. complex_is_finite(value)) then
            scaled = value
            return
        end if

        log_huge = log(huge(1.0_dp))
        if (exponent <= log_huge) then
            factor = exp(exponent)
            scaled = factor * value
            return
        end if

        infinity_value = positive_infinity(0.0_dp)
        xr = real(value, dp)
        yi = aimag(value)
        if (xr == 0.0_dp) then
            xr = 0.0_dp
        else
            xr = sign(infinity_value, xr)
        end if
        if (yi == 0.0_dp) then
            yi = 0.0_dp
        else
            yi = sign(infinity_value, yi)
        end if
        scaled = cmplx(xr, yi, kind=dp)
    end function scale_complex_real

    pure elemental function multiply_complex_exp(value, exponent) result(scaled)
        complex(dp), intent(in) :: value !! finite complex value to multiply
        complex(dp), intent(in) :: exponent !! finite complex exponent
        complex(dp) :: scaled

        complex(dp) :: rotated
        real(dp) :: real_exponent

        if (complex_is_nan(value) .or. complex_is_nan(exponent)) then
            scaled = complex_nan()
            return
        end if
        if (.not. complex_is_finite(value) .or. .not. complex_is_finite(exponent)) then
            scaled = complex_nan()
            return
        end if

        rotated = value * unit_phase(aimag(exponent))
        real_exponent = real(exponent, dp)
        if (real_exponent >= 0.0_dp) then
            scaled = scale_complex_real(rotated, real_exponent)
        else if (real_exponent < log(tiny(1.0_dp))) then
            scaled = (0.0_dp, 0.0_dp)
        else
            scaled = exp(real_exponent) * rotated
        end if
    end function multiply_complex_exp

    pure elemental logical function invalid_inputs(v, z) result(invalid)
        real(dp), intent(in) :: v !! Bessel order to classify
        complex(dp), intent(in) :: z !! complex Bessel argument to classify

        invalid = ieee_is_nan(v) .or. .not. ieee_is_finite(v) .or. &
            .not. complex_is_finite(z)
    end function invalid_inputs

    pure elemental logical function complex_is_zero(z) result(is_zero)
        complex(dp), intent(in) :: z !! complex value to compare with zero

        is_zero = real(z, dp) == 0.0_dp .and. aimag(z) == 0.0_dp
    end function complex_is_zero

    pure elemental logical function complex_is_nan(z) result(is_nan)
        complex(dp), intent(in) :: z !! complex value to classify

        is_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
    end function complex_is_nan

    pure elemental logical function complex_is_finite(z) result(is_finite)
        complex(dp), intent(in) :: z !! complex value to classify

        is_finite = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function complex_is_finite

    pure elemental logical function is_integer_value(v) result(is_integer)
        real(dp), intent(in) :: v !! finite real order to classify

        is_integer = v == anint(v)
    end function is_integer_value

    pure elemental logical function odd_integer_value(v) result(is_odd)
        real(dp), intent(in) :: v !! finite integer-valued real order

        is_odd = modulo(abs(v), 2.0_dp) == 1.0_dp
    end function odd_integer_value

    pure elemental function complex_nan() result(value)
        complex(dp) :: value
        real(dp) :: nan_value

        nan_value = quiet_nan(0.0_dp)
        value = cmplx(nan_value, nan_value, kind=dp)
    end function complex_nan

end module scifort_bessel_complex
