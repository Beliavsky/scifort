! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real Struve and modified Struve functions plus the historical order-zero
! integral helpers. The implementation is independent of SciPy source code and
! follows the defining series, Bessel expansion, recurrences, integral identity,
! and large-argument asymptotics in NIST DLMF Chapter 11.
module scifort_struve
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_jy, only : jv, yv
    use scifort_bessel_order, only : iv
    use scifort_constants, only : scifort_pi
    use scifort_gamma_family, only : gammasgn
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_special_elementary, only : gammaln
    implicit none
    private

    real(dp), parameter :: half_pi = 0.5_dp * scifort_pi
    real(dp), parameter :: integration_anchor = 8.0_dp
    real(dp), parameter :: integration_panel_width = 4.0_dp
    integer, parameter :: max_series_terms = 10000

    ! Eight-point Gauss-Legendre nodes and weights generated from the defining
    ! Legendre polynomial and verified by exact moments through degree 15.
    real(dp), parameter :: gl_nodes(4) = [ &
        0.183434642495649780_dp, 0.525532409916328991_dp, &
        0.796666477413626728_dp, 0.960289856497536176_dp]
    real(dp), parameter :: gl_weights(4) = [ &
        0.362683783378361657_dp, 0.313706645877886881_dp, &
        0.222381034453374427_dp, 0.101228536290377064_dp]

    public :: it2struve0
    public :: itmodstruve0
    public :: itstruve0
    public :: modstruve
    public :: struve

contains

    pure elemental function struve(v, x) result(value)
        real(dp), intent(in) :: v !! real Struve order
        real(dp), intent(in) :: x !! real argument; negative values require integer order
        real(dp) :: value

        value = struve_core(v, x, .false.)
    end function struve

    pure elemental function modstruve(v, x) result(value)
        real(dp), intent(in) :: v !! real modified-Struve order
        real(dp), intent(in) :: x !! real argument; negative values require integer order
        real(dp) :: value

        value = struve_core(v, x, .true.)
    end function modstruve

    pure elemental function itstruve0(x) result(value)
        real(dp), intent(in) :: x !! upper integration limit in integral_0^x H_0(t) dt
        real(dp) :: value
        real(dp) :: ax

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(x)
            return
        end if

        ax = abs(x)
        if (ax <= integration_anchor) then
            value = h0_integral_series(ax, .false.)
        else
            value = h0_integral_series(integration_anchor, .false.) + &
                composite_gauss_h0(integration_anchor, ax, .false.)
        end if
        ! H_0 is odd, so its integral from zero is even.
    end function itstruve0

    pure elemental function it2struve0(x) result(value)
        real(dp), intent(in) :: x !! lower limit in integral_x^infinity H_0(t)/t dt
        real(dp) :: value
        real(dp) :: ax
        real(dp) :: positive_value

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(x)
            return
        end if

        ax = abs(x)
        if (ax <= integration_anchor) then
            positive_value = half_pi - h0_integral_series(ax, .true.)
        else
            positive_value = half_pi - h0_integral_series(integration_anchor, .true.) - &
                composite_gauss_h0(integration_anchor, ax, .true.)
        end if

        if (x < 0.0_dp) then
            ! H_0(t)/t is even and its full-line integral is pi.
            value = scifort_pi - positive_value
        else
            value = positive_value
        end if
    end function it2struve0

    pure elemental function itmodstruve0(x) result(value)
        real(dp), intent(in) :: x !! upper integration limit in integral_0^x L_0(t) dt
        real(dp) :: value
        real(dp) :: ax

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(x)
            return
        end if

        ax = abs(x)
        value = l0_integral_series(ax)
        ! L_0 is odd, so its integral from zero is even.
    end function itmodstruve0

    pure function struve_core(v, x, modified) result(value)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! real argument
        logical, intent(in) :: modified !! true for modified Struve L_v
        real(dp) :: value

        real(dp) :: ax
        real(dp) :: bessel_limit
        real(dp) :: series_limit
        integer :: parity

        if (ieee_is_nan(v) .or. ieee_is_nan(x)) then
            value = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(v) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp .and. .not. is_integer_value(v)) then
            value = quiet_nan(x)
            return
        end if
        if (x >= 0.0_dp .and. is_negative_half_integer(v)) then
            if (modified) then
                value = iv(-v, x)
            else
                value = jv(-v, x)
                if (modulo(int(-v - 0.5_dp), 2) /= 0) value = -value
            end if
            return
        end if
        if (x == 0.0_dp) then
            value = struve_at_zero(v, modified, x)
            return
        end if

        ax = abs(x)
        series_limit = max(8.0_dp, 0.5_dp * abs(v) + 6.0_dp)
        bessel_limit = 24.0_dp + 0.4_dp * abs(v)

        if (modified) then
            if (ax <= bessel_limit) then
                value = struve_power_series(v, ax, .true.)
            else
                value = modified_struve_asymptotic(v, ax)
            end if
        else if (ax <= series_limit) then
            value = struve_power_series(v, ax, .false.)
        else if (ax <= bessel_limit) then
            value = struve_bessel_series(v, ax)
        else
            value = struve_asymptotic(v, ax)
        end if

        if (x < 0.0_dp) then
            parity = int(modulo(abs(anint(v)), 2.0_dp))
            if (parity == 0) then
                value = -value
            end if
        end if
    end function struve_core

    pure function struve_at_zero(v, modified, kind_arg) result(value)
        real(dp), intent(in) :: v !! finite real order
        logical, intent(in) :: modified !! true for modified Struve L_v
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE values
        real(dp) :: value

        integer :: k0
        real(dp) :: exponent
        real(dp) :: log_abs
        real(dp) :: sign_value

        k0 = first_nonzero_series_index(v)
        exponent = 2.0_dp * real(k0, dp) + v + 1.0_dp
        if (exponent > 0.0_dp) then
            value = 0.0_dp
            return
        end if

        sign_value = gammasgn(real(k0, dp) + v + 1.5_dp)
        if (.not. modified .and. modulo(k0, 2) /= 0) sign_value = -sign_value
        if (exponent < 0.0_dp) then
            value = sign(positive_infinity(kind_arg), sign_value)
            return
        end if

        log_abs = -gammaln(real(k0, dp) + 1.5_dp) - &
            gammaln(real(k0, dp) + v + 1.5_dp)
        value = sign(exp(log_abs), sign_value)
    end function struve_at_zero

    pure function struve_power_series(v, x, modified) result(value)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! positive finite argument
        logical, intent(in) :: modified !! true for modified Struve L_v
        real(dp) :: value

        real(dp) :: denominator
        real(dp) :: log_abs
        real(dp) :: ratio
        real(dp) :: sign_value
        real(dp) :: sum_value
        real(dp) :: term
        integer :: k
        integer :: k0

        k0 = first_nonzero_series_index(v)
        k = k0
        log_abs = (2.0_dp * real(k, dp) + v + 1.0_dp) * log(0.5_dp * x) - &
            gammaln(real(k, dp) + 1.5_dp) - gammaln(real(k, dp) + v + 1.5_dp)
        sign_value = gammasgn(real(k, dp) + v + 1.5_dp)
        if (.not. modified .and. modulo(k, 2) /= 0) sign_value = -sign_value
        term = signed_exp_from_log(log_abs, sign_value, x)
        sum_value = term

        do while (k - k0 < max_series_terms)
            denominator = (real(k, dp) + 1.5_dp) * (real(k, dp) + v + 1.5_dp)
            ratio = 0.25_dp * x * x / denominator
            if (.not. modified) ratio = -ratio
            term = term * ratio
            k = k + 1
            sum_value = sum_value + term
            if (k >= k0 + 8 .and. abs(term) <= 4.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(sum_value))) exit
            if (.not. ieee_is_finite(sum_value)) exit
            if (abs(term) > huge(1.0_dp) / 8.0_dp) exit
        end do
        value = sum_value
    end function struve_power_series

    pure function struve_bessel_series(v, x) result(value)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: value

        real(dp) :: coefficient
        real(dp) :: sum_value
        real(dp) :: term
        integer :: k

        coefficient = 2.0_dp
        sum_value = coefficient * jv(v + 0.5_dp, x)
        do k = 0, 499
            coefficient = coefficient * (0.5_dp * x) / real(k + 1, dp) * &
                (real(k, dp) + 0.5_dp) / (real(k, dp) + 1.5_dp)
            term = coefficient * jv(v + real(k, dp) + 1.5_dp, x)
            sum_value = sum_value + term
            if (k >= 10 .and. abs(term) <= 4.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(sum_value))) exit
        end do
        value = sqrt(x / (2.0_dp * scifort_pi)) * sum_value
    end function struve_bessel_series

    pure function struve_asymptotic(v, x) result(value)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! positive large argument
        real(dp) :: value

        value = yv(v, x) + struve_difference_asymptotic(v, x, .false.)
    end function struve_asymptotic

    pure function modified_struve_asymptotic(v, x) result(value)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! positive large argument
        real(dp) :: value

        value = iv(v, x) + struve_difference_asymptotic(v, x, .true.)
    end function modified_struve_asymptotic

    pure function struve_difference_asymptotic(v, x, modified) result(value)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! positive large argument
        logical, intent(in) :: modified !! true for M_v=L_v-I_v; false for K_v=H_v-Y_v
        real(dp) :: value

        real(dp) :: factor
        real(dp) :: log_abs
        real(dp) :: sign_value
        real(dp) :: sum_value
        real(dp) :: term
        real(dp) :: next_term
        integer :: k

        log_abs = -0.5_dp * log(scifort_pi) + (v - 1.0_dp) * log(0.5_dp * x) - &
            gammaln(v + 0.5_dp)
        sign_value = gammasgn(v + 0.5_dp)
        if (modified) sign_value = -sign_value
        if (v + 0.5_dp <= 0.0_dp .and. v + 0.5_dp == anint(v + 0.5_dp)) then
            term = 0.0_dp
        else
            term = signed_exp_from_log(log_abs, sign_value, x)
        end if
        sum_value = term

        do k = 0, 99
            factor = 4.0_dp * (real(k, dp) + 0.5_dp) * &
                (v - 0.5_dp - real(k, dp)) / (x * x)
            if (modified) factor = -factor
            next_term = term * factor
            if (k >= 1 .and. abs(next_term) > abs(term)) exit
            sum_value = sum_value + next_term
            term = next_term
            if (abs(term) <= 2.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(sum_value))) exit
        end do
        value = sum_value
    end function struve_difference_asymptotic

    pure function first_nonzero_series_index(v) result(k0)
        real(dp), intent(in) :: v !! finite real order
        integer :: k0

        real(dp) :: a

        a = v + 1.5_dp
        if (a <= 0.0_dp .and. a == anint(a) .and. abs(a) < real(huge(k0) - 2, dp)) then
            k0 = int(-a) + 1
        else
            k0 = 0
        end if
    end function first_nonzero_series_index

    pure elemental function is_negative_half_integer(x) result(is_half_integer)
        real(dp), intent(in) :: x !! finite real value to test
        logical :: is_half_integer

        is_half_integer = ieee_is_finite(x) .and. x < 0.0_dp .and. &
            abs(x + 0.5_dp) < 2.0_dp ** 52 .and. x + 0.5_dp == anint(x + 0.5_dp)
    end function is_negative_half_integer

    pure elemental function is_integer_value(x) result(is_integer)
        real(dp), intent(in) :: x !! finite real value
        logical :: is_integer

        is_integer = ieee_is_finite(x) .and. abs(x) < 2.0_dp ** 52 .and. x == anint(x)
    end function is_integer_value

    pure function h0_integral_series(x, divide_by_t) result(value)
        real(dp), intent(in) :: x !! nonnegative argument, no larger than integration anchor
        logical, intent(in) :: divide_by_t !! integrate H_0(t)/t when true
        real(dp) :: value

        real(dp) :: add
        real(dp) :: sum_value
        real(dp) :: term
        integer :: k

        if (x == 0.0_dp) then
            value = 0.0_dp
            return
        end if

        term = (0.5_dp * x) * exp(-2.0_dp * gammaln(1.5_dp))
        k = 0
        if (divide_by_t) then
            sum_value = term
        else
            sum_value = 0.5_dp * x * term
        end if

        do while (k < max_series_terms)
            term = -term * 0.25_dp * x * x / (real(k, dp) + 1.5_dp) ** 2
            k = k + 1
            if (divide_by_t) then
                add = term / real(2 * k + 1, dp)
            else
                add = x * term / real(2 * k + 2, dp)
            end if
            sum_value = sum_value + add
            if (k >= 8 .and. abs(add) <= 4.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(sum_value))) exit
        end do
        value = sum_value
    end function h0_integral_series

    pure function l0_integral_series(x) result(value)
        real(dp), intent(in) :: x !! nonnegative finite argument
        real(dp) :: value

        real(dp) :: add
        real(dp) :: sum_value
        real(dp) :: term
        integer :: k

        if (x == 0.0_dp) then
            value = 0.0_dp
            return
        end if

        term = (0.5_dp * x) * exp(-2.0_dp * gammaln(1.5_dp))
        k = 0
        sum_value = 0.5_dp * x * term
        do while (k < max_series_terms)
            term = term * 0.25_dp * x * x / (real(k, dp) + 1.5_dp) ** 2
            k = k + 1
            add = x * term / real(2 * k + 2, dp)
            if (.not. ieee_is_finite(add) .or. .not. ieee_is_finite(sum_value)) then
                value = positive_infinity(x)
                return
            end if
            sum_value = sum_value + add
            if (k >= 8 .and. abs(add) <= 4.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(sum_value))) exit
        end do
        value = sum_value
    end function l0_integral_series

    pure function composite_gauss_h0(a, b, divide_by_t) result(value)
        real(dp), intent(in) :: a !! lower positive integration limit
        real(dp), intent(in) :: b !! upper integration limit, >= a
        logical, intent(in) :: divide_by_t !! integrate H_0(t)/t when true
        real(dp) :: value

        real(dp) :: half_width
        real(dp) :: left
        real(dp) :: midpoint
        real(dp) :: panel_width
        real(dp) :: right
        real(dp) :: x1
        real(dp) :: x2
        real(dp) :: f1
        real(dp) :: f2
        integer :: i
        integer :: p
        integer :: panels

        if (b <= a) then
            value = 0.0_dp
            return
        end if

        panels = max(1, ceiling((b - a) / integration_panel_width))
        panel_width = (b - a) / real(panels, dp)
        value = 0.0_dp
        do p = 0, panels - 1
            left = a + real(p, dp) * panel_width
            right = left + panel_width
            midpoint = 0.5_dp * (left + right)
            half_width = 0.5_dp * panel_width
            do i = 1, 4
                x1 = midpoint - half_width * gl_nodes(i)
                x2 = midpoint + half_width * gl_nodes(i)
                f1 = struve_core(0.0_dp, x1, .false.)
                f2 = struve_core(0.0_dp, x2, .false.)
                if (divide_by_t) then
                    f1 = f1 / x1
                    f2 = f2 / x2
                end if
                value = value + half_width * gl_weights(i) * (f1 + f2)
            end do
        end do
    end function composite_gauss_h0

    pure function signed_exp_from_log(log_abs, sign_value, kind_arg) result(value)
        real(dp), intent(in) :: log_abs !! logarithm of absolute magnitude
        real(dp), intent(in) :: sign_value !! sign of result
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE values
        real(dp) :: value

        real(dp), parameter :: log_huge_dp = log(huge(1.0_dp))
        real(dp), parameter :: log_tiny_dp = log(tiny(1.0_dp))

        if (sign_value == 0.0_dp) then
            value = 0.0_dp
        else if (log_abs > log_huge_dp) then
            value = sign(positive_infinity(kind_arg), sign_value)
        else if (log_abs < log_tiny_dp - 40.0_dp) then
            value = sign(0.0_dp, sign_value)
        else
            value = sign(exp(log_abs), sign_value)
        end if
    end function signed_exp_from_log

end module scifort_struve
