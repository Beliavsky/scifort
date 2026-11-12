! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real ordered cylindrical and modified Bessel compatibility wrappers.
! Integer-order Y uses the Fortran standard intrinsic. Modified I/K reuse
! SciFort's existing scaled/logarithmic kernels and standard order identities.

module scifort_bessel_order
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_i, only : log_besseli_scaled
    use scifort_bessel_k, only : log_besselk, log_besselk_scaled
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: iv
    public :: ive
    public :: kn
    public :: kv
    public :: kve
    public :: yn

contains

    pure elemental function yn(n, x) result(y)
        integer, intent(in) :: n !! integer order of the real Bessel Y function
        real(dp), intent(in) :: x !! nonnegative real argument
        real(dp) :: y

        integer :: n_abs
        real(dp) :: sign_order

        intrinsic :: bessel_yn

        if (ieee_is_nan(x)) then
            y = x
            return
        end if
        if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
            return
        end if
        if (n < -huge(n)) then
            y = quiet_nan(x)
            return
        end if

        n_abs = abs(n)
        sign_order = 1.0_dp
        if (n < 0 .and. mod(n_abs, 2) /= 0) sign_order = -1.0_dp

        if (x == 0.0_dp) then
            if (sign_order > 0.0_dp) then
                y = negative_infinity(x)
            else
                y = positive_infinity(x)
            end if
        else
            y = sign_order * bessel_yn(n_abs, x)
        end if
    end function yn

    pure elemental function kn(n, x) result(y)
        integer, intent(in) :: n !! integer order of the modified Bessel K function
        real(dp), intent(in) :: x !! nonnegative real argument
        real(dp) :: y
        real(dp) :: order

        if (ieee_is_nan(x)) then
            y = x
        else if (x < 0.0_dp) then
            y = quiet_nan(x)
        else if (x == 0.0_dp) then
            y = positive_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            order = abs(real(n, dp))
            y = exp_from_log(log_besselk(order, x), x)
        end if
    end function kn

    pure elemental function kv(v, x) result(y)
        real(dp), intent(in) :: v !! real order of the modified Bessel K function
        real(dp), intent(in) :: x !! nonnegative real argument
        real(dp) :: y

        if (ieee_is_nan(v) .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(v) .or. x < 0.0_dp) then
            y = quiet_nan(x)
        else if (x == 0.0_dp) then
            y = positive_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            y = exp_from_log(log_besselk(v, x), x)
        end if
    end function kv

    pure elemental function kve(v, x) result(y)
        real(dp), intent(in) :: v !! real order of exp(x)*K_v(x)
        real(dp), intent(in) :: x !! nonnegative finite real argument
        real(dp) :: y

        if (ieee_is_nan(v) .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(v) .or. .not. ieee_is_finite(x) .or. x < 0.0_dp) then
            y = quiet_nan(x)
        else if (x == 0.0_dp) then
            y = positive_infinity(x)
        else
            y = exp_from_log(log_besselk_scaled(v, x), x)
        end if
    end function kve

    pure elemental function iv(v, x) result(y)
        real(dp), intent(in) :: v !! real order of the modified Bessel I function
        real(dp), intent(in) :: x !! real argument; negative values require integer order
        real(dp) :: y

        y = besseli_ordered(v, x, .false.)
    end function iv

    pure elemental function ive(v, x) result(y)
        real(dp), intent(in) :: v !! real order of exp(-abs(x))*I_v(x)
        real(dp), intent(in) :: x !! real argument; negative values require integer order
        real(dp) :: y

        y = besseli_ordered(v, x, .true.)
    end function ive

    pure elemental function besseli_ordered(v, x, scaled) result(y)
        real(dp), intent(in) :: v !! real modified-Bessel-I order
        real(dp), intent(in) :: x !! real modified-Bessel-I argument
        logical, intent(in) :: scaled !! true for exp(-abs(x)) scaling
        real(dp) :: y

        logical :: integer_v
        real(dp) :: ax
        real(dp) :: coefficient
        real(dp) :: log_i
        real(dp) :: log_k_term
        real(dp) :: nu
        real(dp) :: sign_x

        if (ieee_is_nan(v) .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(v) .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
            return
        end if

        integer_v = is_integer_value(v)
        if (x < 0.0_dp .and. .not. integer_v) then
            y = quiet_nan(x)
            return
        end if

        if (x == 0.0_dp) then
            if (v == 0.0_dp) then
                y = 1.0_dp
            else if (v > 0.0_dp .or. integer_v) then
                y = 0.0_dp
            else if (scaled) then
                y = quiet_nan(x)
            else
                y = positive_infinity(x)
            end if
            return
        end if

        ax = abs(x)
        sign_x = 1.0_dp
        if (x < 0.0_dp .and. odd_integer_value(v)) sign_x = -1.0_dp

        if (v >= 0.0_dp .or. integer_v) then
            nu = abs(v)
            log_i = log_besseli_scaled(nu, ax)
            if (.not. scaled) log_i = log_i + ax
            y = sign_x * exp_from_log(log_i, x)
            return
        end if

        ! DLMF 10.27.3 rearranged:
        ! I_{-nu}(x) = I_nu(x) + (2/pi) sin(pi*nu) K_nu(x).
        nu = -v
        log_i = log_besseli_scaled(nu, ax)
        if (.not. scaled) log_i = log_i + ax
        coefficient = (2.0_dp / scifort_pi) * sinpi_nonnegative(nu)
        if (coefficient == 0.0_dp) then
            y = exp_from_log(log_i, x)
            return
        end if

        log_k_term = log(abs(coefficient)) + log_besselk(nu, ax)
        if (scaled) log_k_term = log_k_term - ax
        y = signed_log_sum(log_i, log_k_term, coefficient > 0.0_dp, x)
    end function besseli_ordered

    pure elemental logical function is_integer_value(v) result(is_integer)
        real(dp), intent(in) :: v !! finite real order to classify

        is_integer = v == anint(v)
    end function is_integer_value

    pure elemental logical function odd_integer_value(v) result(is_odd)
        real(dp), intent(in) :: v !! finite integer-valued real order

        is_odd = modulo(abs(v), 2.0_dp) == 1.0_dp
    end function odd_integer_value

    pure elemental function sinpi_nonnegative(v) result(y)
        real(dp), intent(in) :: v !! finite nonnegative argument of sin(pi*v)
        real(dp) :: y
        real(dp) :: r

        r = modulo(v, 2.0_dp)
        if (r <= 0.5_dp) then
            y = sin(scifort_pi * r)
        else if (r <= 1.5_dp) then
            y = sin(scifort_pi * (1.0_dp - r))
        else
            y = sin(scifort_pi * (r - 2.0_dp))
        end if
    end function sinpi_nonnegative

    pure elemental function signed_log_sum(log_a, log_b, add_b, kind_arg) result(y)
        real(dp), intent(in) :: log_a !! logarithm of the first positive term
        real(dp), intent(in) :: log_b !! logarithm of the second positive magnitude
        logical, intent(in) :: add_b !! true to add the second term, false to subtract it
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y
        real(dp) :: log_abs

        if (add_b) then
            if (log_a >= log_b) then
                log_abs = log_a + log1p_safe(exp(log_b - log_a))
            else
                log_abs = log_b + log1p_safe(exp(log_a - log_b))
            end if
            y = exp_from_log(log_abs, kind_arg)
            return
        end if

        if (log_a == log_b) then
            y = 0.0_dp
        else if (log_a > log_b) then
            log_abs = log_a + log1p_safe(-exp(log_b - log_a))
            y = exp_from_log(log_abs, kind_arg)
        else
            log_abs = log_b + log1p_safe(-exp(log_a - log_b))
            y = -exp_from_log(log_abs, kind_arg)
        end if
    end function signed_log_sum

    pure elemental function exp_from_log(log_value, kind_arg) result(y)
        real(dp), intent(in) :: log_value !! logarithm of a nonnegative magnitude
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y

        if (ieee_is_nan(log_value)) then
            y = log_value
        else if (log_value > log(huge(1.0_dp))) then
            y = positive_infinity(kind_arg)
        else
            y = exp(log_value)
        end if
    end function exp_from_log

end module scifort_bessel_order
