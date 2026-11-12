! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Fast real scalar Bessel functions matching the common scipy.special names.
! J0/J1/Y0/Y1 wrap the Fortran standard intrinsics. Modified-Bessel I/K
! functions reuse SciFort's existing scaled-I and log-K kernels so the public
! unscaled forms avoid unnecessary intermediate overflow or underflow.
module scifort_bessel_classic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_k, only : log_besselk, log_besselk_scaled
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_special_elementary, only : i0e, i1e
    implicit none
    private

    public :: i0
    public :: i1
    public :: j0
    public :: j1
    public :: k0
    public :: k0e
    public :: k1
    public :: k1e
    public :: y0
    public :: y1

contains

    pure elemental function j0(x) result(y)
        real(dp), intent(in) :: x !! real argument of the order-zero Bessel J function
        real(dp) :: y

        intrinsic :: bessel_j0

        if (ieee_is_nan(x)) then
            y = x
        else if (.not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else
            y = bessel_j0(x)
        end if
    end function j0

    pure elemental function j1(x) result(y)
        real(dp), intent(in) :: x !! real argument of the order-one Bessel J function
        real(dp) :: y

        intrinsic :: bessel_j1

        if (ieee_is_nan(x)) then
            y = x
        else if (.not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else
            y = bessel_j1(x)
        end if
    end function j1

    pure elemental function y0(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real argument of the order-zero Bessel Y function
        real(dp) :: y

        intrinsic :: bessel_y0

        if (ieee_is_nan(x)) then
            y = x
        else if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else if (x == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = bessel_y0(x)
        end if
    end function y0

    pure elemental function y1(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real argument of the order-one Bessel Y function
        real(dp) :: y

        intrinsic :: bessel_y1

        if (ieee_is_nan(x)) then
            y = x
        else if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else if (x == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = bessel_y1(x)
        end if
    end function y1

    pure elemental function i0(x) result(y)
        real(dp), intent(in) :: x !! finite real argument for I0
        real(dp) :: y
        real(dp) :: scaled

        if (ieee_is_nan(x)) then
            y = x
        else if (.not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else
            scaled = i0e(x)
            y = restore_scaled_i(abs(scaled), abs(x), x)
        end if
    end function i0

    pure elemental function i1(x) result(y)
        real(dp), intent(in) :: x !! finite real argument for I1
        real(dp) :: y
        real(dp) :: magnitude, scaled

        if (ieee_is_nan(x)) then
            y = x
        else if (.not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else
            scaled = i1e(x)
            if (scaled == 0.0_dp) then
                y = 0.0_dp
            else
                magnitude = restore_scaled_i(abs(scaled), abs(x), x)
                if (scaled < 0.0_dp) then
                    y = -magnitude
                else
                    y = magnitude
                end if
            end if
        end if
    end function i1

    pure elemental function k0(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real argument for K0
        real(dp) :: y

        y = positive_exp_from_log(log_besselk(0.0_dp, x), x)
    end function k0

    pure elemental function k0e(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real argument; returns exp(x)*K0(x)
        real(dp) :: y

        y = positive_exp_from_log(log_besselk_scaled(0.0_dp, x), x)
    end function k0e

    pure elemental function k1(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real argument for K1
        real(dp) :: y

        y = positive_exp_from_log(log_besselk(1.0_dp, x), x)
    end function k1

    pure elemental function k1e(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real argument; returns exp(x)*K1(x)
        real(dp) :: y

        y = positive_exp_from_log(log_besselk_scaled(1.0_dp, x), x)
    end function k1e

    pure elemental function restore_scaled_i(scaled, ax, kind_arg) result(y)
        real(dp), intent(in) :: scaled !! nonnegative exp(-ax)-scaled modified-Bessel value
        real(dp), intent(in) :: ax !! absolute Bessel argument, finite and nonnegative
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y
        real(dp) :: log_value

        if (scaled == 0.0_dp) then
            y = 0.0_dp
            return
        end if
        log_value = ax + log(scaled)
        y = positive_exp_from_log(log_value, kind_arg)
    end function restore_scaled_i

    pure elemental function positive_exp_from_log(log_value, kind_arg) result(y)
        real(dp), intent(in) :: log_value !! logarithm of a nonnegative result
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y

        if (ieee_is_nan(log_value)) then
            y = log_value
        else if (log_value > log(huge(1.0_dp))) then
            y = positive_infinity(kind_arg)
        else
            y = exp(log_value)
        end if
    end function positive_exp_from_log

end module scifort_bessel_classic
