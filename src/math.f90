! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_math
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_negative_inf, ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    implicit none
    private

    public :: expm1_safe
    public :: log1p_safe
    public :: log1pexp
    public :: negative_infinity
    public :: positive_infinity
    public :: quiet_nan
    public :: valid_loc_scale

contains

    pure elemental function quiet_nan(x) result(y)
        real(dp), intent(in) :: x !! argument whose real kind selects the kind of the result
        real(dp) :: y

        y = ieee_value(x, ieee_quiet_nan)
    end function quiet_nan

    pure elemental function positive_infinity(x) result(y)
        real(dp), intent(in) :: x !! argument whose real kind selects the kind of the result
        real(dp) :: y

        y = ieee_value(x, ieee_positive_inf)
    end function positive_infinity

    pure elemental function negative_infinity(x) result(y)
        real(dp), intent(in) :: x !! argument whose real kind selects the kind of the result
        real(dp) :: y

        y = ieee_value(x, ieee_negative_inf)
    end function negative_infinity

    pure elemental logical function valid_loc_scale(loc, scale) result(valid)
        real(dp), intent(in) :: loc !! location to check
        real(dp), intent(in) :: scale !! scale to check

        valid = ieee_is_finite(loc) .and. ieee_is_finite(scale) .and. scale > 0.0_dp
    end function valid_loc_scale

    ! log(1 + x). For -0.5 <= x <= 1 the argument s = x / (2 + x) satisfies
    ! |s| <= 1/3 and log(1 + x) = 2 * atanh(s) is summed from DLMF 4.6.4.
    ! Outside that interval 1 + x is either exact or large enough that its
    ! rounding error is not amplified. The formulation does not rely on the
    ! exactness of (1 + x) - 1, so value-unsafe compiler optimizations cannot
    ! remove the correction.
    pure elemental function log1p_safe(x) result(y)
        real(dp), intent(in) :: x !! argument of log(1 + x), >= -1
        real(dp) :: y

        integer :: k
        real(dp) :: s
        real(dp) :: s2
        real(dp) :: term

        if (ieee_is_nan(x)) then
            y = x
        else if (x < -1.0_dp) then
            y = quiet_nan(x)
        else if (x <= -1.0_dp) then
            y = negative_infinity(x)
        else if (x < -0.5_dp .or. x > 1.0_dp) then
            y = log(1.0_dp + x)
        else
            s = x / (2.0_dp + x)
            s2 = s * s
            term = s
            y = s
            do k = 1, 40
                term = term * s2
                if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * abs(y)) exit
                y = y + term / real(2 * k + 1, dp)
            end do
            y = 2.0_dp * y
        end if
    end function log1p_safe

    ! exp(x) - 1. For |x| < 0.5 the Maclaurin series of DLMF 4.2.19 without
    ! its leading term is summed in nested form. Elsewhere the subtraction
    ! loses at most a small constant factor of relative accuracy.
    pure elemental function expm1_safe(x) result(y)
        real(dp), intent(in) :: x !! exponent in exp(x) - 1
        real(dp) :: y

        integer, parameter :: n_terms = 20
        integer :: k

        if (ieee_is_nan(x)) then
            y = x
        else if (abs(x) >= 0.5_dp) then
            y = exp(x) - 1.0_dp
        else
            y = 1.0_dp
            do k = n_terms, 2, -1
                y = 1.0_dp + y * x / real(k, dp)
            end do
            y = x * y
        end if
    end function expm1_safe

    pure elemental function log1pexp(x) result(y)
        real(dp), intent(in) :: x !! exponent in log(1 + exp(x))
        real(dp) :: y

        if (x > 0.0_dp) then
            y = x + log1p_safe(exp(-x))
        else
            y = log1p_safe(exp(x))
        end if
    end function log1pexp

end module scifort_math
