! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_math
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_negative_inf, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
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
        real(dp), intent(in) :: x
        real(dp) :: y

        y = ieee_value(x, ieee_quiet_nan)
    end function quiet_nan

    pure elemental function positive_infinity(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        y = ieee_value(x, ieee_positive_inf)
    end function positive_infinity

    pure elemental function negative_infinity(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        y = ieee_value(x, ieee_negative_inf)
    end function negative_infinity

    pure elemental logical function valid_loc_scale(loc, scale) result(valid)
        real(dp), intent(in) :: loc
        real(dp), intent(in) :: scale

        valid = ieee_is_finite(loc) .and. ieee_is_finite(scale) .and. scale > 0.0_dp
    end function valid_loc_scale

    pure elemental function log1p_safe(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        integer :: k
        real(dp) :: term

        if (x < -1.0_dp) then
            y = quiet_nan(x)
        else if (x <= -1.0_dp) then
            y = negative_infinity(x)
        else if (abs(x) > 1.0e-4_dp) then
            y = log(1.0_dp + x)
        else
            y = 0.0_dp
            term = x
            do k = 1, 12
                y = y + term / real(k, dp)
                term = -term * x
            end do
        end if
    end function log1p_safe

    pure elemental function expm1_safe(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        integer :: k
        real(dp) :: term

        if (abs(x) > 1.0e-5_dp) then
            y = exp(x) - 1.0_dp
        else
            y = 0.0_dp
            term = 1.0_dp
            do k = 1, 12
                term = term * x / real(k, dp)
                y = y + term
            end do
        end if
    end function expm1_safe

    pure elemental function log1pexp(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        if (x > 0.0_dp) then
            y = x + log1p_safe(exp(-x))
        else
            y = log1p_safe(exp(x))
        end if
    end function log1pexp

end module scifort_math
