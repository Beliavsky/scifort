! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
module scifort_logrank_bridge
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_empirical_distribution, only : censored_data, logrank, logrank_result
    implicit none
    private
    public :: survival_logrank
contains
    subroutine survival_logrank(x, x_right, y, y_right, alternative, statistic, pvalue, status)
        !! Compare ordinary or right-censored groups without exposing native descriptors.
        !! @dimension nx = size(x)
        !! @dimension rx = size(x_right)
        !! @dimension ny = size(y)
        !! @dimension ry = size(y_right)
        !! @require nx >= 0
        !! @require rx >= 0
        !! @require ny >= 0
        !! @require ry >= 0
        !! All four lengths are independent; their sum must fit the native integer kind.
        !! Empty groups and zero variance return NaNs with status 0.
        real(dp), intent(in) :: x(:) !! shape (nx); finite event times in group X
        real(dp), intent(in) :: x_right(:) !! shape (rx); finite right-censoring times in group X
        real(dp), intent(in) :: y(:) !! shape (ny); finite event times in group Y
        real(dp), intent(in) :: y_right(:) !! shape (ry); finite right-censoring times in group Y
        integer, intent(in) :: alternative !! 0 two-sided; 1 less; 2 greater
        real(dp), intent(out) :: statistic !! signed normal statistic; NaN when undefined or invalid
        real(dp), intent(out) :: pvalue !! tail probability; NaN when undefined or invalid
        integer, intent(out) :: status !! 0 completed, possibly undefined; 1 count overflow; 2 invalid input
        type(censored_data) :: a, b
        type(logrank_result) :: result
        character(len=9), parameter :: alternatives(0:2) = ['two-sided','less     ','greater  ']
        integer(int64) :: total

        statistic = quiet_nan(0.0_dp)
        pvalue = statistic
        total = size(x,kind=int64) + size(x_right,kind=int64) + &
            size(y,kind=int64) + size(y_right,kind=int64)
        status = 1
        if (total > int(huge(0),int64)) return
        status = 2
        if (alternative < 0 .or. alternative > 2) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(x_right))) return
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(y_right))) return
        a%uncensored = x
        a%right = x_right
        b%uncensored = y
        b%right = y_right
        result = logrank(a, b, alternative=alternatives(alternative))
        statistic = result%statistic
        pvalue = result%pvalue
        status = 0
    end subroutine survival_logrank
end module scifort_logrank_bridge
