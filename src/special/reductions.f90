! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special_reductions
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan
    implicit none
    private

    public :: log_softmax
    public :: logsumexp
    public :: softmax

contains

    pure function logsumexp(x) result(value)
        real(dp), intent(in) :: x(:) !! one-dimensional values to reduce
        real(dp) :: value

        real(dp) :: maximum_value
        real(dp) :: total

        if (size(x) == 0) then
            value = negative_infinity(0.0_dp)
            return
        end if
        if (any(ieee_is_nan(x))) then
            value = quiet_nan(0.0_dp)
            return
        end if

        maximum_value = maxval(x)
        if (.not. ieee_is_finite(maximum_value)) then
            value = maximum_value
            return
        end if

        total = sum(exp(x - maximum_value))
        value = maximum_value + log(total)
    end function logsumexp

    pure function softmax(x) result(values)
        real(dp), intent(in) :: x(:) !! one-dimensional values to normalize
        real(dp), allocatable :: values(:)

        real(dp) :: maximum_value
        real(dp) :: total

        allocate(values(size(x)))
        if (size(x) == 0) return
        if (any(ieee_is_nan(x))) then
            values = quiet_nan(0.0_dp)
            return
        end if

        maximum_value = maxval(x)
        if (.not. ieee_is_finite(maximum_value)) then
            values = quiet_nan(0.0_dp)
            return
        end if

        values = exp(x - maximum_value)
        total = sum(values)
        values = values / total
    end function softmax

    pure function log_softmax(x) result(values)
        real(dp), intent(in) :: x(:) !! one-dimensional values to log-normalize
        real(dp), allocatable :: values(:)

        integer :: i
        real(dp) :: maximum_value
        real(dp) :: normalizer

        allocate(values(size(x)))
        if (size(x) == 0) return
        if (any(ieee_is_nan(x))) then
            values = quiet_nan(0.0_dp)
            return
        end if

        maximum_value = maxval(x)
        if (.not. ieee_is_finite(maximum_value)) then
            if (maximum_value > 0.0_dp) then
                do i = 1, size(x)
                    if (x(i) > 0.0_dp .and. .not. ieee_is_finite(x(i))) then
                        values(i) = quiet_nan(0.0_dp)
                    else
                        values(i) = negative_infinity(0.0_dp)
                    end if
                end do
            else
                values = quiet_nan(0.0_dp)
            end if
            return
        end if

        values = x - maximum_value
        normalizer = log(sum(exp(values)))
        values = values - normalizer
    end function log_softmax

end module scifort_special_reductions
