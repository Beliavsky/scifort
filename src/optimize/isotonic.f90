! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Weighted monotone least squares by the pool-adjacent-violators algorithm.
! The implementation is independent of the SciPy and Busing (2022) sources.
!
module scifort_optimize_isotonic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: isotonic_status_success = 0
    integer, parameter, public :: isotonic_status_invalid_input = 1
    integer, parameter, public :: isotonic_status_numerical_failure = 2

    type, public :: isotonic_result
        real(dp), allocatable :: x(:) !! monotone fitted values, same size as input
        real(dp), allocatable :: weights(:) !! total original case weight in each fitted block
        integer, allocatable :: blocks(:) !! one-based start indices, ending at size(x)+1
        integer :: status = isotonic_status_invalid_input !! zero on success
    end type isotonic_result

    public :: isotonic_regression

contains

    pure function isotonic_regression(y, weights, increasing) result(res)
        real(dp), intent(in) :: y(:) !! finite observations in their given order
        real(dp), intent(in), optional :: weights(:) !! finite strictly positive case weights
        logical, intent(in), optional :: increasing !! nondecreasing if true (default true)
        type(isotonic_result) :: res

        integer :: i
        integer :: k
        integer :: n
        integer, allocatable :: first(:)
        logical :: ascending
        real(dp) :: total
        real(dp) :: w
        real(dp), allocatable :: means(:)
        real(dp), allocatable :: pool_weights(:)

        allocate(res%x(0), res%weights(0), res%blocks(0))
        n = size(y)
        if (present(weights)) then
            if (size(weights) /= n) return
            if (any(.not. ieee_is_finite(weights))) return
            if (any(weights <= 0.0_dp)) return
        end if
        if (any(.not. ieee_is_finite(y))) return

        ascending = .true.
        if (present(increasing)) ascending = increasing
        allocate(means(n), pool_weights(n), first(n))
        k = 0
        do i = 1, n
            w = 1.0_dp
            if (present(weights)) w = weights(i)
            k = k + 1
            means(k) = y(i)
            pool_weights(k) = w
            first(k) = i
            do while (k > 1)
                if (ascending) then
                    if (means(k - 1) < means(k)) exit
                else
                    if (means(k - 1) > means(k)) exit
                end if
                if (pool_weights(k - 1) > huge(total) - pool_weights(k)) then
                    res%status = isotonic_status_numerical_failure
                    return
                end if
                total = pool_weights(k - 1) + pool_weights(k)
                ! A convex combination avoids overflow in sums such as
                ! w1*y1 + w2*y2 when observations have opposite signs.
                means(k - 1) = (pool_weights(k - 1) / total) * means(k - 1) + &
                    (pool_weights(k) / total) * means(k)
                pool_weights(k - 1) = total
                k = k - 1
            end do
        end do

        deallocate(res%x, res%weights, res%blocks)
        allocate(res%x(n), res%weights(k), res%blocks(k + 1))
        res%weights = pool_weights(1:k)
        res%blocks(1:k) = first(1:k)
        res%blocks(k + 1) = n + 1
        do i = 1, k
            res%x(first(i):res%blocks(i + 1) - 1) = means(i)
        end do
        res%status = isotonic_status_success
    end function isotonic_regression

end module scifort_optimize_isotonic
