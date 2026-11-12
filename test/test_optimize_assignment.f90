! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference objectives and index conventions checked against SciPy 1.17.0.
program test_optimize_assignment
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_assignment, only : linear_sum_assignment, &
        assignment_status_success, assignment_status_invalid_input, assignment_status_infeasible
    implicit none

    integer :: i
    integer :: j
    integer :: nr
    integer :: nc
    integer :: status
    integer, allocatable :: rows(:)
    integer, allocatable :: cols(:)
    real(dp) :: positive_inf
    real(dp) :: negative_inf
    real(dp) :: nan_value
    real(dp) :: cost(3, 3)
    real(dp) :: small(2, 3)
    real(dp) :: tall(4, 2)
    real(dp), allocatable :: a(:, :)

    positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
    negative_inf = ieee_value(0.0_dp, ieee_negative_inf)
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)

    cost(1, :) = [4.0_dp, 1.0_dp, 3.0_dp]
    cost(2, :) = [2.0_dp, 0.0_dp, 5.0_dp]
    cost(3, :) = [3.0_dp, 2.0_dp, 2.0_dp]
    call check_assignment(cost, .false., 5.0_dp)

    small(1, :) = [8.0_dp, 2.0_dp, 4.0_dp]
    small(2, :) = [7.0_dp, 3.0_dp, 6.0_dp]
    call check_assignment(small, .false., 7.0_dp)
    call check_assignment(small, .true., 14.0_dp)

    tall(1, :) = [7.0_dp, 10.0_dp]
    tall(2, :) = [5.0_dp, 9.0_dp]
    tall(3, :) = [9.0_dp, 1.0_dp]
    tall(4, :) = [2.0_dp, 8.0_dp]
    call check_assignment(tall, .false., 3.0_dp)
    call check_assignment(tall, .true., 19.0_dp)

    do i = 1, 6
        select case (i)
        case (1)
            nr = 3
            nc = 3
        case (2)
            nr = 5
            nc = 7
        case (3)
            nr = 7
            nc = 5
        case (4)
            nr = 8
            nc = 8
        case (5)
            nr = 12
            nc = 16
        case default
            nr = 16
            nc = 12
        end select
        allocate(a(nr, nc))
        do j = 1, nc
            ! Deterministic integer-valued matrices with both cost signs.
            block
                integer :: k
                do k = 1, nr
                    a(k, j) = real(mod(37 * k + 91 * j + 13 * k * j * (k + 1), 257) &
                        - 128, dp)
                end do
            end block
        end do
        select case (i)
        case (1)
            call check_assignment(a, .false., -202.0_dp)
            call check_assignment(a, .true., 81.0_dp)
        case (2)
            call check_assignment(a, .false., -495.0_dp)
            call check_assignment(a, .true., 414.0_dp)
        case (3)
            call check_assignment(a, .false., -479.0_dp)
            call check_assignment(a, .true., 477.0_dp)
        case (4)
            call check_assignment(a, .false., -756.0_dp)
            call check_assignment(a, .true., 504.0_dp)
        case (5)
            call check_assignment(a, .false., -1324.0_dp)
            call check_assignment(a, .true., 1163.0_dp)
        case default
            call check_assignment(a, .false., -1266.0_dp)
            call check_assignment(a, .true., 1308.0_dp)
        end select
        deallocate(a)
    end do

    small(1, :) = [positive_inf, 1.0_dp, 2.0_dp]
    small(2, :) = [3.0_dp, positive_inf, 4.0_dp]
    call check_assignment(small, .false., 4.0_dp)
    small(1, :) = [negative_inf, 5.0_dp, 1.0_dp]
    small(2, :) = [6.0_dp, negative_inf, 4.0_dp]
    call check_assignment(small, .true., 11.0_dp)

    call linear_sum_assignment(reshape([positive_inf, positive_inf, 1.0_dp, 2.0_dp], &
        [2, 2]), rows, cols, status)
    call require(status == assignment_status_infeasible, 'infeasible matrix status')
    call require(size(rows) == 0, 'infeasible outputs cleared')

    cost(1, 1) = nan_value
    call linear_sum_assignment(cost, rows, cols, status)
    call require(status == assignment_status_invalid_input, 'NaN is invalid')
    cost(1, 1) = negative_inf
    call linear_sum_assignment(cost, rows, cols, status)
    call require(status == assignment_status_invalid_input, 'negative infinity for min')
    cost(1, 1) = positive_inf
    call linear_sum_assignment(cost, rows, cols, status, maximize=.true.)
    call require(status == assignment_status_invalid_input, 'positive infinity for max')

    allocate(a(0, 4))
    call linear_sum_assignment(a, rows, cols, status)
    call require(status == assignment_status_success .and. size(rows) == 0, 'empty rows')
    deallocate(a)
    allocate(a(4, 0))
    call linear_sum_assignment(a, rows, cols, status)
    call require(status == assignment_status_success .and. size(cols) == 0, 'empty columns')
    deallocate(a)

    small(1, :) = [1.0e300_dp, 3.0e300_dp, 2.0e300_dp]
    small(2, :) = [4.0e300_dp, 1.0e300_dp, 3.0e300_dp]
    call linear_sum_assignment(small, rows, cols, status)
    call require(status == assignment_status_success .and. all(cols == [1, 2]), &
        'extreme finite values')

    print '(a)', 'test_optimize_assignment: PASS'

contains

    subroutine check_assignment(a, maximum, objective)
        real(dp), intent(in) :: a(:, :) !! reference cost matrix
        logical, intent(in) :: maximum !! compare maximizing rather than minimizing
        real(dp), intent(in) :: objective !! SciPy 1.17 optimum
        integer, allocatable :: r(:)
        integer, allocatable :: c(:)
        integer :: k
        integer :: ierr
        real(dp) :: actual

        call linear_sum_assignment(a, r, c, ierr, maximize=maximum)
        call require(ierr == assignment_status_success, 'assignment status')
        call require(size(r) == min(size(a, 1), size(a, 2)), 'assignment length')
        call require(all(r >= 1) .and. all(r <= size(a, 1)), 'row bounds')
        call require(all(c >= 1) .and. all(c <= size(a, 2)), 'column bounds')
        do k = 2, size(r)
            call require(r(k) > r(k-1), 'row indices sorted')
        end do
        do k = 1, size(c)
            call require(count(c == c(k)) == 1, 'one-to-one assignment')
        end do
        actual = 0.0_dp
        do k = 1, size(r)
            actual = actual + a(r(k), c(k))
        end do
        if (abs(actual - objective) > 1.0e-11_dp * max(1.0_dp, abs(objective))) then
            print '(a,2es24.15)', 'assignment objective mismatch: ', actual, objective
            error stop 1
        end if
    end subroutine check_assignment

    subroutine require(condition, message)
        logical, intent(in) :: condition !! assertion condition
        character(*), intent(in) :: message !! description on failure
        if (.not. condition) then
            print '(a,a)', 'test_optimize_assignment: ', message
            error stop 1
        end if
    end subroutine require

end program test_optimize_assignment
