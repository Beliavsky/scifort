! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rectangular linear assignment by a shortest-augmenting-path primal-dual method.
! This is an independent implementation of the Hungarian assignment algorithm,
! not a translation of SciPy's Jonker-Volgenant implementation.
!
module scifort_optimize_assignment
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: assignment_status_success = 0
    integer, parameter, public :: assignment_status_invalid_input = 1
    integer, parameter, public :: assignment_status_infeasible = 2

    public :: linear_sum_assignment

contains

    pure subroutine linear_sum_assignment(cost_matrix, row_ind, col_ind, status, maximize)
        real(dp), intent(in) :: cost_matrix(:, :) !! finite costs; +Inf forbids a minimum-cost edge
        integer, allocatable, intent(out) :: row_ind(:) !! sorted one-based row indices of pairs
        integer, allocatable, intent(out) :: col_ind(:) !! corresponding one-based column indices
        integer, intent(out) :: status !! zero if feasible, otherwise invalid or infeasible
        logical, intent(in), optional :: maximize !! maximize weights rather than minimize (default false)

        integer :: i
        integer :: i0
        integer :: j
        integer :: j0
        integer :: j1
        integer :: k
        integer :: m
        integer :: n
        integer :: ncol
        integer :: nrow
        integer, allocatable :: match(:)
        integer, allocatable :: predecessor(:)
        integer, allocatable :: partner(:)
        logical :: reversed
        logical :: want_maximum
        logical, allocatable :: visited(:)
        real(dp) :: current
        real(dp) :: delta
        real(dp) :: infinity
        real(dp) :: scale
        real(dp) :: value
        real(dp), allocatable :: min_reduced(:)
        real(dp), allocatable :: potential_col(:)
        real(dp), allocatable :: potential_row(:)

        nrow = size(cost_matrix, 1)
        ncol = size(cost_matrix, 2)
        n = min(nrow, ncol)
        m = max(nrow, ncol)
        reversed = nrow > ncol
        want_maximum = .false.
        if (present(maximize)) want_maximum = maximize
        status = assignment_status_invalid_input
        allocate(row_ind(0), col_ind(0))

        ! A plus infinity is a forbidden edge for minimization; a minus
        ! infinity is forbidden for maximization. The opposite sign and
        ! NaNs are rejected, matching scipy.optimize semantics.
        scale = 1.0_dp
        do j = 1, ncol
            do i = 1, nrow
                value = cost_matrix(i, j)
                if (ieee_is_nan(value)) return
                if (.not. ieee_is_finite(value)) then
                    if (want_maximum .and. value > 0.0_dp) return
                    if (.not. want_maximum .and. value < 0.0_dp) return
                else
                    scale = max(scale, abs(value))
                end if
            end do
        end do
        status = assignment_status_success
        if (n == 0) return

        ! The shorter matrix dimension represents the rows of the working
        ! assignment problem. Keep only O(n+m) auxiliary storage.
        infinity = ieee_value(0.0_dp, ieee_positive_inf)
        allocate(potential_row(0:n), potential_col(0:m), &
            partner(0:m), predecessor(0:m), min_reduced(0:m), visited(0:m))
        potential_row = 0.0_dp
        potential_col = 0.0_dp
        partner = 0

        do i = 1, n
            partner(0) = i
            j0 = 0
            min_reduced = infinity
            predecessor = 0
            visited = .false.
            do
                visited(j0) = .true.
                i0 = partner(j0)
                delta = infinity
                j1 = 0
                do j = 1, m
                    if (visited(j)) cycle
                    if (reversed) then
                        value = cost_matrix(j, i0)
                    else
                        value = cost_matrix(i0, j)
                    end if
                    if (ieee_is_finite(value)) then
                        if (want_maximum) value = -value
                        ! Scaling the costs prevents overflow in the dual
                        ! updates even for input values close to huge(dp).
                        current = (value / scale) - potential_row(i0) - potential_col(j)
                        if (current < min_reduced(j)) then
                            min_reduced(j) = current
                            predecessor(j) = j0
                        end if
                    end if
                    if (min_reduced(j) < delta) then
                        delta = min_reduced(j)
                        j1 = j
                    end if
                end do
                if (j1 == 0) then
                    status = assignment_status_infeasible
                    return
                end if
                do j = 0, m
                    if (visited(j)) then
                        potential_row(partner(j)) = potential_row(partner(j)) + delta
                        potential_col(j) = potential_col(j) - delta
                    else
                        min_reduced(j) = min_reduced(j) - delta
                    end if
                end do
                j0 = j1
                if (partner(j0) == 0) exit
            end do
            do
                j1 = predecessor(j0)
                partner(j0) = partner(j1)
                j0 = j1
                if (j0 == 0) exit
            end do
        end do

        deallocate(row_ind, col_ind)
        allocate(row_ind(n), col_ind(n))
        if (reversed) then
            ! Working columns are the original rows; ascending scan sorts
            ! row_ind as promised by SciPy.
            k = 0
            do j = 1, m
                if (partner(j) == 0) cycle
                k = k + 1
                row_ind(k) = j
                col_ind(k) = partner(j)
            end do
        else
            allocate(match(n))
            match = 0
            do j = 1, m
                if (partner(j) /= 0) match(partner(j)) = j
            end do
            do i = 1, n
                row_ind(i) = i
                col_ind(i) = match(i)
            end do
        end if
    end subroutine linear_sum_assignment

end module scifort_optimize_assignment
