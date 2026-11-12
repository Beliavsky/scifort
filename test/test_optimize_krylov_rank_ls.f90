! SPDX-License-Identifier: BSD-3-Clause
! Rank-revealing LGMRES least-squares regression, including ill conditioning.
program test_optimize_krylov_rank_ls
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_rank_ls, only : krylov_rank_least_squares
    implicit none

    real(dp) :: a(3, 3)
    real(dp) :: b(3)
    real(dp) :: x(3)
    real(dp) :: nearly(2, 2)
    real(dp) :: two_b(2)
    real(dp) :: two_x(2)
    real(dp) :: rotated(4, 4)
    real(dp) :: left_rotation(4, 4)
    real(dp) :: right_rotation(4, 4)
    real(dp) :: singular_values(4, 4)
    real(dp) :: rotated_rhs(4)
    real(dp) :: rotated_x(4)
    real(dp) :: rotated_expected(4)
    real(dp) :: rectangular(5, 3)
    real(dp) :: rectangular_b(5)
    real(dp) :: exact_x(3)
    real(dp), parameter :: delta = 1.0e-10_dp
    real(dp) :: scale
    integer :: rank
    integer :: failures
    integer :: j
    logical :: ok

    failures = 0
    a = 0.0_dp
    a(1, 1) = 1.0_dp
    a(1, 2) = 1.0_dp
    a(3, 3) = 2.0_dp
    b = [2.0_dp, 1.0_dp, 4.0_dp]
    call krylov_rank_least_squares(a, b, x, rank, ok)
    call check(ok .and. rank == 2, 'exact rank deficiency detected')
    call check(maxval(abs(x - [1.0_dp, 1.0_dp, 2.0_dp])) < 1.0e-13_dp, &
        'rank-deficient minimum-norm solution')

    ! Normal equations lose the 1e-10 singular direction, since its squared
    ! singular value is near 1e-20. The unsquared Jacobi system retains it.
    nearly(:, 1) = [1.0_dp, 0.0_dp]
    nearly(:, 2) = [1.0_dp, delta]
    two_b = [0.0_dp, delta]
    call krylov_rank_least_squares(nearly, two_b, two_x, rank, ok)
    call check(ok .and. rank == 2, 'near-dependent columns retain rank')
    call check(maxval(abs(two_x - [-1.0_dp, 1.0_dp])) < 3.0e-6_dp, &
        'small singular direction is resolved')
    call check(maxval(abs(matmul(nearly, two_x) - two_b)) < 1.0e-14_dp, &
        'near-dependent least-squares residual')

    ! Columns smaller than the LAPACK-like relative rank cutoff are dropped.
    nearly(2, 2) = 1.0e-17_dp
    two_b = [2.0_dp, 0.0_dp]
    call krylov_rank_least_squares(nearly, two_b, two_x, rank, ok)
    call check(ok .and. rank == 1 .and. maxval(abs(two_x - 1.0_dp)) < 1.0e-13_dp, &
        'sub-roundoff singular direction is truncated')

    ! Large and small common scales must not change rank or coefficients.
    nearly(2, 2) = delta
    two_b = [0.0_dp, delta]
    do j = 1, 2
        if (j == 1) then
            scale = 1.0e-150_dp
        else
            scale = 1.0e150_dp
        end if
        call krylov_rank_least_squares(scale * nearly, scale * two_b, two_x, rank, ok)
        call check(ok .and. rank == 2, 'scaled near-dependent rank')
        call check(maxval(abs(two_x - [-1.0_dp, 1.0_dp])) < 3.0e-6_dp, &
            'scaled near-dependent coefficients')
    end do

    rectangular = 0.0_dp
    rectangular(:, 1) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    rectangular(:, 2) = [0.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    rectangular(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    exact_x = [2.0_dp, -3.0_dp, 4.0_dp]
    rectangular_b = matmul(rectangular, exact_x)
    call krylov_rank_least_squares(rectangular, rectangular_b, x, rank, ok)
    call check(ok .and. rank == 3, 'rectangular full-column-rank matrix')
    call check(maxval(abs(x - exact_x)) < 2.0e-13_dp, 'rectangular solution')

    ! Dense singular vectors make small singular modes invisible in the
    ! elementwise matrix entries. Compare the recovered truncated-SVD
    ! solution with its independently constructed orthogonal factors.
    left_rotation = 0.0_dp
    right_rotation = 0.0_dp
    do j = 1, 4
        left_rotation(j, j) = 1.0_dp
        right_rotation(j, j) = 1.0_dp
    end do
    left_rotation(1, 1) = cos(0.3_dp)
    left_rotation(3, 3) = cos(0.3_dp)
    left_rotation(1, 3) = -sin(0.3_dp)
    left_rotation(3, 1) = sin(0.3_dp)
    left_rotation(2, 2) = cos(0.5_dp)
    left_rotation(4, 4) = cos(0.5_dp)
    left_rotation(2, 4) = -sin(0.5_dp)
    left_rotation(4, 2) = sin(0.5_dp)
    right_rotation(1, 1) = cos(0.7_dp)
    right_rotation(2, 2) = cos(0.7_dp)
    right_rotation(1, 2) = -sin(0.7_dp)
    right_rotation(2, 1) = sin(0.7_dp)
    right_rotation(3, 3) = cos(0.4_dp)
    right_rotation(4, 4) = cos(0.4_dp)
    right_rotation(3, 4) = -sin(0.4_dp)
    right_rotation(4, 3) = sin(0.4_dp)
    singular_values = 0.0_dp
    singular_values(1, 1) = 3.0_dp
    singular_values(2, 2) = 1.0e-6_dp
    singular_values(3, 3) = 1.0e-10_dp
    rotated = matmul(left_rotation, matmul(singular_values, transpose(right_rotation)))
    rotated_rhs = matmul(left_rotation, [3.0_dp, 1.0e-6_dp, 1.0e-10_dp, 0.0_dp])
    rotated_expected = matmul(right_rotation, [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp])
    call krylov_rank_least_squares(rotated, rotated_rhs, rotated_x, rank, ok)
    call check(ok .and. rank == 3, 'rotated near-deficient dense matrix rank')
    call check(maxval(abs(rotated_x - rotated_expected)) < 1.0e-5_dp, &
        'rotated near-deficient dense matrix coefficients')
    call check(maxval(abs(matmul(rotated, rotated_x) - rotated_rhs)) < 1.0e-12_dp, &
        'rotated near-deficient dense matrix residual')

    a = 0.0_dp
    b = [1.0_dp, 2.0_dp, 3.0_dp]
    call krylov_rank_least_squares(a, b, x, rank, ok)
    call check(ok .and. rank == 0 .and. all(x == 0.0_dp), &
        'zero operator has zero minimum-norm solution')

    ! Invalid shapes do not read out of bounds.
    call krylov_rank_least_squares(a, b(1:2), x, rank, ok)
    call check(.not. ok, 'invalid right-hand side dimensions')
    call krylov_rank_least_squares(a, b, x(1:2), rank, ok)
    call check(.not. ok, 'invalid solution dimensions')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_krylov_rank_ls: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_krylov_rank_ls: PASS'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected result
        character(len=*), intent(in) :: label !! regression case name

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

end program test_optimize_krylov_rank_ls
