! SPDX-License-Identifier: MIT
! TFQMR using the mathematically consistent right-preconditioned operator A*M.
program test_optimize_krylov_right_preconditioner
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_iterative, only : krylov_iterative_solve
    implicit none

    integer, parameter :: n = 64
    real(dp) :: a(n,n), rhs(n), target(n), solution(n), reference(n)
    real(dp) :: relative_error, default_error
    integer :: i, failures, calls, precond_calls, max_precond_calls
    integer :: default_calls, right_calls
    logical :: ok

    failures = 0
    calls = 0
    precond_calls = 0
    max_precond_calls = huge(1)

    a = 0.0_dp
    do i = 1, n
        target(i) = sin(0.3_dp * real(i, dp))
        a(i,i) = 3.0_dp + 0.1_dp * real(i, dp)
    end do
    do i = 1, n-1
        a(i,i+1) = 0.22_dp
        a(i+1,i) = 0.37_dp
    end do
    do i = 1, n-2
        a(i,i+2) = 0.10_dp
    end do
    rhs = matmul(a, target)

    ! SciPy's original M*A recurrence may stop with an inaccurate solution
    ! for an approximate inverse diagonal, even with a tiny bound estimate.
    call krylov_iterative_solve('tfqmr', matvec, inverse_diagonal, &
        rhs, 1.0e-10_dp, 120, solution, ok)
    default_error = norm2(rhs - matmul(a, solution)) / norm2(rhs)
    default_calls = calls
    call check(ok, 'SciPy-compatible TFQMR returns a finite iterate')
    call check(default_error > 1.0e-5_dp, 'left recurrence reproduces optimistic convergence')

    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, inverse_diagonal, &
        rhs, 1.0e-10_dp, 120, solution, ok, &
        verify_residual=.true., right_preconditioned=.true.)
    relative_error = norm2(rhs - matmul(a, solution)) / norm2(rhs)
    right_calls = calls
    call check(ok, 'right-preconditioned diagonal TFQMR verifies convergence')
    call check(relative_error <= 1.0e-10_dp, 'right-preconditioned diagonal true residual')
    call check(maxval(abs(solution - target)) < 1.0e-8_dp, &
        'right-preconditioned diagonal recovers target')
    call check(right_calls < 120, 'right-preconditioned diagonal terminates before budget')

    ! SciPy 1.17.0 TFQMR on B=A*M without an internal preconditioner uses
    ! ten matvecs; the optional original-residual verification uses one more.
    call check(right_calls == 11, 'SciPy right-transformed diagonal work count plus verification')
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, inverse_diagonal, &
        rhs, 1.0e-10_dp, 120, solution, ok, right_preconditioned=.true.)
    call check(ok .and. calls == 10, 'unverified right TFQMR matches SciPy transformed work count')
    call check(norm2(rhs - matmul(a, solution)) <= 1.0e-10_dp * norm2(rhs), &
        'unverified right preconditioned TFQMR has accurate true residual')

    ! A non-diagonal inverse also tests that the two products are performed
    ! in the correct order; in general M*A and A*M do not commute.
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, approximate_triangular_inverse, &
        rhs, 1.0e-10_dp, 120, solution, ok, &
        verify_residual=.true., right_preconditioned=.true.)
    call check(ok, 'right-preconditioned triangular TFQMR converges')
    call check(calls == 11, 'SciPy right-transformed triangular work count plus verification')
    call check(norm2(rhs - matmul(a, solution)) <= 1.0e-10_dp * norm2(rhs), &
        'triangular preconditioner gives true solution')

    ! Setting right preconditioning with the identity is a no-op and must
    ! leave SciPy-compatible iteration counts and solutions unchanged.
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, identity, &
        rhs, 1.0e-11_dp, 120, reference, ok)
    default_calls = calls
    call check(ok, 'identity default succeeds')
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, identity, &
        rhs, 1.0e-11_dp, 120, solution, ok, right_preconditioned=.true.)
    call check(ok .and. calls == default_calls, 'identity unchanged TFQMR call count')
    call check(maxval(abs(solution - reference)) < 1.0e-14_dp, &
        'identity unchanged TFQMR solution')

    ! Ensure both operator errors and preconditioner errors are propagated,
    ! including errors before the first Jacobian-vector product.
    calls = 0
    precond_calls = 0
    max_precond_calls = 0
    call krylov_iterative_solve('tfqmr', matvec, failing_inverse, &
        rhs, 1.0e-10_dp, 120, solution, ok, right_preconditioned=.true.)
    call check(.not. ok .and. calls == 0, 'right mode rejects failed preconditioner before A')
    max_precond_calls = 1
    calls = 0
    precond_calls = 0
    call krylov_iterative_solve('tfqmr', matvec, failing_inverse, &
        rhs, 1.0e-10_dp, 120, solution, ok, right_preconditioned=.true.)
    call check(.not. ok .and. calls >= 1, 'right mode propagates delayed preconditioner failure')
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, nan_inverse, &
        rhs, 1.0e-10_dp, 120, solution, ok, right_preconditioned=.true.)
    call check(.not. ok .and. calls == 0, 'right mode rejects nonfinite preconditioner before A')
    max_precond_calls = huge(1)

    call krylov_iterative_solve('cgs', matvec, identity, &
        rhs, 1.0e-10_dp, 120, solution, ok, right_preconditioned=.true.)
    call check(.not. ok, 'right flag rejected for non-TFQMR method')
    call krylov_iterative_solve('tfqmr', matvec, identity, &
        0.0_dp*rhs, 1.0e-10_dp, 120, solution, ok, right_preconditioned=.true.)
    call check(ok .and. all(abs(solution) < tiny(1.0_dp)), &
        'right TFQMR zero right-hand side')

    if (failures /= 0) then
        write(*,'(a,i0)') 'test_optimize_krylov_right_preconditioner: FAIL ',failures
        stop 1
    end if
    write(*,'(a)') 'test_optimize_krylov_right_preconditioner: PASS'

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition !! expected condition
        character(len=*), intent(in) :: description !! assertion description

        if (.not. condition) then
            failures = failures + 1
            write(*,'(a)') 'FAIL: '//description
        end if
    end subroutine check

    subroutine matvec(v, product, valid)
        real(dp), intent(in) :: v(:) !! vector operand
        real(dp), intent(out) :: product(:) !! A*v
        logical, intent(out) :: valid !! operator evaluation status

        calls = calls + 1
        product = matmul(a, v)
        valid = .true.
    end subroutine matvec

    subroutine identity(v, product, valid)
        real(dp), intent(in) :: v(:) !! vector operand
        real(dp), intent(out) :: product(:) !! unchanged vector
        logical, intent(out) :: valid !! operation status

        product = v
        valid = .true.
    end subroutine identity

    subroutine inverse_diagonal(v, product, valid)
        real(dp), intent(in) :: v(:) !! vector operand
        real(dp), intent(out) :: product(:) !! diagonal inverse product
        logical, intent(out) :: valid !! operation status
        integer :: j

        do j = 1, size(v)
            product(j) = v(j) / a(j,j)
        end do
        valid = .true.
    end subroutine inverse_diagonal

    subroutine approximate_triangular_inverse(v, product, valid)
        real(dp), intent(in) :: v(:) !! vector operand
        real(dp), intent(out) :: product(:) !! approximate triangular inverse action
        logical, intent(out) :: valid !! operation status
        integer :: j

        do j = 1, size(v)
            product(j) = v(j) / a(j,j)
            if (j < size(v)) product(j) = product(j) - 0.06_dp * v(j+1) / a(j+1,j+1)
        end do
        valid = .true.
    end subroutine approximate_triangular_inverse

    subroutine failing_inverse(v, product, valid)
        real(dp), intent(in) :: v(:) !! vector operand
        real(dp), intent(out) :: product(:) !! inverse operator result
        logical, intent(out) :: valid !! status

        precond_calls = precond_calls + 1
        valid = precond_calls <= max_precond_calls
        if (valid) product = v
    end subroutine failing_inverse

    subroutine nan_inverse(v, product, valid)
        real(dp), intent(in) :: v(:) !! vector operand
        real(dp), intent(out) :: product(:) !! deliberately invalid result
        logical, intent(out) :: valid !! claims an otherwise valid operation

        product = ieee_value(v, ieee_quiet_nan)
        valid = .true.
    end subroutine nan_inverse

end program test_optimize_krylov_right_preconditioner
