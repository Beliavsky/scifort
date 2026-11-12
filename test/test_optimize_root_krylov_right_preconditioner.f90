! SPDX-License-Identifier: MIT
! Public TFQMR right-preconditioning, solver forwarding and input validation.
program test_optimize_root_krylov_right_preconditioner
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 16
    real(dp) :: a(n,n), rhs(n), target(n), start(n)
    type(vector_root_result) :: res
    integer :: failures, inverse_calls, i, solver_calls

    failures = 0
    inverse_calls = 0
    solver_calls = 0
    start = 0.0_dp
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

    ! A moderate rdiff avoids cancellation in small finite-difference
    ! Jacobian products near the solution of this linear problem.
    ! The public root wrapper must pass the extra option all the way through
    ! the inexact-Newton controller to the matrix-free TFQMR recurrence.
    res = root(residual, start, method='krylov', inner_method='tfqmr', &
        inner_m=inverse_diagonal, inner_verify_residual=.true., &
        inner_tfqmr_right_preconditioned=.true., inner_maxiter=60, &
        line_search='none', fatol=1.0e-9_dp, rdiff=1.0e-4_dp)
    call check(res%success, 'public root with right-preconditioned TFQMR converges')
    call check(maxval(abs(res%x-target)) <= 1.0e-7_dp, &
        'public right-preconditioned root accuracy')
    call check(maxval(abs(res%fun)) <= 1.0e-9_dp, &
        'root reports accurate original residual')
    call check(inverse_calls > 0, 'public right preconditioner invoked')

    inverse_calls = 0
    res = newton_krylov(residual, start, inner_method='tfqmr', &
        inner_m=inverse_diagonal, inner_verify_residual=.true., &
        inner_tfqmr_right_preconditioned=.true., inner_maxiter=60, &
        line_search='none', f_tol=1.0e-9_dp, rdiff=1.0e-4_dp)
    call check(res%success .and. maxval(abs(res%fun)) <= 1.0e-9_dp, &
        'newton_krylov forwards right TFQMR mode')
    call check(inverse_calls > 0, 'direct newton_krylov calls inverse operator')

    ! Fortran optional logicals set to false are harmless for all methods.
    res = root(residual, start, method='krylov', inner_method='bicgstab', &
        inner_tfqmr_right_preconditioned=.false., inner_maxiter=60, &
        line_search='none', fatol=1.0e-9_dp, rdiff=1.0e-6_dp)
    call check(res%success .and. maxval(abs(res%fun)) <= 1.0e-9_dp, &
        'false flag leaves other Krylov methods available')

    res = root(residual, start, method='krylov', inner_method='cgs', &
        inner_tfqmr_right_preconditioned=.true.)
    call check(.not. res%success .and. index(res%message, 'requires the built-in tfqmr') > 0, &
        'reject true right mode with other Krylov method')
    call check(res%nfev == 0, 'invalid Krylov method/flag combination does not evaluate residual')

    res = root(residual, start, method='broyden1', &
        inner_tfqmr_right_preconditioned=.true.)
    call check(.not. res%success .and. index(res%message, 'requires method=krylov') > 0, &
        'reject true right mode with non-Krylov root')
    call check(res%nfev == 0, 'non-Krylov misuse detected before function call')

    solver_calls = 0
    res = newton_krylov(residual, start, inner_method='tfqmr', &
        inner_solver=custom_solver, inner_tfqmr_right_preconditioned=.true.)
    call check(.not. res%success .and. index(res%message, 'requires the built-in tfqmr') > 0, &
        'reject right mode with custom inner solver')
    call check(solver_calls == 0, 'invalid override rejected before custom inner solver')

    if (failures /= 0) then
        write(*,'(a,i0)') 'test_optimize_root_krylov_right_preconditioner: FAIL ',failures
        stop 1
    end if
    write(*,'(a)') 'test_optimize_root_krylov_right_preconditioner: PASS'

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition !! expected assertion
        character(len=*), intent(in) :: description !! assertion label

        if (.not. condition) then
            failures = failures + 1
            write(*,'(a)') 'FAIL: '//description
        end if
    end subroutine check

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(out) :: fvec(:) !! original residual

        fvec = matmul(a,x) - rhs
    end subroutine residual

    subroutine inverse_diagonal(v, product, ok)
        real(dp), intent(in) :: v(:) !! preconditioner input
        real(dp), intent(out) :: product(:) !! approximate inverse-Jacobian action
        logical, intent(out) :: ok !! preconditioner status
        integer :: j

        inverse_calls = inverse_calls + 1
        do j = 1, size(v)
            product(j) = v(j) / a(j,j)
        end do
        ok = .true.
    end subroutine inverse_diagonal

    subroutine custom_solver(matvec, preconditioner, b, rtol, maxiter, sol, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian-vector action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse action
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(in) :: rtol !! inner relative tolerance
        integer, intent(in) :: maxiter !! iteration budget
        real(dp), intent(out) :: sol(:) !! approximate solution
        logical, intent(out) :: ok !! status

        solver_calls = solver_calls + 1
        sol = 0.0_dp
        ok = .false.
    end subroutine custom_solver

end program test_optimize_root_krylov_right_preconditioner
