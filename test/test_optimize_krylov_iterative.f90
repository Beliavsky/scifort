! SPDX-License-Identifier: MIT
! Independent matrix-free Krylov inner-solver regression tests.
program test_optimize_krylov_iterative
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_iterative, only : krylov_iterative_solve
    implicit none

    integer, parameter :: n = 4
    real(dp), parameter :: target(n) = [1.0_dp, -2.0_dp, 3.0_dp, -1.0_dp]
    real(dp) :: a(n, n), rhs(n), solution(n)
    integer :: failures, matvec_calls, preconditioner_calls, i
    integer, parameter :: scipy_nmatvec(3) = [7, 8, 7]
    logical :: ok
    character(len=8), parameter :: methods(3) = [character(len=8) :: &
        'bicgstab', 'cgs', 'tfqmr']

    failures = 0
    matvec_calls = 0
    preconditioner_calls = 0

    a = 0.0_dp
    a(1, :) = [4.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    a(2, :) = [0.0_dp, 3.0_dp, 1.0_dp, 0.0_dp]
    a(3, :) = [0.0_dp, 0.0_dp, 2.0_dp, 1.0_dp]
    a(4, :) = [1.0_dp, 0.0_dp, 0.0_dp, 5.0_dp]
    rhs = matmul(a, target)

    do i = 1, size(methods)
        matvec_calls = 0
        call krylov_iterative_solve(trim(methods(i)), matrix_action, identity, &
            rhs, 1.0e-10_dp, 60, solution, ok)
        call check(ok, trim(methods(i))//' finite solve')
        call check(maxval(abs(solution - target)) < 1.0e-8_dp, &
            trim(methods(i))//' nonsymmetric system solution')
        call check(matvec_calls == scipy_nmatvec(i), &
            trim(methods(i))//' matrix-vector calls match SciPy 1.17')
    end do

    ! Exact diagonal preconditioning for all nonsymmetric inner methods.
    a = 0.0_dp
    a(1, 1) = 2.0_dp
    a(2, 2) = 3.0_dp
    a(3, 3) = 5.0_dp
    a(4, 4) = 7.0_dp
    rhs = matmul(a, target)
    do i = 1, size(methods)
        preconditioner_calls = 0
        call krylov_iterative_solve(trim(methods(i)), matrix_action, &
            diagonal_inverse, rhs, 1.0e-10_dp, 20, solution, ok)
        call check(ok .and. maxval(abs(solution - target)) < 1.0e-9_dp, &
            trim(methods(i))//' exact diagonal preconditioner')
        call check(preconditioner_calls > 0, &
            trim(methods(i))//' preconditioner invoked')
    end do

    ! Symmetric indefinite matrices are the appropriate MINRES domain.
    a = 0.0_dp
    a(1, :) = [4.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    a(2, :) = [1.0_dp, -3.0_dp, 1.0_dp, 0.0_dp]
    a(3, :) = [0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp]
    a(4, :) = [0.0_dp, 0.0_dp, 1.0_dp, -4.0_dp]
    rhs = matmul(a, target)
    matvec_calls = 0
    call krylov_iterative_solve('minres', matrix_action, identity, rhs, &
        1.0e-10_dp, 60, solution, ok)
    call check(ok, 'MINRES symmetric indefinite solve')
    call check(maxval(abs(solution - target)) < 1.0e-8_dp, &
        'MINRES symmetric indefinite solution')
    call check(matvec_calls == 4, 'MINRES Lanczos iteration count matches SciPy')

    ! A scaled preconditioner is symmetric positive definite even when A is indefinite.
    preconditioner_calls = 0
    call krylov_iterative_solve('minres', matrix_action, scaled_identity, rhs, &
        1.0e-10_dp, 60, solution, ok)
    call check(ok .and. maxval(abs(solution - target)) < 1.0e-8_dp, &
        'preconditioned MINRES solution')
    call check(preconditioner_calls > 0, 'MINRES uses preconditioner')

    call krylov_iterative_solve('minres', matrix_action, negative_identity, rhs, &
        1.0e-10_dp, 60, solution, ok)
    call check(.not. ok, 'MINRES rejects negative-definite preconditioner')

    ! The chosen iterate is returned after a nonconverged, finite iteration budget.
    a(1, :) = [5.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    rhs = matmul(a, target)
    call krylov_iterative_solve('bicgstab', matrix_action, identity, rhs, &
        1.0e-12_dp, 1, solution, ok)
    call check(ok .and. all(abs(solution) < huge(1.0_dp)), &
        'BiCGSTAB returns finite incomplete iterates')
    call check(norm2(solution) > 0.0_dp, 'one BiCGSTAB step makes progress')

    do i = 1, size(methods)
        call krylov_iterative_solve(trim(methods(i)), matrix_action, identity, &
            [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, 3, solution, ok)
        call check(ok .and. all(solution == 0.0_dp), &
            trim(methods(i))//' handles zero right-hand side')
    end do
    call krylov_iterative_solve('minres', matrix_action, identity, &
        [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, 3, solution, ok)
    call check(ok .and. all(solution == 0.0_dp), 'MINRES zero right-hand side')

    do i = 1, size(methods)
        call krylov_iterative_solve(trim(methods(i)), failing_action, identity, &
            rhs, 1.0e-8_dp, 3, solution, ok)
        call check(.not. ok, trim(methods(i))//' matrix-free operator failure')
        call krylov_iterative_solve(trim(methods(i)), nan_action, identity, &
            rhs, 1.0e-8_dp, 3, solution, ok)
        call check(.not. ok, trim(methods(i))//' nonfinite product rejection')
    end do
    call krylov_iterative_solve('minres', failing_action, identity, rhs, &
        1.0e-8_dp, 3, solution, ok)
    call check(.not. ok, 'MINRES operator failure')
    call krylov_iterative_solve('minres', nan_action, identity, rhs, &
        1.0e-8_dp, 3, solution, ok)
    call check(.not. ok, 'MINRES nonfinite product rejection')
    call krylov_iterative_solve('nonsense', matrix_action, identity, rhs, &
        1.0e-8_dp, 3, solution, ok)
    call check(.not. ok, 'invalid method rejected')
    call krylov_iterative_solve('cgs', matrix_action, identity, rhs, &
        1.0e-8_dp, 0, solution, ok)
    call check(.not. ok, 'zero iteration budget rejected')

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_krylov_iterative: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_krylov_iterative: PASS'

contains

    subroutine check(condition, name)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: name !! assertion description

        if (.not. condition) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: '//name
        end if
    end subroutine check

    subroutine matrix_action(v, w, valid)
        real(dp), intent(in) :: v(:) !! matrix operand
        real(dp), intent(out) :: w(:) !! matrix product
        logical, intent(out) :: valid !! application status

        matvec_calls = matvec_calls + 1
        w = matmul(a, v)
        valid = .true.
    end subroutine matrix_action

    subroutine identity(v, w, valid)
        real(dp), intent(in) :: v(:) !! input vector
        real(dp), intent(out) :: w(:) !! unchanged vector
        logical, intent(out) :: valid !! application status

        w = v
        valid = .true.
    end subroutine identity

    subroutine scaled_identity(v, w, valid)
        real(dp), intent(in) :: v(:) !! input vector
        real(dp), intent(out) :: w(:) !! scaled vector
        logical, intent(out) :: valid !! application status

        preconditioner_calls = preconditioner_calls + 1
        w = 0.5_dp * v
        valid = .true.
    end subroutine scaled_identity

    subroutine diagonal_inverse(v, w, valid)
        real(dp), intent(in) :: v(:) !! input vector
        real(dp), intent(out) :: w(:) !! exact inverse-diagonal action
        logical, intent(out) :: valid !! application succeeded
        integer :: k

        preconditioner_calls = preconditioner_calls + 1
        do k = 1, size(v)
            w(k) = v(k) / a(k, k)
        end do
        valid = .true.
    end subroutine diagonal_inverse

    subroutine negative_identity(v, w, valid)
        real(dp), intent(in) :: v(:) !! input vector
        real(dp), intent(out) :: w(:) !! scaled vector
        logical, intent(out) :: valid !! application status

        w = -v
        valid = .true.
    end subroutine negative_identity

    subroutine failing_action(v, w, valid)
        real(dp), intent(in) :: v(:) !! vector operand
        real(dp), intent(out) :: w(:) !! invalid product
        logical, intent(out) :: valid !! failed status

        w = 0.0_dp * v
        valid = .false.
    end subroutine failing_action

    subroutine nan_action(v, w, valid)
        real(dp), intent(in) :: v(:) !! vector operand
        real(dp), intent(out) :: w(:) !! nonfinite product
        logical, intent(out) :: valid !! nominally succeeded status

        w = v
        w(1) = ieee_value(0.0_dp, ieee_quiet_nan)
        valid = .true.
    end subroutine nan_action

end program test_optimize_krylov_iterative
