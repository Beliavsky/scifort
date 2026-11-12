! SPDX-License-Identifier: MIT
! 64-variable deterministic nonsymmetric and symmetric-indefinite systems.
program test_optimize_krylov_large_system
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_iterative, only : krylov_iterative_solve
    implicit none

    integer, parameter :: n = 64
    character(len=8), parameter :: methods(3) = [character(len=8) :: &
        'bicgstab', 'cgs', 'tfqmr']
    real(dp) :: a(n, n), rhs(n), target(n), solution(n)
    integer :: calls, failures, i, k
    logical :: ok

    failures = 0
    calls = 0
    do i = 1, n
        target(i) = cos(0.21_dp * real(i, dp))
    end do

    a = 0.0_dp
    do i = 1, n
        a(i, i) = 2.5_dp + 0.08_dp * real(mod(i, 13), dp)
        if (i < n) then
            a(i, i + 1) = 0.35_dp
            a(i + 1, i) = -0.24_dp
        end if
        if (i + 3 <= n) a(i, i + 3) = 0.07_dp
    end do
    rhs = matmul(a, target)

    do k = 1, size(methods)
        calls = 0
        if (trim(methods(k)) == 'tfqmr') then
            call krylov_iterative_solve('tfqmr', matvec, identity, &
                rhs, 1.0e-11_dp, 160, solution, ok, verify_residual=.true.)
        else
            call krylov_iterative_solve(trim(methods(k)), matvec, inverse_diag, &
                rhs, 1.0e-11_dp, 160, solution, ok)
        end if
        call check(ok, trim(methods(k))//' 64-variable solve reports success')
        call check(norm2(rhs - matmul(a, solution)) <= 1.0e-9_dp * norm2(rhs), &
            trim(methods(k))//' 64-variable true residual')
        call check(calls > 0 .and. calls < 160, &
            trim(methods(k))//' 64-variable matrix-free calls')
    end do

    ! Check 64-variable MINRES with alternating positive and negative
    ! eigenvalues, while keeping its inverse preconditioner positive definite.
    a = 0.0_dp
    do i = 1, n
        a(i, i) = (-1.0_dp)**i * (3.0_dp + 0.02_dp * real(i, dp))
        if (i < n) then
            a(i, i + 1) = 0.18_dp
            a(i + 1, i) = 0.18_dp
        end if
    end do
    rhs = matmul(a, target)
    calls = 0
    call krylov_iterative_solve('minres', matvec, inverse_abs_diag, &
        rhs, 1.0e-11_dp, 160, solution, ok)
    call check(ok, 'MINRES indefinite 64-variable solve succeeds')
    call check(norm2(rhs - matmul(a, solution)) <= 1.0e-9_dp * norm2(rhs), &
        'MINRES indefinite 64-variable true residual')
    call check(calls > 0 .and. calls < 160, &
        'MINRES 64-variable matrix-free evaluations')

    ! A callback that returns success together with NaNs must be rejected,
    ! independently of failures from the matrix-vector operator itself.
    call krylov_iterative_solve('minres', matvec, nan_preconditioner, &
        rhs, 1.0e-9_dp, 30, solution, ok)
    call check(.not. ok, 'MINRES rejects nonfinite preconditioner output')
    call krylov_iterative_solve('tfqmr', matvec, nan_preconditioner, &
        rhs, 1.0e-9_dp, 30, solution, ok, verify_residual=.true.)
    call check(.not. ok, 'TFQMR rejects nonfinite preconditioner output')

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_krylov_large_system: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_krylov_large_system: PASS'

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition !! assertion result
        character(len=*), intent(in) :: description !! assertion name

        if (.not. condition) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: '//description
        end if
    end subroutine check

    subroutine matvec(vec, product, valid)
        real(dp), intent(in) :: vec(:) !! operator operand
        real(dp), intent(out) :: product(:) !! product with test matrix
        logical, intent(out) :: valid !! application status

        calls = calls + 1
        product = matmul(a, vec)
        valid = .true.
    end subroutine matvec

    subroutine identity(vec, product, valid)
        real(dp), intent(in) :: vec(:) !! preconditioner operand
        real(dp), intent(out) :: product(:) !! unchanged vector
        logical, intent(out) :: valid !! application status

        product = vec
        valid = .true.
    end subroutine identity

    subroutine inverse_diag(vec, product, valid)
        real(dp), intent(in) :: vec(:) !! preconditioner operand
        real(dp), intent(out) :: product(:) !! inverse-diagonal result
        logical, intent(out) :: valid !! application status
        integer :: j

        do j = 1, size(vec)
            product(j) = vec(j) / a(j, j)
        end do
        valid = .true.
    end subroutine inverse_diag

    subroutine inverse_abs_diag(vec, product, valid)
        real(dp), intent(in) :: vec(:) !! preconditioner operand
        real(dp), intent(out) :: product(:) !! positive inverse diagonal
        logical, intent(out) :: valid !! application status
        integer :: j

        do j = 1, size(vec)
            product(j) = vec(j) / abs(a(j, j))
        end do
        valid = .true.
    end subroutine inverse_abs_diag

    subroutine nan_preconditioner(vec, product, valid)
        real(dp), intent(in) :: vec(:) !! preconditioner operand
        real(dp), intent(out) :: product(:) !! deliberately nonfinite output
        logical, intent(out) :: valid !! reported application status

        product = vec
        product(1) = ieee_value(0.0_dp, ieee_quiet_nan)
        valid = .true.
    end subroutine nan_preconditioner

end program test_optimize_krylov_large_system
