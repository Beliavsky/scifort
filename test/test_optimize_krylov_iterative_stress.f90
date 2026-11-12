! SPDX-License-Identifier: MIT
! Larger deterministic nonsymmetric and symmetric-indefinite Krylov problems.
program test_optimize_krylov_iterative_stress
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_iterative, only : krylov_iterative_solve
    implicit none

    integer, parameter :: n = 16
    real(dp) :: a(n, n), rhs(n), target(n), solution(n), residual(n)
    integer :: i, failures, call_count, max_calls
    logical :: ok
    character(len=8), parameter :: methods(3) = [character(len=8) :: &
        'bicgstab', 'cgs', 'tfqmr']

    failures = 0
    call_count = 0
    max_calls = huge(1)
    do i = 1, n
        target(i) = sin(real(i, dp) * 0.3_dp)
    end do
    a = 0.0_dp
    do i = 1, n
        a(i, i) = 3.0_dp + 0.1_dp * real(i, dp)
        if (i < n) then
            a(i, i + 1) = 0.22_dp
            a(i + 1, i) = 0.37_dp
        end if
        if (i + 2 <= n) a(i, i + 2) = 0.10_dp
    end do
    rhs = matmul(a, target)
    do i = 1, size(methods)
        call_count = 0
        if (trim(methods(i)) == 'tfqmr') then
            ! SciPy TFQMR's left-preconditioned residual estimator may
            ! terminate early for general approximate diagonal inverses.
            call krylov_iterative_solve(trim(methods(i)), matvec, identity, &
                rhs, 1.0e-10_dp, 100, solution, ok)
        else
            call krylov_iterative_solve(trim(methods(i)), matvec, diagonal_inverse, &
                rhs, 1.0e-10_dp, 100, solution, ok)
        end if
        residual = rhs - matmul(a, solution)
        call check(ok, trim(methods(i))//' finite 16-dimensional solution')
        call check(norm2(residual) <= 1.0e-8_dp * norm2(rhs), &
            trim(methods(i))//' 16-dimensional residual accuracy')
        call check(call_count > 0, trim(methods(i))//' matrix-free evaluations')
    end do

    a = 0.0_dp
    do i = 1, n
        a(i, i) = (-1.0_dp)**i * (3.0_dp + 0.1_dp * real(i, dp))
        if (i < n) then
            a(i, i + 1) = 0.45_dp
            a(i + 1, i) = 0.45_dp
        end if
    end do
    rhs = matmul(a, target)
    call krylov_iterative_solve('minres', matvec, positive_inverse, &
        rhs, 1.0e-10_dp, 100, solution, ok)
    residual = rhs - matmul(a, solution)
    call check(ok, 'MINRES finite indefinite 16-dimensional solution')
    call check(norm2(residual) <= 1.0e-8_dp * norm2(rhs), &
        'MINRES 16-dimensional residual accuracy')

    ! A delayed failure must not leave the previous successful callback's
    ! status set to true when a later matrix-vector product fails.
    max_calls = 2
    do i = 1, size(methods)
        a = 0.0_dp
        call set_diagonal(a)
        rhs = matmul(a, target)
        call_count = 0
        call krylov_iterative_solve(trim(methods(i)), matvec, identity, &
            rhs, 1.0e-14_dp, 100, solution, ok)
        call check(.not. ok, trim(methods(i))//' delayed operator failure')
    end do
    call_count = 0
    call krylov_iterative_solve('minres', matvec, identity, rhs, &
        1.0e-14_dp, 100, solution, ok)
    call check(.not. ok, 'MINRES delayed operator failure')

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_krylov_iterative_stress: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_krylov_iterative_stress: PASS'

contains

    subroutine check(condition, message)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: message !! assertion label

        if (.not. condition) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: '//message
        end if
    end subroutine check

    subroutine set_diagonal(matrix)
        real(dp), intent(out) :: matrix(:, :) !! replacement diagonal operator
        integer :: k

        matrix = 0.0_dp
        do k = 1, min(size(matrix, 1), size(matrix, 2))
            matrix(k, k) = 1.0_dp + real(k, dp) / 3.0_dp
        end do
    end subroutine set_diagonal

    subroutine matvec(v, w, valid)
        real(dp), intent(in) :: v(:) !! matrix operand
        real(dp), intent(out) :: w(:) !! product
        logical, intent(out) :: valid !! false when budget exceeded

        call_count = call_count + 1
        valid = call_count <= max_calls
        if (valid) w = matmul(a, v)
    end subroutine matvec

    subroutine diagonal_inverse(v, w, valid)
        real(dp), intent(in) :: v(:) !! preconditioner operand
        real(dp), intent(out) :: w(:) !! inverse diagonal result
        logical, intent(out) :: valid !! preconditioner status
        integer :: k

        do k = 1, size(v)
            w(k) = v(k) / a(k, k)
        end do
        valid = .true.
    end subroutine diagonal_inverse

    subroutine positive_inverse(v, w, valid)
        real(dp), intent(in) :: v(:) !! preconditioner operand
        real(dp), intent(out) :: w(:) !! positive inverse diagonal result
        logical, intent(out) :: valid !! preconditioner status
        integer :: k

        do k = 1, size(v)
            w(k) = v(k) / abs(a(k, k))
        end do
        valid = .true.
    end subroutine positive_inverse

    subroutine identity(v, w, valid)
        real(dp), intent(in) :: v(:) !! preconditioner operand
        real(dp), intent(out) :: w(:) !! identity preconditioner result
        logical, intent(out) :: valid !! preconditioner status

        w = v
        valid = .true.
    end subroutine identity

end program test_optimize_krylov_iterative_stress
