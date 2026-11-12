! SPDX-License-Identifier: MIT
! Opt-in TFQMR true-residual verification, with matrix-free error injection.
program test_optimize_krylov_true_residual
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_iterative, only : krylov_iterative_solve
    implicit none

    integer, parameter :: n = 16
    real(dp) :: a(n, n), rhs(n), target(n), solution(n), missing(n)
    real(dp) :: poor_residual
    integer :: i, failures, calls, max_calls, default_calls, preconditioned_calls
    logical :: ok

    failures = 0
    calls = 0
    max_calls = huge(1)
    a = 0.0_dp
    do i = 1, n
        target(i) = sin(0.3_dp * real(i, dp))
        a(i, i) = 3.0_dp + 0.1_dp * real(i, dp)
        if (i < n) then
            a(i, i + 1) = 0.22_dp
            a(i + 1, i) = 0.37_dp
        end if
        if (i + 2 <= n) a(i, i + 2) = 0.10_dp
    end do
    rhs = matmul(a, target)

    ! SciPy's estimate-based TFQMR may report success with an optimistic
    ! residual bound when using an approximate diagonal preconditioner.
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, approximate_inverse, &
        rhs, 1.0e-10_dp, 60, solution, ok)
    poor_residual = norm2(rhs - matmul(a, solution)) / norm2(rhs)
    preconditioned_calls = calls
    call check(ok, 'SciPy-compatible finite iterate is accepted')
    call check(poor_residual > 1.0e-5_dp, &
        'reference preconditioner exposes optimistic TFQMR estimate')

    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, approximate_inverse, &
        rhs, 1.0e-10_dp, 60, solution, ok, verify_residual=.true.)
    call check(.not. ok, 'strict mode rejects false convergence')
    call check(calls > preconditioned_calls, &
        'strict mode attempts further steps after optimistic estimate')
    call check(norm2(rhs - matmul(a, solution)) > 1.0e-8_dp * norm2(rhs), &
        'strict rejection is due to an unmet original residual')

    ! For identity preconditioning, verification changes neither the
    ! computed answer nor its mathematical accuracy, only its validation.
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, identity, &
        rhs, 1.0e-10_dp, 60, solution, ok)
    default_calls = calls
    call check(ok, 'unpreconditioned TFQMR reference solve succeeds')
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, identity, &
        rhs, 1.0e-10_dp, 60, solution, ok, verify_residual=.true.)
    call check(ok, 'strict TFQMR succeeds when original residual is small')
    call check(norm2(rhs - matmul(a, solution)) <= 1.0e-10_dp * norm2(rhs), &
        'strict solve satisfies original relative residual tolerance')
    call check(calls > default_calls, 'strict solve performs a verification product')

    ! The last accepted SciPy-compatible iterate can be finite but not yet
    ! converged; strict verification must not silently accept that iterate.
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, identity, &
        rhs, 1.0e-12_dp, 2, solution, ok)
    call check(ok, 'default accepts finite truncated iterate')
    call krylov_iterative_solve('tfqmr', matvec, identity, &
        rhs, 1.0e-12_dp, 2, solution, ok, verify_residual=.true.)
    call check(.not. ok, 'strict rejects inaccurate truncated iterate')
    call check(norm2(solution) > 0.0_dp, 'strict returns partial iterate for diagnosis')

    ! Force a failure exactly when the first would-be converged iterate is
    ! being checked, after all the legacy matrix-vector products succeeded.
    calls = 0
    max_calls = default_calls
    call krylov_iterative_solve('tfqmr', matvec, identity, &
        rhs, 1.0e-10_dp, 60, solution, ok, verify_residual=.true.)
    call check(.not. ok, 'failure during verification rejects the solve')
    call check(calls == default_calls + 1, &
        'late matrix-free operator failure occurs during verification')
    max_calls = huge(1)

    ! Zero right-hand sides must not cause a spurious operator call.
    calls = 0
    call krylov_iterative_solve('tfqmr', matvec, identity, &
        0.0_dp * rhs, 1.0e-10_dp, 60, solution, ok, verify_residual=.true.)
    call check(ok .and. all(solution == 0.0_dp), 'strict zero RHS succeeds')
    call check(calls == 0, 'zero RHS requires no matrix-vector products')

    call krylov_iterative_solve('cgs', matvec, identity, rhs, &
        1.0e-10_dp, 60, solution, ok, verify_residual=.true.)
    call check(.not. ok, 'reject strict verification for unsupported method')
    call krylov_iterative_solve('bicgstab', matvec, identity, rhs, &
        1.0e-10_dp, 60, solution, ok, verify_residual=.false.)
    call check(ok, 'explicit false flag leaves other methods unchanged')
    call krylov_iterative_solve('tfqmr', matvec, identity, rhs, &
        1.0e-10_dp, 60, missing(1:n-1), ok, verify_residual=.true.)
    call check(.not. ok, 'reject mismatched solution shape')
    missing = rhs
    missing(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call krylov_iterative_solve('tfqmr', matvec, identity, missing, &
        1.0e-10_dp, 60, solution, ok, verify_residual=.true.)
    call check(.not. ok, 'reject nonfinite RHS before operator evaluation')

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_krylov_true_residual: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_krylov_true_residual: PASS'

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: description !! assertion label

        if (.not. condition) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: '//description
        end if
    end subroutine check

    subroutine matvec(vec, product, valid)
        real(dp), intent(in) :: vec(:) !! input for matrix-free action
        real(dp), intent(out) :: product(:) !! matrix-vector product
        logical, intent(out) :: valid !! false when call budget exhausted

        calls = calls + 1
        valid = calls <= max_calls
        if (valid) product = matmul(a, vec)
    end subroutine matvec

    subroutine identity(vec, product, valid)
        real(dp), intent(in) :: vec(:) !! preconditioner operand
        real(dp), intent(out) :: product(:) !! unchanged operand
        logical, intent(out) :: valid !! status

        product = vec
        valid = .true.
    end subroutine identity

    subroutine approximate_inverse(vec, product, valid)
        real(dp), intent(in) :: vec(:) !! preconditioner operand
        real(dp), intent(out) :: product(:) !! inverse-diagonal product
        logical, intent(out) :: valid !! status
        integer :: j

        do j = 1, size(vec)
            product(j) = vec(j) / a(j, j)
        end do
        valid = .true.
    end subroutine approximate_inverse

end program test_optimize_krylov_true_residual
