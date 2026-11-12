! SPDX-License-Identifier: MIT
! Seeded diagonal probes and optional exact-coordinate verification.
program test_optimize_root_krylov_random_diagonal
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 64
    real(dp) :: diagonal(n), x0(n), observed(n), saved(n)
    real(dp) :: walsh_error, random_error, distant
    type(vector_root_result) :: res
    integer :: jv_calls, f_calls, solver_calls, failures, i, ns
    integer, allocatable :: random_seed_values(:)
    real(dp) :: rand_before, rand_after, expected_before, expected_after
    logical :: bad_jvp

    failures = 0
    bad_jvp = .false.
    do i = 1, n
        diagonal(i) = 1.0_dp + 0.01_dp * real(i, dp)
    end do
    x0 = 0.25_dp
    distant = 0.5_dp

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, jacobian_product=jvp, inner_solver=inspect)
    walsh_error = sum(abs(observed * diagonal - 1.0_dp)) / real(n, dp)
    call check(jv_calls == 8 .and. solver_calls == 1, 'eight Walsh probes')
    call check(walsh_error > 0.2_dp, 'structured Walsh probes alias stride-eight coupling')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, jacobian_product=jvp, inner_solver=inspect)
    random_error = sum(abs(observed * diagonal - 1.0_dp)) / real(n, dp)
    call check(jv_calls == 8 .and. solver_calls == 1, 'eight seeded probes')
    call check(random_error < walsh_error, 'seeded signs reduce mean aliasing on this matrix')
    saved = observed

    call reset()
    res = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, jacobian_product=jvp, inner_solver=inspect)
    call check(all(saved == observed), 'seed repeatable through both public entry points')
    call check(jv_calls == 8, 'newton_krylov seeded probe accounting')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=159, jacobian_product=jvp, inner_solver=inspect)
    call check(any(saved /= observed), 'different seeds produce different probes')

    ! A single evenly spaced spot check catches the systematic Walsh
    ! aliasing at row one and triggers exact coordinate estimation.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_verify_samples=1, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 8 + 1 + n, 'detected alias triggers exact fallback')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'exact fallback restores true diagonal')
    call check(f_calls == 1 .and. res%nfev == 1, 'analytic checks have no residual calls')

    ! A configurable verification threshold trades extra exact work for
    ! sensitivity to smaller systematic biases.
    ! On a dense/long-range coupling, even a few spot checks may be
    ! lucky. The two-batch Rademacher consistency guard catches the high
    ! variance without consuming another Jv, then falls back to exact.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.15_dp, jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 8 + n, 'batch inconsistency triggers exact fallback without spot calls')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'batch inconsistency fallback recovers true diagonal')

    distant = 0.3_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.25_dp, jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 9, 'default 25 percent tolerance accepts moderate aliasing')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.10_dp, jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 9 + n, 'strict verification catches moderate aliasing')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'strict verification recovers exact diagonal')

    distant = 0.0_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_verify_samples=4, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 12, 'accurate estimate avoids exact fallback')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'verification preserves accurate approximation')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=3, &
        inner_m_diag_verify_rtol=0.1_dp, jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 11, 'consistent independent batches retain the cheap estimate')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'diagonal exact under either seeded probe batch')

    call reset()
    res = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=-123, inner_solver=inspect)
    call check(f_calls == 9 .and. res%nfev == 9, &
        'seeded finite differences consume one residual per probe')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        jacobian_diff_scheme='3-point', inner_solver=inspect)
    call check(f_calls == 17 .and. res%nfev == 17, &
        'centered seeded finite differences consume two residuals per probe')

    ! Ensure random sign generation never modifies Fortran's RNG state.
    call random_seed(size=ns)
    allocate(random_seed_values(ns))
    random_seed_values = 12345
    call random_seed(put=random_seed_values)
    call random_number(rand_before)
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=55, jacobian_product=jvp, inner_solver=inspect)
    call random_number(rand_after)
    call random_seed(put=random_seed_values)
    call random_number(expected_before)
    call random_number(expected_after)
    call check(rand_before == expected_before .and. rand_after == expected_after, &
        'local seeded probing does not change random_number state')

    bad_jvp = .true.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        jacobian_product=jvp, inner_solver=inspect)
    call check(.not. res%success .and. solver_calls == 0, 'failed probe aborts preconditioner')
    bad_jvp = .false.

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_probe_scheme='rademacher')
    call rejected(res, 'probe scheme requires auto')
    call reset()
    res = root(residual, x0, method='hybr', inner_m_diag_probe_scheme='rademacher')
    call rejected(res, 'only method krylov accepts probe scheme')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probe_scheme='unknown')
    call rejected(res, 'unknown scheme rejected')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probe_seed=5)
    call rejected(res, 'seed without rademacher rejected')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_verify_samples=-1)
    call rejected(res, 'negative verification count rejected')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_verify_samples=n+1)
    call rejected(res, 'oversized verification count rejected')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_verify_samples=1)
    call rejected(res, 'verification requires auto')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_verify_rtol=0.1_dp)
    call rejected(res, 'verification tolerance requires checked coordinates')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_verify_samples=1, inner_m_diag_verify_rtol=0.0_dp)
    call rejected(res, 'nonpositive verification tolerance rejected')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_verify_samples=1, inner_m_diag_verify_rtol=1.01_dp)
    call rejected(res, 'oversized verification tolerance rejected')

    if (failures /= 0) then
        print '(a,i0)', 'FAIL randomized diagonal probes: ', failures
        error stop 1
    end if
    print '(a)', 'PASS randomized diagonal probes and verification'
contains
    subroutine reset()
        f_calls = 0
        jv_calls = 0
        solver_calls = 0
    end subroutine reset

    subroutine check(condition, name)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: name
        if (.not. condition) then
            print '(a,a)', 'FAIL: ', name
            failures = failures + 1
        end if
    end subroutine check

    subroutine rejected(result, name)
        type(vector_root_result), intent(in) :: result
        character(len=*), intent(in) :: name
        call check(result%nfev == 0 .and. f_calls == 0 .and. &
            .not. result%success, name)
    end subroutine rejected

    subroutine apply_jacobian(v, product)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: product(:)
        integer :: k
        product = diagonal * v
        do k = 1, n - 8
            product(k) = product(k) + distant * diagonal(k) * v(k + 8)
        end do
    end subroutine apply_jacobian

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        f_calls = f_calls + 1
        call apply_jacobian(x, f)
    end subroutine residual

    subroutine jvp(x, v, product, ok)
        real(dp), intent(in) :: x(:), v(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        jv_calls = jv_calls + 1
        product = 0.0_dp
        ok = .not. bad_jvp
        if (ok) call apply_jacobian(v, product)
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec, preconditioner
        real(dp), intent(in) :: rhs(:), rtol
        integer, intent(in) :: maxiter
        real(dp), intent(out) :: solution(:)
        logical, intent(out) :: ok
        logical :: valid
        solver_calls = solver_calls + 1
        call preconditioner(spread(1.0_dp, 1, size(rhs)), observed, valid)
        call check(valid, 'finite preconditioner output')
        solution = 0.0_dp
        ok = .false.
    end subroutine inspect
end program test_optimize_root_krylov_random_diagonal
