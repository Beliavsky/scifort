! SPDX-License-Identifier: MIT
! Deterministic Walsh probing for automatic Krylov diagonal preconditioning.
program test_optimize_root_krylov_diagonal_probes
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 64
    real(dp) :: diagonal(n), target(n), x0(n), observed(n)
    type(vector_root_result) :: res
    integer :: i, f_calls, jv_calls, solver_calls, failures
    logical :: fail_jvp
    real(dp) :: neighbor, distant

    failures = 0
    fail_jvp = .false.
    neighbor = 0.2_dp
    distant = 0.0_dp
    do i = 1, n
        diagonal(i) = 10.0_dp**(-3.0_dp * real(i - 1, dp) / real(n - 1, dp))
        target(i) = sin(real(i, dp) / 9.0_dp)
    end do
    x0 = target + 0.1_dp

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=4, jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(.not. res%success .and. solver_calls == 1, 'inspect probe diagonal')
    call check(jv_calls == 4, 'four probes cost four Jv products')
    call check(res%nfev == 1 .and. f_calls == 1, 'analytic probing needs no residual evaluations')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 3.0e-14_dp, &
        'four Walsh probes recover the tridiagonal diagonal exactly')

    call reset()
    res = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(solver_calls == 1 .and. jv_calls == 8, 'newton_krylov forwards eight probes')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 3.0e-14_dp, &
        'newton_krylov estimates the expected diagonal')

    distant = 0.3_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=4, jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(maxval(abs(observed * diagonal - 1.0_dp)) > 0.1_dp, &
        'four probes expose deterministic aliasing of stride-four coupling')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 3.0e-14_dp, &
        'eight Walsh probes resolve stride-four coupling')

    ! The historical exact path is preserved, including when the requested
    ! probe count is at least the system dimension.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(jv_calls == n, 'omitting probe option retains n exact Jv products')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 3.0e-14_dp, 'exact baseline')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=n, jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(jv_calls == n, 'n probes selects exact coordinate estimates')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 3.0e-14_dp, 'exact count baseline')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=0, jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(jv_calls == n, 'explicit zero probes selects exact coordinate estimation')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=n+4, jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(jv_calls == n, 'probe count exceeding dimension uses exact method')

    distant = 0.0_dp
    neighbor = 0.0_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=1, jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(jv_calls == 1, 'single dense probe is allowed')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 3.0e-14_dp, &
        'single probe recovers a diagonal Jacobian')
    neighbor = 0.2_dp

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_solver=inspect_preconditioner)
    call check(res%nfev == 9 .and. f_calls == 9, &
        'forward finite differences consume eight residual evaluations')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 1.0e-6_dp, &
        'forward finite-difference diagonal estimate')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, jacobian_diff_scheme='3-point', &
        inner_solver=inspect_preconditioner)
    call check(res%nfev == 17 .and. f_calls == 17, &
        'centered finite differences consume sixteen residual evaluations')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 1.0e-7_dp, &
        'centered finite-difference diagonal estimate')

    call reset()
    res = newton_krylov(residual, x0, iter=3, line_search='none', &
        inner_m_diag_auto=.true., inner_m_diag_probes=8, &
        inner_m_diag_update_interval=1, jacobian_product=jvp, &
        inner_solver=diagonal_step)
    call check(jv_calls == 24 .and. solver_calls == 3, 'three refreshes cost 24 probes')
    call reset()
    res = newton_krylov(residual, x0, iter=3, line_search='none', &
        inner_m_diag_auto=.true., inner_m_diag_probes=8, &
        inner_m_diag_update_interval=0, jacobian_product=jvp, &
        inner_solver=diagonal_step)
    call check(jv_calls == 8 .and. solver_calls == 3, 'initial-only refresh costs eight probes')

    fail_jvp = .true.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, jacobian_product=jvp, inner_solver=inspect_preconditioner)
    call check(.not. res%success .and. solver_calls == 0 .and. res%nfev == 1, &
        'failed probing aborts before the inner solve')
    fail_jvp = .false.

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_probes=8)
    call rejected(res, 'probe count requires auto preconditioning')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.false., &
        inner_m_diag_probes=8)
    call rejected(res, 'probe count rejects explicit false auto preconditioning')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=-1)
    call rejected(res, 'negative probe count rejected')
    call reset()
    res = root(residual, x0, method='hybr', inner_m_diag_probes=8)
    call rejected(res, 'probe count accepted only by krylov')

    if (failures /= 0) then
        print '(a,i0)', 'FAIL diagonal probe tests: ', failures
        error stop 1
    end if
    print '(a)', 'PASS automatic Krylov diagonal Walsh probing'

contains
    subroutine reset()
        f_calls = 0
        jv_calls = 0
        solver_calls = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, name)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: name
        if (.not. condition) then
            print '(a,a)', 'FAIL: ', name
            failures = failures + 1
        end if
    end subroutine check

    subroutine rejected(value, name)
        type(vector_root_result), intent(in) :: value
        character(len=*), intent(in) :: name
        call check(.not. value%success .and. value%nfev == 0 .and. f_calls == 0, name)
    end subroutine rejected

    subroutine apply_jacobian(v, jv)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: jv(:)
        integer :: k, nn
        nn = size(v)
        jv = diagonal(1:nn) * v
        do k = 1, nn-1
            jv(k) = jv(k) + neighbor * diagonal(k) * v(k+1)
        end do
        do k = 1, nn-4
            jv(k) = jv(k) + distant * diagonal(k) * v(k+4)
        end do
    end subroutine apply_jacobian

    subroutine residual(x, fx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx(:)
        f_calls = f_calls + 1
        call apply_jacobian(x - target(1:size(x)), fx)
    end subroutine residual

    subroutine jvp(x, v, jv, ok)
        real(dp), intent(in) :: x(:), v(:)
        real(dp), intent(out) :: jv(:)
        logical, intent(out) :: ok
        jv_calls = jv_calls + 1
        ok = .not. fail_jvp
        jv = 0.0_dp
        if (.not. ok) return
        call apply_jacobian(v, jv)
    end subroutine jvp

    subroutine inspect_preconditioner(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec, preconditioner
        real(dp), intent(in) :: rhs(:), rtol
        integer, intent(in) :: maxiter
        real(dp), intent(out) :: solution(:)
        logical, intent(out) :: ok
        logical :: operation_ok
        solver_calls = solver_calls + 1
        call preconditioner(spread(1.0_dp, 1, size(rhs)), observed, operation_ok)
        solution = 0.0_dp
        ok = .false.
        call check(operation_ok, 'valid preconditioner output')
    end subroutine inspect_preconditioner

    subroutine diagonal_step(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec, preconditioner
        real(dp), intent(in) :: rhs(:), rtol
        integer, intent(in) :: maxiter
        real(dp), intent(out) :: solution(:)
        logical, intent(out) :: ok
        solver_calls = solver_calls + 1
        call preconditioner(rhs, solution, ok)
    end subroutine diagonal_step
end program test_optimize_root_krylov_diagonal_probes
