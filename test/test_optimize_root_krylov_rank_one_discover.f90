! SPDX-License-Identifier: MIT
! Recover an unknown weighted rank-one right factor from exact coordinate Jv products.
program test_optimize_root_krylov_rank_one_discover
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 12
    real(dp) :: a(n), u(n), v(n), target(n), probe(n), x0(n), observed(n)
    real(dp) :: expected(n), diagonal(n), original_a(n), original_u(n)
    integer :: i, failures, f_calls, jv_calls, solver_calls, fail_after
    logical :: incompatible, sparse_coupling, singular_model
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        a(i) = 0.8_dp + 0.18_dp*real(i,dp)
        u(i) = 0.08_dp + 0.018_dp*real(i,dp)
        v(i) = cos(0.49_dp*real(i,dp))
        target(i) = 0.06_dp*sin(0.6_dp*real(i,dp))
        probe(i) = 0.1_dp + 0.039_dp*real(i,dp)
    end do
    v(4) = 0.0_dp
    v(9) = 0.0_dp
    x0 = 0.0_dp
    incompatible = .false.
    sparse_coupling = .false.
    singular_model = .false.
    call reference()

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, &
        line_search='none', maxiter=2, fatol=1.0e-11_dp)
    call check(fit%success .and. maxval(abs(fit%x-target)) < 4.0e-12_dp, &
        'unknown weighted rank-one model solves a linear system')
    call check(maxval(abs(observed-expected)) < 3.0e-12_dp, &
        'discovered right factor yields exact Sherman-Morrison action')
    call check(jv_calls == n .and. solver_calls == 1 .and. fit%nfev == 2, &
        'all coordinate Jacobian products are reused without extra probes')

    call reset()
    fit = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, &
        line_search='none', maxiter=2, f_tol=1.0e-11_dp)
    call check(fit%success .and. maxval(abs(observed-expected)) < 3.0e-12_dp, &
        'newton_krylov forwards rank-one discovery')

    call reset()
    fit = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=1, inner_m_diag_auto=.true., inner_m_diag_rank_one=.true., &
        inner_m_diag_rank_one_discover=.true., jacobian_product=jvp, &
        line_search='none', maxiter=3, fatol=1.0e-11_dp)
    call check(fit%success .and. maxval(abs(fit%x-target)) < 5.0e-11_dp, &
        'discovered inverse supports one-dimensional GMRES convergence')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 5.0e-6_dp, &
        'two-point finite differences identify the unknown right factor')
    call check(fit%nfev == n+2 .and. fit%nfev == f_calls, &
        'two-point discovery uses n finite-difference products')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        inner_solver=inspect, jacobian_diff_scheme='3-point', &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 5.0e-6_dp, &
        'three-point finite differences identify unknown right factors')
    call check(fit%nfev == 2*n+2 .and. fit%nfev == f_calls, &
        'three-point discovery uses 2n residual perturbations')

    original_a = a
    original_u = u
    do i = 1, n
        a(i) = 10.0_dp**(-4.0_dp + 8.0_dp*real(i-1,dp)/real(n-1,dp))
        u(i) = 0.06_dp*a(i)
    end do
    call reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)/max(1.0_dp,abs(expected))) < 2.0e-9_dp, &
        'discovered inverse handles eight orders of diagonal scaling')
    a = original_a
    u = original_u
    call reference()

    incompatible = .true.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 3.0e-13_dp, &
        'incompatible rank-two coupling falls back to exact diagonal')
    call check(jv_calls == n .and. solver_calls == 1, &
        'rank-two fallback needs no repeated Jacobian products')
    incompatible = .false.

    sparse_coupling = .true.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/a)) < 3.0e-13_dp, &
        'nonidentifiable sparse rank-one structure retains diagonal inverse')
    sparse_coupling = .false.

    singular_model = .true.
    diagonal = a*(1.0_dp-v*v/sum(v*v))
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 3.0e-13_dp, &
        'singular discovered Sherman-Morrison system retains diagonal inverse')
    singular_model = .false.
    call reference()

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_rank_one_discover=.true.)
    call check(.not. fit%success .and. fit%nfev == 0 .and. &
        index(fit%message,'requires inner_m_diag_rank_one') > 0, &
        'discovery requires rank-one option')
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        inner_m_diag_rank_one_right=v)
    call check(.not. fit%success .and. fit%nfev == 0 .and. &
        index(fit%message,'conflict') > 0, &
        'cannot supply and discover the same right vector')
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        inner_m_diag_probes=3)
    call check(.not. fit%success .and. fit%nfev == 0, &
        'discovery requires exact coordinate probes')
    fit = root(residual, x0, method='hybr', inner_m_diag_rank_one_discover=.true.)
    call check(.not. fit%success .and. fit%nfev == 0, &
        'non-Krylov methods reject discovery options')

    call reset()
    fail_after = 5
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect)
    call check(.not. fit%success .and. solver_calls == 0 .and. jv_calls == 5, &
        'Jacobian callback failure aborts before the Krylov inner solver')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_rank_one_discover: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_rank_one_discover: PASS'

contains

    subroutine reset()
        jv_calls = 0
        f_calls = 0
        solver_calls = 0
        fail_after = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, title)
        logical, intent(in) :: condition !! expected condition
        character(len=*), intent(in) :: title !! diagnostic description
        if (.not. condition) then
            failures = failures+1
            print '(a)', 'FAIL: '//title
        end if
    end subroutine check

    subroutine reference()
        real(dp) :: w(n), denom
        w = u/a
        denom = 1.0_dp + sum(v*w)
        expected = probe/a - w*(sum(v*probe/a)/denom)
        diagonal = a + u*v
    end subroutine reference

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(out) :: f(:) !! weighted coupled residual
        real(dp) :: delta(n), weights(n), coupling(n)
        f_calls = f_calls + 1
        delta = x-target
        weights = v
        coupling = u
        if (sparse_coupling) then
            weights = 0.0_dp
            weights(2) = 1.0_dp
            coupling = 0.0_dp
            coupling(1) = 0.4_dp
        end if
        if (singular_model) coupling = -a*v/sum(v*v)
        f = a*delta + coupling*sum(weights*delta)
        if (incompatible) f(2) = f(2) + 0.3_dp*delta(5)
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current point, unused for linear residual
        real(dp), intent(in) :: direction(:) !! Jacobian trial direction
        real(dp), intent(out) :: product(:) !! Jacobian times direction
        logical, intent(out) :: ok !! false on injected failure
        real(dp) :: weights(n), coupling(n)
        jv_calls = jv_calls + 1
        ok = .false.
        if (fail_after > 0 .and. jv_calls >= fail_after) return
        weights = v
        coupling = u
        if (sparse_coupling) then
            weights = 0.0_dp
            weights(2) = 1.0_dp
            coupling = 0.0_dp
            coupling(1) = 0.4_dp
        end if
        if (singular_model) coupling = -a*v/sum(v*v)
        product = a*direction + coupling*sum(weights*direction)
        if (incompatible) product(2) = product(2) + 0.3_dp*direction(5)
        ok = .true.
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! current Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: rhs(:) !! Newton right-hand side
        real(dp), intent(in) :: rtol !! requested inner tolerance
        integer, intent(in) :: maxiter !! available Krylov basis budget
        real(dp), intent(out) :: solution(:) !! preconditioned right-hand side
        logical, intent(out) :: ok !! false on inverse failure

        solver_calls = solver_calls + 1
        call preconditioner(probe, observed, ok)
        if (.not. ok) return
        call preconditioner(rhs, solution, ok)
    end subroutine inspect
end program test_optimize_root_krylov_rank_one_discover
