! SPDX-License-Identifier: MIT
! Checked rank-one mean-field preconditioning for matrix-free Krylov roots.
program test_optimize_root_krylov_rank_one
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 12
    real(dp) :: a(n), u(n), target(n), probe(n), observed(n), expected(n), x0(n)
    type(vector_root_result) :: res
    integer :: i, failures, jv_calls, f_calls, solver_calls, fail_after
    logical :: nonsymmetric_extra, singular_base, singular_model

    failures = 0
    do i = 1, n
        a(i) = 0.2_dp + 0.4_dp * real(i, dp)
        u(i) = 0.7_dp + 0.1_dp * real(i, dp)
        target(i) = (-1.0_dp)**i * (0.2_dp + 0.03_dp * real(i, dp))
        probe(i) = 0.05_dp + 0.08_dp * real(i, dp)
    end do
    x0 = 0.0_dp
    nonsymmetric_extra = .false.
    singular_base = .false.
    singular_model = .false.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=2, fatol=1.0e-11_dp)
    expected = probe / a - (u / a) * (sum(probe / a) / (1.0_dp + sum(u / a)))
    call check(res%success .and. maxval(abs(res%x - target)) < 5.0e-12_dp, &
        'rank-one model yields one-step nonlinear solution')
    call check(maxval(abs(observed - expected)) < 1.0e-12_dp, &
        'rank-one preconditioner matches Sherman-Morrison inverse')
    call check(jv_calls == n + 3 .and. solver_calls == 1 .and. res%nfev == 2, &
        'analytic Jacobian uses n coordinate and three structural probes')

    call reset()
    res = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., jacobian_product=jvp, inner_solver=inspect, &
        line_search='none', maxiter=2, f_tol=1.0e-11_dp)
    call check(res%success .and. maxval(abs(observed - expected)) < 1.0e-12_dp, &
        'newton_krylov forwards the checked rank-one option')

    call reset()
    res = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=1, inner_m_diag_auto=.true., inner_m_diag_rank_one=.true., &
        jacobian_product=jvp, line_search='none', maxiter=3, fatol=1.0e-11_dp)
    call check(res%success .and. maxval(abs(res%x - target)) < 1.0e-11_dp, &
        'rank-one correction gives a one-dimensional GMRES Newton step')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed - probe/(a+u))) < 1.0e-13_dp, &
        'omitting option retains original diagonal preconditioner')
    call check(jv_calls == n, 'disabled rank-one option avoids structural probes')

    nonsymmetric_extra = .true.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed - probe/(a+u))) < 1.0e-13_dp, &
        'coupling inconsistent with row-uniform model falls back to diagonal')
    call check(jv_calls >= n + 2 .and. solver_calls == 1, &
        'rejected model still reaches inner solver safely')
    nonsymmetric_extra = .false.

    ! The raw Jacobian diagonal is essential in a strongly ill-scaled
    ! system: the standard inverse diagonal floors tiny entries, while
    ! the row-coupled inverse uses the actual uncoupled entries.
    block
        real(dp) :: original_a(n), original_u(n)
        original_a = a
        original_u = u
        do i = 1, n
            a(i) = 10.0_dp**(-6.0_dp + 12.0_dp * real(i-1, dp) / real(n-1, dp))
            u(i) = 0.3_dp * a(i)
        end do
        call reset()
        res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
            inner_m_diag_rank_one=.true., jacobian_product=jvp, &
            inner_solver=inspect, line_search='none', maxiter=1)
        expected = probe/a - (u/a)*(sum(probe/a)/(1.0_dp+sum(u/a)))
        call check(maxval(abs(observed-expected)/max(1.0_dp,abs(expected))) < 1.0e-10_dp, &
            'row-coupled action uses unfloored diagonal on ill-scaled systems')
        a = original_a
        u = original_u
    end block

    singular_base = .true.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    expected = probe / (a + u)
    expected(3) = probe(3) / u(3)
    call check(maxval(abs(observed - expected)) < 1.0e-13_dp, &
        'singular fitted diagonal safely uses original diagonal inverse')
    singular_base = .false.

    singular_model = .true.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    expected = probe / (a * (1.0_dp - 1.0_dp / real(n, dp)))
    call check(maxval(abs(observed - expected)) < 1.0e-12_dp, &
        'singular Sherman-Morrison denominator safely uses diagonal inverse')
    singular_model = .false.

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_solver=inspect, &
        line_search='none', maxiter=1, fatol=1.0e-8_dp)
    expected = probe/a - (u/a)*(sum(probe/a)/(1.0_dp+sum(u/a)))
    call check(maxval(abs(observed - expected)) < 1.0e-6_dp, &
        'forward differences recognize linear row-coupled Jacobian')
    call check(res%nfev == n + 5 .and. f_calls == res%nfev, &
        'two-point finite differences count three additional products')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_solver=inspect, &
        jacobian_diff_scheme='3-point', line_search='none', maxiter=1, fatol=1.0e-8_dp)
    call check(maxval(abs(observed - expected)) < 1.0e-6_dp, &
        'centered differences recognize linear row-coupled Jacobian')
    call check(res%nfev == 2*n + 8 .and. f_calls == res%nfev, &
        'three-point finite differences count three additional products')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_rank_one=.true.)
    call check(.not. res%success .and. res%nfev == 0, &
        'rank-one option requires automatic diagonal estimation')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_probes=4)
    call check(.not. res%success .and. res%nfev == 0, &
        'rank-one option requires an exact unclipped Jacobian diagonal')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_method='minres')
    call check(.not. res%success .and. res%nfev == 0, &
        'nonsymmetric Woodbury preconditioner rejected for minres')
    res = root(residual, x0, method='hybr', inner_m_diag_rank_one=.true.)
    call check(.not. res%success .and. res%nfev == 0, &
        'rank-one option rejected for non-Krylov methods')

    ! In dimension one the rank-one model is unnecessary; use exact Jii.
    call reset()
    res = root(scalar_residual, [0.0_dp], method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_rank_one=.true., &
        jacobian_product=scalar_jvp, inner_method='gmres', &
        line_search='none', maxiter=3, fatol=1.0e-12_dp)
    call check(res%success .and. abs(res%x(1)-2.0_dp) < 1.0e-12_dp, &
        'dimension one gracefully keeps exact diagonal preconditioning')

    call reset()
    fail_after = n + 1
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., jacobian_product=jvp, inner_solver=inspect)
    call check(.not. res%success .and. solver_calls == 0 .and. jv_calls == n + 1, &
        'failure in row-coupling product aborts before inner solver')
    fail_after = 0

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_rank_one: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_rank_one: PASS'

contains

    subroutine reset()
        jv_calls = 0
        f_calls = 0
        solver_calls = 0
        fail_after = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, title)
        logical, intent(in) :: condition !! expected property
        character(len=*), intent(in) :: title !! diagnostic text

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//title
        end if
    end subroutine check

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! current nonlinear variables
        real(dp), intent(out) :: f(:) !! mean-field residual at x
        real(dp) :: delta(n), base(n), coupling(n)

        f_calls = f_calls + 1
        delta = x - target
        base = a
        coupling = u
        if (singular_base) base(3) = 0.0_dp
        if (singular_model) coupling = -base / real(n, dp)
        f = base * delta + coupling * sum(delta)
        if (nonsymmetric_extra) f(2) = f(2) + 0.35_dp * delta(4)
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear point (constant Jacobian here)
        real(dp), intent(in) :: direction(:) !! Jacobian multiplication direction
        real(dp), intent(out) :: product(:) !! exact Jacobian-vector product
        logical, intent(out) :: ok !! true for valid product
        real(dp) :: base(n), coupling(n)

        jv_calls = jv_calls + 1
        ok = .false.
        if (fail_after > 0 .and. jv_calls >= fail_after) return
        base = a
        coupling = u
        if (singular_base) base(3) = 0.0_dp
        if (singular_model) coupling = -base / real(n, dp)
        product = base * direction + coupling * sum(direction)
        if (nonsymmetric_extra) product(2) = product(2) + 0.35_dp * direction(4)
        ok = .true.
    end subroutine jvp

    subroutine scalar_residual(x, f)
        real(dp), intent(in) :: x(:) !! scalar variable
        real(dp), intent(out) :: f(:) !! scalar residual
        f(1) = 3.0_dp * (x(1)-2.0_dp)
    end subroutine scalar_residual

    subroutine scalar_jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! scalar point
        real(dp), intent(in) :: direction(:) !! scalar direction
        real(dp), intent(out) :: product(:) !! scalar directional product
        logical, intent(out) :: ok !! product valid
        product(1) = 3.0_dp * direction(1)
        ok = .true.
    end subroutine scalar_jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: rhs(:) !! Newton right-hand side
        real(dp), intent(in) :: rtol !! relative inner tolerance
        integer, intent(in) :: maxiter !! available Krylov iterations
        real(dp), intent(out) :: solution(:) !! preconditioned right-hand side
        logical, intent(out) :: ok !! false if preconditioner rejects input

        solver_calls = solver_calls + 1
        call preconditioner(probe, observed, ok)
        if (.not. ok) return
        call preconditioner(rhs, solution, ok)
    end subroutine inspect

end program test_optimize_root_krylov_rank_one
