! SPDX-License-Identifier: MIT
! Checked weighted rank-one correction for matrix-free Krylov roots.
program test_optimize_root_krylov_weighted_rank_one
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 9
    real(dp) :: a(n), u(n), v(n), target(n), probe(n), observed(n), expected(n)
    real(dp) :: x0(n), bad(n), sparse(n), checked_directions(n,3)
    type(vector_root_result) :: res
    integer :: i, failures, jv_calls, f_calls, solver_calls, fail_after
    logical :: extra_coupling

    failures = 0
    v = [1.5_dp, -0.5_dp, 2.0_dp, -1.2_dp, 0.3_dp, 0.0_dp, &
        1.7_dp, -0.75_dp, 0.8_dp]
    do i = 1, n
        a(i) = 0.5_dp + 0.25_dp * real(i, dp)
        u(i) = 0.15_dp + 0.025_dp * real(i, dp)
        target(i) = (-1.0_dp)**i * (0.2_dp + 0.03_dp * real(i, dp))
        probe(i) = 0.05_dp + 0.08_dp * real(i, dp)
    end do
    x0 = 0.0_dp
    extra_coupling = .false.
    expected = probe/a - (u/a)*(dot_product(v, probe/a) / &
        (1.0_dp + dot_product(v, u/a)))

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', &
        maxiter=2, fatol=1.0e-11_dp)
    call check(res%success .and. maxval(abs(res%x-target)) < 2.0e-12_dp, &
        'weighted correction solves linear coupled nonlinear equations')
    call check(maxval(abs(observed-expected)) < 1.0e-12_dp, &
        'weighted Sherman-Morrison inverse action')
    call check(jv_calls == n+3 .and. solver_calls == 1 .and. res%nfev == 2, &
        'analytic weighted fit uses n diagonal and three structural probes')

    call reset()
    res = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', &
        maxiter=2, f_tol=1.0e-11_dp)
    call check(res%success .and. maxval(abs(observed-expected)) < 1.0e-12_dp, &
        'newton_krylov forwards weighted right vector')

    call reset()
    res = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=1, inner_m_diag_auto=.true., inner_m_diag_rank_one=.true., &
        inner_m_diag_rank_one_right=v, jacobian_product=jvp, &
        line_search='none', maxiter=3, fatol=1.0e-11_dp)
    call check(res%success .and. maxval(abs(res%x-target)) < 1.0e-11_dp, &
        'weighted correction gives one-dimensional GMRES step')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=1.0e210_dp*v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 1.0e-12_dp, &
        'large finite weights are normalized before fitting')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=1.0e-210_dp*v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 1.0e-12_dp, &
        'small finite weights are normalized before fitting')

    ! Alternating right-factor signs match the first legacy verification
    ! pattern. Verification must still test two distinct new directions.
    block
        real(dp) :: original_v(n)
        original_v = v
        do i = 1, n
            v(i) = (-1.0_dp)**(i+1) * (0.5_dp + 0.1_dp*real(i,dp))
        end do
        expected = probe/a - (u/a)*(dot_product(v, probe/a) / &
            (1.0_dp + dot_product(v, u/a)))
        call reset()
        res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
            inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=v, &
            jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
        call check(maxval(abs(observed-expected)) < 1.0e-12_dp, &
            'alternating weights retain the correct inverse')
        call check(any(checked_directions(:,2) /= checked_directions(:,1)) .and. &
            any(checked_directions(:,2) /= -checked_directions(:,1)) .and. &
            any(checked_directions(:,3) /= checked_directions(:,1)) .and. &
            any(checked_directions(:,3) /= -checked_directions(:,1)) .and. &
            any(checked_directions(:,3) /= checked_directions(:,2)) .and. &
            any(checked_directions(:,3) /= -checked_directions(:,2)), &
            'two structural checks differ from fitting direction and each other')
        v = original_v
        expected = probe/a - (u/a)*(dot_product(v, probe/a) / &
            (1.0_dp + dot_product(v, u/a)))
    end block

    extra_coupling = .true.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/(a+u*v))) < 1.0e-13_dp, &
        'incompatible weighted coupling falls back to diagonal')
    extra_coupling = .false.

    ! A one-column coupling does not provide enough off-diagonal
    ! information to identify its diagonal row; retain the diagonal model.
    sparse = 0.0_dp
    sparse(4) = 1.0_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=sparse, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/(a+u*v))) < 1.0e-13_dp, &
        'degenerate supplied weight pattern retains diagonal scaling')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=v, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 2.0e-6_dp, &
        'two-point finite differences fit weighted rank-one model')
    call check(res%nfev == n+5 .and. f_calls == res%nfev, &
        'weighted forward differences count three additional products')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=v, &
        inner_solver=inspect, jacobian_diff_scheme='3-point', &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 2.0e-6_dp, &
        'centered finite differences fit weighted rank-one model')
    call check(res%nfev == 2*n+8 .and. f_calls == res%nfev, &
        'weighted centered differences count six additional evaluations')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_rank_one_right=v)
    call check(.not. res%success .and. res%nfev == 0, &
        'right weights require rank-one flag')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=v(:n-1))
    call check(.not. res%success .and. res%nfev == 0, &
        'wrong right-vector shape rejected before evaluation')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=0.0_dp*v)
    call check(.not. res%success .and. res%nfev == 0, &
        'zero right-vector weights rejected before evaluation')
    bad = v
    bad(1) = ieee_value(1.0_dp, ieee_quiet_nan)
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=bad)
    call check(.not. res%success .and. res%nfev == 0, &
        'nonfinite right-vector weights rejected before evaluation')
    res = root(residual, x0, method='hybr', inner_m_diag_rank_one_right=v)
    call check(.not. res%success .and. res%nfev == 0, &
        'right weights rejected for non-Krylov methods')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=v, &
        inner_method='minres')
    call check(.not. res%success .and. res%nfev == 0, &
        'weighted nonsymmetric inverse rejected for minres')

    call reset()
    fail_after = n + 1
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_one_right=v, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(.not. res%success .and. solver_calls == 0 .and. jv_calls == n+1, &
        'weighted coupling Jv callback failure aborts before inner solver')
    fail_after = 0

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_weighted_rank_one: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_weighted_rank_one: PASS'

contains

    subroutine reset()
        jv_calls = 0
        f_calls = 0
        solver_calls = 0
        fail_after = 0
        observed = 0.0_dp
        checked_directions = 0.0_dp
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
        real(dp), intent(out) :: f(:) !! weighted rank-one residual
        real(dp) :: delta(n)
        f_calls = f_calls + 1
        delta = x - target
        f = a*delta + u*dot_product(v, delta)
        if (extra_coupling) f(2) = f(2) + 0.4_dp*delta(4)
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(in) :: direction(:) !! Jacobian-vector direction
        real(dp), intent(out) :: product(:) !! exact Jacobian-vector product
        logical, intent(out) :: ok !! valid Jacobian-vector product
        jv_calls = jv_calls + 1
        if (jv_calls > n .and. jv_calls <= n+3) then
            checked_directions(:, jv_calls-n) = direction
        end if
        ok = .false.
        if (fail_after > 0 .and. jv_calls >= fail_after) return
        product = a*direction + u*dot_product(v, direction)
        if (extra_coupling) product(2) = product(2) + 0.4_dp*direction(4)
        ok = .true.
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: rhs(:) !! Newton right-hand side
        real(dp), intent(in) :: rtol !! inner tolerance
        integer, intent(in) :: maxiter !! number of inner iterations available
        real(dp), intent(out) :: solution(:) !! approximate Newton correction
        logical, intent(out) :: ok !! whether the inner solve succeeded
        solver_calls = solver_calls + 1
        call preconditioner(probe, observed, ok)
        if (.not. ok) return
        call preconditioner(rhs, solution, ok)
    end subroutine inspect
end program test_optimize_root_krylov_weighted_rank_one
