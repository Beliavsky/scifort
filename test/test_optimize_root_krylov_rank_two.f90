! SPDX-License-Identifier: MIT
! Checked rank-two Woodbury preconditioner for matrix-free Krylov roots.
program test_optimize_root_krylov_rank_two
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 9
    real(dp) :: a(n), u(n,2), v(n,2), target(n), probe(n), observed(n)
    real(dp) :: expected(n), diag(n), x0(n), bad(n,2), sparse(n,2)
    real(dp) :: saved_u(n,2), saved_v(n,2), work(n,2), small(2,2), projected(2)
    type(vector_root_result) :: res
    integer :: i, failures, jv_calls, f_calls, solver_calls, fail_after
    logical :: extra_coupling

    failures = 0
    v(:,1) = [1.5_dp, -0.5_dp, 2.0_dp, -1.2_dp, 0.3_dp, 0.0_dp, &
        1.7_dp, -0.75_dp, 0.8_dp]
    v(:,2) = [0.2_dp, 1.1_dp, -0.4_dp, 0.8_dp, -0.9_dp, 1.2_dp, &
        0.6_dp, -0.3_dp, 0.75_dp]
    do i = 1, n
        a(i) = 0.5_dp + 0.25_dp*real(i,dp)
        u(i,1) = 0.15_dp + 0.025_dp*real(i,dp)
        u(i,2) = 0.12_dp - 0.009_dp*real(i,dp)
        target(i) = (-1.0_dp)**i*(0.2_dp+0.03_dp*real(i,dp))
        probe(i) = 0.05_dp + 0.08_dp*real(i,dp)
    end do
    extra_coupling = .false.
    x0 = 0.0_dp
    call reference()

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', &
        maxiter=2, fatol=1.0e-11_dp)
    call check(res%success .and. maxval(abs(res%x-target)) < 2.0e-12_dp, &
        'rank-two inverse solves coupled linear roots')
    call check(maxval(abs(observed-expected)) < 2.0e-12_dp, &
        'rank-two Woodbury inverse action')
    call check(jv_calls == n+4 .and. solver_calls == 1 .and. res%nfev == 2, &
        'analytic model uses n diagonal and four structural products')

    call reset()
    res = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', &
        maxiter=2, f_tol=1.0e-11_dp)
    call check(res%success .and. maxval(abs(observed-expected)) < 2.0e-12_dp, &
        'newton_krylov forwards rank-two right matrix')

    call reset()
    res = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=1, inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_right=v, jacobian_product=jvp, &
        line_search='none', maxiter=3, fatol=1.0e-11_dp)
    call check(res%success .and. maxval(abs(res%x-target)) < 1.0e-11_dp, &
        'one-dimensional GMRES with the exact rank-two inverse')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=1.0e210_dp*v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 2.0e-12_dp, &
        'large finite weights are normalized independently')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=1.0e-210_dp*v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 2.0e-12_dp, &
        'small finite weights are normalized independently')

    extra_coupling = .true.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diag)) < 2.0e-13_dp, &
        'unmodeled coupling causes diagonal fallback')
    extra_coupling = .false.

    bad(:,1) = v(:,1)
    bad(:,2) = 2.0_dp*v(:,1)
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=bad, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diag)) < 2.0e-13_dp, &
        'dependent right columns safely use diagonal inverse')

    sparse = 0.0_dp
    sparse(3,1) = 1.0_dp
    sparse(4,2) = 1.0_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=sparse, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diag)) < 2.0e-13_dp, &
        'unidentifiable sparse columns safely use diagonal inverse')

    ! A very small or very large coupling diagonal should not cause the
    ! rank-two fit to use the floored inverse Jacobian diagonal instead of D.
    saved_u = u
    do i = 1, n
        a(i) = 10.0_dp**(-6.0_dp + 12.0_dp*real(i-1,dp)/real(n-1,dp))
        u(i,:) = [0.3_dp*a(i),-0.2_dp*a(i)]
    end do
    call reference()
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)/max(1.0_dp,abs(expected))) < 1.0e-9_dp, &
        'ill-scaled rank-two fit uses unfloored diagonal')
    do i = 1, n
        a(i) = 0.5_dp + 0.25_dp*real(i,dp)
    end do
    u = saved_u
    call reference()

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 2.0e-6_dp, &
        'forward differences fit the rank-two linear model')
    call check(res%nfev == n+6 .and. f_calls == res%nfev, &
        'two-point differences count four extra Jacobian probes')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_solver=inspect, jacobian_diff_scheme='3-point', &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 2.0e-6_dp, &
        'centered differences fit the rank-two linear model')
    call check(res%nfev == 2*n+10 .and. f_calls == res%nfev, &
        'centered differences count eight extra evaluations')

    res = root(residual, x0, method='krylov', inner_m_diag_rank_two=.true.)
    call check(.not. res%success .and. res%nfev == 0, &
        'rank-two option requires right matrix')
    res = root(residual, x0, method='krylov', inner_m_diag_rank_two_right=v)
    call check(.not. res%success .and. res%nfev == 0, &
        'right matrix requires enabled rank-two option')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v(:n-1,:))
    call check(.not. res%success .and. res%nfev == 0, 'right matrix rows validated')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v(:,1:1))
    call check(.not. res%success .and. res%nfev == 0, 'right matrix columns validated')
    bad = v
    bad(3,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=bad)
    call check(.not. res%success .and. res%nfev == 0, 'nonfinite weights rejected')
    bad = v
    bad(:,2) = 0.0_dp
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=bad)
    call check(.not. res%success .and. res%nfev == 0, 'zero weight column rejected')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_one=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_right=v)
    call check(.not. res%success .and. res%nfev == 0, 'rank options are exclusive')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_m_diag_probes=4)
    call check(.not. res%success .and. res%nfev == 0, 'exact diagonal required')
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_method='minres')
    call check(.not. res%success .and. res%nfev == 0, 'rank-two minres rejected')
    res = root(residual, x0, method='hybr', inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_right=v)
    call check(.not. res%success .and. res%nfev == 0, 'rank-two non-Krylov rejected')

    ! An exact singular Woodbury matrix must keep the ordinary diagonal
    ! preconditioner instead of attempting to invert the two-by-two factor.
    saved_u = u
    saved_v = v
    v(:,2) = v(:,2) - dot_product(v(:,1),v(:,2)) / &
        dot_product(v(:,1),v(:,1))*v(:,1)
    u(:,1) = -a*v(:,1)/dot_product(v(:,1),v(:,1))
    u(:,2) = 0.0_dp
    diag = a + u(:,1)*v(:,1)
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diag)) < 2.0e-12_dp, &
        'singular rank-two Woodbury system retains diagonal scaling')
    v = saved_v
    u = saved_u
    call reference()

    ! Four independent directions are needed for two fitting and two
    ! verification products. Smaller systems use the exact diagonal.
    call reset()
    res = root(residual_small, x0(:3), method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v(:3,:), &
        jacobian_product=jvp_small, inner_solver=inspect_small, &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed(:3)-probe(:3)/diag(:3))) < 2.0e-13_dp, &
        'small rank-two system uses exact diagonal preconditioner')
    call check(jv_calls == 3, 'small system avoids unnecessary probes')

    call reset()
    fail_after = n+1
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(.not. res%success .and. solver_calls == 0 .and. jv_calls == n+1, &
        'failed Jacobian callback stops before the inner solver')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_rank_two: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_rank_two: PASS'

contains

    subroutine reference()
        real(dp) :: determinant, coeff(2)
        work(:,1) = u(:,1)/a
        work(:,2) = u(:,2)/a
        small = matmul(transpose(v),work)
        small(1,1) = small(1,1)+1.0_dp
        small(2,2) = small(2,2)+1.0_dp
        determinant = small(1,1)*small(2,2)-small(1,2)*small(2,1)
        projected = matmul(transpose(v),probe/a)
        coeff(1) = (small(2,2)*projected(1)-small(1,2)*projected(2))/determinant
        coeff(2) = (small(1,1)*projected(2)-small(2,1)*projected(1))/determinant
        expected = probe/a - matmul(work,coeff)
        diag = a + u(:,1)*v(:,1) + u(:,2)*v(:,2)
    end subroutine reference

    subroutine reset()
        jv_calls = 0
        f_calls = 0
        solver_calls = 0
        fail_after = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition !! property to verify
        character(len=*), intent(in) :: description !! failure label
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine check

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! candidate nonlinear parameters
        real(dp), intent(out) :: f(:) !! model residuals
        real(dp) :: delta(n)
        f_calls = f_calls+1
        delta = x-target
        f = a*delta + matmul(u,matmul(transpose(v),delta))
        if (extra_coupling) f(3) = f(3) + 0.4_dp*delta(6)
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear parameters
        real(dp), intent(in) :: direction(:) !! Jacobian-vector direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! action was evaluated successfully
        jv_calls = jv_calls + 1
        ok = .false.
        if (fail_after > 0 .and. jv_calls >= fail_after) return
        product = a*direction + matmul(u,matmul(transpose(v),direction))
        if (extra_coupling) product(3) = product(3) + 0.4_dp*direction(6)
        ok = .true.
    end subroutine jvp

    subroutine residual_small(x, f)
        real(dp), intent(in) :: x(:) !! small-system parameters
        real(dp), intent(out) :: f(:) !! small-system residuals
        real(dp) :: delta(3)
        f_calls = f_calls+1
        delta = x-target(:3)
        f = a(:3)*delta + matmul(u(:3,:),matmul(transpose(v(:3,:)),delta))
    end subroutine residual_small

    subroutine jvp_small(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! small-system point
        real(dp), intent(in) :: direction(:) !! small Jacobian direction
        real(dp), intent(out) :: product(:) !! small Jacobian product
        logical, intent(out) :: ok !! whether the product is finite
        jv_calls = jv_calls+1
        product = a(:3)*direction + matmul(u(:3,:), &
            matmul(transpose(v(:3,:)),direction))
        ok = .true.
    end subroutine jvp_small

    subroutine inspect_small(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! matrix-free small Jacobian
        procedure(vector_root_linear_operator) :: preconditioner !! approximate small inverse
        real(dp), intent(in) :: rhs(:) !! Newton right-hand side
        real(dp), intent(in) :: rtol !! relative linear tolerance
        integer, intent(in) :: maxiter !! available Krylov iterations
        real(dp), intent(out) :: solution(:) !! computed correction
        logical, intent(out) :: ok !! preconditioner succeeded
        solver_calls = solver_calls+1
        call preconditioner(probe(:3),observed(:3),ok)
        if (.not. ok) return
        call preconditioner(rhs,solution,ok)
    end subroutine inspect_small

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! matrix-free Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! current inverse approximation
        real(dp), intent(in) :: rhs(:) !! Newton right-hand side
        real(dp), intent(in) :: rtol !! linear solve tolerance
        integer, intent(in) :: maxiter !! allowed inner iterations
        real(dp), intent(out) :: solution(:) !! Newton correction
        logical, intent(out) :: ok !! preconditioner application succeeded
        solver_calls = solver_calls+1
        call preconditioner(probe,observed,ok)
        if (.not. ok) return
        call preconditioner(rhs,solution,ok)
    end subroutine inspect
end program test_optimize_root_krylov_rank_two
