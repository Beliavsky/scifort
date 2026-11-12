! SPDX-License-Identifier: MIT
! Opt-in three-probe diagonal recovery for a checked rank-two Woodbury inverse.
program test_optimize_root_krylov_rank_two_sketch
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 12
    real(dp) :: a(n), u(n,2), v(n,2), target(n), probe(n), expected(n)
    real(dp) :: x0(n), observed(n), diagonal(n), bad(n,2)
    real(dp) :: w(n,2), small(2,2), projected(2), coefficients(2)
    real(dp) :: determinant
    real(dp) :: saved_a(n), saved_u(n,2)
    integer :: i, failures, f_calls, jv_calls, solver_calls, fail_after
    logical :: incompatible
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        a(i) = 0.5_dp + 0.25_dp*real(i,dp)
        u(i,1) = 0.12_dp + 0.022_dp*real(i,dp)
        u(i,2) = 0.19_dp - 0.015_dp*real(i,dp)
        v(i,1) = 0.4_dp + sin(0.37_dp*real(i,dp))
        v(i,2) = cos(0.51_dp*real(i,dp))
        probe(i) = 0.1_dp + 0.031_dp*real(i,dp)
        target(i) = 0.14_dp*cos(real(i,dp))
    end do
    x0 = 0.0_dp
    incompatible = .false.
    call reference()

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_m_diag_rank_two_sketch=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=2, fatol=1.0e-11_dp)
    call check(fit%success .and. maxval(abs(fit%x-target)) < 1.0e-12_dp, &
        'sketched rank-two preconditioner converges on structured system')
    call check(maxval(abs(observed-expected)) < 2.0e-12_dp, &
        'sketched model reconstructs Woodbury inverse action')
    call check(jv_calls == 7 .and. solver_calls == 1 .and. fit%nfev == f_calls, &
        'three fit and four structural Jv products replace n+4 products')

    call reset()
    fit = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_m_diag_rank_two_sketch=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=2, f_tol=1.0e-11_dp)
    call check(fit%success .and. jv_calls == 7 .and. &
        maxval(abs(observed-expected)) < 2.0e-12_dp, &
        'newton_krylov forwards the sketch option')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_m_diag_rank_two_sketch=.true., inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 1.0e-6_dp, &
        'forward-difference sketch reproduces inverse action')
    call check(fit%nfev == 9 .and. f_calls == fit%nfev, &
        'seven forward Jv probes cost seven residual evaluations')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_m_diag_rank_two_sketch=.true., inner_solver=inspect, &
        jacobian_diff_scheme='3-point', line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 1.0e-6_dp, &
        'centered-difference sketch reproduces inverse action')
    call check(fit%nfev == 16 .and. f_calls == fit%nfev, &
        'seven centered Jv probes cost fourteen residual evaluations')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=1.0e210_dp*v, &
        inner_m_diag_rank_two_sketch=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jv_calls == 7 .and. maxval(abs(observed-expected)) < 2.0e-12_dp, &
        'large finite right weights are normalized for the sketch')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=1.0e-210_dp*v, &
        inner_m_diag_rank_two_sketch=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jv_calls == 7 .and. maxval(abs(observed-expected)) < 2.0e-12_dp, &
        'small finite right weights are normalized for the sketch')

    saved_a = a
    saved_u = u
    do i = 1, n
        a(i) = 10.0_dp**(-5.0_dp + 10.0_dp*real(i-1,dp)/real(n-1,dp))
        u(i,1) = 0.25_dp*a(i)
        u(i,2) = -0.18_dp*a(i)
    end do
    call reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_m_diag_rank_two_sketch=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jv_calls == 7 .and. &
        maxval(abs(observed-expected)/max(1.0_dp,abs(expected))) < 1.0e-9_dp, &
        'sketch retains correct inverse over ten orders of diagonal scale')
    a = saved_a
    u = saved_u
    call reference()

    ! A rank-two model cannot represent an extra off-diagonal coupling.
    ! Rejected sketch estimates must recover the exact diagonal fallback.
    incompatible = .true.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_m_diag_rank_two_sketch=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 3.0e-13_dp, &
        'incompatible coupling falls back to exact diagonal')
    call check(jv_calls >= n .and. solver_calls == 1, &
        'incompatible sketch invokes exact diagonal recovery')
    incompatible = .false.

    bad(:,1) = v(:,1)
    bad(:,2) = 2.0_dp*v(:,1)
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=bad, &
        inner_m_diag_rank_two_sketch=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 3.0e-13_dp, &
        'dependent right columns use the exact diagonal fallback')
    call check(jv_calls >= n, 'degenerate right factors do not enable a sketch')

    ! Fewer than five variables cannot support three independent fitting
    ! directions and the later two independent model checks.
    call reset()
    fit = root(residual_small, x0(:3), method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_right=v(:3,:), &
        inner_m_diag_rank_two_sketch=.true., jacobian_product=jvp_small, &
        inner_solver=inspect_small, line_search='none', maxiter=1)
    call check(jv_calls == 3 .and. &
        maxval(abs(observed(:3)-probe(:3)/diagonal(:3))) < 2.0e-13_dp, &
        'small systems skip sketch products and use exact diagonal')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_rank_two_sketch=.true.)
    call check(.not. fit%success .and. fit%nfev == 0 .and. &
        index(fit%message,'requires inner_m_diag_rank_two') > 0, &
        'sketch requires enabled rank-two option')
    call reset()
    fit = root(residual, x0, method='hybr', inner_m_diag_rank_two_sketch=.true.)
    call check(.not. fit%success .and. fit%nfev == 0, &
        'sketch cannot be requested for non-Krylov methods')

    call reset()
    fail_after = 2
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=v, &
        inner_m_diag_rank_two_sketch=.true., jacobian_product=jvp, &
        inner_solver=inspect)
    call check(.not. fit%success .and. jv_calls == 2 .and. solver_calls == 0, &
        'failed sketch Jacobian callback aborts before applying preconditioner')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_rank_two_sketch: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_rank_two_sketch: PASS'

contains
    subroutine reference()
        w(:,1) = u(:,1)/a
        w(:,2) = u(:,2)/a
        small = matmul(transpose(v),w)
        small(1,1) = small(1,1)+1.0_dp
        small(2,2) = small(2,2)+1.0_dp
        determinant = small(1,1)*small(2,2)-small(1,2)*small(2,1)
        projected = matmul(transpose(v),probe/a)
        coefficients(1) = (small(2,2)*projected(1)-small(1,2)*projected(2))/determinant
        coefficients(2) = (small(1,1)*projected(2)-small(2,1)*projected(1))/determinant
        expected = probe/a-matmul(w,coefficients)
        diagonal = a+u(:,1)*v(:,1)+u(:,2)*v(:,2)
    end subroutine reference

    subroutine reset()
        f_calls = 0
        jv_calls = 0
        solver_calls = 0
        fail_after = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition !! assertion result
        character(len=*), intent(in) :: description !! failed-test description
        if (.not. condition) then
            failures = failures+1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine check

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! current nonlinear point
        real(dp), intent(out) :: f(:) !! coupled linear residuals
        real(dp) :: delta(n)
        f_calls = f_calls+1
        delta = x-target
        f = a*delta+matmul(u,matmul(transpose(v),delta))
        if (incompatible) f(3) = f(3)+0.4_dp*delta(6)
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! Jacobian linearization point
        real(dp), intent(in) :: direction(:) !! Jacobian-vector input
        real(dp), intent(out) :: product(:) !! exact Jacobian-vector output
        logical, intent(out) :: ok !! callback completed successfully
        jv_calls = jv_calls+1
        ok = .false.
        if (fail_after > 0 .and. jv_calls >= fail_after) return
        product = a*direction+matmul(u,matmul(transpose(v),direction))
        if (incompatible) product(3) = product(3)+0.4_dp*direction(6)
        ok = .true.
    end subroutine jvp

    subroutine residual_small(x, f)
        real(dp), intent(in) :: x(:) !! short vector iterate
        real(dp), intent(out) :: f(:) !! first three coupled residuals
        real(dp) :: delta(3)
        f_calls = f_calls+1
        delta = x-target(:3)
        f = a(:3)*delta+matmul(u(:3,:),matmul(transpose(v(:3,:)),delta))
    end subroutine residual_small

    subroutine jvp_small(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! short Jacobian point
        real(dp), intent(in) :: direction(:) !! three-element direction
        real(dp), intent(out) :: product(:) !! short Jacobian product
        logical, intent(out) :: ok !! true on successful evaluation
        jv_calls = jv_calls+1
        product = a(:3)*direction+matmul(u(:3,:), &
            matmul(transpose(v(:3,:)),direction))
        ok = .true.
    end subroutine jvp_small

    subroutine inspect_small(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! short Jacobian operator
        procedure(vector_root_linear_operator) :: preconditioner !! short inverse action
        real(dp), intent(in) :: rhs(:) !! short Newton right-hand side
        real(dp), intent(in) :: rtol !! requested inner tolerance
        integer, intent(in) :: maxiter !! inner iteration budget
        real(dp), intent(out) :: solution(:) !! short Newton correction
        logical, intent(out) :: ok !! whether inversion succeeded
        solver_calls = solver_calls+1
        call preconditioner(probe(:3),observed(:3),ok)
        if (.not. ok) return
        call preconditioner(rhs,solution,ok)
    end subroutine inspect_small

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! linearized Jacobian
        procedure(vector_root_linear_operator) :: preconditioner !! fitted inverse action
        real(dp), intent(in) :: rhs(:) !! Newton right-hand side
        real(dp), intent(in) :: rtol !! desired inner linear tolerance
        integer, intent(in) :: maxiter !! inner Krylov budget
        real(dp), intent(out) :: solution(:) !! Newton correction
        logical, intent(out) :: ok !! preconditioner action succeeded
        solver_calls = solver_calls+1
        call preconditioner(probe,observed,ok)
        if (.not. ok) return
        call preconditioner(rhs,solution,ok)
    end subroutine inspect
end program test_optimize_root_krylov_rank_two_sketch
