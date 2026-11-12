! SPDX-License-Identifier: MIT
! Nonsymmetric rank-two discovery using Jv and transpose-Jv products.
program test_optimize_root_krylov_rank_two_nonsymmetric
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 24
    real(dp) :: d(n), u(n,2), v(n,2), w(n), z(n), target(n)
    real(dp) :: probe(n), q(n), diagonal(n), observed(n), x0(n), bad_d(n), saved_v1(n)
    real(dp) :: sm11, sm12, sm21, sm22, alpha, coefficient
    integer :: i, failures, forward_calls, transpose_calls, f_calls, solve_calls
    integer :: fail_transpose_at, fail_forward_at
    logical :: incompatible, wrong_transpose
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        d(i) = 0.9_dp + 0.035_dp*real(i,dp)
        u(i,1) = 0.10_dp + 0.065_dp*sin(0.63_dp*real(i,dp))
        u(i,2) = -0.07_dp + 0.08_dp*cos(0.38_dp*real(i,dp))
        v(i,1) = 0.13_dp*cos(0.31_dp*real(i,dp)) + 0.03_dp
        v(i,2) = 0.14_dp*sin(0.51_dp*real(i,dp)) - 0.02_dp
        w(i) = 0.17_dp*cos(0.97_dp*real(i,dp))
        z(i) = 0.15_dp*sin(0.79_dp*real(i,dp))
        target(i) = 0.03_dp*cos(0.46_dp*real(i,dp))
        probe(i) = 0.1_dp + sin(0.19_dp*real(i,dp))
    end do
    x0 = 0.0_dp
    incompatible = .false.
    wrong_transpose = .false.
    call update_reference()

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(solve_calls == 1 .and. forward_calls == 4 .and. transpose_calls == 2, &
        'nonsymmetric discovery uses four Jv and two transpose Jv calls')
    call check(maxval(abs(observed-probe)) < 1.0e-11_dp, &
        'reconstructed nonsymmetric Woodbury inverse is exact')
    call check(fit%nfev == 2 .and. fit%nfev == f_calls, &
        'analytic transpose and forward products are not counted as residual evaluations')

    call reset()
    fit = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(forward_calls == 4 .and. transpose_calls == 2 .and. &
        maxval(abs(observed-probe)) < 1.0e-11_dp, &
        'newton_krylov forwards nonsymmetric options')

    call reset()
    fit = root(residual, x0, method='krylov', inner_method='gmres', inner_maxiter=1, &
        inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, maxiter=4, fatol=1.0e-11_dp, &
        line_search='none')
    call check(fit%success .and. maxval(abs(fit%x-target)) < 1.0e-10_dp, &
        'one-vector GMRES converges with learned nonsymmetric rank-two inverse')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        jacobian_transpose_product=transpose_jvp, inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe)) < 1.0e-5_dp, &
        'forward finite differences reconstruct nonsymmetric rank-two coupling')
    call check(fit%nfev == 6 .and. fit%nfev == f_calls .and. transpose_calls == 2, &
        'forward finite differences need four perturbation evaluations')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        jacobian_transpose_product=transpose_jvp, inner_solver=inspect, &
        jacobian_diff_scheme='3-point', line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe)) < 1.0e-5_dp, &
        'centered finite differences reconstruct nonsymmetric coupling')
    call check(fit%nfev == 10 .and. fit%nfev == f_calls, &
        'centered finite differences need eight perturbation evaluations')

    incompatible = .true.
    call update_reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-q/diagonal)) < 1.0e-11_dp, &
        'rank-three coupling falls back to exact diagonal')
    call check(forward_calls >= n+3 .and. transpose_calls == 2, &
        'rejected rank-three structure pays for exact coordinate diagonal')
    incompatible = .false.
    call update_reference()

    wrong_transpose = .true.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-q/diagonal)) < 1.0e-11_dp, &
        'inconsistent transpose products trigger diagonal fallback')
    wrong_transpose = .false.

    ! Keep the coupling rank two, but make its Woodbury 2-by-2 matrix
    ! singular by shifting the first right factor along the first left.
    saved_v1 = v(:,1)
    sm11 = 1.0_dp + dot_product(v(:,1),u(:,1)/d)
    sm12 = dot_product(v(:,1),u(:,2)/d)
    sm21 = dot_product(v(:,2),u(:,1)/d)
    sm22 = 1.0_dp + dot_product(v(:,2),u(:,2)/d)
    coefficient = dot_product(u(:,1),u(:,1)/d)*sm22 - &
        dot_product(u(:,1),u(:,2)/d)*sm21
    alpha = -(sm11*sm22-sm12*sm21)/coefficient
    v(:,1) = v(:,1) + alpha*u(:,1)
    call update_reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-q/diagonal)) < 1.0e-11_dp, &
        'singular rank-two Woodbury matrix reverts to exact diagonal scaling')
    call check(forward_calls >= n+4 .and. transpose_calls == 2, &
        'singular Woodbury fallback completes model fit before exact diagonal')
    v(:,1) = saved_v1
    call update_reference()

    call reset()
    fail_transpose_at = 2
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(fit%status == 0 .and. .not. fit%success .and. solve_calls == 0, &
        'failed transpose callback prevents inner solver invocation')

    call reset()
    fail_forward_at = 3
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(fit%status == 0 .and. .not. fit%success .and. solve_calls == 0, &
        'failed forward verification prevents inner solver invocation')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'missing transpose callback rejected before residual evaluation')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        jacobian_product=jvp, jacobian_transpose_product=transpose_jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'transpose callback without base rejected')

    call reset()
    fit = root(residual, x0, method='hybr', &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        jacobian_transpose_product=transpose_jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'nonsymmetric options rejected for non-Krylov methods')

    bad_d = d
    bad_d(3) = ieee_value(1.0_dp,ieee_quiet_nan)
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=bad_d, &
        jacobian_transpose_product=transpose_jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'nonfinite base diagonal rejected')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d(:n-1), &
        jacobian_transpose_product=transpose_jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'incorrect base diagonal length rejected')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        jacobian_transpose_product=transpose_jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'symmetric and nonsymmetric modes are exclusive')

    if (failures /= 0) error stop 1
    print '(a)', 'nonsymmetric rank-two low-product tests passed'

contains
    subroutine reset()
        forward_calls = 0
        transpose_calls = 0
        f_calls = 0
        solve_calls = 0
        fail_forward_at = 0
        fail_transpose_at = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, name)
        logical, intent(in) :: condition !! condition expected by the test
        character(len=*), intent(in) :: name !! assertion label
        if (.not. condition) then
            failures = failures+1
            print '(a)', 'FAILED: '//name
        end if
    end subroutine check

    subroutine update_reference()
        integer :: k
        call action(probe,q)
        do k = 1, n
            diagonal(k) = d(k) + u(k,1)*v(k,1) + u(k,2)*v(k,2)
            if (incompatible) diagonal(k) = diagonal(k) + w(k)*z(k)
        end do
    end subroutine update_reference

    subroutine action(vec, product)
        real(dp), intent(in) :: vec(:) !! Jacobian right direction
        real(dp), intent(out) :: product(:) !! Jacobian product
        product = d*vec + u(:,1)*dot_product(v(:,1),vec) + &
            u(:,2)*dot_product(v(:,2),vec)
        if (incompatible) product = product + w*dot_product(z,vec)
    end subroutine action

    subroutine transpose_action(vec, product)
        real(dp), intent(in) :: vec(:) !! Jacobian left direction
        real(dp), intent(out) :: product(:) !! transpose product
        product = d*vec + v(:,1)*dot_product(u(:,1),vec) + &
            v(:,2)*dot_product(u(:,2),vec)
        if (incompatible) product = product + z*dot_product(w,vec)
    end subroutine transpose_action

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current nonlinear iterate
        real(dp), intent(out) :: fvec(:) !! nonlinear residual
        call action(x-target,fvec)
        f_calls = f_calls+1
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(in) :: vec(:) !! Jacobian direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector action
        logical, intent(out) :: ok !! success flag
        forward_calls = forward_calls+1
        ok = fail_forward_at == 0 .or. forward_calls < fail_forward_at
        product = 0.0_dp
        if (ok) call action(vec,product)
    end subroutine jvp

    subroutine transpose_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(in) :: vec(:) !! transpose direction
        real(dp), intent(out) :: product(:) !! transpose-Jacobian action
        logical, intent(out) :: ok !! success flag
        transpose_calls = transpose_calls+1
        ok = fail_transpose_at == 0 .or. transpose_calls < fail_transpose_at
        product = 0.0_dp
        if (.not. ok) return
        if (wrong_transpose) then
            call action(vec,product)
        else
            call transpose_action(vec,product)
        end if
    end subroutine transpose_jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse action
        real(dp), intent(in) :: rhs(:) !! nonlinear solve right-hand side
        real(dp), intent(in) :: rtol !! relative inner tolerance
        integer, intent(in) :: maxiter !! Krylov iteration limit
        real(dp), intent(out) :: solution(:) !! proposed Newton correction
        logical, intent(out) :: ok !! action status
        solve_calls = solve_calls+1
        call preconditioner(q,observed,ok)
        solution = rhs
    end subroutine inspect
end program test_optimize_root_krylov_rank_two_nonsymmetric
