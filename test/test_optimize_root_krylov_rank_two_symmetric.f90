! SPDX-License-Identifier: MIT
! Low-product symmetric rank-two discovery with a known uncoupled diagonal.
program test_optimize_root_krylov_rank_two_symmetric
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 16
    real(dp) :: d(n), u(n,2), w(n), target(n), probe(n), q(n)
    real(dp) :: observed(n), diagonal(n), x0(n), nan_d(n)
    integer :: i, failures, jv_calls, f_calls, solve_calls, fail_after
    logical :: incompatible
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        d(i) = 0.8_dp + 0.07_dp*real(i,dp)
        u(i,1) = 0.12_dp + 0.08_dp*sin(0.67_dp*real(i,dp))
        u(i,2) = 0.13_dp*cos(0.39_dp*real(i,dp))
        w(i) = 0.18_dp*sin(0.95_dp*real(i,dp))
        target(i) = 0.03_dp*cos(0.4_dp*real(i,dp))
        probe(i) = sin(0.17_dp*real(i,dp)) + 0.2_dp*real(i,dp)/real(n,dp)
    end do
    x0 = 0.0_dp
    incompatible = .false.
    call update_reference()

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(solve_calls == 1 .and. jv_calls == 4, &
        'symmetric rank-two discovery needs only four products')
    call check(maxval(abs(observed-probe)) < 1.0e-11_dp, &
        'learned Woodbury inverse is exact for symmetric rank-two coupling')

    call reset()
    fit = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jv_calls == 4 .and. maxval(abs(observed-probe)) < 1.0e-11_dp, &
        'newton_krylov forwards symmetric rank-two base')

    call reset()
    fit = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=1, inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, jacobian_product=jvp, &
        line_search='none', maxiter=3, fatol=1.0e-11_dp)
    call check(fit%success .and. maxval(abs(fit%x-target)) < 1.0e-10_dp, &
        'one-vector GMRES converges using the learned rank-two inverse')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, inner_solver=inspect, &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe)) < 1.0e-5_dp, &
        'four forward finite-difference products reconstruct the coupling')
    call check(fit%nfev == 6 .and. fit%nfev == f_calls, &
        'two-point difference count: four perturbations plus two residuals')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, inner_solver=inspect, &
        jacobian_diff_scheme='3-point', line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe)) < 1.0e-5_dp, &
        'centered finite differences reconstruct symmetric coupling')
    call check(fit%nfev == 10 .and. fit%nfev == f_calls, &
        'three-point difference count: eight perturbations plus two residuals')

    incompatible = .true.
    call update_reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-q/diagonal)) < 1.0e-11_dp, &
        'rank-three contamination falls back to the exact Jacobian diagonal')
    call check(jv_calls >= n+3 .and. jv_calls <= n+4, &
        'failed structural check uses n exact diagonal products')
    incompatible = .false.
    call update_reference()

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d(:n-1), jacobian_product=jvp, &
        inner_solver=inspect)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'invalid symmetric base length rejected before residual evaluation')

    nan_d = d
    nan_d(4) = ieee_value(1.0_dp,ieee_quiet_nan)
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=nan_d, jacobian_product=jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'nonfinite supplied diagonal rejected before residual evaluation')

    nan_d = d
    nan_d(4) = 0.0_dp
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=nan_d, jacobian_product=jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'singular supplied diagonal rejected before residual evaluation')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_stream=.true., &
        inner_m_diag_rank_two_symmetric_base=d, jacobian_product=jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'stream and symmetric discovery cannot both be enabled')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, jacobian_product=jvp, &
        inner_solver=inspect, maxiter=1, line_search='none')
    call check(jv_calls == 4, 'all successful analytic products are counted')

    call reset()
    fail_after = 3
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, jacobian_product=jvp, &
        inner_solver=inspect, maxiter=1, line_search='none')
    call check(fit%status == 0 .and. .not. fit%success .and. solve_calls == 0, &
        'failed verification callback propagates without using a preconditioner')
    fail_after = 0

    if (failures /= 0) error stop 1
    print '(a)', 'symmetric rank-two low-product tests passed'

contains

    subroutine reset()
        jv_calls = 0
        f_calls = 0
        solve_calls = 0
        fail_after = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, name)
        logical, intent(in) :: condition !! true for a passing assertion
        character(len=*), intent(in) :: name !! assertion description
        if (.not. condition) then
            failures = failures+1
            print '(a)', 'FAILED: '//name
        end if
    end subroutine check

    subroutine update_reference()
        integer :: k
        call action(probe,q)
        do k = 1, n
            diagonal(k) = d(k) + u(k,1)**2 - 0.65_dp*u(k,2)**2
            if (incompatible) diagonal(k) = diagonal(k) + w(k)**2
        end do
    end subroutine update_reference

    subroutine action(vec, product)
        real(dp), intent(in) :: vec(:) !! Jacobian direction
        real(dp), intent(out) :: product(:) !! exact Jacobian action
        product = d*vec + u(:,1)*dot_product(u(:,1),vec) - &
            0.65_dp*u(:,2)*dot_product(u(:,2),vec)
        if (incompatible) product = product + w*dot_product(w,vec)
    end subroutine action

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(out) :: fvec(:) !! linear system residual
        call action(x-target,fvec)
        f_calls = f_calls + 1
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! Jacobian evaluation point
        real(dp), intent(in) :: vec(:) !! multiplication direction
        real(dp), intent(out) :: product(:) !! Jacobian action
        logical, intent(out) :: ok !! successful evaluation
        jv_calls = jv_calls+1
        ok = fail_after == 0 .or. jv_calls < fail_after
        product = 0.0_dp
        if (ok) call action(vec,product)
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian operator
        procedure(vector_root_linear_operator) :: preconditioner !! inverse action
        real(dp), intent(in) :: rhs(:) !! nonlinear system right side
        real(dp), intent(in) :: rtol !! requested tolerance
        integer, intent(in) :: maxiter !! linear iteration budget
        real(dp), intent(out) :: solution(:) !! proposed Newton solution
        logical, intent(out) :: ok !! successful preconditioner application
        solve_calls = solve_calls+1
        call preconditioner(q,observed,ok)
        solution = rhs
    end subroutine inspect
end program test_optimize_root_krylov_rank_two_symmetric
