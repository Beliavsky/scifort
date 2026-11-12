! SPDX-License-Identifier: MIT
! Adaptive refresh scheduling from accepted nonlinear residual progress.
program test_optimize_root_krylov_adaptive_refresh
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none
    integer, parameter :: n = 5
    integer :: failures, jv_calls, function_calls, solver_calls
    logical :: transition
    real(dp) :: x0(n), relaxation
    type(vector_root_result) :: fit

    failures = 0
    x0 = 1.0_dp
    relaxation = 0.1_dp
    transition = .false.

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=6, fatol=1.0e-20_dp, inner_solver=damped_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=4, &
        jacobian_product=jvp)
    call check(fit%nit == 6 .and. jv_calls == 2*n, &
        'fixed four-step schedule refreshes at iterations 1 and 5')

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=6, fatol=1.0e-20_dp, inner_solver=damped_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=4, &
        inner_m_diag_update_adaptive=.true., jacobian_product=jvp)
    call check(fit%nit == 6 .and. jv_calls == 6*n, &
        'stalled progress triggers early refresh on each accepted step')

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=6, fatol=1.0e-20_dp, inner_solver=damped_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=4, &
        inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=0.95_dp, jacobian_product=jvp)
    call check(fit%nit == 6 .and. jv_calls == 2*n, &
        'high stall threshold preserves maximum interval')

    relaxation = 0.5_dp
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=6, fatol=1.0e-20_dp, inner_solver=damped_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=4, &
        inner_m_diag_update_adaptive=.true., jacobian_product=jvp)
    call check(fit%nit == 6 .and. jv_calls == 2*n, &
        'rapid progress does not trigger extra refreshes')

    ! A single weak step brings a refresh forward.  Later strong steps
    ! restart the maximum-age clock at that early refresh, not at step 1.
    relaxation = 0.1_dp
    transition = .true.
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=6, fatol=1.0e-20_dp, inner_solver=damped_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=4, &
        inner_m_diag_update_adaptive=.true., jacobian_product=jvp)
    call check(fit%nit == 6 .and. jv_calls == 3*n, &
        'an early refresh resets the maximum-age clock')
    transition = .false.

    relaxation = 0.1_dp
    call reset()
    fit = newton_krylov(residual, x0, line_search='none', &
        maxiter=6, f_tol=1.0e-20_dp, inner_solver=damped_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=4, &
        inner_m_diag_update_adaptive=.true., jacobian_product=jvp)
    call check(fit%nit == 6 .and. jv_calls == 6*n, &
        'newton_krylov forwards adaptive refresh settings')

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=6, fatol=1.0e-20_dp, inner_solver=damped_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_adaptive=.true., &
        jacobian_product=jvp)
    call check(fit%nit == 6 .and. jv_calls == 6*n, &
        'adaptive mode defaults to four maximum iterations between updates')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_update_adaptive=.true.)
    call invalid('adaptive refresh requires auto diagonal estimation')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_interval=0)
    call invalid('adaptive refresh rejects interval zero')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_stall_ratio=0.8_dp)
    call invalid('stall threshold requires adaptive mode')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_stall_ratio=0.0_dp)
    call invalid('stall threshold rejects zero')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_stall_ratio=1.1_dp)
    call invalid('stall threshold rejects values over one')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=ieee_value(0.0_dp,ieee_quiet_nan))
    call invalid('stall threshold rejects NaN')

    if (failures > 0) error stop 1
    print '(a)', 'adaptive diagonal refresh tests passed'
contains
    subroutine reset()
        jv_calls = 0
        function_calls = 0
        solver_calls = 0
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition !! expected test condition
        character(len=*), intent(in) :: description !! explanation of check
        if (.not. condition) then
            write (*,'(a)') 'FAIL: '//description
            failures = failures + 1
        end if
    end subroutine check

    subroutine invalid(description)
        character(len=*), intent(in) :: description !! invalid API input scenario
        call check(fit%status == 0 .and. .not. fit%success .and. &
            fit%nfev == 0 .and. function_calls == 0, description)
    end subroutine invalid

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(out) :: fvec(:) !! linear residual
        function_calls = function_calls + 1
        fvec = x
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: direction(:) !! Jacobian direction
        real(dp), intent(out) :: product(:) !! exact identity action
        logical, intent(out) :: ok !! true on successful application
        jv_calls = jv_calls + 1
        product = direction
        ok = .true.
    end subroutine jvp

    subroutine damped_solver(matvec, preconditioner, rhs, rtol, &
            maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse action
        real(dp), intent(in) :: rhs(:) !! current residual
        real(dp), intent(in) :: rtol !! requested relative error
        integer, intent(in) :: maxiter !! inner work budget
        real(dp), intent(out) :: solution(:) !! deliberately damped exact solve
        logical, intent(out) :: ok !! success indicator
        solver_calls = solver_calls + 1
        if (transition .and. solver_calls > 1) then
            solution = 0.5_dp*rhs
        else
            solution = relaxation*rhs
        end if
        ok = .true.
    end subroutine damped_solver
end program test_optimize_root_krylov_adaptive_refresh
