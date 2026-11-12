! SPDX-License-Identifier: MIT
! Adaptive Krylov preconditioner refresh based on inner work and rebuild cost.
program test_optimize_root_krylov_work_adaptive
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 5
    integer :: failures, jv_calls, function_calls, solver_calls
    integer :: work_mode
    real(dp) :: x0(n)
    type(vector_root_result) :: fit

    failures = 0
    x0 = 1.0_dp
    work_mode = 1

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=5, fatol=1.0e-20_dp, inner_solver=work_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., jacobian_product=jvp)
    call check(fit%nit == 5 .and. jv_calls == n, &
        'default adaptive mode ignores increased inner work during strong progress')

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=5, fatol=1.0e-20_dp, inner_solver=work_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_rebuild_cost=1, jacobian_product=jvp)
    call check(fit%nit == 5 .and. jv_calls == 2*n, &
        'growing inner work over rebuild cost triggers one early refresh')

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=5, fatol=1.0e-20_dp, inner_solver=work_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_rebuild_cost=4, jacobian_product=jvp)
    call check(fit%nit == 5 .and. jv_calls == n, &
        'high estimated rebuild cost suppresses extra refresh')

    call reset()
    fit = newton_krylov(residual, x0, line_search='none', &
        maxiter=5, f_tol=1.0e-20_dp, inner_solver=work_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_rebuild_cost=1, jacobian_product=jvp)
    call check(fit%nit == 5 .and. jv_calls == 2*n, &
        'newton_krylov forwards work-aware refresh options')

    work_mode = 2
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=5, fatol=1.0e-20_dp, inner_solver=work_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
        jacobian_product=jvp)
    call check(fit%nit == 5 .and. jv_calls == n, &
        'zero baseline work does not spuriously trigger adaptive refresh')

    work_mode = 1
    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_update_work_growth=2.0_dp)
    call invalid('growth threshold requires adaptive refresh')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_rebuild_cost=1)
    call invalid('rebuild cost requires work growth')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=1.0_dp)
    call invalid('growth factor must exceed one')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_rebuild_cost=-1)
    call invalid('rebuild cost must be nonnegative')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_work_growth=ieee_value(0.0_dp, ieee_quiet_nan))
    call invalid('growth factor must be finite')

    call reset()
    fit = root(residual, x0, method='hybr', &
        inner_m_diag_update_work_growth=2.0_dp)
    call invalid('non-Krylov methods reject the new option')

    if (failures /= 0) error stop 1
    print '(a)', 'work-aware Krylov refresh tests passed'

contains

    subroutine reset()
        jv_calls = 0
        function_calls = 0
        solver_calls = 0
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition !! true for a successful test
        character(len=*), intent(in) :: description !! test description
        if (.not. condition) then
            print '(a)', 'FAIL: '//description
            failures = failures + 1
        end if
    end subroutine check

    subroutine invalid(description)
        character(len=*), intent(in) :: description !! rejected option case
        call check(fit%status == 0 .and. .not. fit%success .and. &
            fit%nfev == 0 .and. function_calls == 0, description)
    end subroutine invalid

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(out) :: fvec(:) !! identity residual
        function_calls = function_calls + 1
        fvec = x
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(in) :: direction(:) !! probe direction
        real(dp), intent(out) :: product(:) !! exact identity action
        logical, intent(out) :: ok !! whether evaluation succeeded
        jv_calls = jv_calls + 1
        product = direction
        ok = .true.
    end subroutine jvp

    subroutine work_solver(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! approximate inverse
        real(dp), intent(in) :: rhs(:) !! residual
        real(dp), intent(in) :: rtol !! requested relative linear tolerance
        integer, intent(in) :: maxiter !! allowed Krylov basis size
        real(dp), intent(out) :: solution(:) !! damped search direction
        logical, intent(out) :: ok !! whether preconditioner succeeded
        real(dp) :: tmp(size(rhs))
        integer :: j, calls

        solver_calls = solver_calls + 1
        calls = 1
        if (solver_calls == 2) calls = 5
        if (work_mode == 2 .and. solver_calls == 1) calls = 0
        do j = 1, calls
            call preconditioner(rhs, tmp, ok)
            if (.not. ok) return
        end do
        solution = 0.5_dp * rhs
        ok = .true.
    end subroutine work_solver
end program test_optimize_root_krylov_work_adaptive
