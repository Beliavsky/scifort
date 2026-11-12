! SPDX-License-Identifier: MIT
! Optional CPU-time cost model for adaptive Krylov preconditioner updates.
program test_optimize_root_krylov_timed_refresh
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 5
    integer :: failures, jv_calls, solver_calls, residual_calls
    real(dp) :: x0(n), start_time, end_time, sink
    logical :: clock_available
    type(vector_root_result) :: fit

    failures = 0
    sink = 0.0_dp
    x0 = 1.0_dp
    call cpu_time(start_time)
    clock_available = start_time >= 0.0_dp

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', &
        maxiter=5, fatol=1.0e-20_dp, inner_solver=slow_solver, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
        jacobian_product=jvp)
    call check(fit%nit == 5 .and. jv_calls == n, &
        'existing work-count policy is unchanged')

    if (clock_available) then
        call reset()
        fit = root(residual, x0, method='krylov', line_search='none', &
            maxiter=5, fatol=1.0e-20_dp, inner_solver=slow_solver, &
            inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
            inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
            inner_m_diag_update_timed_cost=.true., &
            inner_m_diag_update_timed_min_seconds=0.0005_dp, jacobian_product=jvp)
        call check(fit%nit == 5 .and. jv_calls == 2*n, &
            'increase in inner CPU time triggers exactly one early refresh')

        call reset()
        fit = newton_krylov(residual, x0, line_search='none', &
            maxiter=5, f_tol=1.0e-20_dp, inner_solver=slow_solver, &
            inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
            inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
            inner_m_diag_update_timed_cost=.true., &
            inner_m_diag_update_timed_min_seconds=0.0005_dp, jacobian_product=jvp)
        call check(fit%nit == 5 .and. jv_calls == 2*n, &
            'newton_krylov forwards CPU-timed refresh options')

        call reset()
        fit = root(residual, x0, method='krylov', line_search='none', &
            maxiter=5, fatol=1.0e-20_dp, inner_solver=slow_solver, &
            inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
            inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
            inner_m_diag_update_timed_cost=.true., &
            inner_m_diag_update_timed_min_seconds=0.5_dp, jacobian_product=jvp)
        call check(fit%nit == 5 .and. jv_calls == n, &
            'coarse or too-short timing intervals do not cause early refresh')
    end if

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_update_timed_cost=.true.)
    call invalid('timed cost requires work-growth and adaptive refresh')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_update_timed_cost=.true.)
    call invalid('timed and product-count cost models conflict')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_rebuild_cost=0, inner_m_diag_update_timed_cost=.true.)
    call invalid('timed and manually specified costs conflict')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_timed_min_seconds=0.0_dp)
    call invalid('minimum timing interval requires timed mode')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_timed_cost=.true., &
        inner_m_diag_update_timed_min_seconds=-1.0_dp)
    call invalid('negative minimum CPU seconds rejected')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_timed_cost=.true., &
        inner_m_diag_update_timed_min_seconds=ieee_value(0.0_dp, ieee_quiet_nan))
    call invalid('nonfinite minimum CPU seconds rejected')

    call reset()
    fit = root(residual, x0, method='hybr', &
        inner_m_diag_update_timed_cost=.true.)
    call invalid('non-Krylov method rejects timed options')

    call cpu_time(end_time)
    if (failures /= 0) error stop 1
    print '(a)', 'CPU-timed Krylov refresh tests passed'

contains

    subroutine reset()
        jv_calls = 0
        solver_calls = 0
        residual_calls = 0
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition !! expected behavior
        character(len=*), intent(in) :: description !! test case
        if (.not. condition) then
            print '(a)', 'FAIL: '//description
            failures = failures + 1
        end if
    end subroutine check

    subroutine invalid(description)
        character(len=*), intent(in) :: description !! invalid input case
        call check(fit%status == 0 .and. .not. fit%success .and. &
            fit%nfev == 0 .and. residual_calls == 0, description)
    end subroutine invalid

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(out) :: fvec(:) !! linear residual
        residual_calls = residual_calls + 1
        fvec = x
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! state at which J acts
        real(dp), intent(in) :: direction(:) !! Jacobian multiplication direction
        real(dp), intent(out) :: product(:) !! identity Jacobian action
        logical, intent(out) :: ok !! action validity
        jv_calls = jv_calls + 1
        product = direction
        ok = .true.
    end subroutine jvp

    subroutine slow_solver(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: rhs(:) !! linear residual
        real(dp), intent(in) :: rtol !! requested relative tolerance
        integer, intent(in) :: maxiter !! allowed linear iterations
        real(dp), intent(out) :: solution(:) !! deliberately damped linear update
        logical, intent(out) :: ok !! solver success
        real(dp) :: start, now, target
        real(dp) :: product(size(rhs))
        integer :: k

        solver_calls = solver_calls + 1
        call preconditioner(rhs, product, ok)
        if (.not. ok) return
        ! Windows CPU_TIME can advance in roughly 16 ms ticks. Use a baseline
        ! above one tick and a tenfold increase, with enough iterations to
        ! reach both targets on fast builds. Retain a finite fallback cap.
        target = 0.020_dp
        if (solver_calls == 2) target = 0.200_dp
        call cpu_time(start)
        do k = 1, 1000000000
            sink = sink + 1.0e-12_dp
            if (mod(k, 500) == 0) then
                call cpu_time(now)
                if (start < 0.0_dp .or. now < start) exit
                if (start >= 0.0_dp .and. now - start >= target) exit
            end if
        end do
        if (clock_available) then
            call check(now - start >= target, 'synthetic CPU workload reaches its target')
        end if
        solution = 0.5_dp * rhs
        ok = .true.
    end subroutine slow_solver
end program test_optimize_root_krylov_timed_refresh
