! SPDX-License-Identifier: MIT
! Opt-in system_clock wall-time cost model for adaptive Krylov refresh.
program test_optimize_root_krylov_wall_clock_refresh
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 5
    integer :: failures, jv_calls, solver_calls, residual_calls
    real(dp) :: x0(n), sink
    integer(int64) :: clock_tick, clock_rate
    logical :: clock_available
    type(vector_root_result) :: fit

    failures = 0
    sink = 0.0_dp
    x0 = 1.0_dp
    call system_clock(count=clock_tick, count_rate=clock_rate)
    clock_available = clock_tick >= 0_int64 .and. clock_rate > 0_int64

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
            inner_m_diag_update_wall_clock=.true., &
            inner_m_diag_update_timed_min_seconds=0.0005_dp, jacobian_product=jvp)
        call check(fit%nit == 5 .and. jv_calls == 2*n, &
            'increase in inner wall time triggers exactly one early refresh')

        call reset()
        fit = newton_krylov(residual, x0, line_search='none', &
            maxiter=5, f_tol=1.0e-20_dp, inner_solver=slow_solver, &
            inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
            inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
            inner_m_diag_update_timed_cost=.true., &
            inner_m_diag_update_wall_clock=.true., &
            inner_m_diag_update_timed_min_seconds=0.0005_dp, jacobian_product=jvp)
        call check(fit%nit == 5 .and. jv_calls == 2*n, &
            'newton_krylov forwards wall-clock refresh option')

        call reset()
        fit = root(residual, x0, method='krylov', line_search='none', &
            maxiter=5, fatol=1.0e-20_dp, inner_solver=slow_solver, &
            inner_m_diag_auto=.true., inner_m_diag_update_interval=6, &
            inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
            inner_m_diag_update_timed_cost=.true., &
            inner_m_diag_update_wall_clock=.true., &
            inner_m_diag_update_timed_min_seconds=0.5_dp, jacobian_product=jvp)
        call check(fit%nit == 5 .and. jv_calls == n, &
            'minimum wall-time threshold prevents early refresh')
    end if

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_update_wall_clock=.true.)
    call invalid('wall clock requires CPU-timed refresh mode')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_wall_clock=.true.)
    call invalid('wall-clock option requires timed cost even with adaptive policy')

    call reset()
    fit = root(residual, x0, method='hybr', &
        inner_m_diag_update_wall_clock=.true.)
    call invalid('non-Krylov method rejects wall-clock option')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_timed_cost=.true., &
        inner_m_diag_update_wall_clock=.false., &
        inner_m_diag_update_timed_min_seconds=-1.0_dp)
    call invalid('explicit false wall option retains CPU-time validation')

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

    if (failures /= 0) error stop 1
    print '(a)', 'Wall-clock Krylov refresh tests passed'

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
        real(dp) :: elapsed, target
        real(dp) :: product(size(rhs))
        integer(int64) :: start, now, rate
        integer :: k

        solver_calls = solver_calls + 1
        call preconditioner(rhs, product, ok)
        if (.not. ok) return
        ! Drive this workload with the same clock as the wall-time cost model.
        ! CPU_TIME can have a much coarser resolution, masking the intended
        ! growth ratio. Keep a wide margin over scheduling and clock jitter.
        target = 0.020_dp
        if (solver_calls == 2) target = 0.200_dp
        call system_clock(count=start, count_rate=rate)
        do k = 1, 1000000000
            sink = sink + 1.0e-12_dp
            if (mod(k, 500) == 0) then
                call system_clock(count=now)
                if (rate <= 0_int64 .or. start < 0_int64 .or. now < start) exit
                elapsed = real(now - start, dp) / real(rate, dp)
                if (elapsed >= target) exit
            end if
        end do
        solution = 0.5_dp * rhs
        ok = .true.
    end subroutine slow_solver
end program test_optimize_root_krylov_wall_clock_refresh
