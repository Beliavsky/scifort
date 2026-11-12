! SPDX-License-Identifier: MIT
program test_optimize_dual_annealing_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_dual_annealing
    implicit none
    type(dual_annealing_result) :: res, repeated
    real(dp) :: bounds(2,2), one(1,2)
    integer :: checks, callback_calls

    checks = 0
    callback_calls = 0
    bounds(1,:) = [-5.0_dp, 5.0_dp]
    bounds(2,:) = [-5.0_dp, 5.0_dp]
    one(1,:) = [-3.0_dp, 3.0_dp]

    res = dual_annealing_core(quadratic, bounds, maxiter=0, seed=13, x0=[2.0_dp,-2.0_dp])
    call must(res%status == da_maxiter .and. res%success, 'zero iterations succeeds')
    call must(res%nit == 0 .and. res%nfev == 1, 'zero iteration eval count')
    call must(all(res%x == [2.0_dp,-2.0_dp]), 'x0 retained in zero iterations')
    call must(abs(res%fun - quadratic(res%x)) < 1.0e-15_dp, 'initial objective reported')
    call must(res%nlocal == 0, 'no implicit local minimizer')
    res = dual_annealing_core(quadratic, bounds, maxiter=40, seed=456, x0=[3.0_dp,-4.0_dp])
    repeated = dual_annealing_core(quadratic, bounds, maxiter=40, seed=456, x0=[3.0_dp,-4.0_dp])
    call must(res%status == da_maxiter .and. res%success, 'normal iteration termination')
    call must(all(res%x == repeated%x), 'exact seed reproducibility x')
    call must(res%fun == repeated%fun, 'exact seed reproducibility f')
    call must(res%naccepted == repeated%naccepted, 'exact seed reproducibility count')
    call must(res%nfev == 161 .and. res%nit == 40, 'objective evaluation counts')
    call must(res%fun <= quadratic([3.0_dp,-4.0_dp]), 'global best improves')
    call must(all(res%x >= bounds(:,1)) .and. all(res%x <= bounds(:,2)), 'feasible solution')
    call must(res%naccepted >= 1 .and. res%naccepted <= 160, 'accepted count bounds')
    call must(res%nbest_updates > 0, 'global best updates')

    res = dual_annealing_core(quadratic, bounds, maxiter=30, maxfun=5, seed=42)
    call must(res%status == da_maxfun .and. .not. res%success, 'maxfun termination')
    call must(res%nfev == 5, 'no additional objective calls after budget')

    res = dual_annealing_core(quadratic, bounds, maxiter=9, &
        restart_temp_ratio=0.99_dp, seed=18)
    call must(res%nrestarts > 0, 'reannealing path exercised')
    call must(res%status == da_maxiter, 'reannealing completes')
    call must(abs(res%fun-quadratic(res%x)) < 1.0e-12_dp, 'best after restart evaluated correctly')

    res = dual_annealing_core(quadratic, bounds, maxiter=100, seed=15, x0=[4.0_dp,4.0_dp], &
        callback=stop_on_improvement)
    call must(res%status == da_callback_stop .and. .not. res%success, 'early stop callback')
    call must(callback_calls == 1, 'callback called exactly once before stop')
    call must(res%nbest_updates == 1, 'callback only on new best')
    call must(abs(res%fun - quadratic(res%x)) < 1.0e-12_dp, 'callback returned true best')

    res = dual_annealing_core(quadratic, bounds, maxiter=10, seed=15, &
        minimizer=exact_local, local_interval=1)
    call must(res%status == da_maxiter .and. res%nlocal >= 1, 'user local minimizer used')
    call must(res%fun < 1.0e-25_dp, 'local optimizer finds known optimum')
    call must(res%nfev >= res%nlocal*2+1, 'local objective evaluations included')
    res = dual_annealing_core(quadratic, bounds, maxiter=10, seed=15, &
        minimizer=exact_local, no_local_search=.true.)
    call must(res%nlocal == 0, 'disable local minimizer')
    res = dual_annealing_core(quadratic, bounds, no_local_search=.false.)
    call must(res%status == da_invalid, 'local callback required if enabled')

    res = dual_annealing_core(quadratic, bounds, maxiter=15, seed=8, visit=1.2_dp)
    call must(res%status == da_maxiter, 'low visiting parameter')
    res = dual_annealing_core(quadratic, bounds, maxiter=15, seed=8, visit=2.95_dp)
    call must(res%status == da_maxiter, 'high visiting parameter')
    res = dual_annealing_core(quadratic, bounds, maxiter=15, seed=8, accept=-10000.0_dp)
    call must(res%status == da_maxiter, 'acceptance parameter endpoint')
    res = dual_annealing_core(negative_one, one, maxiter=40, seed=2)
    call must(res%status == da_maxiter .and. res%fun == -1.0_dp, 'constant objective')

    res = dual_annealing_core(quadratic, bounds(:,1:1))
    call must(res%status == da_invalid, 'reject bounds column shape')
    res = dual_annealing_core(quadratic, bounds(1:0,:))
    call must(res%status == da_invalid, 'reject no variables')
    res = dual_annealing_core(quadratic, bounds, x0=[0.0_dp])
    call must(res%status == da_invalid, 'reject x0 length')
    res = dual_annealing_core(quadratic, bounds, x0=[0.0_dp,20.0_dp])
    call must(res%status == da_invalid, 'reject infeasible x0')
    bounds(1,2) = bounds(1,1)
    res = dual_annealing_core(quadratic, bounds)
    call must(res%status == da_invalid, 'reject zero-width bounds')
    bounds(1,:) = [-5.0_dp, 5.0_dp]
    bounds(1,2) = ieee_value(1.0_dp,ieee_positive_inf)
    res = dual_annealing_core(quadratic, bounds)
    call must(res%status == da_invalid, 'reject infinite bound')
    bounds(1,:) = [-5.0_dp, 5.0_dp]
    res = dual_annealing_core(quadratic, bounds, maxiter=-1)
    call must(res%status == da_invalid, 'reject negative maxiter')
    res = dual_annealing_core(always_infinite, bounds, maxfun=3, seed=10)
    call must(res%status == da_no_initial .and. res%nfev == 3, &
        'nonfinite initialization exhausts evaluation budget')
    res = dual_annealing_core(quadratic, bounds, maxfun=0)
    call must(res%status == da_invalid, 'reject zero maxfun')
    res = dual_annealing_core(quadratic, bounds, local_interval=0)
    call must(res%status == da_invalid, 'reject invalid local interval')
    res = dual_annealing_core(quadratic, bounds, visit=1.0_dp)
    call must(res%status == da_invalid, 'reject visiting parameter at 1')
    res = dual_annealing_core(quadratic, bounds, visit=3.0_dp)
    call must(res%status == da_invalid, 'reject unsupported singular visit endpoint')
    res = dual_annealing_core(quadratic, bounds, initial_temp=0.01_dp)
    call must(res%status == da_invalid, 'reject too-small initial temperature')
    res = dual_annealing_core(quadratic, bounds, restart_temp_ratio=1.0_dp)
    call must(res%status == da_invalid, 'reject invalid restart ratio')
    res = dual_annealing_core(quadratic, bounds, accept=0.0_dp)
    call must(res%status == da_invalid, 'reject invalid accept parameter')
    bounds(1,:) = [-huge(1.0_dp),huge(1.0_dp)]
    res = dual_annealing_core(quadratic, bounds)
    call must(res%status == da_invalid, 'reject overflowed bound span')

    print '(a,1x,i0)', 'DUAL_ANNEALING_NATIVE_CHECKS', checks
contains
    function quadratic(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = (x(1) - 0.75_dp)**2 + (x(2) + 1.25_dp)**2
    end function quadratic

    function always_infinite(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = ieee_value(1.0_dp, ieee_positive_inf)
    end function always_infinite

    function negative_one(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = -1.0_dp + 0.0_dp*sum(x)
    end function negative_one

    subroutine exact_local(fun, lower, upper, start, x, value, nfev, success)
        procedure(da_objective) :: fun
        real(dp), intent(in) :: lower(:)
        real(dp), intent(in) :: upper(:)
        real(dp), intent(in) :: start(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: value
        integer, intent(out) :: nfev
        logical, intent(out) :: success
        real(dp) :: ignored
        ignored = fun(start)
        x = min(upper,max(lower,[0.75_dp,-1.25_dp]))
        value = fun(x)
        nfev = 2
        success = .true.
    end subroutine exact_local

    subroutine stop_on_improvement(x, value, context, stop)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: value
        integer, intent(in) :: context
        logical, intent(out) :: stop
        callback_calls = callback_calls + 1
        stop = .true.
    end subroutine stop_on_improvement

    subroutine must(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
        checks = checks + 1
    end subroutine must
end program test_optimize_dual_annealing_reference
