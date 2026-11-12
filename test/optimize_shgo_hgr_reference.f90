! SPDX-License-Identifier: MIT
! Homology-birth-count stabilization stopping and boundary cases.
program optimize_shgo_hgr_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_shgo
    implicit none
    type(shgo_result) :: result, again
    real(dp) :: interval(1,2), box(2,2)
    integer :: checks, callback_calls, objective_calls
    logical :: stop_at_third

    checks = 0
    callback_calls = 0
    objective_calls = 0
    stop_at_third = .false.
    interval(1,:) = [-1.0_dp,1.0_dp]
    box(:,1) = -1.0_dp
    box(:,2) = 1.0_dp

    ! HGR stopping is explicitly opt-in and does not change earlier defaults.
    result = shgo_core(constant,interval,iters=4)
    call check(result%status == shgo_completed .and. result%nit == 4, &
        'HGR disabled by default')
    call check(result%nfev == 17 .and. result%hgr_streak == 0, &
        'disabled HGR preserves full sampling')

    ! Three rounds means two genuine comparisons; the initial observation
    ! is never miscounted as evidence of stability.
    result = shgo_core(constant,interval,iters=5,hgr_stability_iters=2)
    call check(result%status == shgo_hgr_stable .and. result%success, &
        'stable HGR stops with distinct success status')
    call check(result%nit == 3 .and. result%nfev == 9, &
        'minimum rounds and two equal-rank comparisons')
    call check(result%hgr_estimate == 1 .and. result%hgr_streak == 2, &
        'HGR rank and streak diagnostics')
    again = shgo_core(constant,interval,iters=5,hgr_stability_iters=2)
    call check(again%status == result%status .and. again%nfev == result%nfev, &
        'HGR stopping is reproducible')
    result = shgo_core(constant,interval,iters=5,hgr_stability_iters=1)
    call check(result%status == shgo_hgr_stable .and. result%nit == 3, &
        'default three-round minimum applies even to one comparison')
    result = shgo_core(constant,interval,iters=5,hgr_stability_iters=1, &
        hgr_min_iters=4)
    call check(result%status == shgo_hgr_stable .and. result%nit == 4 .and. &
        result%nfev == 17, 'user-specified minimum round count')
    result = shgo_core(constant,interval,iters=4,hgr_stability_iters=4)
    call check(result%status == shgo_completed .and. result%hgr_streak == 3, &
        'insufficient consecutive comparisons does not stop')

    ! New basins reset the counter, rather than allowing a stale rank to
    ! trigger early termination. Here a second well appears on level 2.
    result = shgo_core(double_well,interval,iters=5,hgr_stability_iters=2, &
        hgr_f_tol=0.0_dp)
    call check(result%status == shgo_hgr_stable .and. result%nit == 4, &
        'basin discovery postpones rank stabilization')
    call check(result%hgr_estimate == 2 .and. result%hgr_streak == 2, &
        'two independent graph minima detected')

    ! An improving incumbent must reset stability even if the graph rank
    ! remains unchanged. No local solver is used in this experiment.
    result = shgo_core(shifted,interval,iters=3,hgr_stability_iters=1, &
        hgr_f_tol=0.0_dp)
    call check(result%status == shgo_completed .and. result%nit == 3, &
        'material sampling improvement prevents premature HGR stop')
    call check(result%hgr_estimate == 1 .and. result%hgr_streak == 0, &
        'improvement resets HGR streak')
    result = shgo_core(shifted,interval,iters=5,hgr_stability_iters=1, &
        hgr_f_tol=0.0_dp)
    call check(result%status == shgo_hgr_stable .and. result%nit == 4, &
        'stable minimum after improvement permits later HGR stop')

    ! Adaptive mesh has its own nonconforming graph but the same semantics.
    result = shgo_core(constant,box,iters=5,adaptive=.true., &
        hgr_stability_iters=2,local_iter=0)
    call check(result%status == shgo_hgr_stable .and. result%nit == 3, &
        'adaptive mesh supports HGR stopping')
    call check(result%hgr_estimate == 1 .and. result%hgr_streak == 2, &
        'adaptive graph rank remains stable')
    call check(result%nmesh == result%nfev .and. result%nrefined > 0, &
        'adaptive mesh diagnostics remain consistent')
    result = shgo_core(two_plateaus,box,iters=5,adaptive=.true., &
        adaptive_fraction=1.0_dp,hgr_stability_iters=2,local_iter=5)
    call check(result%status == shgo_hgr_stable .and. result%nit == 3, &
        'disjoint equal-valued plateaus remain distinct')
    call check(result%hgr_estimate == 2 .and. result%nminima == 2, &
        'connected equal-valued vertices form one minimum per plateau')
    result = shgo_core(double_well,interval,iters=5,adaptive=.true., &
        adaptive_fraction=1.0_dp,hgr_stability_iters=2,hgr_f_tol=0.0_dp)
    call check(result%status == shgo_hgr_stable .and. result%nit == 4 .and. &
        result%hgr_estimate == 2, 'fully refined adaptive basin-rank agreement')

    ! If polishing is deferred, run it at a tentative stopping level.
    objective_calls = 0
    result = shgo_core(constant,interval,iters=5,hgr_stability_iters=2, &
        minimize_every_iter=.false.,minimizer=local_check)
    call check(result%status == shgo_hgr_stable .and. result%nit == 3, &
        'deferred local search is not skipped at HGR termination')
    call check(result%nlocal == 1 .and. result%nlfev == 1, &
        'deferred polishing runs once at early-stop round')
    call check(result%nfev == objective_calls, 'deferred HGR local objective accounting')

    ! A local improvement at a tentative stopping level cancels termination.
    result = shgo_core(shifted,interval,iters=6,hgr_stability_iters=1, &
        minimize_every_iter=.false.,minimizer=local_exact,hgr_f_tol=0.0_dp)
    call check(result%status == shgo_hgr_stable .and. result%nit == 5, &
        'deferred local improvement postpones HGR stop')
    call check(result%fun <= 1.0e-25_dp .and. result%nlocal >= 1, &
        'deferred local polishing improves sampled minimum')

    ! Callback and evaluation budgets have precedence over HGR exit.
    callback_calls = 0
    stop_at_third = .true.
    result = shgo_core(constant,interval,iters=5,hgr_stability_iters=2, &
        callback=callback)
    call check(result%status == shgo_callback_stop .and. result%nit == 3, &
        'callback overrides same-round HGR stabilization')
    call check(callback_calls == 3, 'callback invoked on all completed rounds')
    stop_at_third = .false.
    result = shgo_core(constant,interval,iters=5,maxfun=9,hgr_stability_iters=2)
    call check(result%status == shgo_hgr_stable .and. result%nfev == 9, &
        'HGR can stop exactly at evaluation budget')
    result = shgo_core(constant,interval,iters=5,maxfun=5,hgr_stability_iters=2)
    call check(result%status == shgo_maxfun .and. result%nit == 2, &
        'insufficient budget takes precedence before HGR stabilization')
    result = shgo_core(nan_value,interval,iters=4,hgr_stability_iters=1)
    call check(result%status == shgo_numerical .and. .not. result%success, &
        'all nonfinite values never produce a false stable rank')
    call check(result%hgr_estimate == 0 .and. result%hgr_streak == 0, &
        'empty graph never increments stability counter')
    result = shgo_core(constant,interval,f_min=0.0_dp,hgr_stability_iters=1)
    call check(result%status == shgo_target .and. result%nfev == 1, &
        'known target wins at first sample')

    call invalid(shgo_core(constant,interval,hgr_stability_iters=-1), &
        'negative stabilization length')
    call invalid(shgo_core(constant,interval,hgr_min_iters=1), &
        'invalid minimum iteration count')
    call invalid(shgo_core(constant,interval,hgr_f_tol=-1.0_dp), &
        'negative improvement tolerance')
    call invalid(shgo_core(constant,interval, &
        hgr_f_tol=ieee_value(0.0_dp,ieee_quiet_nan)), 'NaN improvement tolerance')

    print '(a,i0,a)', 'SHGO HGR regression: ', checks, ' checks passed'
contains
    function constant(x) result(f)
        real(dp), intent(in) :: x(:) !! box coordinate
        real(dp) :: f !! constant scalar objective
        objective_calls = objective_calls + 1
        f = 0.0_dp + 0.0_dp*sum(x)
    end function constant

    function double_well(x) result(f)
        real(dp), intent(in) :: x(:) !! one-dimensional point
        real(dp) :: f !! two minima at -0.5 and +0.5
        f = (x(1)**2-0.25_dp)**2
    end function double_well

    function shifted(x) result(f)
        real(dp), intent(in) :: x(:) !! one-dimensional point
        real(dp) :: f !! shifted convex objective
        f = (x(1)-0.3_dp)**2
    end function shifted

    function two_plateaus(x) result(f)
        real(dp), intent(in) :: x(:) !! two-dimensional point
        real(dp) :: f !! two disjoint flat minima in the first coordinate
        f = 1.0_dp
        if (abs(x(1)) > 0.25_dp) f = 0.0_dp
    end function two_plateaus

    function nan_value(x) result(f)
        real(dp), intent(in) :: x(:) !! one-dimensional point
        real(dp) :: f !! intentionally nonfinite objective
        f = ieee_value(0.0_dp,ieee_quiet_nan)
        if (size(x) /= 1) f = 0.0_dp
    end function nan_value

    subroutine local_check(fun,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: fun !! budget-guarded function
        real(dp), intent(in) :: lower(:) !! search lower bound
        real(dp), intent(in) :: upper(:) !! search upper bound
        real(dp), intent(in) :: start(:) !! current candidate
        integer, intent(in) :: maxfev !! available evaluation calls
        real(dp), intent(out) :: x(:) !! returned local point
        real(dp), intent(out) :: value !! returned objective
        integer, intent(out) :: nfev !! reported evaluation count
        logical, intent(out) :: success !! accepted solution flag
        x = min(upper,max(lower,start))
        value = fun(x)
        nfev = 1
        success = maxfev >= 1
    end subroutine local_check

    subroutine local_exact(fun,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: fun !! budget-guarded function
        real(dp), intent(in) :: lower(:) !! search lower bound
        real(dp), intent(in) :: upper(:) !! search upper bound
        real(dp), intent(in) :: start(:) !! current candidate
        integer, intent(in) :: maxfev !! available evaluation calls
        real(dp), intent(out) :: x(:) !! returned local point
        real(dp), intent(out) :: value !! returned objective
        integer, intent(out) :: nfev !! reported evaluation count
        logical, intent(out) :: success !! accepted solution flag
        x = min(upper,max(lower,0.3_dp))
        value = fun(x)
        nfev = 1
        success = maxfev >= 1 .and. size(start) == size(x)
    end subroutine local_exact

    subroutine callback(x,value,stop)
        real(dp), intent(in) :: x(:) !! best coordinate so far
        real(dp), intent(in) :: value !! best function value so far
        logical, intent(out) :: stop !! halt request at third completed round
        callback_calls = callback_calls + 1
        stop = stop_at_third .and. callback_calls == 3
        if (size(x) /= 1 .or. value >= huge(1.0_dp)) error stop 'bad callback'
    end subroutine callback

    subroutine invalid(res,description)
        type(shgo_result), intent(in) :: res !! invalid result
        character(len=*), intent(in) :: description !! assertion name
        call check(res%status == shgo_invalid .and. res%nfev == 0,description)
    end subroutine invalid

    subroutine check(condition,description)
        logical, intent(in) :: condition !! checked property
        character(len=*), intent(in) :: description !! assertion name
        if (.not. condition) then
            print '(a,a)', 'SHGO HGR FAILED: ', description
            error stop 1
        end if
        checks = checks + 1
    end subroutine check
end program optimize_shgo_hgr_reference
