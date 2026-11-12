! SPDX-License-Identifier: MIT
! Conforming 2D adaptive SHGO topology, budgets, plateaus and stopping.
program optimize_shgo_conforming_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_shgo
    implicit none
    type(shgo_result) :: r, previous
    real(dp) :: box(2,2), box4(4,2), interval(1,2)
    integer :: checks, objective_calls, callback_calls
    logical :: request_stop

    checks = 0
    objective_calls = 0
    callback_calls = 0
    request_stop = .false.
    box(:,1) = -1.0_dp
    box(:,2) = 1.0_dp
    box4(:,1) = -1.0_dp
    box4(:,2) = 1.0_dp
    interval(1,:) = [-1.0_dp,1.0_dp]

    r = shgo_core(quadratic,box,iters=1,local_iter=0,adaptive=.true.,conforming=.true.)
    call check(r%success .and. r%status == shgo_completed,'initial mesh success')
    call check(r%nfev == 13 .and. r%nmesh == 13,'initial 9 corners + 4 centers')
    call check(r%nleaf == 4 .and. r%ntriangles == 16,'initial four fans')
    call check(r%conforming_valid .and. r%hgr_estimate == 1,'validated initial complex')

    objective_calls = 0
    r = shgo_core(quadratic,box,iters=2,local_iter=0,adaptive=.true., &
        adaptive_fraction=0.25_dp,conforming=.true.)
    call check(r%status == shgo_completed .and. r%nit == 2,'one selective split')
    call check(r%nfev == 21 .and. r%nmesh == 21,'refinement caches old center')
    call check(r%nleaf == 7 .and. r%nrefined == 1,'four children replace one cell')
    call check(r%ntriangles == 30 .and. r%conforming_valid,'coarse-fine face triangles')
    call check(objective_calls == r%nfev,'no repeated evaluations')
    call check(r%nminima == 1 .and. r%fun < 0.01_dp,'graph minimum detected')
    previous = shgo_core(quadratic,box,iters=2,local_iter=0,adaptive=.true., &
        adaptive_fraction=0.25_dp,conforming=.true.)
    call check(previous%fun == r%fun .and. all(previous%x == r%x), &
        'deterministic coordinates and objective')
    call check(previous%ntriangles == r%ntriangles .and. previous%nfev == r%nfev, &
        'deterministic topology and evaluation count')

    r = shgo_core(quadratic,box,iters=4,local_iter=0,adaptive=.true., &
        adaptive_fraction=0.25_dp,conforming=.true.)
    call check(r%conforming_valid .and. r%nleaf == 25 .and. r%nrefined == 7, &
        'multiple asymmetric refinement rounds remain conforming')
    call check(r%ntriangles == 115 .and. r%nfev == 64, &
        'multiple hanging face vertices included')
    call check(r%fun < 0.001_dp,'selected refinement improves basin sample')
    r = shgo_core(quadratic,box,iters=7,local_iter=0,adaptive=.true., &
        adaptive_fraction=0.01_dp,conforming=.true.)
    call check(r%status == shgo_completed .and. r%conforming_valid, &
        'deeply uneven 2D mesh has conforming faces')
    call check(r%nrefined == 6 .and. r%nleaf == 22, &
        'small adaptive fraction splits only one leaf per round')

    r = shgo_core(quadratic,box,iters=2,local_iter=0,adaptive=.true., &
        adaptive_fraction=1.0_dp,conforming=.true.)
    call check(r%conforming_valid .and. r%nleaf == 16,'fully refined cells')
    call check(r%nmesh == 41 .and. r%ntriangles == 64, &
        'uniform grid plus 16 face-center fan vertices')
    r = shgo_core(quadratic,box,iters=3,local_iter=0,adaptive=.true., &
        adaptive_fraction=1.0_dp,conforming=.true.)
    call check(r%conforming_valid .and. r%nmesh == 145 .and. r%ntriangles == 256, &
        'fully refined second level remains conforming')
    call check(r%fun < 0.01_dp,'full refinement finds shifted minimum')

    r = shgo_core(quadratic,box,iters=3,local_iter=0,adaptive=.true., &
        maxfun=12,conforming=.true.)
    call check(r%status == shgo_maxfun .and. r%nmesh == 0, &
        'initial budget includes cell centers')
    call check(r%nfev == 1 .and. r%nit == 0,'initial budget preflight atomic')
    r = shgo_core(quadratic,box,iters=3,local_iter=0,adaptive=.true., &
        max_mesh_points=12,conforming=.true.)
    call check(r%status == shgo_mesh_limit .and. r%nmesh == 0, &
        'initial vertex cap includes centers')
    r = shgo_core(quadratic,box,iters=3,local_iter=0,adaptive=.true., &
        maxfun=20,conforming=.true.)
    call check(r%status == shgo_maxfun .and. r%nfev == 13 .and. r%nit == 1, &
        'incomplete split rejected without any extra objective calls')
    call check(r%nrefined == 0 .and. r%nmesh == 13,'no partially sampled split')
    r = shgo_core(quadratic,box,iters=3,local_iter=0,adaptive=.true., &
        maxfun=21,conforming=.true.)
    call check(r%status == shgo_maxfun .and. r%nfev == 21 .and. r%nit == 2, &
        'split commits exactly at evaluation budget')
    call check(r%nmesh == 21 .and. r%nrefined == 1,'last complete split retained')
    r = shgo_core(quadratic,box,iters=3,local_iter=0,adaptive=.true., &
        adaptive_fraction=1.0_dp,maxfun=21,conforming=.true.)
    call check(r%status == shgo_maxfun .and. r%nrefined == 1 .and. r%nit == 1, &
        'partial selection round retains committed cells without incrementing nit')
    call check(.not. r%conforming_valid .and. r%ntriangles == 0, &
        'incomplete round correctly marks topology as not revalidated')
    r = shgo_core(quadratic,box,iters=3,local_iter=0,adaptive=.true., &
        max_mesh_points=20,conforming=.true.)
    call check(r%status == shgo_mesh_limit .and. r%nmesh == 13, &
        'cell refinement enforces vertex memory cap')

    r = shgo_core(flat,box,iters=5,local_iter=0,adaptive=.true., &
        conforming=.true.,hgr_stability_iters=2)
    call check(r%status == shgo_hgr_stable .and. r%nit == 3, &
        'flat function triggers stable graph rank only after three rounds')
    call check(r%hgr_estimate == 1 .and. r%hgr_streak == 2, &
        'connected equal-value plateau gives one birth')
    call check(r%nminima == 1 .and. r%conforming_valid, &
        'flat objective produces one local candidate')
    r = shgo_core(two_plateaus,box,iters=3,local_iter=3,adaptive=.true., &
        conforming=.true.,adaptive_fraction=1.0_dp)
    call check(r%conforming_valid .and. r%hgr_estimate == 2, &
        'two disconnected flat basins remain distinct')
    call check(r%nminima >= 2 .and. all(r%funl == 0.0_dp), &
        'both flat basins represented in sorted result pool')
    r = shgo_core(nonfinite,box,iters=2,local_iter=0,adaptive=.true., &
        conforming=.true.)
    call check(r%status == shgo_numerical .and. .not. r%success, &
        'nonfinite sampled function produces no false minimum')
    call check(r%conforming_valid .and. r%hgr_estimate == 0, &
        'geometric validity independent of objective feasibility')

    r = shgo_core(quadratic,box,iters=3,maxfun=400,adaptive=.true., &
        conforming=.true.,minimizer=polish_exact)
    call check(r%success .and. r%fun <= 1.0e-26_dp,'local polishing reaches exact optimum')
    call check(r%nlfev == r%nlocal .and. r%nfev == r%nmesh+r%nlfev, &
        'mesh/local objective budgets counted independently')
    request_stop = .true.
    callback_calls = 0
    r = shgo_core(quadratic,box,iters=5,adaptive=.true., &
        conforming=.true.,local_iter=0,callback=iteration_callback)
    call check(r%status == shgo_callback_stop .and. r%nit == 1, &
        'iteration callback terminates completed conforming mesh')
    call check(callback_calls == 1 .and. r%conforming_valid, &
        'callback invoked on completed and validated triangulation')
    request_stop = .false.
    r = shgo_core(flat,box,iters=3,adaptive=.true.,conforming=.true.,f_min=0.0_dp)
    call check(r%status == shgo_target .and. r%nfev == 1 .and. r%nmesh == 0, &
        'known target stops before mesh construction')

    call invalid(shgo_core(quadratic,box,conforming=.true.), &
        'conforming without adaptive is invalid')
    call invalid(shgo_core(quadratic,interval,adaptive=.true.,conforming=.true.), &
        'conforming requires 2D or 3D bounds')
    call invalid(shgo_core(quadratic,box4,adaptive=.true.,conforming=.true.), &
        'conforming rejects 4D without silently falling back')
    call invalid(shgo_core(quadratic,box,adaptive=.true.,conforming=.true.,iters=21), &
        'maximum dyadic refinement level guarded')
    print '(a,i0)', 'native conforming SHGO assertions passed: ', checks
contains
    function quadratic(x) result(value)
        real(dp), intent(in) :: x(:) !! point in the search bounds
        real(dp) :: value !! squared distance to known minimum
        objective_calls = objective_calls+1
        value = sum((x-[0.3_dp,-0.2_dp])**2)
    end function quadratic
    function flat(x) result(value)
        real(dp), intent(in) :: x(:) !! feasible point
        real(dp) :: value !! constant objective
        value = 0.0_dp*sum(x)
    end function flat
    function two_plateaus(x) result(value)
        real(dp), intent(in) :: x(:) !! point in the plane
        real(dp) :: value !! flat minima at left and right
        value = 1.0_dp
        if (abs(x(1)) > 0.25_dp) value = 0.0_dp
    end function two_plateaus
    function nonfinite(x) result(value)
        real(dp), intent(in) :: x(:) !! point in the plane
        real(dp) :: value !! deliberately nonfinite everywhere
        value = ieee_value(0.0_dp,ieee_quiet_nan)
        if (size(x) /= 2) value = 0.0_dp
    end function nonfinite
    subroutine polish_exact(objective,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: objective !! budget-guarded objective
        real(dp), intent(in) :: lower(:) !! lower bounds for search
        real(dp), intent(in) :: upper(:) !! upper bounds for search
        real(dp), intent(in) :: start(:) !! input local candidate
        integer, intent(in) :: maxfev !! maximum remaining calls
        real(dp), intent(out) :: x(:) !! final point
        real(dp), intent(out) :: value !! final objective
        integer, intent(out) :: nfev !! objective calls performed
        logical, intent(out) :: success !! solver reported convergence
        x = min(upper,max(lower,[0.3_dp,-0.2_dp]))
        value = objective(x)
        nfev = 1
        success = maxfev >= 1 .and. size(start) == size(x)
    end subroutine polish_exact
    subroutine iteration_callback(x,value,stop)
        real(dp), intent(in) :: x(:) !! incumbent best point
        real(dp), intent(in) :: value !! incumbent objective
        logical, intent(out) :: stop !! callback requests termination
        callback_calls = callback_calls+1
        stop = request_stop
        if (size(x) /= 2 .or. value > 1.0e5_dp) stop = .true.
    end subroutine iteration_callback
    subroutine invalid(result,description)
        type(shgo_result), intent(in) :: result !! expected invalid result
        character(len=*), intent(in) :: description !! tested parameter rule
        call check(result%status == shgo_invalid .and. result%nfev == 0, description)
    end subroutine invalid
    subroutine check(condition,description)
        logical, intent(in) :: condition !! success condition
        character(len=*), intent(in) :: description !! name of the assertion
        if (.not. condition) then
            print '(a,a)', 'CONFORMING SHGO FAILED: ',description
            error stop 1
        end if
        checks = checks+1
    end subroutine check
end program optimize_shgo_conforming_reference
