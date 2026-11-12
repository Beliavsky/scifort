! SPDX-License-Identifier: MIT
! Native regression tests for conforming 3D adaptive SHGO tetrahedral meshes.
program optimize_shgo_conforming_3d_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_shgo
    implicit none
    type(shgo_result) :: result, previous
    real(dp) :: box(3,2), box4(4,2)
    integer :: checks, objective_calls, callback_calls
    logical :: stop_now

    checks = 0
    objective_calls = 0
    callback_calls = 0
    stop_now = .false.
    box(:,1) = -1.0_dp
    box(:,2) = 1.0_dp
    box4(:,1) = -1.0_dp
    box4(:,2) = 1.0_dp

    result = shgo_core(quadratic,box,iters=1,local_iter=0, &
        adaptive=.true.,conforming=.true.)
    call check(result%status == shgo_completed .and. result%success, &
        'initial tetrahedral mesh completes')
    call check(result%nmesh == 71 .and. result%nfev == 71, &
        'initial 27 corners plus 36 shared face centers and 8 cell centers')
    call check(result%nleaf == 8 .and. result%ntetrahedra == 192, &
        'each initial cube has 24 oriented tetrahedra')
    call check(result%conforming_valid .and. result%ntriangles == 0, &
        'tetrahedral mesh validated, planar diagnostic unused')
    call check(result%hgr_estimate == 1 .and. result%nminima >= 1, &
        'one sampled sphere minimum')

    objective_calls = 0
    result = shgo_core(quadratic,box,iters=2,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.01_dp)
    call check(result%status == shgo_completed .and. result%nit == 2, &
        'single adaptive split completes')
    call check(result%nmesh == 127 .and. result%nfev == 127, &
        'split counts unique face centers exactly')
    call check(result%nleaf == 15 .and. result%nrefined == 1, &
        'seven net new leaf cells in one octree split')
    call check(result%ntetrahedra == 414 .and. result%conforming_valid, &
        'coarse-fine internal faces share triangles')
    call check(objective_calls == result%nfev, &
        'no duplicate objective calls for shared face/edge vertices')
    call check(result%fun < 0.02_dp, 'adaptive split improves incumbent')
    previous = shgo_core(quadratic,box,iters=2,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.01_dp)
    call check(previous%fun == result%fun .and. all(previous%x == result%x), &
        'fully reproducible sampled best point')
    call check(previous%nmesh == result%nmesh .and. &
        previous%ntetrahedra == result%ntetrahedra, &
        'fully reproducible tetrahedron/vertex counts')

    result = shgo_core(quadratic,box,iters=4,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.25_dp)
    call check(result%status == shgo_completed .and. result%conforming_valid, &
        'uneven multi-level octree interfaces validated')
    call check(result%nmesh == 1145 .and. result%nleaf == 176, &
        'unique vertices and leaf topology across four rounds')
    call check(result%ntetrahedra == 5208 .and. result%nrefined == 24, &
        'multi-level tetrahedralization count')
    call check(result%fun < 0.001_dp, 'multi-level refinement improves sample')

    result = shgo_core(quadratic,box,iters=7,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.01_dp)
    call check(result%status == shgo_completed .and. result%conforming_valid, &
        'deeply unbalanced three-dimensional mesh still conforming')
    call check(result%nrefined == 6 .and. result%nleaf == 50, &
        'deeply unbalanced octree selects only one cube per round')

    result = shgo_core(quadratic,box,iters=2,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=1.0_dp)
    call check(result%status == shgo_completed .and. result%conforming_valid, &
        'whole octree refinement validates')
    call check(result%nmesh == 429 .and. result%nleaf == 64, &
        'fully refined unique vertex count')
    call check(result%ntetrahedra == 1536 .and. result%nrefined == 8, &
        '64 uniform cubes triangulated with 24 tets each')

    result = shgo_core(quadratic,box,iters=3,local_iter=0, &
        adaptive=.true.,conforming=.true.,maxfun=70)
    call check(result%status == shgo_maxfun .and. result%nfev == 1, &
        'initial complete mesh budget includes face centers')
    call check(result%nmesh == 0 .and. result%nit == 0, &
        'initial incomplete mesh not exposed')
    result = shgo_core(quadratic,box,iters=3,local_iter=0, &
        adaptive=.true.,conforming=.true.,max_mesh_points=70)
    call check(result%status == shgo_mesh_limit .and. result%nfev == 1, &
        'initial exact vertex memory preflight')
    result = shgo_core(quadratic,box,iters=2,local_iter=0, &
        adaptive=.true.,conforming=.true.,maxfun=126,adaptive_fraction=0.01_dp)
    call check(result%status == shgo_maxfun .and. result%nfev == 71, &
        'split is rejected atomically one call short')
    call check(result%nrefined == 0 .and. result%ntetrahedra == 192 .and. &
        result%conforming_valid,'old validated mesh retained on rejected split')
    result = shgo_core(quadratic,box,iters=2,local_iter=0, &
        adaptive=.true.,conforming=.true.,maxfun=127,adaptive_fraction=0.01_dp)
    call check(result%status == shgo_completed .and. result%nfev == 127, &
        'split may commit with exactly sufficient budget')
    call check(result%conforming_valid .and. result%ntetrahedra == 414, &
        'exact-budget split is triangulated and validated')
    result = shgo_core(quadratic,box,iters=2,local_iter=0, &
        adaptive=.true.,conforming=.true.,maxfun=127,adaptive_fraction=1.0_dp)
    call check(result%status == shgo_maxfun .and. result%nrefined == 1 .and. &
        result%nit == 1,'partially selected batch does not increment completed rounds')
    call check(.not. result%conforming_valid .and. result%ntetrahedra == 0, &
        'partially committed round invalidates prior triangulation report')
    result = shgo_core(quadratic,box,iters=2,local_iter=0, &
        adaptive=.true.,conforming=.true.,max_mesh_points=126)
    call check(result%status == shgo_mesh_limit .and. result%nfev == 71, &
        'cell vertex memory cap rejects oversized split atomically')

    result = shgo_core(flat,box,iters=5,local_iter=0, &
        adaptive=.true.,conforming=.true.,hgr_stability_iters=2)
    call check(result%status == shgo_hgr_stable .and. result%nit == 3, &
        'flat mesh stabilization uses three complete rounds')
    call check(result%hgr_estimate == 1 .and. result%hgr_streak == 2, &
        'connected flat tetrahedral complex is one plateau')
    call check(result%nminima == 1 .and. result%conforming_valid, &
        'one flat plateau is one representative minimum')
    result = shgo_core(two_basins,box,iters=3,local_iter=4, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=1.0_dp)
    call check(result%conforming_valid .and. result%hgr_estimate == 2, &
        'two flat disconnected basin components survive tetrahedral graph')
    call check(result%nminima >= 2 .and. all(result%funl == 0.0_dp), &
        'two separate flat local minimum plateaus retained')
    result = shgo_core(nonfinite,box,iters=2,local_iter=0, &
        adaptive=.true.,conforming=.true.)
    call check(result%status == shgo_numerical .and. .not. result%success, &
        'all-nonfinite objective cannot produce false minimizer')
    call check(result%conforming_valid .and. result%hgr_estimate == 0, &
        'geometry validity independent of objective feasibility')

    result = shgo_core(quadratic,box,iters=2,adaptive=.true., &
        conforming=.true.,minimizer=polish_exact,maxfun=400)
    call check(result%success .and. result%fun <= 1.0e-26_dp, &
        'bounded custom polishing reaches exact minimum')
    call check(result%nlfev == result%nlocal .and. &
        result%nfev == result%nmesh+result%nlfev, &
        'polishing objective budget reported separately')
    stop_now = .true.
    callback_calls = 0
    result = shgo_core(quadratic,box,iters=5,adaptive=.true., &
        conforming=.true.,local_iter=0,callback=on_iteration)
    call check(result%status == shgo_callback_stop .and. result%nit == 1, &
        'callback can stop completed tetrahedral round')
    call check(callback_calls == 1 .and. result%conforming_valid, &
        'callback sees validated mesh')
    stop_now = .false.
    result = shgo_core(flat,box,iters=3,adaptive=.true., &
        conforming=.true.,f_min=0.0_dp)
    call check(result%status == shgo_target .and. result%nfev == 1 .and. &
        result%nmesh == 0,'known objective target stops before first mesh')

    call invalid(shgo_core(quadratic,box,conforming=.true.), &
        'conforming only available in adaptive mode')
    call invalid(shgo_core(quadratic,box4,adaptive=.true.,conforming=.true.), &
        'four-dimensional conforming mode is rejected')
    call invalid(shgo_core(quadratic,box,adaptive=.true.,conforming=.true.,iters=21), &
        'dyadic level overflow refused')

    print '(a,i0)', 'native conforming 3D SHGO assertions passed: ',checks
contains
    function quadratic(x) result(value)
        real(dp), intent(in) :: x(:) !! three-dimensional bounded point
        real(dp) :: value !! shifted strictly convex quadratic
        objective_calls = objective_calls+1
        value = sum((x-[0.3_dp,-0.2_dp,0.15_dp])**2)
    end function quadratic
    function flat(x) result(value)
        real(dp), intent(in) :: x(:) !! bounded point
        real(dp) :: value !! constant objective
        value = 0.0_dp*sum(x)
    end function flat
    function two_basins(x) result(value)
        real(dp), intent(in) :: x(:) !! bounded point
        real(dp) :: value !! disconnected flat low-value regions
        value = 1.0_dp
        if (abs(x(1)) >= 0.5_dp) value = 0.0_dp
    end function two_basins
    function nonfinite(x) result(value)
        real(dp), intent(in) :: x(:) !! bounded point
        real(dp) :: value !! deliberately invalid function everywhere
        value = ieee_value(0.0_dp,ieee_quiet_nan)
        if (size(x) /= 3) value = 0.0_dp
    end function nonfinite
    subroutine polish_exact(objective,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: objective !! bounded budget-guarded objective
        real(dp), intent(in) :: lower(:) !! lower bounds
        real(dp), intent(in) :: upper(:) !! upper bounds
        real(dp), intent(in) :: start(:) !! proposed initial candidate
        integer, intent(in) :: maxfev !! remaining objective calls
        real(dp), intent(out) :: x(:) !! optimized point
        real(dp), intent(out) :: value !! objective at optimized point
        integer, intent(out) :: nfev !! callback evaluation count
        logical, intent(out) :: success !! local optimization succeeded
        x = min(upper,max(lower,[0.3_dp,-0.2_dp,0.15_dp]))
        value = objective(x)
        nfev = 1
        success = maxfev >= 1 .and. size(start) == size(x)
    end subroutine polish_exact
    subroutine on_iteration(x,value,stop)
        real(dp), intent(in) :: x(:) !! best current point
        real(dp), intent(in) :: value !! best objective value
        logical, intent(out) :: stop !! early stop requested
        callback_calls = callback_calls+1
        stop = stop_now
        if (size(x) /= 3 .or. value > 1.0e5_dp) stop = .true.
    end subroutine on_iteration
    subroutine invalid(r,description)
        type(shgo_result), intent(in) :: r !! expected invalid result
        character(len=*), intent(in) :: description !! test case
        call check(r%status == shgo_invalid .and. r%nfev == 0,description)
    end subroutine invalid
    subroutine check(condition,description)
        logical, intent(in) :: condition !! asserted property
        character(len=*), intent(in) :: description !! assertion description
        if (.not. condition) then
            print '(a,a)', 'CONFORMING 3D SHGO FAILED: ',description
            error stop 1
        end if
        checks = checks+1
    end subroutine check
end program optimize_shgo_conforming_3d_reference
