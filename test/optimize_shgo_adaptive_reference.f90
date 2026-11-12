! SPDX-License-Identifier: MIT
! Adaptive dyadic SHGO graph, strict caps, vertex caching, and local search.
program optimize_shgo_adaptive_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_shgo
    implicit none
    type(shgo_result) :: result, again, uniform
    real(dp) :: box(2,2), bx1(1,2), box3(3,2), shift(2)
    integer :: checks, objective_calls, callback_calls
    logical :: force_stop

    checks = 0
    objective_calls = 0
    callback_calls = 0
    force_stop = .false.
    box(:,1) = -2.0_dp
    box(:,2) = 2.0_dp
    bx1(1,:) = [-2.0_dp,2.0_dp]
    box3(:,1) = -1.0_dp
    box3(:,2) = 1.0_dp
    shift = [0.625_dp, -0.425_dp]

    ! Splitting all leaves reproduces a complete nested uniform grid.
    result = shgo_core(double_well,bx1,iters=2,maxfun=200, &
        adaptive=.true.,adaptive_fraction=1.0_dp)
    call check(result%success .and. result%status == shgo_completed, &
        'full adaptive 1D search completes')
    call check(result%nfev == 5 .and. result%nmesh == 5, '1D shared vertices')
    call check(result%nleaf == 4 .and. result%nrefined == 2,'1D leaf accounting')
    call check(result%nminima == 2 .and. result%funl(1) < result%funl(2), &
        '1D double-well minima')
    uniform = shgo_core(double_well,bx1,iters=2)
    call check(all(abs(result%funl-uniform%funl)<1.0e-12_dp), &
        '1D graph agreement with uniform grid')

    objective_calls = 0
    result = shgo_core(quadratic,box,iters=3,maxfun=500,local_iter=0, &
        adaptive=.true.,adaptive_fraction=1.0_dp)
    call check(result%nfev == 81 .and. result%nmesh == 81, &
        '2D uniform-equivalent vertex count')
    call check(result%nleaf == 64 .and. result%nrefined == 20, &
        '2D uniform-equivalent leaf count')
    call check(objective_calls == result%nfev, '2D objective caching')
    uniform = shgo_core(quadratic,box,iters=3,maxfun=500,local_iter=0)
    call check(abs(result%fun-uniform%fun) < 1.0e-12_dp, &
        '2D uniform-equivalent best objective')
    call check(maxval(abs(result%x-uniform%x)) < 1.0e-12_dp, &
        '2D uniform-equivalent best point')

    ! A quarter of leaves per round explores far fewer unique points.
    result = shgo_core(quadratic,box,iters=3,maxfun=500,local_iter=0,adaptive=.true.)
    call check(result%nit == 3 .and. result%status == shgo_completed, &
        'selective adaptive iterations complete')
    call check(result%nfev >= 14 .and. result%nfev < 81, &
        'selective sampling reduces evaluations')
    call check(result%nmesh == result%nfev .and. result%nleaf == 13 .and. &
        result%nrefined == 3, 'selective mesh and leaf diagnostics')
    call check(result%fun < 0.03_dp, 'refinement finds shifted basin')
    again = shgo_core(quadratic,box,iters=3,maxfun=500,local_iter=0, &
        adaptive=.true.,adaptive_fraction=0.25_dp)
    call check(again%nfev == result%nfev .and. again%nleaf == result%nleaf, &
        'default adaptive fraction is 0.25')
    call check(all(again%x == result%x) .and. again%fun == result%fun, &
        'adaptive trajectory is reproducible')

    ! Exercise a higher-dimensional adaptive mesh and uniform equivalence.
    result = shgo_core(sphere,box3,iters=2,maxfun=600, &
        adaptive=.true.,adaptive_fraction=0.25_dp)
    call check(result%status == shgo_completed .and. result%nit == 2, &
        '3D adaptive mesh completes')
    call check(result%nmesh > 27 .and. result%nmesh < 125, &
        '3D selective vertex count')
    call check(result%nleaf == 22 .and. result%nrefined == 2, &
        '3D adaptive leaf accounting')
    result = shgo_core(sphere,box3,iters=2,maxfun=600, &
        adaptive=.true.,adaptive_fraction=1.0_dp)
    call check(result%nfev == 125 .and. result%nleaf == 64, &
        '3D all-cell refinement matches uniform mesh')
    call check(result%fun < 1.0e-25_dp,'3D origin retained')

    ! No sampling or cell is partially committed when a split cannot fit.
    objective_calls = 0
    result = shgo_core(quadratic,box,iters=3,maxfun=8,adaptive=.true.)
    call check(result%status == shgo_maxfun .and. result%nfev == 1, &
        'insufficient initial mesh budget')
    call check(result%nmesh == 0 .and. result%nit == 0, &
        'no incomplete initial mesh reported')
    result = shgo_core(quadratic,box,iters=3,maxfun=14,adaptive=.true.)
    call check(result%status == shgo_maxfun .and. result%nfev == 14, &
        'exhausted budget after complete adaptive split')
    call check(result%nmesh == 14 .and. result%nleaf == 7 .and. &
        result%nit == 2 .and. result%nrefined == 1, &
        'last complete adaptive split preserved')
    result = shgo_core(quadratic,box,iters=3,max_mesh_points=9, &
        adaptive=.true.)
    call check(result%status == shgo_mesh_limit .and. result%nfev == 9, &
        'adaptive mesh vertex cap honored')
    call check(result%nmesh == 9 .and. result%nrefined == 0, &
        'mesh cap stops before split')
    result = shgo_core(quadratic,box,iters=3,n=20,max_mesh_points=20, &
        adaptive=.true.)
    call check(result%status == shgo_mesh_limit .and. result%nfev == 1, &
        'initial adaptive mesh checked before allocation')
    call check(result%nit == 0 .and. result%nmesh == 0, &
        'initial grid cap leaves no mesh')

    result = shgo_core(quadratic,box,iters=3,maxfun=20,adaptive=.true.)
    call check(result%status == shgo_maxfun .and. result%nfev <= 20, &
        'partial refinement round has strict budget')
    call check(result%nit == 2 .and. result%nmesh == result%nfev .and. &
        result%nrefined >= 1, 'partial round retains only committed cells')

    ! Local callback sees guarded objective, and duplicate minima are merged.
    objective_calls = 0
    result = shgo_core(quadratic,box,iters=3,maxfun=200,adaptive=.true., &
        minimizer=exact_local)
    call check(result%fun < 1.0e-25_dp .and. result%nminima == 1, &
        'adaptive local minimizer finds global optimum')
    call check(result%nlfev == result%nlocal .and. &
        result%nfev == result%nmesh+result%nlfev, &
        'adaptive mesh and local evaluation accounting')
    call check(result%nfev == objective_calls, 'custom minimizer call accounting')

    force_stop = .true.
    callback_calls = 0
    result = shgo_core(quadratic,box,iters=4,adaptive=.true., &
        callback=hook,local_iter=0)
    call check(result%status == shgo_callback_stop .and. result%nit == 1, &
        'adaptive callback terminates after first full mesh')
    call check(callback_calls == 1, 'adaptive callback invoked exactly once')
    force_stop = .false.
    result = shgo_core(quadratic,box,iters=3,adaptive=.true., &
        callback=hook,local_iter=0)
    call check(callback_calls == 4 .and. result%nit == 3, &
        'adaptive callback runs after each full round')

    result = shgo_core(constant_zero,box,adaptive=.true.,f_min=0.0_dp)
    call check(result%status == shgo_target .and. result%nfev == 1 .and. &
        result%nmesh == 0, 'adaptive target stops at initial center')
    result = shgo_core(infinite_objective,box,iters=3,adaptive=.true.)
    call check(result%status == shgo_numerical .and. .not. result%success, &
        'all infeasible adaptive objective values')
    call check(result%nfev == result%nmesh, 'nonfinite samples cached')
    result = shgo_core(partly_nonfinite,box,iters=3,adaptive=.true.)
    call check(result%success .and. result%fun <= 1.0e-12_dp, &
        'partly nonfinite adaptive search')

    call invalid(shgo_core(quadratic,box,adaptive=.true., &
        adaptive_fraction=0.0_dp),'zero adaptive fraction')
    call invalid(shgo_core(quadratic,box,adaptive=.true., &
        adaptive_fraction=-1.0_dp),'negative adaptive fraction')
    call invalid(shgo_core(quadratic,box,adaptive=.true., &
        adaptive_fraction=1.5_dp),'fraction greater than one')
    call invalid(shgo_core(quadratic,box,adaptive=.true., &
        adaptive_fraction=ieee_value(0.0_dp,ieee_quiet_nan)),'NaN adaptive fraction')
    call invalid(shgo_core(quadratic,box,adaptive=.true.,iters=21), &
        'refinement depth guard')

    print '(a,i0)', 'native adaptive SHGO assertions passed: ',checks
contains
    function double_well(x) result(f)
        real(dp), intent(in) :: x(:) !! single coordinate
        real(dp) :: f !! tilted bistable potential
        f = (x(1)**2-1.0_dp)**2+0.1_dp*x(1)
    end function double_well
    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! two-dimensional position
        real(dp) :: f !! shifted quadratic
        objective_calls = objective_calls+1
        f = sum((x-shift)**2)
    end function quadratic
    function sphere(x) result(f)
        real(dp), intent(in) :: x(:) !! any dimension
        real(dp) :: f !! nonnegative paraboloid
        f = sum(x**2)
    end function sphere
    function constant_zero(x) result(f)
        real(dp), intent(in) :: x(:) !! fixed search point
        real(dp) :: f !! constant zero
        f = 0.0_dp
        if (size(x) /= 2) error stop 'unexpected objective dimension'
    end function constant_zero
    function infinite_objective(x) result(f)
        real(dp), intent(in) :: x(:) !! point with no finite score
        real(dp) :: f !! NaN to mimic infeasible sampling
        f = ieee_value(0.0_dp,ieee_quiet_nan)
        if (size(x) /= 2) error stop 'unexpected objective dimension'
    end function infinite_objective
    function partly_nonfinite(x) result(f)
        real(dp), intent(in) :: x(:) !! two-dimensional point
        real(dp) :: f !! domain-limited scalar objective
        if (x(1) < 0.0_dp) then
            f = ieee_value(0.0_dp,ieee_quiet_nan)
        else
            f = sum(x**2)
        end if
    end function partly_nonfinite
    subroutine exact_local(fun,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: fun !! budget guarded objective
        real(dp), intent(in) :: lower(:) !! inclusive lower domain
        real(dp), intent(in) :: upper(:) !! inclusive upper domain
        real(dp), intent(in) :: start(:) !! original candidate
        integer, intent(in) :: maxfev !! residual budget
        real(dp), intent(out) :: x(:) !! local solution
        real(dp), intent(out) :: value !! objective at local solution
        integer, intent(out) :: nfev !! local objective calls
        logical, intent(out) :: success !! local solver success
        x = max(lower,min(upper,shift))
        value = fun(x)
        nfev = 1
        success = maxfev >= 1 .and. size(start) == size(x)
    end subroutine exact_local
    subroutine hook(x,value,stop)
        real(dp), intent(in) :: x(:) !! current minimum coordinate
        real(dp), intent(in) :: value !! objective value
        logical, intent(out) :: stop !! early-stop flag
        callback_calls = callback_calls+1
        stop = force_stop
        if (size(x) /= 2 .or. value >= huge(1.0_dp)) error stop 'invalid callback'
    end subroutine hook
    subroutine invalid(r,description)
        type(shgo_result), intent(in) :: r !! failed invocation
        character(len=*), intent(in) :: description !! assertion name
        call check(r%status == shgo_invalid .and. r%nfev == 0,description)
    end subroutine invalid
    subroutine check(ok,description)
        logical, intent(in) :: ok !! assertion predicate
        character(len=*), intent(in) :: description !! assertion label
        if (.not. ok) then
            print *, 'ADAPTIVE SHGO FAILED: ',description
            error stop 1
        end if
        checks = checks+1
    end subroutine check
end program optimize_shgo_adaptive_reference
