! SPDX-License-Identifier: MIT
! Stand-alone tests of simplicial vertex graph, budgets, callbacks and local adapter.
program optimize_shgo_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_shgo
    implicit none
    type(shgo_result) :: result, again
    real(dp) :: box(2,2), bx1(1,2), bad(2,2), shift(2), v
    integer :: checks, objective_calls, callback_calls, j
    logical :: force_stop

    checks = 0
    objective_calls = 0
    callback_calls = 0
    force_stop = .false.
    box(:,1) = -2.0_dp
    box(:,2) = 2.0_dp
    bx1(1,:) = [-2.0_dp, 2.0_dp]

    ! A one-dimensional double-well has two H_0 births in the five-point mesh.
    result = shgo_core(double_well,bx1,iters=2,maxfun=200)
    call check(result%success .and. result%nit == 2,'dyadic mesh completed')
    call check(result%nmesh == 5 .and. result%nfev == 5,'nested vertex caching')
    call check(result%nminima == 2, 'two local minima in sublevel complex')
    call check(result%funl(1) < result%funl(2), 'sorted minima')
    call check(abs(result%xl(1,1)+1.0_dp) < 1.0e-12_dp, 'negative basin')
    call check(abs(result%xl(1,2)-1.0_dp) < 1.0e-12_dp, 'positive basin')
    call check(result%fun < -0.09_dp, 'best of two basins')
    again = shgo_core(double_well,bx1,iters=2,maxfun=200)
    call check(all(result%xl == again%xl) .and. all(result%funl == again%funl), &
        'fully reproducible without random state')
    call check(result%nlfev == 0 .and. result%nlocal == 0, 'core no-local mode')

    ! Test the topology and correct value ordering of a two-dimensional mesh.
    result = shgo_core(two_bowls,box,iters=2,maxfun=1000)
    call check(result%nfev == 25 .and. result%nmesh == 25, '2D caching')
    call check(result%nminima >= 2,'distinct 2D basin candidates')
    call check(result%funl(1) <= result%funl(2),'2D minima sorting')
    call check(result%fun <= 0.5_dp, '2D sampled best')

    shift = [0.4_dp,-0.7_dp]
    objective_calls = 0
    result = shgo_core(quadratic,box,iters=2,maxfun=2000,minimizer=exact_local)
    call check(result%success .and. result%nit == 2, 'quadratic refinement')
    call check(result%fun < 1.0e-25_dp, 'quadratic local polishing')
    call check(maxval(abs(result%x-shift)) < 1.0e-12_dp, 'quadratic minimizer location')
    call check(result%nminima == 1, 'duplicate local solutions merged')
    call check(result%nlfev == result%nlocal, 'local objective-call accounting')
    call check(objective_calls == result%nfev, 'all objective calls accounted')
    call check(result%nfev == 25 + result%nlfev, 'local and cached sample costs')

    ! Find a known target in the very first center evaluation.
    objective_calls = 0
    result = shgo_core(constant_zero,box,f_min=0.0_dp)
    call check(result%status == shgo_target .and. result%success, 'known optimum')
    call check(result%nfev == 1 .and. result%nit == 0, 'target early stop')
    call check(result%nminima == 1, 'early-stop minima arrays allocated')

    ! No complete initial mesh fits into these budgets, but center is retained.
    objective_calls = 0
    result = shgo_core(quadratic,box,maxfun=1)
    call check(result%status == shgo_maxfun .and. result%success, 'small evaluation budget')
    call check(result%nfev == 1 .and. objective_calls == 1,'budget=1 strict calls')
    call check(result%nit == 0 .and. result%nmesh == 0,'incomplete grid not reported')
    result = shgo_core(quadratic,box,maxfun=8)
    call check(result%nfev == 1 .and. result%status == shgo_maxfun, &
        'mesh is committed atomically')
    result = shgo_core(quadratic,box,maxfun=9)
    call check(result%nit == 1 .and. result%nfev == 9,'exact mesh budget')
    result = shgo_core(quadratic,box,max_mesh_points=8)
    call check(result%status == shgo_mesh_limit .and. result%nfev == 1, &
        'mesh cap before allocation')
    result = shgo_core(quadratic,box,n=20,maxfun=25)
    call check(result%nit == 1 .and. result%nmesh == 25, 'n chooses larger dyadic grid')

    ! A malicious local routine must not break the global objective budget.
    objective_calls = 0
    result = shgo_core(quadratic,box,maxfun=11,minimizer=overcalling_local)
    call check(result%nfev == 11 .and. objective_calls == 11,'local eval hard cap')
    call check(result%nlfev == 2 .and. result%nfev <= 11,'local budget accounting')
    call check(result%success .and. result%status == shgo_maxfun, 'budget exhaustion status')
    call check(result%fun < huge(1.0_dp), 'best sampled point retained')

    ! Other local routines may fail or produce invalid points without corrupting results.
    result = shgo_core(quadratic,box,maxfun=100,minimizer=bad_local)
    call check(result%success .and. result%nit == 1, 'reject nonfinite local minimizer')
    call check(all(result%x >= box(:,1)) .and. all(result%x <= box(:,2)), &
        'reject infeasible local return')
    result = shgo_core(quadratic,box,maxfun=100,minimizer=failed_but_improved)
    call check(result%success .and. result%fun < 1.0e-25_dp, &
        'track best actual evaluation from a failing local solver')
    call check(result%funl(1) == result%fun, 'global result matches sorted minima')

    force_stop = .true.
    callback_calls = 0
    result = shgo_core(quadratic,box,iters=3,callback=hook)
    call check(result%status == shgo_callback_stop .and. result%success, 'callback stop')
    call check(callback_calls == 1 .and. result%nit == 1,'callback per completed mesh')
    force_stop = .false.
    result = shgo_core(quadratic,box,iters=3,callback=hook,local_iter=0)
    call check(result%nit == 3 .and. callback_calls == 4,'callbacks over all iterations')
    call check(result%nfev == 81,'3-level 2D grid and caching')
    result = shgo_core(quadratic,box,iters=3,minimize_every_iter=.false.)
    call check(result%nfev == 81 .and. result%nit == 3, 'delay local searches')

    ! NaNs and infinities are treated as unusable objective samples, not fatal.
    result = shgo_core(infinite_objective,box,maxfun=100)
    call check(result%status == shgo_numerical .and. .not. result%success, &
        'all nonfinite samples')
    call check(result%nfev == 9 .and. result%nminima == 0,'nonfinite accounting')
    result = shgo_core(partly_nonfinite,box,iters=2)
    call check(result%success .and. result%fun < 0.01_dp, 'partly nonfinite objective')

    bad = box
    bad(2,2) = bad(2,1)
    call invalid(shgo_core(quadratic,bad),'zero width')
    call invalid(shgo_core(quadratic,box(:,1:1)),'bad bounds shape')
    bad = box
    bad(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call invalid(shgo_core(quadratic,bad),'NaN bound')
    bad = box
    bad(1,2) = ieee_value(0.0_dp,ieee_positive_inf)
    call invalid(shgo_core(quadratic,bad),'infinite bound')
    call invalid(shgo_core(quadratic,box,iters=0),'invalid iterations')
    call invalid(shgo_core(quadratic,box,maxfun=0),'invalid maxfun')
    call invalid(shgo_core(quadratic,box,n=-1),'invalid n')
    call invalid(shgo_core(quadratic,box,local_iter=-1),'invalid local_iter')
    call invalid(shgo_core(quadratic,box,merge_tol=-1.0_dp),'invalid merge_tol')
    call invalid(shgo_core(quadratic,box,f_tol=-1.0_dp),'invalid f_tol')
    call invalid(shgo_core(quadratic,box, &
        f_tol=ieee_value(0.0_dp,ieee_quiet_nan)),'NaN f_tol')
    result = shgo_core(sphere,spread(box(1,:),1,9))
    call invalid(result,'dimension limit')
    result = shgo_core(sphere,spread(box(1,:),1,8),n=huge(1))
    call check(result%status == shgo_mesh_limit .and. result%nfev == 1, &
        'large requested grid does not overflow 64-bit size')
    bad = box
    bad(:,1) = -huge(1.0_dp)
    bad(:,2) = huge(1.0_dp)
    call invalid(shgo_core(quadratic,bad),'overflowing box widths')
    ! Deliberately huge level request must not overflow 64-bit vertex counts.
    result = shgo_core(quadratic,box,iters=19,maxfun=200)
    call check(result%status == shgo_maxfun .and. result%nit <= 3, &
        'large iteration request capped before unsafe mesh construction')

    print '(a,i0)', 'native SHGO assertions passed: ', checks
contains
    function double_well(x) result(f)
        real(dp), intent(in) :: x(:) !! scalar coordinate
        real(dp) :: f !! tilted double-well potential
        f = (x(1)**2-1.0_dp)**2+0.1_dp*x(1)
    end function double_well
    function two_bowls(x) result(f)
        real(dp), intent(in) :: x(:) !! two coordinates
        real(dp) :: f !! smaller of two offset paraboloids
        f = min((x(1)-1.0_dp)**2+(x(2)-1.0_dp)**2, &
            (x(1)+1.0_dp)**2+(x(2)+1.0_dp)**2+0.25_dp)
    end function two_bowls
    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! two-dimensional point
        real(dp) :: f !! shifted quadratic bowl
        objective_calls = objective_calls + 1
        f = sum((x-shift)**2)
    end function quadratic
    function sphere(x) result(f)
        real(dp), intent(in) :: x(:) !! any-dimensional coordinate
        real(dp) :: f !! squared Euclidean length
        f = sum(x**2)
    end function sphere
    function constant_zero(x) result(f)
        real(dp), intent(in) :: x(:) !! unused test coordinate
        real(dp) :: f !! identically zero objective
        f = 0.0_dp
        if (size(x) /= 2) error stop 'bad test shape'
    end function constant_zero
    function infinite_objective(x) result(f)
        real(dp), intent(in) :: x(:) !! unused test coordinate
        real(dp) :: f !! never finite
        f = ieee_value(0.0_dp,ieee_positive_inf)
        if (size(x) /= 2) error stop 'bad test shape'
    end function infinite_objective
    function partly_nonfinite(x) result(f)
        real(dp), intent(in) :: x(:) !! two coordinates
        real(dp) :: f !! finite only in nonnegative first coordinate
        if (x(1) < 0.0_dp) then
            f = ieee_value(0.0_dp,ieee_quiet_nan)
        else
            f = x(1)**2+x(2)**2
        end if
    end function partly_nonfinite
    subroutine exact_local(fun,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: fun !! objective to be called once
        real(dp), intent(in) :: lower(:) !! lower bounds
        real(dp), intent(in) :: upper(:) !! upper bounds
        real(dp), intent(in) :: start(:) !! input seed
        integer, intent(in) :: maxfev !! permitted evaluations
        real(dp), intent(out) :: x(:) !! analytic minimizer of test quadratic
        real(dp), intent(out) :: value !! evaluated objective
        integer, intent(out) :: nfev !! reported evaluation count
        logical, intent(out) :: success !! true for test quadratic
        x = max(lower,min(upper,shift))
        value = fun(x)
        nfev = 1
        success = maxfev >= 1 .and. size(start) == size(x)
    end subroutine exact_local
    subroutine overcalling_local(fun,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: fun !! deliberately called beyond maxfev
        real(dp), intent(in) :: lower(:) !! lower bounds
        real(dp), intent(in) :: upper(:) !! upper bounds
        real(dp), intent(in) :: start(:) !! fixed test point
        integer, intent(in) :: maxfev !! available calls
        real(dp), intent(out) :: x(:) !! unchanged point
        real(dp), intent(out) :: value !! first actual evaluation
        integer, intent(out) :: nfev !! attempted number of calls
        logical, intent(out) :: success !! return not suitable due budget misuse
        x = max(lower,min(upper,start))
        value = fun(x)
        do j = 1, maxfev+3
            v = fun(x)
        end do
        nfev = maxfev+4
        success = .false.
    end subroutine overcalling_local
    subroutine bad_local(fun,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: fun !! unused objective
        real(dp), intent(in) :: lower(:) !! lower bounds
        real(dp), intent(in) :: upper(:) !! upper bounds
        real(dp), intent(in) :: start(:) !! starting candidate
        integer, intent(in) :: maxfev !! unused evaluation limit
        real(dp), intent(out) :: x(:) !! out of box on purpose
        real(dp), intent(out) :: value !! invalid value
        integer, intent(out) :: nfev !! zero evaluations
        logical, intent(out) :: success !! simulate erroneous success
        x = upper + 1.0_dp
        value = ieee_value(0.0_dp,ieee_quiet_nan)
        nfev = 0
        success = .true.
        if (size(lower)+size(start)+maxfev < 0) error stop 'unreachable'
    end subroutine bad_local
    subroutine failed_but_improved(fun,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: fun !! budget-guarded objective
        real(dp), intent(in) :: lower(:) !! box lower bounds
        real(dp), intent(in) :: upper(:) !! box upper bounds
        real(dp), intent(in) :: start(:) !! initial local point
        integer, intent(in) :: maxfev !! available objective calls
        real(dp), intent(out) :: x(:) !! point deliberately left unchanged
        real(dp), intent(out) :: value !! failed minimizer value
        integer, intent(out) :: nfev !! actual calls made
        logical, intent(out) :: success !! deliberately false
        real(dp) :: best
        x = start
        best = fun(max(lower,min(upper,shift)))
        value = huge(1.0_dp)
        nfev = 1
        success = .false.
        if (best < 0.0_dp .or. maxfev < 1) error stop 'unexpected test arguments'
    end subroutine failed_but_improved

    subroutine hook(x,value,stop)
        real(dp), intent(in) :: x(:) !! best current point
        real(dp), intent(in) :: value !! best objective
        logical, intent(out) :: stop !! caller-requested exit
        callback_calls = callback_calls + 1
        stop = force_stop
        if (size(x) /= 2 .or. value >= huge(1.0_dp)) error stop 'bad hook'
    end subroutine hook
    subroutine invalid(r,description)
        type(shgo_result), intent(in) :: r !! expected invalid result
        character(len=*), intent(in) :: description !! assertion name
        call check(r%status == shgo_invalid .and. r%nfev == 0,description)
    end subroutine invalid
    subroutine check(ok,description)
        logical, intent(in) :: ok !! test outcome
        character(len=*), intent(in) :: description !! failing assertion label
        if (.not. ok) then
            print *, 'SHGO FAILED: ', description
            error stop 1
        end if
        checks = checks+1
    end subroutine check
end program optimize_shgo_reference
