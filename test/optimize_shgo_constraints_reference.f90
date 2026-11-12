! SPDX-License-Identifier: MIT
! Standalone checked/optimized SHGO constraint regression tests.
program optimize_shgo_constraints_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result, shgo_completed, &
        shgo_infeasible, shgo_invalid, shgo_maxfun
    implicit none
    real(dp) :: bounds(2,2), bounds3(3,2), a(1,2), two_constraints(2,2)
    type(shgo_result) :: res
    integer :: checks

    checks = 0
    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp
    a(1,:) = [1.0_dp,1.0_dp]

    res = shgo_core(quadratic,bounds,iters=1,maxfun=100,local_iter=0, &
        linear_a=a,linear_lb=[0.9_dp],linear_ub=[huge(1.0_dp)])
    call assert(res%success .and. res%feasible, 'linear lower inequality feasible')
    call assert(res%status == shgo_completed, 'linear lower inequality completes')
    call assert(res%maxcv <= 1.0e-8_dp, 'linear inequality violation')
    call assert(res%nfev < res%nmesh, 'infeasible samples skip objective')
    call assert(res%ncev >= res%nmesh, 'count constraint evaluations')
    call assert(res%fun == quadratic(res%x), 'raw untransformed objective returned')
    call assert(all(matmul(a,res%x) >= 0.9_dp-1.0e-8_dp), 'returned x satisfies bound')
    call assert(res%nminima >= 1, 'feasible minima pool')
    call assert(all(res%funl <= 1.0_dp), 'raw minima, not graph scores')

    res = shgo_core(quadratic,bounds,iters=1,maxfun=100,local_iter=1, &
        minimizer=project_inequality, linear_a=a,linear_lb=[0.9_dp], &
        linear_ub=[huge(1.0_dp)])
    call assert(res%success, 'constrained local polishing')
    call assert(abs(res%fun - 0.18_dp) < 1.0e-12_dp, 'projected quadratic optimum')
    call assert(maxval(abs(res%x - [0.5_dp,0.4_dp])) < 1.0e-12_dp, 'projected x')
    call assert(res%nlfev >= 1, 'local objective counted')
    call assert(res%nfev <= 100, 'local objective budget respected')
    call assert(res%nlocal >= 1, 'polisher invoked')
    call assert(res%nminima >= 1 .and. abs(res%funl(1)-0.18_dp)<1.0e-12_dp, &
        'local result in sorted minima')

    ! x+y=0.7 contains no first-level mesh vertex. The feasibility-first
    ! graph must nevertheless produce a start for equality restoration.
    res = shgo_core(shifted,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=nonlinear_sum,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp])
    call assert(.not. res%success, 'no feasible equality mesh sample')
    call assert(res%status == shgo_infeasible, 'infeasible status for equality mesh')
    call assert(res%nfev == 0, 'infeasible mesh skips user objective')
    call assert(res%nminima == 0, 'no infeasible minima accepted')
    call assert(res%ncev >= 9, 'constraint checked throughout mesh')
    call assert(res%maxcv > 0.0_dp, 'best infeasible violation reported')

    res = shgo_core(shifted,bounds,iters=1,maxfun=100,local_iter=2, &
        minimizer=project_equality, nonlinear=nonlinear_sum, &
        nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp])
    call assert(res%status == shgo_completed .and. res%success, 'equality local restoration')
    call assert(res%maxcv < 1.0e-12_dp, 'equality residual')
    call assert(res%fun < 1.0e-12_dp, 'nonlinear equality optimum')
    call assert(abs(sum(res%x)-0.7_dp)<1.0e-12_dp, 'nonlinear equality x')
    call assert(res%nfev >= 1 .and. res%nfev <= 100, 'equality objective budget')

    res = shgo_core(shifted,bounds,iters=1,maxfun=100,local_iter=2, &
        minimizer=project_equality, linear_a=a,linear_lb=[0.7_dp], &
        linear_ub=[0.7_dp])
    call assert(res%success .and. res%fun<1.0e-12_dp, 'linear equality restoration')
    call assert(res%maxcv<1.0e-12_dp, 'linear equality residual')

    ! Upper inequality and equality can coexist in the same call.
    two_constraints(1,:) = [1.0_dp,1.0_dp]
    two_constraints(2,:) = [1.0_dp,-1.0_dp]
    res = shgo_core(shifted,bounds,iters=2,maxfun=200,local_iter=1, &
        minimizer=project_equality,linear_a=two_constraints, &
        linear_lb=[0.7_dp,-huge(1.0_dp)], &
        linear_ub=[0.7_dp,0.8_dp],nonlinear=nonlinear_sum, &
        nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp])
    call assert(res%success .and. res%fun<1.0e-12_dp, 'mixed constraints')
    call assert(res%maxcv<1.0e-12_dp, 'mixed constraints feasible')

    res = shgo_core(shifted,bounds,iters=1,maxfun=100,local_iter=1, &
        minimizer=invalid_minimizer, nonlinear=nonlinear_sum, &
        nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp])
    call assert(.not. res%success, 'reject infeasible claimed local solution')
    call assert(res%nminima == 0, 'do not store invalid local solution')

    res = shgo_core(quadratic,bounds,iters=1,maxfun=100,local_iter=1, &
        minimizer=fake_feasible_minimizer,linear_a=a, &
        linear_lb=[0.9_dp],linear_ub=[huge(1.0_dp)])
    call assert(res%success, 'feasible problem despite fake successful minimizer')
    call assert(res%fun >= 0.0_dp, 'fake feasible minimizer value must be re-evaluated')
    call assert(res%fun == quadratic(res%x), 'returned objective always evaluated')
    res = shgo_core(quadratic,bounds,iters=1,maxfun=100,local_iter=1, &
        minimizer=nan_minimizer,linear_a=a, &
        linear_lb=[0.9_dp],linear_ub=[huge(1.0_dp)])
    call assert(res%success, 'nan local result leaves valid mesh result intact')
    call assert(res%fun >= 0.0_dp, 'nan local result cannot pollute objective')

    res = shgo_core(quadratic,bounds,iters=2,maxfun=100,local_iter=0, &
        linear_a=a,linear_lb=[3.0_dp],linear_ub=[huge(1.0_dp)])
    call assert(res%status == shgo_infeasible, 'impossible constraint status')
    call assert(.not. res%success .and. .not. res%feasible, 'impossible constraints fail')
    call assert(res%nminima==0 .and. res%nfev==0, 'no infeasible objective evaluated')
    call assert(abs(res%maxcv-1.0_dp)<1.0e-12_dp, 'smallest violation kept')
    call assert(maxval(abs(res%x-[1.0_dp,1.0_dp]))<1.0e-12_dp, 'least violating point')

    res = shgo_core(shifted,bounds,iters=2,maxfun=400,local_iter=0, &
        adaptive=.true., conforming=.true., &
        nonlinear=nonlinear_sum, nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp])
    call assert(res%status==shgo_infeasible, 'conforming 2d infeasible status')
    call assert(res%conforming_valid .and. res%ntriangles>0, 'conforming mesh valid with constraints')
    call assert(res%maxcv>=0.0_dp, 'conforming infeasible violation')

    bounds3(:,1)=0.0_dp
    bounds3(:,2)=1.0_dp
    res = shgo_core(sphere3,bounds3,iters=1,maxfun=300,local_iter=0, &
        adaptive=.true.,conforming=.true.,nonlinear=nonlinear_sum3, &
        nonlinear_lb=[1.5_dp],nonlinear_ub=[1.5_dp])
    call assert(res%conforming_valid .and. res%ntetrahedra==192, 'conforming 3d constraints')
    call assert(res%success .and. res%feasible, 'conforming 3d equality feasible')
    call assert(res%maxcv<=1.0e-8_dp, 'conforming 3d equality violation')

    res = shgo_core(quadratic,bounds,iters=1,maxfun=1,local_iter=0, &
        linear_a=a,linear_lb=[0.9_dp],linear_ub=[huge(1.0_dp)])
    call assert(res%nfev<=1, 'budget cap constrained')
    call assert(res%status==shgo_maxfun, 'budget status with initial feasible sample')
    call assert(res%success, 'feasible center retained at budget limit')

    res = shgo_core(quadratic,bounds, linear_a=a,linear_lb=[0.7_dp])
    call assert(res%status==shgo_invalid,'missing linear_ub rejected')
    res = shgo_core(quadratic,bounds, nonlinear=nonlinear_sum,nonlinear_lb=[0.7_dp])
    call assert(res%status==shgo_invalid,'missing nonlinear_ub rejected')
    res = shgo_core(quadratic,bounds, linear_a=a, &
        linear_lb=[0.7_dp,0.8_dp],linear_ub=[1.0_dp])
    call assert(res%status==shgo_invalid,'linear length mismatch rejected')
    res = shgo_core(quadratic,bounds, nonlinear=nonlinear_sum, &
        nonlinear_lb=[0.0_dp],nonlinear_ub=[-1.0_dp])
    call assert(res%status==shgo_invalid,'reversed constraint bounds rejected')
    res = shgo_core(quadratic,bounds, nonlinear=nonlinear_sum, &
        nonlinear_lb=[0.0_dp],nonlinear_ub=[1.0_dp],constraint_tol=-1.0_dp)
    call assert(res%status==shgo_invalid,'negative tolerance rejected')
    res = shgo_core(quadratic,bounds, nonlinear=nonlinear_sum, &
        nonlinear_lb=[ieee_value(0.0_dp,ieee_quiet_nan)],nonlinear_ub=[1.0_dp])
    call assert(res%status==shgo_invalid,'nan constraint bounds rejected')

    print '(a,i0)', 'SHGO constrained native assertions passed: ', checks
contains
    subroutine assert(ok,msg)
        logical, intent(in) :: ok !! condition to assert
        character(len=*), intent(in) :: msg !! failure description
        checks=checks+1
        if (.not. ok) then
            print '(a,a)', 'FAIL: ',msg
            error stop 'SHGO constraint test failure'
        end if
    end subroutine assert

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! 2D bounded point
        real(dp) :: f !! shifted quadratic
        f=(x(1)-0.2_dp)**2 + (x(2)-0.1_dp)**2
    end function quadratic

    function shifted(x) result(f)
        real(dp), intent(in) :: x(:) !! 2D bounded point
        real(dp) :: f !! quadratic with exact optimum on x+y=0.7
        f=(x(1)-0.55_dp)**2 + (x(2)-0.15_dp)**2
    end function shifted

    function sphere3(x) result(f)
        real(dp), intent(in) :: x(:) !! 3D point
        real(dp) :: f !! sum of squares
        f=sum((x-0.5_dp)**2)
    end function sphere3

    subroutine nonlinear_sum(x,c)
        real(dp), intent(in) :: x(:) !! 2D point
        real(dp), intent(out) :: c(:) !! sum constraint
        c(1) = sum(x)
    end subroutine nonlinear_sum

    subroutine nonlinear_sum3(x,c)
        real(dp), intent(in) :: x(:) !! 3D point
        real(dp), intent(out) :: c(:) !! sum equality
        c(1)=sum(x)
    end subroutine nonlinear_sum3

    subroutine project_inequality(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! budget guarded objective
        real(dp), intent(in) :: lower(:) !! lower box bounds
        real(dp), intent(in) :: upper(:) !! upper box bounds
        real(dp), intent(in) :: start(:) !! supplied start
        integer, intent(in) :: maxfev !! available objective calls
        real(dp), intent(out) :: x(:) !! polished feasible point
        real(dp), intent(out) :: value !! evaluated objective
        integer, intent(out) :: nfev !! reported calls
        logical, intent(out) :: success !! accepted when budget permits
        x = [0.5_dp,0.4_dp]
        value = huge(1.0_dp)
        nfev=0
        success=maxfev>=1
        if (success) then
            value=fun(x)
            nfev=1
        end if
        if (sum(start)+sum(lower)+sum(upper)<-huge(1.0_dp)) error stop 'unreachable'
    end subroutine project_inequality

    subroutine project_equality(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! budget guarded objective
        real(dp), intent(in) :: lower(:) !! lower bounds
        real(dp), intent(in) :: upper(:) !! upper bounds
        real(dp), intent(in) :: start(:) !! initial point
        integer, intent(in) :: maxfev !! objective budget
        real(dp), intent(out) :: x(:) !! equality point
        real(dp), intent(out) :: value !! objective there
        integer, intent(out) :: nfev !! reported objective count
        logical, intent(out) :: success !! whether a call took place
        x=[0.55_dp,0.15_dp]
        success=maxfev>=1
        nfev=0
        value=huge(1.0_dp)
        if (success) then
            value=fun(x)
            nfev=1
        end if
        if (sum(start)+sum(lower)+sum(upper)<-huge(1.0_dp)) error stop 'unreachable'
    end subroutine project_equality

    subroutine invalid_minimizer(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! ignored by deliberately invalid minimizer
        real(dp), intent(in) :: lower(:) !! bounding box
        real(dp), intent(in) :: upper(:) !! bounding box
        real(dp), intent(in) :: start(:) !! starting point
        integer, intent(in) :: maxfev !! unused budget
        real(dp), intent(out) :: x(:) !! returned infeasible point
        real(dp), intent(out) :: value !! faked objective
        integer, intent(out) :: nfev !! zero reported calls
        logical, intent(out) :: success !! deliberately misleading success
        x=[0.0_dp,0.0_dp]
        value=-1000.0_dp
        nfev=0
        success=.true.
        if (sum(start)+sum(lower)+sum(upper)+real(maxfev,dp)<-huge(1.0_dp)) &
            error stop 'unreachable'
    end subroutine invalid_minimizer
    subroutine fake_feasible_minimizer(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! ignored for adversarial regression
        real(dp), intent(in) :: lower(:) !! box lower bounds
        real(dp), intent(in) :: upper(:) !! box upper bounds
        real(dp), intent(in) :: start(:) !! candidate start
        integer, intent(in) :: maxfev !! evaluation budget
        real(dp), intent(out) :: x(:) !! feasible but not optimized point
        real(dp), intent(out) :: value !! deliberately fabricated value
        integer, intent(out) :: nfev !! falsely reported calls
        logical, intent(out) :: success !! deliberately fabricated success
        x=[1.0_dp,1.0_dp]
        value=-1000.0_dp
        nfev=0
        success=.true.
    end subroutine fake_feasible_minimizer

    subroutine nan_minimizer(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! unused in adversarial regression
        real(dp), intent(in) :: lower(:) !! box lower bounds
        real(dp), intent(in) :: upper(:) !! box upper bounds
        real(dp), intent(in) :: start(:) !! candidate start
        integer, intent(in) :: maxfev !! evaluation budget
        real(dp), intent(out) :: x(:) !! nonfinite solver return
        real(dp), intent(out) :: value !! falsely successful objective
        integer, intent(out) :: nfev !! zero reported objective calls
        logical, intent(out) :: success !! deliberately misleading success
        x=ieee_value(0.0_dp,ieee_quiet_nan)
        value=-1000.0_dp
        nfev=0
        success=.true.
    end subroutine nan_minimizer
end program optimize_shgo_constraints_reference
