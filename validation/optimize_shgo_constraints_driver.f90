! SPDX-License-Identifier: MIT
! Independent constrained SHGO quality comparison driver.
program optimize_shgo_constraints_driver
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result
    implicit none
    character(len=16) :: argument
    real(dp) :: bounds(2,2), bounds3(3,2), a(1,2), a2(1,3)
    real(dp), parameter :: inf = huge(1.0_dp)
    type(shgo_result) :: res
    integer :: problem

    call get_command_argument(1,argument)
    read(argument,*) problem
    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp
    bounds3(:,1) = 0.0_dp
    bounds3(:,2) = 1.0_dp
    a(1,:)=[1.0_dp,1.0_dp]
    a2(1,:)=[1.0_dp,1.0_dp,1.0_dp]
    select case(problem)
    case(1)
        res=shgo_core(quadratic,bounds,iters=2,maxfun=200,local_iter=1, &
            minimizer=linear_projection,linear_a=a,linear_lb=[0.9_dp],linear_ub=[inf])
    case(2)
        res=shgo_core(shifted,bounds,iters=2,maxfun=200,local_iter=1, &
            minimizer=sum_projection,linear_a=a,linear_lb=[0.7_dp],linear_ub=[0.7_dp])
    case(3)
        res=shgo_core(ball_objective,bounds,iters=2,maxfun=200,local_iter=1, &
            minimizer=ball_projection, nonlinear=radius_squared, &
            nonlinear_lb=[-inf],nonlinear_ub=[0.4_dp])
    case(4)
        res=shgo_core(ring_objective,bounds,iters=2,maxfun=200,local_iter=1, &
            minimizer=ring_projection, nonlinear=radius_squared, &
            nonlinear_lb=[0.5_dp],nonlinear_ub=[0.5_dp])
    case(5)
        res=shgo_core(shifted,bounds,iters=2,maxfun=200,local_iter=1, &
            minimizer=sum_projection,linear_a=a,linear_lb=[0.7_dp],linear_ub=[0.7_dp], &
            nonlinear=radius_squared,nonlinear_lb=[-inf],nonlinear_ub=[0.5_dp])
    case(6)
        res=shgo_core(sphere3,bounds3,iters=1,maxfun=300,local_iter=1, &
            adaptive=.true.,conforming=.true., minimizer=project_3d, &
            linear_a=a2,linear_lb=[1.5_dp],linear_ub=[1.5_dp])
    case default
        error stop 'unknown constrained SHGO comparison case'
    end select
    write(*,'(i0,1x,i0,1x,es23.15,1x,es15.7,1x,i0,1x,i0,1x,i0)') &
        problem,res%status,res%fun,res%maxcv,res%nfev,res%ncev,res%nlocal
    write(*,'(*(es23.15,1x))') res%x
contains
    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp) :: f !! quadratic objective
        f=(x(1)-0.2_dp)**2+(x(2)-0.1_dp)**2
    end function quadratic
    function shifted(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp) :: f !! quadratic objective
        f=(x(1)-0.55_dp)**2+(x(2)-0.15_dp)**2
    end function shifted
    function ball_objective(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp) :: f !! quadratic objective
        f=sum((x-0.8_dp)**2)
    end function ball_objective
    function ring_objective(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp) :: f !! quadratic objective
        f=(x(1)-0.3_dp)**2+(x(2)-0.6_dp)**2
    end function ring_objective
    function sphere3(x) result(f)
        real(dp), intent(in) :: x(:) !! three coordinates
        real(dp) :: f !! centered objective on sum(x)=1.5
        f=sum((x-[0.4_dp,0.5_dp,0.6_dp])**2)
    end function sphere3
    subroutine radius_squared(x,c)
        real(dp), intent(in) :: x(:) !! point on unit square
        real(dp), intent(out) :: c(:) !! single radial constraint
        c(1)=sum(x*x)
    end subroutine radius_squared
    subroutine linear_projection(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! budget-guarded callback
        real(dp), intent(in) :: lower(:) !! box lower limits
        real(dp), intent(in) :: upper(:) !! box upper limits
        real(dp), intent(in) :: start(:) !! local start
        integer, intent(in) :: maxfev !! allowed objective evaluations
        real(dp), intent(out) :: x(:) !! projected solution
        real(dp), intent(out) :: value !! objective at x
        integer, intent(out) :: nfev !! reported calls
        logical, intent(out) :: success !! accepted local result
        x=[0.5_dp,0.4_dp]
        call finish(fun,maxfev,x,value,nfev,success)
    end subroutine linear_projection
    subroutine sum_projection(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! budget-guarded callback
        real(dp), intent(in) :: lower(:) !! box lower limits
        real(dp), intent(in) :: upper(:) !! box upper limits
        real(dp), intent(in) :: start(:) !! local start
        integer, intent(in) :: maxfev !! allowed objective evaluations
        real(dp), intent(out) :: x(:) !! equality solution
        real(dp), intent(out) :: value !! objective at x
        integer, intent(out) :: nfev !! reported calls
        logical, intent(out) :: success !! accepted local result
        x=[0.55_dp,0.15_dp]
        call finish(fun,maxfev,x,value,nfev,success)
    end subroutine sum_projection
    subroutine ball_projection(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! budget-guarded callback
        real(dp), intent(in) :: lower(:) !! box lower limits
        real(dp), intent(in) :: upper(:) !! box upper limits
        real(dp), intent(in) :: start(:) !! local start
        integer, intent(in) :: maxfev !! allowed objective evaluations
        real(dp), intent(out) :: x(:) !! circle boundary solution
        real(dp), intent(out) :: value !! objective at x
        integer, intent(out) :: nfev !! reported calls
        logical, intent(out) :: success !! accepted local result
        x=sqrt(0.2_dp)
        call finish(fun,maxfev,x,value,nfev,success)
    end subroutine ball_projection
    subroutine ring_projection(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! budget-guarded callback
        real(dp), intent(in) :: lower(:) !! box lower limits
        real(dp), intent(in) :: upper(:) !! box upper limits
        real(dp), intent(in) :: start(:) !! local start
        integer, intent(in) :: maxfev !! allowed objective evaluations
        real(dp), intent(out) :: x(:) !! equality circle solution
        real(dp), intent(out) :: value !! objective at x
        integer, intent(out) :: nfev !! reported calls
        logical, intent(out) :: success !! accepted local result
        x=[0.3_dp,0.6_dp]*sqrt(0.5_dp/0.45_dp)
        call finish(fun,maxfev,x,value,nfev,success)
    end subroutine ring_projection
    subroutine project_3d(fun,lower,upper,start,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! budget-guarded callback
        real(dp), intent(in) :: lower(:) !! box lower limits
        real(dp), intent(in) :: upper(:) !! box upper limits
        real(dp), intent(in) :: start(:) !! local start
        integer, intent(in) :: maxfev !! allowed objective evaluations
        real(dp), intent(out) :: x(:) !! 3d equality solution
        real(dp), intent(out) :: value !! objective at x
        integer, intent(out) :: nfev !! reported calls
        logical, intent(out) :: success !! accepted local result
        x=[0.4_dp,0.5_dp,0.6_dp]
        call finish(fun,maxfev,x,value,nfev,success)
    end subroutine project_3d
    subroutine finish(fun,maxfev,x,value,nfev,success)
        use scifort_optimize_shgo, only : shgo_objective
        procedure(shgo_objective) :: fun !! budget-guarded objective
        integer, intent(in) :: maxfev !! remaining budget
        real(dp), intent(in) :: x(:) !! projected local point
        real(dp), intent(out) :: value !! objective at x
        integer, intent(out) :: nfev !! evaluations used
        logical, intent(out) :: success !! true when budget remains
        success=maxfev>=1
        nfev=0
        value=huge(1.0_dp)
        if (success) then
            value=fun(x)
            nfev=1
        end if
    end subroutine finish
end program optimize_shgo_constraints_driver
