! SPDX-License-Identifier: MIT
! Standalone Fortran driver for independent SciPy constrained comparisons.
program optimize_shgo_restoration_driver
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result
    implicit none
    real(dp) :: bounds(2,2), a(1,2)
    type(shgo_result) :: res
    integer :: which
    character(len=10) :: arg
    call get_command_argument(1,arg)
    read(arg,*) which
    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp
    a(1,:) = [1.0_dp,1.0_dp]
    select case (which)
    case (1)
        res = shgo_core(objective,bounds,iters=1,maxfun=200,local_iter=0, &
            linear_a=a,linear_lb=[0.7_dp],linear_ub=[0.7_dp],restore_feasibility=.true.)
    case (2)
        res = shgo_core(objective,bounds,iters=1,maxfun=200,local_iter=0, &
            nonlinear=sum_con,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
            nonlinear_jac=sum_jac,restore_feasibility=.true.)
    case (3)
        res = shgo_core(objective,bounds,iters=1,maxfun=200,local_iter=0, &
            nonlinear=radius_con,nonlinear_lb=[0.13_dp],nonlinear_ub=[0.13_dp], &
            nonlinear_jac=radius_jac,restore_feasibility=.true.)
    case (4)
        res = shgo_core(objective,bounds,iters=1,maxfun=200,local_iter=0, &
            nonlinear=product_con,nonlinear_lb=[0.11_dp],nonlinear_ub=[0.11_dp], &
            nonlinear_jac=product_jac,restore_feasibility=.true.)
    case (5)
        res = shgo_core(objective,bounds,iters=1,maxfun=200,local_iter=0, &
            nonlinear=radius_con,nonlinear_lb=[0.13_dp],nonlinear_ub=[0.13_dp], &
            restore_feasibility=.true.)
    case (6)
        res = shgo_core(objective,bounds,iters=1,maxfun=200,local_iter=0, &
            linear_a=a,linear_lb=[3.0_dp],linear_ub=[3.0_dp],restore_feasibility=.true.)
    case default
        error stop 'invalid restoration comparison'
    end select
    write(*,'(i0,1x,i0,1x,l1,1x,es25.15e3,1x,es15.7e3,5(1x,i0))') &
        which,res%status,res%success,res%fun,res%maxcv, &
        res%nfev,res%ncev,res%ncjev,res%nrestore,res%nrestored
    write(*,'(*(es25.15e3,1x))') res%x
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded decision vector
        real(dp) :: f !! unconstrained quadratic
        f = sum((x-0.25_dp)**2)
    end function objective
    subroutine sum_con(x,c)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: c(:) !! sum of coordinates
        c(1) = sum(x)
    end subroutine sum_con
    subroutine sum_jac(x,j)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: j(:,:) !! exact sum Jacobian
        j = 1.0_dp
    end subroutine sum_jac
    subroutine radius_con(x,c)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: c(:) !! squared radius
        c(1) = sum(x*x)
    end subroutine radius_con
    subroutine radius_jac(x,j)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: j(:,:) !! squared-radius Jacobian
        j(1,:) = 2.0_dp*x
    end subroutine radius_jac
    subroutine product_con(x,c)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: c(:) !! coordinate product
        c(1) = product(x)
    end subroutine product_con
    subroutine product_jac(x,j)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: j(:,:) !! product Jacobian
        j(1,:) = [x(2),x(1)]
    end subroutine product_jac
end program optimize_shgo_restoration_driver
