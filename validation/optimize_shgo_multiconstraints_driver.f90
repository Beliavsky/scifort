! SPDX-License-Identifier: MIT
! Emit numerical results from multiple independent SHGO nonlinear constraints.
program optimize_shgo_multiconstraints_driver
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result, shgo_nonlinear_constraint
    implicit none
    type(shgo_nonlinear_constraint) :: blocks(3)
    type(shgo_result) :: res
    real(dp) :: bounds(2,2), bounds3(3,2), inf
    integer :: problem, status
    character(len=24) :: arg

    read(*,*) problem
    inf = ieee_value(0.0_dp,ieee_positive_inf)
    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp
    bounds3(:,1) = 0.0_dp
    bounds3(:,2) = 1.0_dp
    blocks(1)%fun => constraint_sum
    blocks(1)%jac => jac_sum
    blocks(1)%lb = [0.7_dp]
    blocks(1)%ub = [0.7_dp]
    select case (problem)
    case (1)
        blocks(2)%fun => constraint_difference
        blocks(2)%lb = [0.1_dp]
        blocks(2)%ub = [0.1_dp]
        res = shgo_core(objective,bounds,iters=1,maxfun=150,local_iter=0, &
            constraints=blocks(:2),restore_feasibility=.true.,restoration_starts=9)
    case (2)
        blocks(2)%fun => constraint_radius
        blocks(2)%lb = [0.25_dp]
        blocks(2)%ub = [0.25_dp]
        res = shgo_core(objective,bounds,iters=1,maxfun=150,local_iter=0, &
            constraints=blocks(:2),restore_feasibility=.true.,restoration_starts=9)
    case (3)
        blocks(2)%fun => constraint_first
        blocks(2)%jac => jac_first
        blocks(2)%lb = [0.4_dp]
        blocks(2)%ub = [0.4_dp]
        res = shgo_core(objective,bounds,iters=1,maxfun=150,local_iter=0, &
            constraints=blocks(:2),restore_feasibility=.true.,restoration_starts=9)
    case (4)
        blocks(2)%fun => constraint_first
        blocks(2)%jac => jac_first
        blocks(2)%lb = [0.2_dp]
        blocks(2)%ub = [0.2_dp]
        blocks(3)%fun => constraint_difference23
        blocks(3)%lb = [0.1_dp]
        blocks(3)%ub = [0.1_dp]
        res = shgo_core(objective,bounds3,iters=1,maxfun=150,local_iter=0, &
            constraints=blocks,restore_feasibility=.true.,restoration_starts=9)
    case (5)
        blocks(2)%fun => constraint_difference
        blocks(2)%lb = [0.1_dp]
        blocks(2)%ub = [0.1_dp]
        blocks(3)%fun => constraint_product_radius
        blocks(3)%lb = [0.05_dp,-inf]
        blocks(3)%ub = [inf,0.35_dp]
        res = shgo_core(objective,bounds,iters=1,maxfun=150,local_iter=0, &
            constraints=blocks,restore_feasibility=.true.,restoration_starts=9)
    case default
        error stop 'unknown comparison problem'
    end select
    write(*,'(i0,1x,l1,1x,es24.16,1x,es14.6,1x,i0,1x,i0,1x,*(es24.16,1x))') &
        problem,res%success,res%fun,res%maxcv,res%nfev,res%ncjev,res%x
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp) :: f !! quadratic objective
        f = sum((x-0.25_dp)**2)
    end function objective
    subroutine constraint_sum(x,c)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: c(:) !! sum
        c(1) = sum(x)
    end subroutine constraint_sum
    subroutine jac_sum(x,a)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: a(:,:) !! Jacobian of sum
        a = 1.0_dp
    end subroutine jac_sum
    subroutine constraint_difference(x,c)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: c(:) !! difference
        c(1) = x(1)-x(2)
    end subroutine constraint_difference
    subroutine constraint_difference23(x,c)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: c(:) !! second-third difference
        c(1) = x(2)-x(3)
    end subroutine constraint_difference23
    subroutine constraint_radius(x,c)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: c(:) !! squared Euclidean radius
        c(1) = sum(x*x)
    end subroutine constraint_radius
    subroutine constraint_first(x,c)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: c(:) !! first coordinate
        c(1) = x(1)
    end subroutine constraint_first
    subroutine jac_first(x,a)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: a(:,:) !! derivative of first coordinate
        a = 0.0_dp
        a(1,1) = 1.0_dp
    end subroutine jac_first
    subroutine constraint_product_radius(x,c)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: c(:) !! product and radius
        c(1) = x(1)*x(2)
        c(2) = sum(x*x)
    end subroutine constraint_product_radius
end program optimize_shgo_multiconstraints_driver
