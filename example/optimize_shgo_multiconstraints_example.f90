! SPDX-License-Identifier: MIT
! Multiple independent nonlinear constraint objects with mixed Jacobians.
program optimize_shgo_multiconstraints_example
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_nonlinear_constraint, shgo_result
    use scifort_optimize_shgo_minimize, only : shgo
    implicit none
    real(dp) :: bounds(2,2), infinity
    type(shgo_nonlinear_constraint) :: blocks(3)
    type(shgo_result) :: res

    infinity = ieee_value(0.0_dp,ieee_positive_inf)
    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp

    ! The first block enforces x(1)+x(2) = 0.7 and supplies an exact Jacobian.
    blocks(1)%fun => sum_constraint
    blocks(1)%jac => sum_jacobian
    blocks(1)%lb = [0.7_dp]
    blocks(1)%ub = [0.7_dp]

    ! The second block has no Jacobian: SHGO estimates only its own rows.
    blocks(2)%fun => difference_constraint
    blocks(2)%lb = [0.1_dp]
    blocks(2)%ub = [0.1_dp]

    ! A third block returns multiple components with independent bounds.
    blocks(3)%fun => product_radius
    blocks(3)%lb = [0.05_dp,-infinity]
    blocks(3)%ub = [infinity,0.35_dp]

    res = shgo(objective,bounds,iters=1,maxfun=150, &
        constraints=blocks,restore_feasibility=.true.,restoration_starts=9)
    print '(a,l1)', 'feasible: ',res%feasible
    print '(a,*(f12.8,1x))', 'best x: ',res%x
    print '(a,es15.7)', 'objective: ',res%fun
    print '(a,es15.7)', 'maximum constraint violation: ',res%maxcv
    print '(a,4(i0,1x))', 'nfev ncev ncjev nrestored: ', &
        res%nfev,res%ncev,res%ncjev,res%nrestored
    if (.not. res%success .or. res%maxcv > 1.0e-8_dp) &
        error stop 'SHGO example failed feasibility check'

contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! two bounded variables
        real(dp) :: f !! objective value
        f = sum((x-0.25_dp)**2)
    end function objective

    subroutine sum_constraint(x,c)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp), intent(out) :: c(:) !! sum equality
        c(1) = sum(x)
    end subroutine sum_constraint

    subroutine sum_jacobian(x,a)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp), intent(out) :: a(:,:) !! one analytic constraint row
        a = 1.0_dp
    end subroutine sum_jacobian

    subroutine difference_constraint(x,c)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp), intent(out) :: c(:) !! difference equality
        c(1) = x(1)-x(2)
    end subroutine difference_constraint

    subroutine product_radius(x,c)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp), intent(out) :: c(:) !! product and squared radius
        c(1) = x(1)*x(2)
        c(2) = sum(x*x)
    end subroutine product_radius
end program optimize_shgo_multiconstraints_example
