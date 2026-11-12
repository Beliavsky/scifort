! SPDX-License-Identifier: MIT
! Public adapter API smoke test using the validation-only local solver stub.
program optimize_shgo_adapter_smoke
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_shgo_minimize, only : shgo
    use scifort_optimize_shgo, only : shgo_result, shgo_completed, shgo_invalid, &
        shgo_hgr_stable, shgo_nonlinear_constraint
    implicit none
    real(dp) :: bounds(2,2), bounds3(3,2), a(1,2)
    type(shgo_result) :: res
    type(shgo_nonlinear_constraint) :: blocks(3), bad(1)
    real(dp) :: infinity
    integer :: jac_calls, vec_jac_calls
    jac_calls = 0
    vec_jac_calls = 0
    infinity = ieee_value(0.0_dp,ieee_positive_inf)
    bounds(:,1)=-2.0_dp
    bounds(:,2)=2.0_dp
    res=shgo(fun,bounds,iters=3,adaptive=.true.,adaptive_fraction=0.25_dp, &
        method='Nelder-Mead',local_iter=1,maxfun=200)
    if (res%status /= shgo_completed .or. res%nrefined /= 3) &
        error stop 'SHGO adaptive public adapter failed'
    if (.not. res%success .or. res%nmesh == 0 .or. res%nleaf /= 13) &
        error stop 'SHGO adapter invalid diagnostics'
    res=shgo(fun,bounds,iters=2,adaptive=.true.,method='powell', &
        local_iter=1,maxfun=200)
    if (res%status /= shgo_completed .or. .not. res%success) &
        error stop 'SHGO Powell adapter failed'
    res=shgo(fun,bounds,adaptive=.true.,adaptive_fraction=-1.0_dp)
    if (res%status /= shgo_invalid) error stop 'SHGO adapter option validation failed'
    res=shgo(flat,bounds,iters=5,hgr_stability_iters=2,local_iter=0)
    if (res%status /= shgo_hgr_stable .or. res%nit /= 3 .or. &
        res%hgr_estimate /= 1) error stop 'SHGO HGR adapter option forwarding failed'
    res=shgo(fun,bounds,iters=2,adaptive=.true.,conforming=.true., &
        local_iter=0,maxfun=100)
    if (res%status /= shgo_completed .or. .not. res%conforming_valid .or. &
        res%ntriangles /= 30 .or. res%nmesh /= 21) &
        error stop 'SHGO conforming adapter forwarding failed'
    bounds3(:,1) = -1.0_dp
    bounds3(:,2) = 1.0_dp
    res=shgo(fun3,bounds3,iters=2,adaptive=.true.,conforming=.true., &
        local_iter=0,maxfun=200)
    if (res%status /= shgo_completed .or. .not. res%conforming_valid .or. &
        res%ntetrahedra /= 596 .or. res%nmesh /= 175) &
        error stop 'SHGO 3D conforming adapter forwarding failed'
    a(1,:) = [1.0_dp,1.0_dp]
    res=shgo(fun,bounds,iters=1,maxfun=200, &
        linear_a=a,linear_lb=[0.5_dp],linear_ub=[huge(1.0_dp)])
    if (res%status == shgo_invalid .or. res%ncev <= res%nmesh .or. &
        res%nlocal <= 0) error stop 'SHGO constrained SLSQP adapter failed'
    res=shgo(fun,bounds,iters=1,maxfun=200, &
        nonlinear=nonlinear_sum,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp])
    if (res%status == shgo_invalid .or. res%ncev <= res%nmesh .or. &
        res%nlocal <= 0) error stop 'SHGO nonlinear equality adapter forwarding failed'
    jac_calls = 0
    res=shgo(fun,bounds,iters=1,maxfun=200, &
        linear_a=a,linear_lb=[0.7_dp],linear_ub=[0.7_dp], &
        nonlinear=nonlinear_sum,nonlinear_lb=[-0.5_dp],nonlinear_ub=[0.8_dp], &
        nonlinear_jac=sum_jac)
    if (res%status == shgo_invalid .or. res%ncjev == 0 .or. &
        res%ncjev /= jac_calls) error stop 'SHGO combined analytic SLSQP Jacobian failed'
    ! The stub independently checks SLSQP's equality-first/inequality-next
    ! Jacobian rows using finite differences of the canonical constraints.
    ! The first and third blocks are analytic; the middle one uses finite
    ! differences. All three must enter the local solver in correct order.
    blocks(1)%fun => nonlinear_sum
    blocks(1)%jac => sum_jac
    blocks(1)%lb = [0.0_dp]
    blocks(1)%ub = [0.0_dp]
    blocks(2)%fun => nonlinear_difference
    blocks(2)%lb = [-1.0_dp]
    blocks(2)%ub = [1.0_dp]
    blocks(3)%fun => nonlinear_product_and_radius
    blocks(3)%jac => product_radius_jac
    blocks(3)%lb = [-infinity,-infinity]
    blocks(3)%ub = [1.0_dp,3.0_dp]
    jac_calls = 0
    vec_jac_calls = 0
    res=shgo(fun,bounds,iters=1,maxfun=250,constraints=blocks)
    if (res%status == shgo_invalid .or. res%nlocal == 0 .or. &
        res%ncjev <= 0 .or. res%ncjev /= jac_calls+vec_jac_calls .or. &
        vec_jac_calls <= 0 .or. res%ncev <= res%nmesh) &
        error stop 'SHGO multiple SLSQP constraints and mixed Jacobians failed'
    res=shgo(fun,bounds,iters=1,maxfun=250,constraints=blocks, &
        nonlinear=nonlinear_sum,nonlinear_lb=[0.0_dp],nonlinear_ub=[0.0_dp])
    if (res%status /= shgo_invalid) error stop 'SHGO adapter accepted mixed constraint APIs'
    res=shgo(fun,bounds,iters=1,maxfun=250,constraints=bad)
    if (res%status /= shgo_invalid) error stop 'SHGO adapter accepted invalid constraint record'
    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp
    res=shgo(fun,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=nonlinear_sum,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        nonlinear_jac=sum_jac,restore_feasibility=.true.)
    if (.not. res%success .or. res%nrestored < 1 .or. res%maxcv > 1.0e-8_dp) &
        error stop 'SHGO adapter feasibility restoration forwarding failed'
    res=shgo(fun,bounds,iters=1,method='powell', &
        linear_a=a,linear_lb=[0.7_dp],linear_ub=[0.7_dp])
    if (res%status /= shgo_invalid) error stop 'nonconstrained local method accepted'
    res=shgo(fun,bounds,iters=1,method='slsqp')
    if (res%status /= shgo_completed) error stop 'unconstrained SLSQP adapter failed'
    res=shgo(fun,bounds,conforming=.true.)
    if (res%status /= shgo_invalid) error stop 'SHGO conforming adapter invalid option failed'
    res=shgo(flat,bounds,hgr_stability_iters=-1)
    if (res%status /= shgo_invalid) error stop 'SHGO HGR adapter invalid option failed'
    print '(a)', 'SHGO adaptive/HGR/conforming/constraint-Jacobian adapter test passed'
contains
    function fun(x) result(f)
        real(dp), intent(in) :: x(:) !! two coordinates
        real(dp) :: f !! shifted convex objective
        f = sum((x-[0.5_dp,-0.5_dp])**2)
    end function fun
    function fun3(x) result(f)
        real(dp), intent(in) :: x(:) !! three-dimensional bounded input
        real(dp) :: f !! smooth shifted quadratic
        f = sum((x-[0.3_dp,-0.2_dp,0.15_dp])**2)
    end function fun3
    function flat(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp) :: f !! constant objective
        f = 0.0_dp*sum(x)
    end function flat
    subroutine nonlinear_sum(x,c)
        real(dp), intent(in) :: x(:) !! box coordinates
        real(dp), intent(out) :: c(:) !! nonlinear equality value
        c(1) = sum(x)
    end subroutine nonlinear_sum
    subroutine sum_jac(x,a)
        real(dp), intent(in) :: x(:) !! SLSQP local coordinates
        real(dp), intent(out) :: a(:,:) !! Jacobian of sum(x)
        a = 1.0_dp
        jac_calls = jac_calls+1
    end subroutine sum_jac
    subroutine nonlinear_difference(x,c)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: c(:) !! first-minus-second
        c(1) = x(1)-x(2)
    end subroutine nonlinear_difference
    subroutine nonlinear_product_and_radius(x,c)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: c(:) !! product and radius
        c(1) = x(1)*x(2)
        c(2) = sum(x*x)
    end subroutine nonlinear_product_and_radius
    subroutine product_radius_jac(x,a)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: a(:,:) !! gradients of product and radius
        a(1,:) = [x(2),x(1)]
        a(2,:) = 2.0_dp*x
        vec_jac_calls = vec_jac_calls + 1
    end subroutine product_radius_jac
end program optimize_shgo_adapter_smoke
