! SPDX-License-Identifier: MIT
! Independent feasibility-restoration regression tests for SHGO.
program optimize_shgo_restoration_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result, shgo_completed, &
        shgo_invalid, shgo_infeasible
    implicit none
    type(shgo_result) :: res
    real(dp) :: bounds(2,2), bounds3(3,2), a(1,2), a2(2,2)
    integer :: checks, jac_calls

    checks = 0
    jac_calls = 0
    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp
    a(1,:) = [1.0_dp,1.0_dp]

    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp])
    call assert(res%status == shgo_infeasible, 'default leaves narrow equality unsampled')
    call assert(res%nrestore == 0 .and. res%nrestored == 0, 'default disables restoration')
    call assert(res%ncjev == 0, 'default does not evaluate any nonlinear Jacobian')

    jac_calls = 0
    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        nonlinear_jac=sum_jac,restore_feasibility=.true.)
    call assert(res%success .and. res%feasible, 'analytic Jacobian restores equality')
    call assert(res%status == shgo_completed, 'successful restoration completes normally')
    call assert(abs(sum(res%x)-0.7_dp) <= 1.0e-8_dp, 'nonlinear equality residual')
    call assert(res%maxcv <= 1.0e-8_dp, 'returned maximum violation')
    call assert(res%nrestore > 0 .and. res%nrestored > 0, 'restoration diagnostics')
    call assert(res%ncjev > 0 .and. res%ncjev == jac_calls, 'analytic Jacobian counter')
    call assert(res%nminima >= 1, 'restored feasible point is in minima pool')
    call assert(abs(res%fun-objective(res%x)) < 1.0e-13_dp, 'true objective evaluated')
    call assert(res%nfev <= 100 .and. res%nfev > 0, 'objective budget not exceeded')

    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        restore_feasibility=.true.)
    call assert(res%success .and. res%feasible, 'finite difference restoration')
    call assert(res%maxcv <= 1.0e-8_dp, 'finite difference equality residual')
    call assert(res%ncjev == 0 .and. res%ncev > 9, 'FD uses constraint-vector calls')
    call assert(res%nrestored > 0, 'FD restoration counter')

    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        linear_a=a,linear_lb=[0.7_dp],linear_ub=[0.7_dp], &
        restore_feasibility=.true.)
    call assert(res%success .and. res%feasible, 'exact linear Jacobian restoration')
    call assert(abs(sum(res%x)-0.7_dp) < 1.0e-8_dp, 'linear equality restored')
    call assert(res%ncjev == 0, 'linear Jacobian has no user callback')
    call assert(res%nrestored > 0, 'linear equality restored count')

    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=radius_constraint,nonlinear_lb=[0.13_dp],nonlinear_ub=[0.13_dp], &
        nonlinear_jac=radius_jac,restore_feasibility=.true.)
    call assert(res%success .and. res%feasible, 'curved equality restored')
    call assert(abs(sum(res%x**2)-0.13_dp) <= 1.0e-8_dp, 'curved residual')
    call assert(res%ncjev > 0, 'curved analytic Jacobian used')

    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=radius_constraint,nonlinear_lb=[0.13_dp],nonlinear_ub=[0.13_dp], &
        restore_feasibility=.true.)
    call assert(res%success, 'curved FD equality restored')
    call assert(abs(sum(res%x**2)-0.13_dp) <= 1.0e-8_dp, 'curved FD equality residual')

    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=product_constraint,nonlinear_lb=[0.11_dp],nonlinear_ub=[0.11_dp], &
        nonlinear_jac=product_jac,restore_feasibility=.true.)
    call assert(res%success .and. res%feasible, 'nonlinear product equality')
    call assert(abs(product(res%x)-0.11_dp) <= 1.0e-8_dp, 'product equality residual')

    ! An unattainable equality must never be reported as feasible.
    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        linear_a=a,linear_lb=[3.0_dp],linear_ub=[3.0_dp], &
        restore_feasibility=.true.)
    call assert(res%status == shgo_infeasible, 'impossible equality status')
    call assert(.not. res%success .and. .not. res%feasible, 'impossible equality not accepted')
    call assert(res%nrestored == 0, 'no fabricated solution')
    call assert(res%nfev == 0, 'do not evaluate objective for impossible equality')
    call assert(res%maxcv >= 1.0_dp, 'least violation for impossible equality')

    a2(1,:) = [1.0_dp,1.0_dp]
    a2(2,:) = [1.0_dp,1.0_dp]
    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        linear_a=a2,linear_lb=[0.7_dp,0.8_dp],linear_ub=[0.7_dp,0.8_dp], &
        restore_feasibility=.true.)
    call assert(.not. res%success, 'contradictory equalities cannot be restored')
    call assert(res%nrestored == 0, 'inconsistent equalities never accepted')

    ! The restoration phase must respect strict objective call accounting.
    res = shgo_core(objective,bounds,iters=1,maxfun=9,local_iter=0, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        nonlinear_jac=sum_jac,restore_feasibility=.true.)
    call assert(res%nfev <= 9, 'strict maxfun on restored solution')
    call assert(res%success, 'restoration works with small objective budget')

    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        nonlinear_jac=sum_jac,restore_feasibility=.true.,restoration_starts=0)
    call assert(.not. res%success .and. res%nrestore == 0, 'zero starts disables projection')

    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        restoration_starts=-1,restore_feasibility=.true.)
    call assert(res%status == shgo_invalid, 'negative projection starts rejected')
    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        restoration_maxiter=0,restore_feasibility=.true.)
    call assert(res%status == shgo_invalid, 'zero projection iterations rejected')
    res = shgo_core(objective,bounds,nonlinear_jac=sum_jac)
    call assert(res%status == shgo_invalid, 'Jacobian without nonlinear function rejected')
    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        nonlinear_jac=invalid_jac,restore_feasibility=.true.)
    call assert(.not. res%success, 'nonfinite analytic Jacobian never accepted')
    call assert(res%nrestored == 0, 'invalid analytic Jacobian has no accepted points')

    ! Exercise the conforming 3D graph with an equality absent from mesh.
    bounds3(:,1) = 0.0_dp
    bounds3(:,2) = 1.0_dp
    res = shgo_core(objective,bounds3,iters=1,maxfun=200,local_iter=0, &
        adaptive=.true.,conforming=.true.,nonlinear=sum_constraint, &
        nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        nonlinear_jac=sum_jac,restore_feasibility=.true.)
    call assert(res%conforming_valid .and. res%ntetrahedra == 192, '3D mesh unaffected')
    call assert(res%success .and. res%maxcv <= 1.0e-8_dp, '3D equality restored')
    call assert(res%nrestored > 0 .and. res%ncjev > 0, '3D diagnostic counters')

    print '(a,i0)', 'SHGO restoration native assertions passed: ', checks
contains
    subroutine assert(ok,msg)
        logical, intent(in) :: ok !! condition to verify
        character(len=*), intent(in) :: msg !! assertion name
        checks = checks+1
        if (.not. ok) then
            print '(a,a)', 'FAIL: ',msg
            error stop 'SHGO restoration failure'
        end if
    end subroutine assert

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded optimization coordinates
        real(dp) :: f !! original scalar objective
        f = sum((x-0.25_dp)**2)
    end function objective

    subroutine sum_constraint(x,c)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: c(:) !! sum constraint output
        c(1) = sum(x)
    end subroutine sum_constraint
    subroutine sum_jac(x,a)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: a(:,:) !! derivative of sum(x)
        a = 1.0_dp
        jac_calls = jac_calls+1
    end subroutine sum_jac
    subroutine radius_constraint(x,c)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: c(:) !! squared-radius constraint output
        c(1) = sum(x*x)
    end subroutine radius_constraint
    subroutine radius_jac(x,a)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: a(:,:) !! derivative of squared radius
        a(1,:) = 2.0_dp*x
    end subroutine radius_jac
    subroutine product_constraint(x,c)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: c(:) !! product constraint output
        c(1) = product(x)
    end subroutine product_constraint
    subroutine product_jac(x,a)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: a(:,:) !! product constraint derivative
        a(1,:) = [x(2),x(1)]
    end subroutine product_jac
    subroutine invalid_jac(x,a)
        real(dp), intent(in) :: x(:) !! bounded coordinates
        real(dp), intent(out) :: a(:,:) !! deliberately invalid derivative
        a = ieee_value(0.0_dp,ieee_quiet_nan)
    end subroutine invalid_jac
end program optimize_shgo_restoration_reference
