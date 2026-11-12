! SPDX-License-Identifier: MIT
! Two-dimensional constrained global search with mixed linear and nonlinear bounds.
program optimize_shgo_constraints_example
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_shgo_minimize, only : shgo
    use scifort_optimize_shgo, only : shgo_result
    implicit none
    real(dp) :: bounds(2,2), a(1,2), neg_inf, pos_inf
    type(shgo_result) :: res

    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp
    a(1,:) = [1.0_dp, 1.0_dp]
    neg_inf = ieee_value(0.0_dp, ieee_negative_inf)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    ! The unconstrained minimizer (0.2,0.1) is excluded by x+y >= 0.9.
    ! Also require x(1)**2+x(2)**2 <= 0.8.
    ! The constrained optimum is (0.5,0.4) with f=0.18.
    res = shgo(objective, bounds, iters=2, maxfun=500, &
        linear_a=a, linear_lb=[0.9_dp], linear_ub=[pos_inf], &
        nonlinear=radius2, nonlinear_lb=[neg_inf], nonlinear_ub=[0.8_dp], &
        constraint_tol=1.0e-7_dp)
    print '(a,i0)', 'status: ', res%status
    print '(a,l1)', 'success: ', res%success
    print '(a,es14.6)', 'objective: ', res%fun
    print '(a,2f14.7)', 'x: ', res%x
    print '(a,es14.6)', 'max constraint violation: ', res%maxcv
    print '(a,2(i0,1x))', 'objective/constraint calls: ', res%nfev, res%ncev
contains
    function objective(x) result(value)
        real(dp), intent(in) :: x(:) !! point in [0,1]^2
        real(dp) :: value !! objective value
        value = (x(1)-0.2_dp)**2 + (x(2)-0.1_dp)**2
    end function objective

    subroutine radius2(x,c)
        real(dp), intent(in) :: x(:) !! point in [0,1]^2
        real(dp), intent(out) :: c(:) !! nonlinear squared-radius constraint
        c(1) = sum(x*x)
    end subroutine radius2
end program optimize_shgo_constraints_example
