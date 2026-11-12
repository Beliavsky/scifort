! SPDX-License-Identifier: MIT
! Feasibility restoration for an equality manifold missed by mesh vertices.
program optimize_shgo_restoration_example
    use scifort_kinds, only : dp
    use scifort_optimize_shgo_minimize, only : shgo
    use scifort_optimize_shgo, only : shgo_result
    implicit none
    type(shgo_result) :: res
    real(dp) :: bounds(2,2)
    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp
    ! The 3-by-3 initial grid contains no point with x(1)+x(2)=0.7.
    ! Restoration is a feasibility projection, not a local optimizer.
    res = shgo(objective,bounds,iters=1,local_iter=0,maxfun=100, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp], &
        nonlinear_jac=sum_jac,restore_feasibility=.true.)
    print '(a,l1)', 'feasible: ',res%feasible
    print '(a,2f15.9)', 'restored x: ',res%x
    print '(a,es15.6)', 'objective: ',res%fun
    print '(a,es15.6)', 'maximum constraint violation: ',res%maxcv
    print '(a,4(i0,1x))', 'nfev ncev ncjev nrestored: ', &
        res%nfev,res%ncev,res%ncjev,res%nrestored
    if (.not. res%success) error stop 'expected feasible equality projection'
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded decision vector
        real(dp) :: f !! quadratic objective
        f = sum((x-0.25_dp)**2)
    end function objective
    subroutine sum_constraint(x,c)
        real(dp), intent(in) :: x(:) !! bounded decision vector
        real(dp), intent(out) :: c(:) !! equality value
        c(1) = sum(x)
    end subroutine sum_constraint
    subroutine sum_jac(x,a)
        real(dp), intent(in) :: x(:) !! bounded decision vector
        real(dp), intent(out) :: a(:,:) !! one constraint derivative per row
        a = 1.0_dp
    end subroutine sum_jac
end program optimize_shgo_restoration_example
