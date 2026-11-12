! Tests the DE-to-L-BFGS-B bridge interface against the test double.
program main
    use scifort_kinds, only: dp
    use scifort_optimize_differential_evolution
    use scifort_optimize_differential_evolution_polish, only:de_lbfgsb_polisher
    implicit none
    type(differential_evolution_result) :: r
    real(dp) :: b(2,2)
    b(:,1)=[-2.0_dp,-2.0_dp]
    b(:,2)=[2.0_dp,2.0_dp]
    r=differential_evolution(f,b,polisher=de_lbfgsb_polisher,maxiter=3,seed=19)
    if (.not.allocated(r%x)) error stop 'no x'
    print *, 'L-BFGS-B BRIDGE INTERFACE TEST PASSED',r%nfev
contains
    real(dp) function f(x)
        real(dp), intent(in) :: x(:)
        f=sum((x-0.75_dp)**2)
    end function
end program
