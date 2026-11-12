! SPDX-License-Identifier: MIT
! Small executable for independent SciPy comparisons of constrained DE.
program reference_de_constraints_driver
    use scifort_kinds, only : dp
    use scifort_optimize_differential_evolution
    implicit none
    character(len=32) :: problem, strategy, seed_arg
    type(differential_evolution_result) :: r
    integer :: seed, status
    real(dp) :: b2(2,2), b3(3,2), a2(1,2), a3(1,3), x(3)
    logical :: mask(3)
    call get_command_argument(1,problem,status=status)
    if (status /= 0) error stop 'missing problem'
    call get_command_argument(2,seed_arg,status=status)
    if (status /= 0) error stop 'missing seed'
    read(seed_arg,*,iostat=status) seed
    if (status /= 0) error stop 'invalid seed'
    call get_command_argument(3,strategy,status=status)
    if (status /= 0) error stop 'missing strategy'

    b2(:,1)=[-1.0_dp,-1.0_dp]
    b2(:,2)=[2.0_dp,2.0_dp]
    b3(:,1)=[-3.6_dp,-3.0_dp,1.2_dp]
    b3(:,2)=[4.4_dp,4.0_dp,5.7_dp]
    mask=[.true.,.false.,.true.]
    a2(1,:)=[1.0_dp,1.0_dp]
    a3(1,:)=[1.0_dp,1.0_dp,0.0_dp]
    select case (trim(problem))
    case ('integer')
        r=differential_evolution(integer_objective,b3,seed=seed,strategy=trim(strategy), &
                                 integrality=mask,maxiter=180,popsize=12, &
                                 tol=0.0_dp,atol=0.0_dp,updating='deferred')
    case ('linear')
        r=differential_evolution(linear_objective,b2,seed=seed,strategy=trim(strategy), &
                                 linear_a=a2,linear_lb=[1.5_dp],linear_ub=[huge(1.0_dp)], &
                                 maxiter=180,popsize=12,tol=0.0_dp,atol=0.0_dp, &
                                 updating='deferred')
    case ('circle')
        r=differential_evolution(circle_objective,b2,seed=seed,strategy=trim(strategy), &
                                 nonlinear=circle,nonlinear_lb=[-huge(1.0_dp)], &
                                 nonlinear_ub=[1.0_dp],maxiter=180,popsize=12, &
                                 tol=0.0_dp,atol=0.0_dp,updating='deferred')
    case ('combined')
        b3(:,1)=[-2.0_dp,-2.0_dp,-2.0_dp]
        b3(:,2)=[2.0_dp,2.0_dp,2.0_dp]
        mask=[.true.,.false.,.false.]
        r=differential_evolution(combined_objective,b3,seed=seed,strategy=trim(strategy), &
                                 integrality=mask,linear_a=a3,linear_lb=[0.5_dp], &
                                 linear_ub=[huge(1.0_dp)],nonlinear=circle_first_two, &
                                 nonlinear_lb=[-huge(1.0_dp)],nonlinear_ub=[0.81_dp], &
                                 maxiter=180,popsize=12,tol=0.0_dp,atol=0.0_dp, &
                                 updating='deferred')
    case ('infeasible')
        b2(:,2)=[1.0_dp,1.0_dp]
        r=differential_evolution(linear_objective,b2,seed=seed,strategy=trim(strategy), &
                                 linear_a=a2,linear_lb=[10.0_dp], &
                                 linear_ub=[huge(1.0_dp)],maxiter=60,popsize=12, &
                                 tol=0.0_dp,atol=0.0_dp,updating='deferred')
    case default
        error stop 'unknown problem'
    end select
    if (.not.allocated(r%x)) error stop 'missing result'
    x = 0.0_dp
    x(1:size(r%x)) = r%x
    write(*,'(a,1x,i0,1x,es27.17e3,1x,es27.17e3,1x,l1,1x,3(es27.17e3,1x),i0,1x,i0,1x,i0)') &
        'RESULT',r%status,r%fun,r%maxcv,r%feasible,x,r%nfev,r%ncev,r%nit
contains
    real(dp) function integer_objective(v) result(z)
        real(dp), intent(in) :: v(:) !! three-dimensional candidate
        z=(v(1)-2.0_dp)**2+(v(2)+0.75_dp)**2+(v(3)-3.0_dp)**2
    end function integer_objective
    real(dp) function linear_objective(v) result(z)
        real(dp), intent(in) :: v(:) !! two-dimensional candidate
        z=(v(1)-0.1_dp)**2+(v(2)-0.2_dp)**2
    end function linear_objective
    real(dp) function circle_objective(v) result(z)
        real(dp), intent(in) :: v(:) !! two-dimensional candidate
        z=(v(1)-1.5_dp)**2+(v(2)-1.3_dp)**2
    end function circle_objective
    real(dp) function combined_objective(v) result(z)
        real(dp), intent(in) :: v(:) !! three-dimensional candidate
        z=(v(1)-1.0_dp)**2+(v(2)-0.7_dp)**2+v(3)**2
    end function combined_objective
    subroutine circle(v,c)
        real(dp), intent(in) :: v(:) !! candidate
        real(dp), intent(out) :: c(:) !! squared radius
        c(1)=sum(v**2)
    end subroutine circle
    subroutine circle_first_two(v,c)
        real(dp), intent(in) :: v(:) !! candidate
        real(dp), intent(out) :: c(:) !! squared radius in x(1:2)
        c(1)=sum(v(1:2)**2)
    end subroutine circle_first_two
end program reference_de_constraints_driver
