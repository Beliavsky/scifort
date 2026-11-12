! SPDX-License-Identifier: MIT
! Native tests for integer and constrained differential evolution.
program test_optimize_de_constraints_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_differential_evolution
    implicit none

    type(differential_evolution_result) :: r, s
    real(dp) :: b(3,2), linear(2,3), cmin(2), cmax(2), fixed(3,2)
    real(dp) :: b2(2,2), a2(1,2), reference
    integer :: checks, polish_calls
    logical :: mask(3)
    checks = 0
    polish_calls = 0
    b(:,1) = [-3.6_dp,-3.0_dp,1.2_dp]
    b(:,2) = [ 4.4_dp, 4.0_dp,5.7_dp]
    mask = [.true.,.false.,.true.]

    r = differential_evolution(mixed_objective,b,integrality=mask,maxiter=220,popsize=12,seed=91)
    call must(r%fun < 1.0e-5_dp, 'mixed-integer objective')
    call must(r%feasible.and.r%maxcv<=0.0_dp, 'mixed-integer feasibility')
    call must(all(abs(r%population(1,:)-anint(r%population(1,:)))<=0.0_dp), 'integer population row 1')
    call must(all(abs(r%population(3,:)-anint(r%population(3,:)))<=0.0_dp), 'integer population row 3')
    call must(all(r%population(1,:) >= -3.0_dp), 'fractional lower integer bound')
    call must(all(r%population(3,:) <= 5.0_dp), 'fractional upper integer bound')
    call must(r%ncev==0, 'no constraint callbacks without constraints')
    checks = checks + 7

    s = differential_evolution(mixed_objective,b,integrality=mask,maxiter=220,popsize=12,seed=91)
    call must(all(abs(s%population-r%population)<=0.0_dp), 'mixed integer RNG reproducible')
    call must(s%nfev==r%nfev, 'mixed integer nfev reproducible')
    checks = checks + 2
    r = differential_evolution(mixed_objective,b,maxiter=0,seed=29,integrality=mask, &
                               x0=[1.8_dp,-0.75_dp,3.1_dp])
    call must(any(all(abs(r%population - spread([2.0_dp,-0.75_dp,3.0_dp],2,size(r%population,2)))<=0.0_dp,dim=1)), &
              'x0 integer projection')
    checks = checks + 1
    r = differential_evolution(mixed_objective,b,maxiter=60,seed=33,integrality=mask, &
                               polisher=mock_polisher)
    call must(polish_calls==1, 'polisher called for continuous coordinates')
    call must(abs(r%x(1)-2.0_dp)<1.0e-12_dp.and.abs(r%x(3)-3.0_dp)<1.0e-12_dp, &
              'polisher freezes integer coordinates')
    call must(abs(r%x(2)+0.75_dp)<1.0e-12_dp, 'polisher improves continuous coordinate')
    call must(r%nfev>=1, 'polishing objective accounting')
    checks = checks + 4

    b2(:,1) = [-1.0_dp,-1.0_dp]
    b2(:,2) = [ 2.0_dp, 2.0_dp]
    a2(1,:) = [1.0_dp,1.0_dp]
    r = differential_evolution(linear_objective,b2,linear_a=a2,linear_lb=[1.5_dp], &
                               linear_ub=[huge(1.0_dp)],maxiter=220,popsize=14,seed=137)
    call must(r%feasible.and.r%maxcv<1.0e-9_dp, 'linear inequality feasibility')
    call must(sum(r%x)>=1.5_dp-1.0e-9_dp, 'linear inequality satisfied')
    call must(r%fun<0.78_dp, 'linear inequality reaches boundary optimum')
    call must(r%ncev>=r%nfev, 'constraint evaluations counted')
    call must(any(r%population_constraint_violation>0.0_dp).eqv..false., &
              'all population members converge to feasible points')
    checks = checks + 5

    r = differential_evolution(outside_objective,b2,nonlinear=circle_constraint, &
                               nonlinear_lb=[-huge(1.0_dp)],nonlinear_ub=[1.0_dp], &
                               maxiter=240,popsize=13,seed=51)
    reference = (sqrt(1.5_dp**2+1.3_dp**2)-1.0_dp)**2
    call must(r%feasible, 'nonlinear inequality feasible')
    call must(sum(r%x**2)<1.0_dp+1.0e-9_dp, 'nonlinear inequality respected')
    call must(abs(r%fun-reference)<0.025_dp, 'nonlinear constrained objective')
    checks = checks + 3

    b(:,1)=[-2.0_dp,-2.0_dp,-2.0_dp]
    b(:,2)=[ 2.0_dp, 2.0_dp, 2.0_dp]
    mask=[.true.,.false.,.false.]
    linear(1,:)=[1.0_dp,1.0_dp,0.0_dp]
    cmin(1)=0.5_dp
    cmax(1)=huge(1.0_dp)
    r = differential_evolution(combined_objective,b,integrality=mask, &
                               linear_a=linear(1:1,:),linear_lb=cmin(1:1),linear_ub=cmax(1:1), &
                               nonlinear=circle_first_two, nonlinear_lb=[-huge(1.0_dp)], &
                               nonlinear_ub=[0.81_dp],maxiter=220,popsize=16,seed=29)
    call must(r%feasible.and.r%maxcv<1.0e-9_dp, 'mixed constraints feasible')
    call must(abs(r%x(1))<=0.0_dp, 'nonlinear constraint forces integer zero')
    call must(r%fun<1.01_dp, 'mixed constrained objective')
    call must(all(abs(r%population(1,:)-anint(r%population(1,:)))<=0.0_dp), 'mixed integer population')
    checks = checks + 4

    b2(:,1)=[-1.0_dp,-1.0_dp]
    b2(:,2)=[1.0_dp,1.0_dp]
    r = differential_evolution(linear_objective,b2,linear_a=a2,linear_lb=[10.0_dp], &
                               linear_ub=[huge(1.0_dp)], maxiter=80,popsize=8,seed=27)
    call must(.not.r%feasible.and.r%status==de_status_infeasible, 'infeasible status')
    call must(r%nfev==0, 'never call objective for infeasible candidates')
    call must(r%maxcv>=8.0_dp-1.0e-9_dp, 'constraint violation reported')
    call must(all(r%population_energies>1.0e20_dp), 'invalid energies reported')
    checks = checks + 4

    fixed(:,1)=[1.0_dp,1.0_dp,1.0_dp]
    fixed(:,2)=fixed(:,1)
    linear(1,:)=[1.0_dp,1.0_dp,1.0_dp]
    r=differential_evolution(mixed_objective,fixed,linear_a=linear(1:1,:), &
                             linear_lb=[8.0_dp],linear_ub=[huge(1.0_dp)])
    call must(r%status==de_status_infeasible.and..not.r%feasible, 'infeasible fixed constraints')
    call must(r%nfev==0.and.r%ncev==1, 'fixed infeasible objective not called')
    checks = checks + 2

    r=differential_evolution(mixed_objective,fixed,nonlinear=circle_first_two, &
                             nonlinear_lb=[0.0_dp],nonlinear_ub=[1.0_dp])
    call must(r%status==de_status_infeasible.and..not.r%feasible, &
              'fixed nonlinear constraint callback')
    call must(r%nfev==0.and.r%ncev==1, 'fixed nonlinear evaluation accounting')
    checks = checks + 2

    b2(:,1)=[0.0_dp,0.0_dp]
    b2(:,2)=[1.0_dp,1.0_dp]
    r=differential_evolution(linear_objective,b2,nonlinear=nan_constraint, &
                             nonlinear_lb=[0.0_dp],nonlinear_ub=[1.0_dp], &
                             maxiter=120,seed=77)
    call must(r%feasible, 'NaN constraint treated as infeasible')
    call must(r%x(1)>=0.5_dp, 'NaN region avoided')
    checks = checks + 2

    b(:,1)=[0.2_dp,0.0_dp,-3.0_dp]
    b(:,2)=[0.8_dp,1.0_dp,2.0_dp]
    r=differential_evolution(mixed_objective,b,integrality=mask)
    call must(r%status==de_status_invalid, 'empty integer interval rejected')
    r=differential_evolution(mixed_objective,b,integrality=[.true.])
    call must(r%status==de_status_invalid, 'wrong integrality shape')
    r=differential_evolution(mixed_objective,b,linear_a=linear)
    call must(r%status==de_status_invalid, 'incomplete linear constraint args')
    r=differential_evolution(mixed_objective,b,nonlinear=circle_first_two)
    call must(r%status==de_status_invalid, 'incomplete nonlinear constraint args')
    r=differential_evolution(mixed_objective,b,constraint_tol=-0.1_dp)
    call must(r%status==de_status_invalid, 'negative constraint_tol rejected')
    r=differential_evolution(mixed_objective,b,linear_a=linear, &
                             linear_lb=[0.0_dp],linear_ub=[1.0_dp])
    call must(r%status==de_status_invalid, 'wrong linear constraint shape')
    r=differential_evolution(mixed_objective,b,nonlinear=circle_first_two, &
                             nonlinear_lb=[0.0_dp],nonlinear_ub=[1.0_dp,2.0_dp])
    call must(r%status==de_status_invalid, 'wrong nonlinear constraint shape')
    checks = checks + 7

    b(:,1)=[0.0_dp,0.0_dp,0.0_dp]
    b(:,2)=[2.0_dp,2.0_dp,2.0_dp]
    mask=[.true.,.true.,.true.]
    polish_calls = 0
    r=differential_evolution(mixed_objective,b,integrality=mask,maxiter=20, &
                             seed=6,polisher=mock_polisher)
    call must(polish_calls==0, 'all-integer solution skips polishing')
    checks = checks + 1

    print '(a,i0)', 'EXTENDED DE NATIVE CHECKS PASSED: ',checks
contains
    real(dp) function mixed_objective(x) result(value)
        real(dp), intent(in) :: x(:) !! three-dimensional test objective
        value = (x(1)-2.0_dp)**2+(x(2)+0.75_dp)**2+(x(3)-3.0_dp)**2
    end function mixed_objective
    real(dp) function linear_objective(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional test objective
        value = (x(1)-0.1_dp)**2+(x(2)-0.2_dp)**2
    end function linear_objective
    real(dp) function outside_objective(x) result(value)
        real(dp), intent(in) :: x(:) !! target outside unit disk
        value = (x(1)-1.5_dp)**2+(x(2)-1.3_dp)**2
    end function outside_objective
    real(dp) function combined_objective(x) result(value)
        real(dp), intent(in) :: x(:) !! three-dimensional test objective
        value = (x(1)-1.0_dp)**2+(x(2)-0.7_dp)**2+x(3)**2
    end function combined_objective
    subroutine circle_constraint(x,c)
        real(dp), intent(in) :: x(:) !! candidate vector
        real(dp), intent(out) :: c(:) !! squared distance from origin
        c(1) = sum(x**2)
    end subroutine circle_constraint
    subroutine circle_first_two(x,c)
        real(dp), intent(in) :: x(:) !! candidate vector
        real(dp), intent(out) :: c(:) !! first two squared coordinates
        c(1) = sum(x(1:2)**2)
    end subroutine circle_first_two
    subroutine nan_constraint(x,c)
        real(dp), intent(in) :: x(:) !! candidate vector
        real(dp), intent(out) :: c(:) !! NaN for excluded region
        c(1) = x(1)
        if (x(1) < 0.5_dp) c(1) = ieee_value(0.0_dp,ieee_quiet_nan)
    end subroutine nan_constraint
    subroutine mock_polisher(objective,x0,lo,hi,x,value,nfev,success)
        procedure(de_objective) :: objective !! objective function
        real(dp), intent(in) :: x0(:) !! starting point
        real(dp), intent(in) :: lo(:) !! lower bounds
        real(dp), intent(in) :: hi(:) !! upper bounds
        real(dp), intent(out) :: x(:) !! polished point
        real(dp), intent(out) :: value !! objective
        integer, intent(out) :: nfev !! evaluation count
        logical, intent(out) :: success !! valid output
        polish_calls = polish_calls + 1
        if (any(lo>hi)) error stop 'bad polish bounds'
        x=x0
        x(2)=-0.75_dp
        value=objective(x)
        nfev=1
        success=.true.
    end subroutine mock_polisher
    subroutine must(ok,message)
        logical, intent(in) :: ok !! condition to check
        character(len=*), intent(in) :: message !! failure diagnostic
        if (.not.ok) then
            write(*,'(a)') 'FAILED: '//message
            error stop 1
        end if
    end subroutine must
end program test_optimize_de_constraints_reference
