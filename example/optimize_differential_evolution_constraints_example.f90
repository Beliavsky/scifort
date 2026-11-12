! SPDX-License-Identifier: MIT
program optimize_differential_evolution_constraints_example
    use scifort_kinds, only : dp
    use scifort_optimize_differential_evolution, only : differential_evolution, &
        differential_evolution_result
    implicit none
    type(differential_evolution_result) :: result
    real(dp) :: bounds(3,2), a(1,3)
    logical :: integer_variables(3)

    bounds(:,1) = [-2.0_dp, -2.0_dp, -2.0_dp]
    bounds(:,2) = [ 2.0_dp,  2.0_dp,  2.0_dp]
    integer_variables = [.true., .false., .false.]
    a(1,:) = [1.0_dp, 1.0_dp, 0.0_dp]

    result = differential_evolution(objective, bounds, seed=2026, maxiter=250, &
        popsize=15, integrality=integer_variables, &
        linear_a=a, linear_lb=[0.5_dp], linear_ub=[huge(1.0_dp)], &
        nonlinear=unit_disk, nonlinear_lb=[-huge(1.0_dp)], nonlinear_ub=[0.81_dp])
    print '(a,3f12.6)', 'Solution: ', result%x
    print '(a,f15.9)', 'Objective: ', result%fun
    print '(a,l1,1x,es12.4)', 'Feasible and maxcv: ', result%feasible, result%maxcv
    print '(a,i0,1x,i0)', 'Objective/constraint evaluations: ', result%nfev, result%ncev
    if (.not. result%feasible .or. result%fun > 1.01_dp) error stop 'unexpected solution'

contains
    real(dp) function objective(x) result(value)
        real(dp), intent(in) :: x(:) !! decision vector, size 3
        value = (x(1)-1.0_dp)**2 + (x(2)-0.7_dp)**2 + x(3)**2
    end function objective

    subroutine unit_disk(x, c)
        real(dp), intent(in) :: x(:) !! decision vector, size 3
        real(dp), intent(out) :: c(:) !! squared disk radius, size 1
        c(1) = x(1)**2 + x(2)**2
    end subroutine unit_disk
end program optimize_differential_evolution_constraints_example
