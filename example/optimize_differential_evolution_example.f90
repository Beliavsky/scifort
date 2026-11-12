! SPDX-License-Identifier: MIT
program optimize_differential_evolution_example
    use scifort_kinds, only : dp
    use scifort_optimize_differential_evolution, only : differential_evolution, differential_evolution_result
    implicit none
    type(differential_evolution_result) :: res
    real(dp) :: bounds(2,2)

    bounds(:,1) = [-5.12_dp,-5.12_dp]
    bounds(:,2) = [ 5.12_dp, 5.12_dp]
    res = differential_evolution(rastrigin, bounds, strategy='best1bin', &
        seed=2026, maxiter=250, popsize=15)
    print '(a,l1)', 'converged: ',res%success
    print '(a,es14.6)', 'fun: ',res%fun
    print '(a,2f12.7)', 'x: ',res%x
    print '(a,i0,a,i0)', 'generations: ',res%nit,', evaluations: ',res%nfev
contains
    function rastrigin(x) result(f)
        real(dp), intent(in) :: x(:) !! candidate vector
        real(dp) :: f !! objective to minimize
        real(dp), parameter :: pi = acos(-1.0_dp)
        f = 10.0_dp*real(size(x),dp) + sum(x**2 - 10.0_dp*cos(2.0_dp*pi*x))
    end function rastrigin
end program optimize_differential_evolution_example
