! SPDX-License-Identifier: MIT
program optimize_dual_annealing_example
    use scifort_kinds, only : dp
    use scifort_optimize_dual_annealing, only : dual_annealing_core, dual_annealing_result
    implicit none
    type(dual_annealing_result) :: res
    real(dp) :: bounds(2,2)

    bounds(:,1) = [-4.0_dp,-4.0_dp]
    bounds(:,2) = [ 4.0_dp, 4.0_dp]
    res = dual_annealing_core(objective, bounds, maxiter=400, seed=49)
    if (.not. res%success) error stop 'dual annealing did not finish'
    print '(a,2(1x,f12.6))', 'best position:', res%x
    print '(a,1x,es16.8)', 'best objective:', res%fun
    print '(a,1x,i0)', 'function evaluations:', res%nfev
    print '(a,1x,i0)', 'accepted moves:', res%naccepted
    print '(a,1x,i0)', 'temperature restarts:', res%nrestarts
    if (res%fun > 1.0e-2_dp) error stop 'example optimization did not improve enough'

contains
    function objective(x) result(value)
        real(dp), intent(in) :: x(:) !! candidate variables, shape (2)
        real(dp) :: value !! shifted quadratic objective
        value = (x(1) - 0.75_dp)**2 + (x(2) + 1.25_dp)**2
    end function objective
end program optimize_dual_annealing_example
