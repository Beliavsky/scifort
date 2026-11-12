! SPDX-License-Identifier: MIT
program optimize_basinhopping_example
    use scifort_kinds, only : dp
    use scifort_optimize, only : basinhopping, basinhopping_result
    implicit none
    type(basinhopping_result) :: res

    res = basinhopping(doublewell, [2.4_dp], niter=80, seed=47, &
        stepsize=1.5_dp, temperature=0.4_dp, method='nelder-mead')
    if (.not.res%success) error stop 'basinhopping failed'
    print '(a,es16.8)', 'best objective: ', res%fun
    print '(a,es16.8)', 'best parameter: ', res%x(1)
    print '(a,i0)', 'steps attempted: ', res%nit
    print '(a,i0)', 'accepted moves: ', res%naccepted
    print '(a,i0)', 'local solves: ', res%nlocal
    if (res%fun > -0.2_dp) error stop 'did not reach deep minimum'
contains
    function doublewell(x) result(fx)
        real(dp), intent(in) :: x(:)
        real(dp) :: fx
        fx = (x(1)**2 - 1.0_dp)**2 + 0.25_dp*x(1)
    end function doublewell
end program optimize_basinhopping_example
