! Example for scipy.optimize.linprog-style dense linear programming.
program optimize_linprog_example
    use scifort_kinds, only : dp
    use scifort_optimize_linprog, only : linprog, linprog_result, linprog_optimal
    implicit none
    type(linprog_result) :: result
    real(dp) :: aub(2,2), b(2), c(2)

    ! Maximize x + 2*y subject to x+y <= 4, x+3*y <= 6, x,y >= 0.
    ! Since linprog minimizes, negate the objective coefficients.
    c = [-1.0_dp, -2.0_dp]
    aub(1,:) = [1.0_dp,1.0_dp]
    aub(2,:) = [1.0_dp,3.0_dp]
    b = [4.0_dp,6.0_dp]
    result = linprog(c,a_ub=aub,b_ub=b)
    if (result%status /= linprog_optimal) error stop 'linprog example failed'
    if (maxval(abs(result%x - [3.0_dp,1.0_dp])) > 1.0e-9_dp) &
        error stop 'unexpected linprog optimum'
    print '(a,2f12.6)', 'optimal x,y: ',result%x
    print '(a,f12.6)', 'maximum revenue: ', -result%fun
    print '(a,i0)', 'simplex pivots: ', result%nit
end program optimize_linprog_example
