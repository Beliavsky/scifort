! SPDX-License-Identifier: MIT
program optimize_milp_example
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_milp, only : milp, milp_result, milp_optimal
    implicit none
    type(milp_result) :: res
    real(dp) :: inf

    inf = ieee_value(1.0_dp, ieee_positive_inf)

    ! Published SciPy example: maximize x(2) using integer decisions.
    ! SciFort minimizes c^T*x, so the second coefficient is negative.
    res = milp([0.0_dp,-1.0_dp], integrality=[1,1], &
        a=reshape([-1.0_dp,3.0_dp,2.0_dp,1.0_dp,2.0_dp,3.0_dp],[3,2]), &
        b_l=[-inf,-inf,-inf], b_u=[1.0_dp,12.0_dp,12.0_dp], &
        mip_rel_gap=0.0_dp)
    if (res%status /= milp_optimal) error stop 'integer example failed'
    print '(a,2(f8.3,1x))', 'integer x: ',res%x
    print '(a,f10.4)', 'minimum objective: ',res%fun
    print '(a,i0)', 'LP relaxations solved: ',res%mip_node_count

    ! A semi-continuous variable may be exactly zero or between 3 and 5.
    ! A separate row requires x >= 2, so the optimum must jump to x=3.
    res = milp([1.0_dp], integrality=[2], lower=[3.0_dp], upper=[5.0_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[2.0_dp], b_u=[inf])
    if (res%status /= milp_optimal) error stop 'semi-continuous example failed'
    print '(a,f8.3)', 'semi-continuous x: ',res%x(1)
end program optimize_milp_example
