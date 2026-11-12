! SPDX-License-Identifier: MIT
program test_optimize_milp
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_milp, only : milp, milp_result, milp_optimal, &
        milp_node_limit, milp_infeasible, milp_unbounded, milp_numerical, milp_invalid
    implicit none
    type(milp_result) :: r
    real(dp) :: inf
    integer :: checks
    checks = 0
    inf = ieee_value(1.0_dp, ieee_positive_inf)

    ! Published SciPy example: maximize y subject to three inequalities.
    r = milp([0.0_dp,-1.0_dp], integrality=[1,1], &
        a=reshape([-1.0_dp,3.0_dp,2.0_dp,1.0_dp,2.0_dp,3.0_dp],[3,2]), &
        b_l=[-inf,-inf,-inf], b_u=[1.0_dp,12.0_dp,12.0_dp], &
        mip_rel_gap=0.0_dp)
    call check(r%status==milp_optimal, 'published SciPy example status')
    call check(abs(r%fun+2.0_dp)<1e-7_dp, 'published SciPy example objective')
    call check(maxval(abs(r%x-[1.0_dp,2.0_dp]))<1e-7_dp, 'published SciPy example permutation')

    r = milp([-10.0_dp,-6.0_dp], integrality=[1,1], &
        a=reshape([5.0_dp,4.0_dp],[1,2]), b_l=[-inf], b_u=[8.0_dp], &
        upper=[1.0_dp,1.0_dp], mip_rel_gap=0.0_dp)
    call check(r%status==milp_optimal, 'binary knapsack optimal')
    call check(r%success, 'binary knapsack success')
    call check(allocated(r%x), 'binary knapsack incumbent')
    call check(abs(r%fun+10.0_dp) < 1e-7_dp, 'binary knapsack objective')
    call check(maxval(abs(r%x-[1.0_dp,0.0_dp])) < 1e-7_dp, 'binary knapsack solution')
    call check(r%mip_gap <= 1e-12_dp, 'certified zero gap')
    call check(r%mip_node_count > 0, 'count solved nodes')

    r = milp([-5.0_dp,-3.0_dp], integrality=[1,0], &
        a=reshape([2.0_dp,1.0_dp],[1,2]), b_l=[-inf], b_u=[4.5_dp], &
        upper=[4.0_dp,5.0_dp], mip_rel_gap=0.0_dp)
    call check(r%status==milp_optimal, 'mixed integer continuous optimal')
    call check(abs(r%fun+13.5_dp) < 1e-7_dp, 'mixed objective')
    call check(abs(r%x(1)) < 1e-7_dp, 'mixed integer coordinate')
    call check(abs(r%x(2)-4.5_dp) < 1e-7_dp, 'mixed continuous coordinate')

    r = milp([2.0_dp,1.0_dp], integrality=[1,1], &
        a=reshape([1.0_dp,1.0_dp],[1,2]), b_l=[3.0_dp], b_u=[3.0_dp], &
        upper=[3.0_dp,3.0_dp], mip_rel_gap=0.0_dp)
    call check(r%status==milp_optimal, 'integer equality optimal')
    call check(abs(r%fun-3.0_dp)<1e-7_dp, 'integer equality objective')
    call check(maxval(abs(r%x-[0.0_dp,3.0_dp]))<1e-7_dp, 'integer equality solution')

    r = milp([1.0_dp], integrality=[2], lower=[3.0_dp], upper=[5.0_dp])
    call check(r%status==milp_optimal, 'semi-continuous zero status')
    call check(abs(r%x(1))<1e-7_dp, 'semi-continuous zero allowed')
    r = milp([1.0_dp], integrality=[2], lower=[3.0_dp], upper=[5.0_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[2.0_dp], b_u=[inf])
    call check(r%status==milp_optimal, 'semi-continuous nonzero status')
    call check(abs(r%fun-3.0_dp)<1e-7_dp, 'semi-continuous gap respected')
    r = milp([1.0_dp], integrality=[3], lower=[2.2_dp], upper=[5.8_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[1.0_dp], b_u=[inf])
    call check(r%status==milp_optimal, 'semi-integer nonzero status')
    call check(abs(r%fun-3.0_dp)<1e-7_dp, 'semi-integer integrality')
    r = milp([-1.0_dp], integrality=[3], lower=[2.2_dp], upper=[5.8_dp])
    call check(r%status==milp_optimal, 'semi-integer maximize status')
    call check(abs(r%x(1)-5.0_dp)<1e-7_dp, 'semi-integer maximize')
    r = milp([1.0_dp], integrality=[2], lower=[-5.0_dp], upper=[-2.0_dp])
    call check(r%status==milp_optimal, 'negative semi-continuous status')
    call check(abs(r%x(1)+5.0_dp)<1e-7_dp, 'negative semi-continuous allowed')
    r = milp([-1.0_dp], integrality=[2], lower=[-5.0_dp], upper=[-2.0_dp])
    call check(r%status==milp_optimal, 'negative semi-continuous zero status')
    call check(abs(r%x(1))<1e-7_dp, 'negative semi-continuous zero')

    r = milp([1.0_dp], integrality=[1], a=reshape([1.0_dp],[1,1]), &
        b_l=[0.5_dp], b_u=[0.5_dp])
    call check(r%status==milp_infeasible, 'fractional integer equality infeasible')
    call check(.not. allocated(r%x), 'no incumbent on infeasible equality')
    r = milp([1.0_dp], integrality=[1], lower=[0.2_dp], upper=[0.8_dp])
    call check(r%status==milp_infeasible, 'empty integer interval infeasible')
    r = milp([1.0_dp], integrality=[3], lower=[2.2_dp], upper=[2.8_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[1.0_dp], b_u=[inf])
    call check(r%status==milp_infeasible, 'empty semi-integer positive interval')
    r = milp([1.0_dp], integrality=[1], upper=[1.0e12_dp+2.0_dp], &
        a=reshape([1.0_dp],[1,1]), &
        b_l=[1.0e12_dp+0.5_dp], b_u=[1.0e12_dp+0.5_dp])
    call check(r%status==milp_infeasible, 'large-scale fractional equality still infeasible')
    r = milp([1.0_dp], integrality=[1], upper=[1.0e12_dp+2.0_dp], &
        a=reshape([1.0_dp],[1,1]), &
        b_l=[1.0e12_dp+1.0_dp], b_u=[1.0e12_dp+1.0_dp])
    call check(r%status==milp_optimal, 'large-scale integer equality optimal')
    call check(abs(r%x(1)-(1.0e12_dp+1.0_dp))<1e-6_dp, 'large-scale exact integer')

    r = milp([0.0_dp], a=reshape([1.0_dp],[1,1]), b_l=[2.0_dp], b_u=[1.0_dp])
    call check(r%status==milp_invalid, 'invalid inconsistent row bounds')
    r = milp([0.0_dp], a=reshape([1.0_dp],[1,1]), b_l=[1.0_dp], b_u=[1.0_dp], &
        upper=[0.0_dp])
    call check(r%status==milp_infeasible, 'continuous infeasible constraints')

    r = milp([-1.0_dp], lower=[0.0_dp], upper=[inf])
    call check(r%status==milp_unbounded, 'continuous unbounded')
    r = milp([-1.0_dp], integrality=[1])
    call check(r%status==milp_numerical, 'integer unbounded relaxation inconclusive')
    r = milp([1.0_dp,2.0_dp], a=reshape([1.0_dp,1.0_dp],[1,2]), &
        b_l=[2.0_dp], b_u=[inf])
    call check(r%status==milp_optimal, 'all-continuous nonempty row lower')
    call check(abs(r%fun-2.0_dp)<1e-7_dp, 'all-continuous objective')
    call check(r%mip_node_count==1, 'continuous one node')
    r = milp([1.0_dp], a=reshape([1.0_dp],[1,1]), b_l=[-inf], b_u=[inf])
    call check(r%status==milp_optimal, 'unbounded row skipped')

    r = milp([-1.0_dp], integrality=[1], upper=[2.0_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[-inf], b_u=[1.5_dp], &
        node_limit=0, mip_rel_gap=0.0_dp)
    call check(r%status==milp_node_limit, 'zero node budget status')
    call check(r%mip_node_count==0, 'zero node budget counted')
    call check(.not. allocated(r%x), 'zero budget no solution')
    r = milp([-1.0_dp], integrality=[1], upper=[2.0_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[-inf], b_u=[1.5_dp], &
        node_limit=1, mip_rel_gap=0.0_dp)
    call check(r%status==milp_node_limit, 'one node limit')
    call check(r%mip_node_count==1, 'one node used')
    call check(abs(r%mip_dual_bound+1.5_dp)<1e-7_dp, 'limit valid dual bound')
    r = milp([-1.0_dp], integrality=[1], upper=[2.0_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[-inf], b_u=[1.5_dp], &
        node_limit=2, mip_rel_gap=0.0_dp)
    call check(r%status==milp_node_limit, 'two node limit')
    call check(allocated(r%x), 'limited search retains incumbent')
    call check(abs(r%fun+1.0_dp)<1e-7_dp, 'limited incumbent objective')
    call check(r%mip_gap > 0.0_dp, 'positive live gap')
    r = milp([-1.0_dp], integrality=[1], upper=[2.0_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[-inf], b_u=[1.5_dp], &
        node_limit=2, mip_rel_gap=0.6_dp)
    call check(r%status==milp_optimal, 'requested relative gap termination')
    call check(r%mip_gap <= 0.6_dp, 'gap diagnostic limited')
    r = milp([-1.0_dp], integrality=[1], upper=[2.0_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[-inf], b_u=[1.5_dp], &
        node_limit=3, mip_rel_gap=0.0_dp)
    call check(r%status==milp_optimal, 'proof completes in three nodes')
    call check(r%mip_gap <= 1e-12_dp, 'proof gap zero')
    r = milp([-1.0_dp], integrality=[1], upper=[2.0_dp], &
        a=reshape([1.0_dp],[1,1]), b_l=[-inf], b_u=[1.5_dp], &
        time_limit=0.0_dp)
    call check(r%status==milp_node_limit, 'zero CPU budget')

    r = milp([1.0_dp], integrality=[-1])
    call check(r%status==milp_invalid, 'negative integrality')
    r = milp([1.0_dp], integrality=[4])
    call check(r%status==milp_invalid, 'unsupported integrality')
    r = milp([1.0_dp], integrality=[1,0])
    call check(r%status==milp_invalid, 'wrong integrality shape')
    r = milp([1.0_dp], lower=[2.0_dp], upper=[1.0_dp])
    call check(r%status==milp_invalid, 'invalid variable bounds')
    r = milp([1.0_dp], lower=[inf])
    call check(r%status==milp_invalid, 'invalid positive infinite lower')
    r = milp([1.0_dp], upper=[-inf])
    call check(r%status==milp_invalid, 'invalid negative infinite upper')
    r = milp([1.0_dp], a=reshape([1.0_dp],[1,1]), b_l=[0.0_dp])
    call check(r%status==milp_invalid, 'missing constraint upper')
    r = milp([1.0_dp], a=reshape([1.0_dp],[1,1]), b_u=[0.0_dp])
    call check(r%status==milp_invalid, 'missing constraint lower')
    r = milp([1.0_dp], a=reshape([1.0_dp],[1,1]), b_l=[inf], b_u=[inf])
    call check(r%status==milp_invalid, 'invalid row positive infinite lower')
    r = milp([1.0_dp], node_limit=-1)
    call check(r%status==milp_invalid, 'negative node limit')
    r = milp([1.0_dp], mip_rel_gap=-0.1_dp)
    call check(r%status==milp_invalid, 'negative gap')
    r = milp([1.0_dp], tol=0.0_dp)
    call check(r%status==milp_invalid, 'zero tolerance')
    r = milp([1.0_dp], lp_maxiter=-1)
    call check(r%status==milp_invalid, 'negative LP iterations')
    r = milp([ieee_value(0.0_dp,ieee_quiet_nan)])
    call check(r%status==milp_invalid, 'nonfinite objective')
    r = milp([1.0_dp], a=reshape([ieee_value(0.0_dp,ieee_quiet_nan)],[1,1]), &
        b_l=[0.0_dp],b_u=[inf])
    call check(r%status==milp_invalid, 'NaN matrix')
    r = milp([1.0_dp], upper=[ieee_value(0.0_dp,ieee_quiet_nan)])
    call check(r%status==milp_invalid, 'NaN upper bound')
    r = milp([real(dp)::])
    call check(r%status==milp_optimal, 'zero variable LP')
    call check(allocated(r%x), 'zero variable x allocated')
    call check(size(r%x)==0, 'zero variable x length')

    print '(a,i0)', 'milp native checks passed: ', checks
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition !! assertion value
        character(len=*), intent(in) :: label !! test label
        if (.not. condition) then
            print '(a,a)', 'FAILED: ',label
            error stop 1
        end if
        checks = checks + 1
    end subroutine check
end program test_optimize_milp
