program test_optimize_linprog
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_linprog, only : linprog, linprog_result, &
        linprog_optimal, linprog_infeasible, linprog_unbounded, linprog_iteration_limit, linprog_invalid
    implicit none
    type(linprog_result) :: res
    real(dp) :: inf
    integer :: checks
    checks = 0
    inf = ieee_value(1.0_dp,ieee_positive_inf)

    ! Bounded LP with active inequality optimum x=(4,0).
    res = linprog([-1.0_dp,4.0_dp], a_ub=reshape([-3.0_dp,1.0_dp,-1.0_dp,2.0_dp], [2,2]), &
        b_ub=[6.0_dp,4.0_dp], lower=[-3.0_dp,0.0_dp], upper=[inf,inf])
    call check(res%status == linprog_optimal, 'simple bounded LP status')
    call check(abs(res%fun+4.0_dp)<1e-7_dp, 'simple bounded LP objective')

    res = linprog([1.0_dp,1.0_dp], a_eq=reshape([1.0_dp,1.0_dp],[1,2]), b_eq=[3.0_dp])
    call check(res%status == linprog_optimal, 'equality status')
    call check(abs(res%fun-3.0_dp)<1e-7_dp, 'equality objective')
    call check(maxval(abs(res%con))<1e-7_dp, 'equality residual')

    res = linprog([1.0_dp, 2.0_dp], lower=[-inf,0.0_dp], upper=[inf,5.0_dp], &
        a_eq=reshape([1.0_dp,1.0_dp],[1,2]), b_eq=[1.0_dp])
    call check(res%status == linprog_optimal, 'free variable status')
    call check(abs(res%fun-1.0_dp)<1e-7_dp, 'free variable objective')

    res = linprog([1.0_dp], lower=[-inf], upper=[3.0_dp])
    call check(res%status==linprog_unbounded, 'unbounded above transformed variable')

    res = linprog([1.0_dp], a_ub=reshape([1.0_dp,-1.0_dp],[2,1]), &
        b_ub=[0.0_dp,-1.0_dp])
    call check(res%status==linprog_infeasible, 'infeasible contradictory inequalities')

    res = linprog([-1.0_dp],lower=[0.0_dp],upper=[inf])
    call check(res%status==linprog_unbounded, 'unbounded objective')

    res = linprog([1.0_dp, 2.0_dp], &
        a_eq=reshape([1.0_dp,2.0_dp,1.0_dp,2.0_dp],[2,2]), b_eq=[3.0_dp,6.0_dp])
    call check(res%status==linprog_optimal, 'redundant equality status')
    call check(abs(res%fun-3.0_dp)<1e-7_dp, 'redundant equality objective')

    res = linprog([1.0_dp], lower=[-5.0_dp],upper=[-2.0_dp])
    call check(res%status==linprog_optimal, 'negative bounds status')
    call check(abs(res%x(1)+5.0_dp)<1e-7_dp, 'negative bounds optimum')

    res = linprog([0.0_dp], a_eq=reshape([0.0_dp],[1,1]), b_eq=[0.0_dp])
    call check(res%status==linprog_optimal, 'redundant zero equality')
    res = linprog([0.0_dp], a_eq=reshape([0.0_dp],[1,1]), b_eq=[1.0_dp])
    call check(res%status==linprog_infeasible, 'inconsistent zero equality')
    res = linprog([1.0_dp],maxiter=0)
    call check(res%status==linprog_optimal, 'already optimal no iterations')
    res = linprog([-1.0_dp], a_ub=reshape([1.0_dp],[1,1]),b_ub=[2.0_dp],maxiter=0)
    call check(res%status==linprog_iteration_limit, 'iteration limit')

    res = linprog([-1.0_dp], a_ub=reshape([1.0_dp],[1,1]), &
        b_ub=[1.0e12_dp])
    call check(res%status==linprog_optimal, 'large RHS status')
    call check(abs(res%x(1)/1.0e12_dp-1.0_dp)<1e-10_dp, 'large RHS pivot')

    res = linprog([0.0_dp], a_ub=reshape([1.0_dp,-1.0_dp],[2,1]), &
        b_ub=[1.0e12_dp,-(1.0e12_dp+1.0_dp)])
    call check(res%status==linprog_infeasible, 'large RHS infeasibility')

    ! Phase I must select the true minimum ratio, not treat order-one
    ! differences around 1e12 as ties due to the simplex pivot tolerance.
    res = linprog([1.0_dp], a_eq=reshape([1.0_dp],[1,1]), &
        b_eq=[1.0e12_dp+1.0_dp], lower=[0.0_dp], upper=[1.0e12_dp+2.0_dp])
    call check(res%status==linprog_optimal, 'large RHS feasible equality')
    call check(abs(res%x(1)-(1.0e12_dp+1.0_dp))<1e-6_dp, 'large RHS equality solution')

    res = linprog([1.0_dp], lower=[1.0_dp],upper=[0.0_dp])
    call check(res%status==linprog_invalid, 'inconsistent bounds')
    res = linprog([1.0_dp], a_ub=reshape([1.0_dp],[1,1]))
    call check(res%status==linprog_invalid, 'missing inequality RHS')
    res = linprog([1.0_dp,2.0_dp],a_ub=reshape([1.0_dp],[1,1]),b_ub=[1.0_dp])
    call check(res%status==linprog_invalid, 'wrong A_ub shape')
    res = linprog([ieee_value(1.0_dp,ieee_quiet_nan)])
    call check(res%status==linprog_invalid, 'nonfinite objective')
    res = linprog([-1.0_dp],a_ub=reshape([1.0_dp],[1,1]),b_ub=[1.0_dp],maxiter=1000)
    call check(res%status==linprog_optimal, 'bounded LP with iteration budget')

    res = linprog([real(dp)::])
    call check(res%status==linprog_optimal, 'zero variable LP')
    print '(a,i0)', 'linprog native checks passed: ',checks
contains
    subroutine check(condition,label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            print '(a,a)', 'FAILED: ',label
            error stop 1
        end if
        checks=checks+1
    end subroutine check
end program test_optimize_linprog
