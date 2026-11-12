! SPDX-License-Identifier: MIT
! Independent multi-object nonlinear constraint tests for SciFort SHGO.
program optimize_shgo_multiconstraints_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result, shgo_nonlinear_constraint, &
        shgo_constraint_blocks_valid, shgo_invalid, shgo_infeasible
    implicit none
    type(shgo_nonlinear_constraint) :: blocks(3), bad(2), empty(0)
    type(shgo_result) :: res
    real(dp) :: bounds(2,2), bounds3(3,2), nan, infinity
    integer :: checks, sum_jac_calls, diff_calls, objective_calls, vector_calls

    checks = 0
    sum_jac_calls = 0
    diff_calls = 0
    objective_calls = 0
    vector_calls = 0
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    infinity = ieee_value(0.0_dp,ieee_positive_inf)
    bounds(:,1) = 0.0_dp
    bounds(:,2) = 1.0_dp
    bounds3(:,1) = 0.0_dp
    bounds3(:,2) = 1.0_dp

    blocks(1)%fun => sum_constraint
    blocks(1)%jac => sum_jacobian
    blocks(1)%lb = [0.7_dp]
    blocks(1)%ub = [0.7_dp]
    blocks(2)%fun => difference_constraint
    ! No Jacobian: restoration should finite-difference only this block.
    blocks(2)%lb = [0.1_dp]
    blocks(2)%ub = [0.1_dp]
    blocks(3)%fun => vector_constraint
    blocks(3)%lb = [0.05_dp,-infinity]
    blocks(3)%ub = [infinity,0.35_dp]

    call assert(shgo_constraint_blocks_valid(blocks), 'valid mixed-derivative blocks')
    call assert(shgo_constraint_blocks_valid(empty), 'empty block sequence accepted')

    res = shgo_core(objective,bounds,iters=1,maxfun=150,local_iter=0, &
        constraints=blocks,restore_feasibility=.true.,restoration_starts=9)
    call assert(res%success .and. res%feasible, 'two nonlinear equalities jointly restored')
    call assert(all(abs(res%x-[0.4_dp,0.3_dp]) < 1.0e-7_dp), 'exact intersection recovered')
    call assert(res%maxcv <= 1.0e-8_dp, 'all four component bounds respected')
    call assert(res%nrestored > 0 .and. res%nrestore > 0, 'restoration attempts accounted')
    call assert(res%ncjev > 0, 'analytic block Jacobian used')
    call assert(res%ncjev == sum_jac_calls, 'only analytic callback counted')
    call assert(diff_calls > 0 .and. vector_calls > 0, 'separate block callbacks called')
    call assert(res%ncev > res%ncjev, 'finite difference evaluations counted')
    call assert(res%nfev == objective_calls, 'actual objective calls correctly counted')
    call assert(res%nfev <= 150, 'objective budget respected')
    call assert(abs(res%fun-objective(res%x)) <= 1.0e-13_dp, 'objective value matches x')
    call assert(res%nminima > 0, 'restored point recorded as local candidate')

    ! A second run with an analytic Jacobian for each block avoids the
    ! finite-difference fallback entirely.
    blocks(2)%jac => difference_jacobian
    blocks(3)%jac => vector_jacobian
    sum_jac_calls = 0
    res = shgo_core(objective,bounds,iters=1,maxfun=150,local_iter=0, &
        constraints=blocks,restore_feasibility=.true.,restoration_starts=9)
    call assert(res%success .and. res%maxcv <= 1.0e-8_dp, 'all-analytic restoration succeeds')
    call assert(all(abs(res%x-[0.4_dp,0.3_dp]) < 1.0e-7_dp), 'analytic Jacobian intersection')
    call assert(res%ncjev >= 3 .and. mod(res%ncjev,3) == 0, 'each analytic block counted')
    call assert(res%ncjev == 3*sum_jac_calls, 'per-callback Jacobian counting')
    nullify(blocks(2)%jac)
    nullify(blocks(3)%jac)

    ! An incompatible bound is detectable even if one block is fully feasible.
    blocks(2)%lb = [2.0_dp]
    blocks(2)%ub = [2.0_dp]
    res = shgo_core(objective,bounds,iters=1,maxfun=150,local_iter=0, &
        constraints=blocks,restore_feasibility=.true.)
    call assert(res%status == shgo_infeasible, 'unreachable equality returns infeasible')
    call assert(.not. res%success .and. res%nrestored == 0, 'no impossible point accepted')
    blocks(2)%lb = [0.1_dp]
    blocks(2)%ub = [0.1_dp]

    ! Fewer objective calls than the initial grid does not permit unbudgeted
    ! restoration or a fabricated feasible minimizer.
    objective_calls = 0
    res = shgo_core(objective,bounds,iters=1,maxfun=3,local_iter=0, &
        constraints=blocks,restore_feasibility=.true.)
    call assert(res%nfev <= 3 .and. res%nfev == objective_calls, 'strict small maxfun')

    ! Constraint collections can be combined with independent linear bounds.
    res = shgo_core(objective,bounds,iters=1,maxfun=150,local_iter=0, &
        constraints=blocks,linear_a=reshape([1.0_dp,1.0_dp],[1,2]), &
        linear_lb=[0.7_dp],linear_ub=[0.7_dp],restore_feasibility=.true.)
    call assert(res%success .and. res%maxcv <= 1.0e-8_dp, 'linear and block constraints coexist')

    ! Compatibility and validation: old single-vector callback is still valid.
    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0, &
        nonlinear=sum_constraint,nonlinear_lb=[0.7_dp], &
        nonlinear_ub=[0.7_dp],nonlinear_jac=sum_jacobian,restore_feasibility=.true.)
    call assert(res%success .and. res%maxcv <= 1.0e-8_dp, 'legacy callback still supported')
    res = shgo_core(objective,bounds,constraints=blocks,nonlinear=sum_constraint, &
        nonlinear_lb=[0.7_dp],nonlinear_ub=[0.7_dp])
    call assert(res%status == shgo_invalid, 'old and new nonlinear API cannot be mixed')
    res = shgo_core(objective,bounds,constraints=blocks,nonlinear_jac=sum_jacobian)
    call assert(res%status == shgo_invalid, 'legacy Jacobian conflicts with blocks')

    res = shgo_core(objective,bounds,iters=1,maxfun=100,local_iter=0,constraints=empty)
    call assert(res%success .and. res%ncev == 0, 'zero blocks behave as unconstrained')

    ! Validate missing procedures, allocation, inconsistent lengths, and NaNs.
    call assert(.not. shgo_constraint_blocks_valid(bad), 'unassociated procedure rejected')
    bad(1)%fun => sum_constraint
    call assert(.not. shgo_constraint_blocks_valid(bad), 'missing bound arrays rejected')
    bad(1)%lb = [0.0_dp]
    bad(1)%ub = [0.0_dp,1.0_dp]
    call assert(.not. shgo_constraint_blocks_valid(bad), 'mismatched bound vectors rejected')
    bad(1)%ub = [-1.0_dp]
    call assert(.not. shgo_constraint_blocks_valid(bad), 'inconsistent bounds rejected')
    bad(1)%ub = [nan]
    call assert(.not. shgo_constraint_blocks_valid(bad), 'NaN bound rejected')
    bad(1)%ub = [infinity]
    bad(1)%lb = [infinity]
    call assert(.not. shgo_constraint_blocks_valid(bad), 'positive infinite lower rejected')
    bad(1)%lb = [0.0_dp]
    bad(2)%fun => difference_constraint
    bad(2)%lb = [0.0_dp]
    bad(2)%ub = [-infinity]
    call assert(.not. shgo_constraint_blocks_valid(bad), 'negative infinite upper rejected')
    bad(2)%ub = [infinity]
    call assert(shgo_constraint_blocks_valid(bad), 'valid block list after repairs')
    res = shgo_core(objective,bounds,constraints=bad(1:1))
    call assert(res%status /= shgo_invalid, 'valid single block accepted')
    nullify(bad(2)%fun)
    res = shgo_core(objective,bounds,constraints=bad)
    call assert(res%status == shgo_invalid, 'invalid block reported without dereference')
    call assert(index(res%message,'constraint block') > 0, 'useful invalid-block message')

    ! In 3D, two independently specified surfaces with different Jacobians.
    blocks(1)%fun => sum_constraint
    blocks(1)%lb = [0.7_dp]
    blocks(1)%ub = [0.7_dp]
    blocks(2)%fun => first_constraint
    blocks(2)%jac => first_jacobian
    blocks(2)%lb = [0.2_dp]
    blocks(2)%ub = [0.2_dp]
    res = shgo_core(objective,bounds3,iters=1,maxfun=150,local_iter=0, &
        constraints=blocks(1:2),adaptive=.true.,conforming=.true., &
        restore_feasibility=.true.,restoration_starts=8)
    call assert(res%conforming_valid .and. res%ntetrahedra == 192, '3D topology preserved')
    call assert(res%success .and. res%maxcv <= 1.0e-8_dp, '3D multi-constraint feasibility')
    call assert(abs(sum(res%x)-0.7_dp) <= 1.0e-8_dp .and. &
        abs(res%x(1)-0.2_dp) <= 1.0e-8_dp, '3D independent equality residuals')

    print '(a,i0)', 'SHGO multiconstraint native assertions passed: ',checks
contains
    subroutine assert(ok,msg)
        logical, intent(in) :: ok !! expected condition
        character(len=*), intent(in) :: msg !! failure identifier
        checks = checks+1
        if (.not. ok) then
            print '(a,a)', 'FAIL: ', msg
            error stop 'SHGO multiconstraint failure'
        end if
    end subroutine assert

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp) :: f !! objective value
        objective_calls = objective_calls+1
        f = sum((x-0.25_dp)**2)
    end function objective

    subroutine sum_constraint(x,c)
        real(dp), intent(in) :: x(:) !! bounded variables
        real(dp), intent(out) :: c(:) !! single sum component
        c(1) = sum(x)
    end subroutine sum_constraint

    subroutine sum_jacobian(x,a)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: a(:,:) !! sum derivative
        a = 1.0_dp
        sum_jac_calls = sum_jac_calls + 1
    end subroutine sum_jacobian

    subroutine difference_constraint(x,c)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: c(:) !! first minus second
        diff_calls = diff_calls+1
        c(1) = x(1)-x(2)
    end subroutine difference_constraint

    subroutine difference_jacobian(x,a)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: a(:,:) !! gradient of difference
        a(1,:) = 0.0_dp
        a(1,1) = 1.0_dp
        a(1,2) = -1.0_dp
    end subroutine difference_jacobian

    subroutine vector_constraint(x,c)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: c(:) !! product and quadratic radius
        vector_calls = vector_calls+1
        c(1) = x(1)*x(2)
        c(2) = sum(x*x)
    end subroutine vector_constraint

    subroutine vector_jacobian(x,a)
        real(dp), intent(in) :: x(:) !! coordinates
        real(dp), intent(out) :: a(:,:) !! Jacobian with two rows
        a(1,:) = [x(2),x(1)]
        a(2,:) = 2.0_dp*x
    end subroutine vector_jacobian

    subroutine first_constraint(x,c)
        real(dp), intent(in) :: x(:) !! variables
        real(dp), intent(out) :: c(:) !! first coordinate
        c(1) = x(1)
    end subroutine first_constraint

    subroutine first_jacobian(x,a)
        real(dp), intent(in) :: x(:) !! variables
        real(dp), intent(out) :: a(:,:) !! gradient of x1
        a = 0.0_dp
        a(1,1) = 1.0_dp
    end subroutine first_jacobian
end program optimize_shgo_multiconstraints_reference
