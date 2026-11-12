! SPDX-License-Identifier: MIT
! Per-block BFGS/SR1 updates for nonlinear constraint Hessians.
program test_optimize_trust_constr_constraint_updates
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : linear_constraint, minimize, &
        multivariate_minimize_result, nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_constraint_updates_reference, only : &
        sr1_equality_x, bfgs_equality_x, sr1_bounded_x
    implicit none

    type(nonlinear_constraint) :: blocks(2)
    type(multivariate_minimize_result) :: r
    type(multivariate_minimize_result) :: sr1_short
    type(multivariate_minimize_result) :: bfgs_short
    real(dp) :: infinity
    integer :: failures
    integer :: exact_hessian_calls

    failures = 0
    exact_hessian_calls = 0
    infinity = ieee_value(0.0_dp, ieee_positive_inf)
    blocks(1)%fun => nonlinear_values
    blocks(1)%jac => nonlinear_jacobian
    blocks(1)%lower = [0.2_dp, -infinity]
    blocks(1)%upper = [0.2_dp, infinity]
    blocks(1)%hess_update = 'sr1'
    blocks(2)%fun => other_values
    blocks(2)%jac => other_jacobian
    blocks(2)%hess => other_hessian
    blocks(2)%lower = [0.1_dp]

    r = trust_constr(objective, [0.85_dp, 0.5_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, &
        maxiter=300, gtol=1.0e-9_dp)
    call check_solution('SR1 equality with exact other block', r, sr1_equality_x)
    call check_hessian_counts(r)
    if (exact_hessian_calls < 1) call fail('other block exact Hessian not used')

    blocks(1)%hess_update = 'BFGS'
    r = minimize(objective, [0.85_dp, 0.5_dp], method='trust-constr', &
        jac=gradient, hess=objective_hessian, constraint_objects=blocks, &
        maxiter=300, gtol=1.0e-9_dp)
    call check_solution('BFGS equality via minimize', r, bfgs_equality_x)
    call check_hessian_counts(r)

    blocks(1)%hess_update = 'sr1'
    blocks(1)%lower = [0.2_dp, 0.8_dp]
    r = trust_constr(objective, [0.85_dp, 0.5_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, &
        lower=[0.0_dp, 0.0_dp], upper=[1.5_dp, 1.5_dp], &
        bounds_keep_feasible=[.true., .true.], &
        maxiter=300, gtol=1.0e-9_dp)
    call check_solution('SR1 barrier with exact other block', r, sr1_bounded_x)
    call check_hessian_counts(r)
    if (size(r%constr_nhev) /= 3) call fail('bounds block count not appended')

    ! No exact objective Hessian is required with an updated constraint Hessian.
    r = minimize(objective, [0.85_dp, 0.5_dp], method='trust-constr', &
        jac=gradient, constraint_objects=blocks, &
        lower=[0.0_dp, 0.0_dp], upper=[1.5_dp, 1.5_dp], &
        maxiter=300, gtol=1.0e-8_dp)
    if (.not. r%success) call fail('SR1 constraint and BFGS objective did not converge')
    if (r%success) then
        if (maxval(abs(r%x - sr1_bounded_x)) > 1.0e-2_dp) &
            call fail('SR1 plus BFGS objective solution differs')
    end if

    ! Numeric constraint Jacobians remain usable with a Hessian update strategy.
    nullify(blocks(1)%jac)
    r = trust_constr(objective, [0.85_dp, 0.5_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, &
        lower=[0.0_dp, 0.0_dp], upper=[1.5_dp, 1.5_dp], &
        maxiter=3)
    if (r%status == -1) call fail('SR1 rejected numerical constraint Jacobian')
    if (allocated(r%constr_nfev)) then
        if (r%constr_nfev(1) < 3) call fail('no numerical Jacobian work recorded')
    end if
    blocks(1)%jac => nonlinear_jacobian

    ! A different update strategy must genuinely change the short-run path.
    blocks(1)%lower = [0.2_dp, -infinity]
    sr1_short = trust_constr(objective, [0.85_dp, 0.5_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=5)
    blocks(1)%hess_update = 'bfgs'
    bfgs_short = trust_constr(objective, [0.85_dp, 0.5_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=5)
    if (maxval(abs(sr1_short%x - bfgs_short%x)) < 1.0e-9_dp) &
        call fail('SR1 and BFGS updates gave identical short-run iterates')

    blocks(1)%hess_update = 'bad'
    r = trust_constr(objective, [0.85_dp, 0.5_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=1)
    if (r%status /= -1) call fail('unknown constraint Hessian strategy accepted')

    blocks(1)%hess_update = 'sr1'
    blocks(1)%hess => other_hessian
    r = trust_constr(objective, [0.85_dp, 0.5_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=1)
    if (r%status /= -1) call fail('SR1 conflicted with exact Hessian')
    nullify(blocks(1)%hess)

    blocks(1)%hess_fd_scheme = '3-point'
    r = trust_constr(objective, [0.85_dp, 0.5_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=1)
    if (r%status /= -1) call fail('SR1 conflicted with numerical Hessian')
    blocks(1)%hess_fd_scheme = ''

    blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[0.2_dp], upper=[0.2_dp])
    blocks(1)%hess_update = 'sr1'
    r = trust_constr(objective, [0.85_dp, 0.5_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=1)
    if (r%status /= -1) call fail('linear block accepted a Hessian update')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_constraint_updates: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_constraint_updates: PASS'

contains

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp) :: f !! scalar objective
        f = (x(1) - 0.55_dp)**2 + 2.0_dp * (x(2) - 0.3_dp)**2
    end function objective

    subroutine gradient(x, g)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: g(:) !! objective gradient
        g = [2.0_dp * (x(1) - 0.55_dp), 4.0_dp * (x(2) - 0.3_dp)]
    end subroutine gradient

    subroutine objective_hessian(x, h)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: h(:, :) !! objective Hessian
        if (size(x) /= 2) error stop 'invalid objective size'
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 4.0_dp
    end subroutine objective_hessian

    subroutine nonlinear_values(x, c)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: c(:) !! two constraint values
        c = [x(1)**2 - x(2)**2, x(1) + x(2)**2]
    end subroutine nonlinear_values

    subroutine nonlinear_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: a(:, :) !! two-row analytic Jacobian
        a(1, :) = [2.0_dp * x(1), -2.0_dp * x(2)]
        a(2, :) = [1.0_dp, 2.0_dp * x(2)]
    end subroutine nonlinear_jacobian

    subroutine other_values(x, c)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: c(:) !! exact-Hessian constraint value
        c = [x(1) * x(2)]
    end subroutine other_values

    subroutine other_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: a(:, :) !! exact-Hessian constraint Jacobian
        a(1, :) = [x(2), x(1)]
    end subroutine other_jacobian

    subroutine other_hessian(x, v, h)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(in) :: v(:) !! Lagrange multipliers
        real(dp), intent(out) :: h(:, :) !! weighted exact Hessian
        if (size(x) /= 2) error stop 'invalid Hessian argument'
        exact_hessian_calls = exact_hessian_calls + 1
        h = 0.0_dp
        h(1, 2) = v(1)
        h(2, 1) = v(1)
    end subroutine other_hessian

    subroutine check_solution(label, result, target)
        character(len=*), intent(in) :: label !! test case name
        type(multivariate_minimize_result), intent(in) :: result !! optimizer result
        real(dp), intent(in) :: target(:) !! SciPy reference minimizer
        if (.not. result%success) then
            print '(a,1x,a,i0,1x,a)', trim(label), 'status=', result%status, trim(result%message)
            call fail('failed to converge')
            return
        end if
        if (maxval(abs(result%x - target)) > 5.0e-3_dp) then
            print '(a,2(1x,f12.7))', trim(label), result%x
            call fail('minimum differs from SciPy')
        end if
        if (result%constr_violation > 1.0e-6_dp) call fail('constraint violation too large')
    end subroutine check_solution

    subroutine check_hessian_counts(result)
        type(multivariate_minimize_result), intent(in) :: result !! optimizer result
        if (.not. allocated(result%constr_nhev)) then
            call fail('no per-block Hessian accounting')
            return
        end if
        if (result%constr_nhev(1) /= 0) call fail('update Hessian was counted as exact')
        if (result%constr_nhev(2) < 1) call fail('exact Hessian call was not accounted')
        if (result%constr_njev(1) <= 1) call fail('update Jacobian work not accounted')
        if (result%constr_njev(1) /= result%constr_nfev(1)) &
            call fail('redundant Jacobian calls during Hessian update')
    end subroutine check_hessian_counts

    subroutine fail(message)
        character(len=*), intent(in) :: message !! failure description
        print '(a)', 'FAIL: ' // message
        failures = failures + 1
    end subroutine fail

end program test_optimize_trust_constr_constraint_updates
