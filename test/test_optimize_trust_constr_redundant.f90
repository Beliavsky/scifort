! SPDX-License-Identifier: MIT
! Overdetermined equality regressions against SciPy 1.17 SVDFactorization.
program test_optimize_trust_constr_redundant
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : linear_constraint, minimize, &
        multivariate_minimize_result, nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_redundant_reference, only : &
        linear_x, linear_nit, linear_status, &
        bounded_x, bounded_nit, bounded_status, &
        nonlinear_x, nonlinear_nit, nonlinear_status, &
        fixed_x, fixed_status, &
        inconsistent_x, inconsistent_cv, inconsistent_status
    implicit none

    type(multivariate_minimize_result) :: r
    type(nonlinear_constraint) :: blocks(2)
    real(dp) :: posinf
    integer :: failures
    integer :: hess_calls

    failures = 0
    hess_calls = 0
    posinf = ieee_value(0.0_dp, ieee_positive_inf)

    ! More equality rows than variables is legal when the extra rows are
    ! compatible.  Gram-pseudoinverse projection replaces the old size error.
    r = trust_constr(objective, [0.0_dp, 0.0_dp], jac=gradient, hess=hessian, &
        constraints=linear_residual, constraint_jac=linear_jacobian, &
        meq=3, maxiter=150)
    call check_solution('overdetermined linear', r, linear_x, linear_status)
    if (r%nit /= linear_nit) call fail('linear iteration count')
    if (trim(r%method) /= 'equality_constrained_sqp') call fail('linear solver method')
    if (size(r%v) /= 3 .or. size(r%constr_jac, 1) /= 3) &
        call fail('redundant linear multiplier dimensions')

    ! Independent constraint objects and variable bounds exercise the
    ! interior-point conversion with more equalities than variables.
    blocks(1) = linear_constraint(reshape([1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp], [2, 2]), &
        lower=[3.0_dp, 6.0_dp], upper=[3.0_dp, 6.0_dp])
    blocks(2) = linear_constraint(reshape([1.0_dp, -1.0_dp], [1, 2]), &
        lower=[1.0_dp], upper=[1.0_dp])
    r = minimize(objective, [0.0_dp, 0.0_dp], method='trust-constr', &
        jac=gradient, hess=hessian, constraint_objects=blocks, &
        lower=[0.0_dp, 0.0_dp], maxiter=150)
    call check_solution('overdetermined with bounds', r, bounded_x, bounded_status)
    if (r%nit /= bounded_nit) call fail('bounded iteration count')
    if (trim(r%method) /= 'tr_interior_point') call fail('bounded solver method')
    if (size(r%constr_nfev) /= 3) then
        call fail('bounded constraint block count')
    else if (any(r%constr_nfev /= 0) .or. any(r%constr_njev /= 0) .or. &
            any(r%constr_nhev /= 0)) then
        call fail('linear/bounds blocks evaluated callbacks')
    end if

    ! The same rank-deficient projection handles redundant nonlinear
    ! equalities and preserves analytic weighted constraint Hessians.
    hess_calls = 0
    r = trust_constr(objective, [1.8_dp, 0.8_dp], jac=gradient, hess=hessian, &
        constraints=nonlinear_residual, constraint_jac=nonlinear_jacobian, &
        constraint_hess=nonlinear_hessian, meq=3, maxiter=150)
    call check_solution('overdetermined nonlinear', r, nonlinear_x, nonlinear_status)
    if (r%nit /= nonlinear_nit) call fail('nonlinear iteration count')
    if (hess_calls < 1) call fail('nonlinear Hessian not invoked')
    if (size(r%constr_nhev) /= 1) then
        call fail('nonlinear Hessian counter shape')
    else if (r%constr_nhev(1) /= hess_calls) then
        call fail('nonlinear Hessian evaluation count')
    end if

    ! A fixed variable is itself another equality, which can make the
    ! canonicalized system even more overdetermined.
    r = trust_constr(objective, [2.0_dp, 1.0_dp], jac=gradient, hess=hessian, &
        constraints=linear_residual, constraint_jac=linear_jacobian, &
        meq=3, lower=[2.0_dp, 0.0_dp], upper=[2.0_dp, posinf], maxiter=150)
    call check_solution('fixed bound and redundant equations', r, fixed_x, fixed_status)

    ! Rank deficiency does not mean feasibility: incompatible equations
    ! must not be declared solved solely because optimality is small.
    r = trust_constr(objective, [0.0_dp, 0.0_dp], jac=gradient, hess=hessian, &
        constraints=inconsistent_residual, constraint_jac=linear_jacobian, &
        meq=3, maxiter=150)
    if (r%success .or. r%status /= inconsistent_status) &
        call fail('inconsistent equality system accepted')
    if (maxval(abs(r%x - inconsistent_x)) > 2.0e-6_dp) &
        call fail('inconsistent equality least-squares point')
    if (abs(r%constr_violation - inconsistent_cv) > 1.0e-6_dp) &
        call fail('inconsistent constraint violation')

    ! BFGS and SR1 objective approximations also use the rank-deficient
    ! projection when exact objective Hessians are not provided.
    r = trust_constr(objective, [0.0_dp, 0.0_dp], jac=gradient, &
        constraints=linear_residual, constraint_jac=linear_jacobian, &
        meq=3, maxiter=150)
    call check_solution('overdetermined BFGS', r, linear_x, linear_status)
    r = trust_constr(objective, [0.0_dp, 0.0_dp], jac=gradient, &
        constraints=linear_residual, constraint_jac=linear_jacobian, &
        meq=3, objective_hess_update='sr1', maxiter=150)
    call check_solution('overdetermined SR1', r, linear_x, linear_status)

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_redundant: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_redundant: PASS'

contains

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp) :: f !! objective value
        f = (x(1) - 0.1_dp)**2 + 2.0_dp * (x(2) - 0.3_dp)**2
    end function objective

    subroutine gradient(x, g)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: g(:) !! objective gradient
        g = [2.0_dp * (x(1) - 0.1_dp), 4.0_dp * (x(2) - 0.3_dp)]
    end subroutine gradient

    subroutine hessian(x, h)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: h(:, :) !! objective Hessian
        if (size(x) /= 2) error stop 'test objective must be two-dimensional'
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 4.0_dp
    end subroutine hessian

    subroutine linear_residual(x, c)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: c(:) !! redundant equality residuals
        c = [sum(x) - 3.0_dp, 2.0_dp * sum(x) - 6.0_dp, x(1) - x(2) - 1.0_dp]
    end subroutine linear_residual

    subroutine inconsistent_residual(x, c)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: c(:) !! incompatible equality residuals
        c = [sum(x) - 3.0_dp, 2.0_dp * sum(x) - 7.0_dp, x(1) - x(2) - 1.0_dp]
    end subroutine inconsistent_residual

    subroutine linear_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: a(:, :) !! equality Jacobian
        if (size(x) /= 2) error stop 'test Jacobian must be two-dimensional'
        a(1, :) = [1.0_dp, 1.0_dp]
        a(2, :) = [2.0_dp, 2.0_dp]
        a(3, :) = [1.0_dp, -1.0_dp]
    end subroutine linear_jacobian

    subroutine nonlinear_residual(x, c)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: c(:) !! redundant nonlinear equalities
        c = [sum(x**2) - 5.0_dp, 2.0_dp * sum(x**2) - 10.0_dp, &
            x(1) - x(2) - 1.0_dp]
    end subroutine nonlinear_residual

    subroutine nonlinear_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: a(:, :) !! nonlinear equality Jacobian
        a(1, :) = 2.0_dp * x
        a(2, :) = 4.0_dp * x
        a(3, :) = [1.0_dp, -1.0_dp]
    end subroutine nonlinear_jacobian

    subroutine nonlinear_hessian(x, v, h)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(in) :: v(:) !! equality multipliers
        real(dp), intent(out) :: h(:, :) !! Hessian of multiplier-weighted constraints
        integer :: i
        if (size(x) /= 2) error stop 'test Hessian must be two-dimensional'
        hess_calls = hess_calls + 1
        h = 0.0_dp
        do i = 1, 2
            h(i, i) = 2.0_dp * (v(1) + 2.0_dp * v(2))
        end do
    end subroutine nonlinear_hessian

    subroutine check_solution(label, r, expected, status)
        character(len=*), intent(in) :: label !! scenario name
        type(multivariate_minimize_result), intent(in) :: r !! solver result
        real(dp), intent(in) :: expected(:) !! SciPy 1.17 solution
        integer, intent(in) :: status !! SciPy termination status
        if (.not. r%success .or. r%status /= status) then
            print '(a,a,1x,i0,1x,a)', 'FAILED: ', label, r%status, trim(r%message)
            failures = failures + 1
        else if (maxval(abs(r%x - expected)) > 3.0e-6_dp .or. &
                r%constr_violation > 3.0e-7_dp) then
            print '(a,a,2(1x,es12.4))', 'FAILED solution: ', label, r%x
            failures = failures + 1
        end if
    end subroutine check_solution

    subroutine fail(what)
        character(len=*), intent(in) :: what !! failure description
        print '(a,a)', 'FAILED: ', what
        failures = failures + 1
    end subroutine fail
end program test_optimize_trust_constr_redundant
