! SPDX-License-Identifier: MIT
! Signed relative steps for numerical nonlinear constraint Jacobians.
program test_optimize_trust_constr_constraint_rel_step
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, multivariate_minimize_result, &
        nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_constraint_rel_step_reference
    implicit none

    type(nonlinear_constraint) :: objects(1)
    type(multivariate_minimize_result) :: r
    integer :: failures
    integer :: observed_calls
    integer :: bound_violations
    real(dp) :: probes(2, 20)
    real(dp) :: infinity
    real(dp) :: quiet_nan
    logical :: monitor_bounds

    failures = 0
    observed_calls = 0
    bound_violations = 0
    monitor_bounds = .false.
    infinity = ieee_value(0.0_dp, ieee_positive_inf)
    quiet_nan = ieee_value(0.0_dp, ieee_quiet_nan)

    objects(1)%fun => constraint_values
    objects(1)%lower = [0.0_dp, 0.0_dp]
    objects(1)%upper = objects(1)%lower

    call check_case('dense signed two point', [2.0_dp, -3.0_dp], &
        [-0.1_dp, 0.2_dp], '2-point', .false., dense_two_jac, dense_two_nfev)
    if (observed_calls >= 3) then
        if (maxval(abs(probes(:, 2) - [1.8_dp, -3.0_dp])) > 1.0e-12_dp) &
            call fail('dense first signed probe')
        if (maxval(abs(probes(:, 3) - [2.0_dp, -3.6_dp])) > 1.0e-12_dp) &
            call fail('dense second signed probe')
    end if

    call check_case('grouped signed two point', [2.0_dp, -3.0_dp], &
        [-0.1_dp, 0.2_dp], '2-point', .true., sparse_two_jac, sparse_two_nfev)
    if (observed_calls >= 2) then
        if (maxval(abs(probes(:, 2) - [1.8_dp, -3.6_dp])) > 1.0e-12_dp) &
            call fail('grouped columns did not share signed probe')
    end if

    call check_case('dense signed three point', [2.0_dp, -3.0_dp], &
        [-0.1_dp, 0.2_dp], '3-point', .false., dense_three_jac, dense_three_nfev)
    call check_case('grouped signed three point', [2.0_dp, -3.0_dp], &
        [-0.1_dp, 0.2_dp], '3-point', .true., sparse_three_jac, sparse_three_nfev)

    ! A two-point step requested outside a retained bound must reverse its
    ! direction, while the other coordinate keeps its negative displacement.
    call check_case('bounded signed two point', [2.0_dp, -3.0_dp], &
        [-0.1_dp, 0.2_dp], '2-point', .false., bounded_two_jac, bounded_two_nfev, &
        lower_bound=[2.0_dp, -infinity], upper_bound=[infinity, -3.0_dp])
    call check_case('bounded grouped two point', [2.0_dp, -3.0_dp], &
        [-0.1_dp, 0.2_dp], '2-point', .true., bounded_sparse_two_jac, &
        bounded_sparse_two_nfev, lower_bound=[2.0_dp, -infinity], &
        upper_bound=[infinity, -3.0_dp])
    call check_case('bounded grouped three point', [2.0_dp, -3.0_dp], &
        [-0.1_dp, 0.2_dp], '3-point', .true., bounded_sparse_three_jac, &
        bounded_sparse_three_nfev, lower_bound=[2.0_dp, -infinity], &
        upper_bound=[infinity, -3.0_dp])

    ! SciPy replaces a requested zero or a sub-ULP step by an automatic step.
    call check_case('zero relative step', [0.0_dp, -3.0_dp], &
        [0.0_dp, 0.2_dp], '2-point', .false., zero_step_two_jac, zero_step_two_nfev)
    if (observed_calls >= 2) then
        if (probes(1, 2) <= 0.0_dp) call fail('zero relative step used a zero probe')
    end if
    call check_case('small coordinates', [0.5_dp, -0.5_dp], &
        [0.2_dp, 0.2_dp], '2-point', .false., small_two_jac, small_two_nfev)
    call check_case('reverse signs', [2.0_dp, -3.0_dp], &
        [0.1_dp, -0.2_dp], '2-point', .false., reverse_two_jac, reverse_two_nfev)
    call check_case('sub-ULP fallback', [2.0_dp, -3.0_dp], &
        [1.0e-300_dp, 1.0e-300_dp], '2-point', .true., [4.0_dp, -6.0_dp], 2)

    ! The generic minimize front end must preserve the signed step vector.
    objects(1)%fd_scheme = '2-point'
    objects(1)%finite_diff_rel_step = [-0.1_dp, 0.2_dp]
    if (allocated(objects(1)%finite_diff_jac_sparsity)) &
        deallocate(objects(1)%finite_diff_jac_sparsity)
    r = minimize(objective, [2.0_dp, -3.0_dp], method='trust-constr', &
        jac=gradient, hess=hessian, constraint_objects=objects, maxiter=1)
    call check_jacobian('minimize forwarding', r, dense_two_jac)

    ! Analytic Jacobians override the finite-difference settings.
    objects(1)%jac => exact_constraint_jacobian
    observed_calls = 0
    r = trust_constr(objective, [2.0_dp, -3.0_dp], jac=gradient, &
        hess=hessian, constraint_objects=objects, maxiter=1)
    call check_jacobian('analytic Jacobian precedence', r, [4.0_dp, -6.0_dp])
    if (allocated(r%constr_njev)) then
        if (r%constr_njev(1) /= 1) call fail('analytic Jacobian call count')
    else
        call fail('analytic Jacobian count not allocated')
    end if
    nullify(objects(1)%jac)

    ! Invalid relative steps must be rejected before the constraint is called.
    objects(1)%finite_diff_rel_step = [0.1_dp]
    observed_calls = 0
    r = trust_constr(objective, [2.0_dp, -3.0_dp], jac=gradient, &
        hess=hessian, constraint_objects=objects, maxiter=1)
    call check_invalid('wrong-sized relative step', r)
    objects(1)%finite_diff_rel_step = [quiet_nan, 0.1_dp]
    observed_calls = 0
    r = trust_constr(objective, [2.0_dp, -3.0_dp], jac=gradient, &
        hess=hessian, constraint_objects=objects, maxiter=1)
    call check_invalid('NaN relative step', r)
    objects(1)%finite_diff_rel_step = [infinity, 0.1_dp]
    observed_calls = 0
    r = trust_constr(objective, [2.0_dp, -3.0_dp], jac=gradient, &
        hess=hessian, constraint_objects=objects, maxiter=1)
    call check_invalid('infinite relative step', r)

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_constraint_rel_step: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_constraint_rel_step: PASS'

contains

    subroutine check_case(label, x0, rel, scheme, colored, expected, calls, &
            lower_bound, upper_bound)
        character(len=*), intent(in) :: label !! test case description
        real(dp), intent(in) :: x0(:) !! initial physical variables
        real(dp), intent(in) :: rel(:) !! custom relative differences
        character(len=*), intent(in) :: scheme !! finite difference scheme
        logical, intent(in) :: colored !! use disjoint-column grouping
        real(dp), intent(in) :: expected(:) !! SciPy reference Jacobian diagonal
        integer, intent(in) :: calls !! SciPy constraint evaluation count
        real(dp), intent(in), optional :: lower_bound(:) !! optional finite lower endpoints
        real(dp), intent(in), optional :: upper_bound(:) !! optional finite upper endpoints

        objects(1)%fd_scheme = scheme
        objects(1)%finite_diff_rel_step = rel
        if (allocated(objects(1)%finite_diff_jac_sparsity)) &
            deallocate(objects(1)%finite_diff_jac_sparsity)
        if (colored) then
            allocate(objects(1)%finite_diff_jac_sparsity(2, 2))
            objects(1)%finite_diff_jac_sparsity = .false.
            objects(1)%finite_diff_jac_sparsity(1, 1) = .true.
            objects(1)%finite_diff_jac_sparsity(2, 2) = .true.
        end if
        observed_calls = 0
        bound_violations = 0
        monitor_bounds = present(lower_bound)
        if (present(lower_bound) .and. present(upper_bound)) then
            r = trust_constr(objective, x0, jac=gradient, hess=hessian, &
                constraint_objects=objects, lower=lower_bound, upper=upper_bound, &
                bounds_keep_feasible=[.true., .true.], maxiter=1)
        else
            r = trust_constr(objective, x0, jac=gradient, hess=hessian, &
                constraint_objects=objects, maxiter=1)
        end if
        monitor_bounds = .false.
        call check_jacobian(label, r, expected)
        if (observed_calls /= calls) then
            print '(a,1x,a,2(1x,i0))', '  count mismatch', trim(label), observed_calls, calls
            call fail('constraint evaluation count differs from SciPy')
        end if
        if (allocated(r%constr_nfev)) then
            if (r%constr_nfev(1) /= observed_calls) &
                call fail('reported constraint count differs from callback count')
        else
            call fail('constraint count unavailable')
        end if
        if (bound_violations /= 0) call fail('constraint probe crossed retained bounds')
    end subroutine check_case

    subroutine check_jacobian(label, result, expected)
        character(len=*), intent(in) :: label !! test case description
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        real(dp), intent(in) :: expected(:) !! SciPy Jacobian diagonal
        real(dp) :: tolerance

        tolerance = 5.0e-7_dp
        if (result%status == -1) then
            print '(a,1x,a,1x,a)', '  invalid result', trim(label), trim(result%message)
            call fail('solver rejected valid relative steps')
        end if
        if (.not. allocated(result%constr_jac)) then
            call fail('constraint Jacobian unavailable')
            return
        end if
        ! With variable bounds the result appends an identity bounds block
        ! after the two user-constraint rows.
        if (size(result%constr_jac, 1) < 2 .or. size(result%constr_jac, 2) /= 2) then
            print '(a,1x,a)', '  wrong shape:', trim(label)
            call fail('incorrect constraint Jacobian shape')
            return
        end if
        if (abs(result%constr_jac(1, 1) - expected(1)) > tolerance .or. &
                abs(result%constr_jac(2, 2) - expected(2)) > tolerance .or. &
                abs(result%constr_jac(1, 2)) > tolerance .or. &
                abs(result%constr_jac(2, 1)) > tolerance) then
            print '(a,1x,a)', '  Jacobian mismatch:', trim(label)
            print '(a,2(1x,es16.7))', '  observed:', result%constr_jac(1, 1), result%constr_jac(2, 2)
            print '(a,2(1x,es16.7))', '  expected:', expected
            call fail('finite-difference Jacobian disagrees with SciPy')
        end if
    end subroutine check_jacobian

    subroutine check_invalid(label, result)
        character(len=*), intent(in) :: label !! test case description
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        if (result%status /= -1 .or. result%success .or. observed_calls /= 0) then
            print '(a,1x,a)', '  not rejected:', trim(label)
            call fail('invalid relative step was accepted')
        end if
    end subroutine check_invalid

    subroutine fail(message)
        character(len=*), intent(in) :: message !! reason for test failure
        print '(a)', 'FAIL: '//message
        failures = failures + 1
    end subroutine fail

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: f !! objective value
        f = sum(x**2)
    end function objective

    subroutine gradient(x, g)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: g(:) !! objective gradient
        g = 2.0_dp * x
    end subroutine gradient

    subroutine hessian(x, h)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: h(:, :) !! objective Hessian
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 2.0_dp
    end subroutine hessian

    subroutine constraint_values(x, c)
        real(dp), intent(in) :: x(:) !! physical variables
        real(dp), intent(out) :: c(:) !! independent quadratic constraints
        observed_calls = observed_calls + 1
        if (observed_calls <= size(probes, 2)) probes(:, observed_calls) = x
        if (monitor_bounds) then
            if (x(1) < 2.0_dp - 1.0e-12_dp .or. x(2) > -3.0_dp + 1.0e-12_dp) &
                bound_violations = bound_violations + 1
        end if
        c = x**2
    end subroutine constraint_values

    subroutine exact_constraint_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! physical variables
        real(dp), intent(out) :: a(:, :) !! exact Jacobian
        a = 0.0_dp
        a(1, 1) = 2.0_dp * x(1)
        a(2, 2) = 2.0_dp * x(2)
    end subroutine exact_constraint_jacobian

end program test_optimize_trust_constr_constraint_rel_step
