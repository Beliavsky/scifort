! SPDX-License-Identifier: MIT
! LinearConstraint-style block regressions against SciPy 1.17.0.
program test_optimize_trust_constr_linear
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : linear_constraint, minimize, &
        multivariate_minimize_result, nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_linear_reference
    implicit none

    type(nonlinear_constraint) :: blocks(2)
    type(nonlinear_constraint) :: invalid(1)
    type(multivariate_minimize_result) :: r
    integer :: failures
    real(dp) :: infinity
    real(dp) :: nan
    real(dp) :: max_weight

    failures = 0
    max_weight = 0.0_dp
    infinity = ieee_value(0.0_dp, ieee_positive_inf)
    nan = ieee_value(0.0_dp, ieee_quiet_nan)

    blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[2.0_dp], upper=[2.0_dp])
    blocks(2) = linear_constraint(reshape([1.0_dp, 0.0_dp], [1, 2]), &
        lower=[0.5_dp], keep_feasible=[.false.])
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=blocks)
    call check_solution('two linear blocks', r, linear_only_x, linear_only_fun, &
        linear_only_nit, linear_only_nfev)
    call check_values('linear constraints', r%constr, linear_only_constr, 1.0e-10_dp)
    call check_counts('linear nfev', r%constr_nfev, linear_only_constr_nfev)
    call check_counts('linear njev', r%constr_njev, linear_only_constr_njev)
    call check_counts('linear nhev', r%constr_nhev, linear_only_constr_nhev)
    call check_values('linear first jac', r%constr_jac(1, :), &
        [1.0_dp, 1.0_dp], 1.0e-14_dp)
    call check_values('linear second jac', r%constr_jac(2, :), &
        [1.0_dp, 0.0_dp], 1.0e-14_dp)

    ! Multiple constant-matrix rows may combine equality and one-sided
    ! inequality components inside a single constraint object.
    blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2]), &
        lower=[2.0_dp, 0.5_dp], upper=[2.0_dp, infinity])
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=blocks(1:1))
    call check_solution('two-row linear block', r, linear_multirow_x, &
        linear_multirow_fun, linear_multirow_nit, linear_multirow_nfev)
    call check_values('two-row linear values', r%constr, linear_multirow_constr, 1.0e-10_dp)
    call check_values('two-row first Jacobian', r%constr_jac(1, :), &
        [1.0_dp, 1.0_dp], 1.0e-14_dp)
    call check_values('two-row second Jacobian', r%constr_jac(2, :), &
        [1.0_dp, 0.0_dp], 1.0e-14_dp)
    call check_counts('two-row linear nfev', r%constr_nfev, linear_multirow_constr_nfev)
    call check_counts('two-row linear njev', r%constr_njev, linear_multirow_constr_njev)
    call check_counts('two-row linear nhev', r%constr_nhev, linear_multirow_constr_nhev)

    blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[2.0_dp], upper=[2.0_dp])

    ! A matrix object also works through the generic minimize dispatcher.
    r = minimize(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=blocks)
    call check_solution('generic linear blocks', r, linear_only_x, linear_only_fun, &
        linear_only_nit, linear_only_nfev)
    call check_counts('generic linear nfev', r%constr_nfev, linear_only_constr_nfev)

    ! A linear equality paired with an analytic nonlinear inequality.
    blocks(2) = nonlinear_constraint()
    blocks(2)%fun => first_value
    blocks(2)%jac => first_jac
    blocks(2)%hess => zero_hess
    blocks(2)%lower = [0.5_dp]
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=blocks)
    call check_solution('linear/analytic nonlinear', r, linear_nonlinear_x, &
        linear_nonlinear_fun, linear_nonlinear_nit, linear_nonlinear_nfev)
    call check_counts('linear/nonlinear nfev', r%constr_nfev, linear_nonlinear_constr_nfev)
    call check_counts('linear/nonlinear njev', r%constr_njev, linear_nonlinear_constr_njev)
    call check_counts('linear/nonlinear nhev', r%constr_nhev, linear_nonlinear_constr_nhev)

    ! A linear block never receives finite-difference function probes.
    nullify(blocks(2)%jac)
    blocks(2)%fd_scheme = '3-point'
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=blocks)
    call check_solution('linear/three-point nonlinear', r, linear_three_point_x, &
        linear_three_point_fun, linear_three_point_nit, linear_three_point_nfev)
    call check_counts('linear/three-point nfev', r%constr_nfev, linear_three_point_constr_nfev)
    call check_counts('linear/three-point njev', r%constr_njev, linear_three_point_constr_njev)
    call check_counts('linear/three-point nhev', r%constr_nhev, linear_three_point_constr_nhev)

    ! Adding variable bounds appends another zero-work block in SciPy order.
    blocks(2)%jac => first_jac
    blocks(2)%fd_scheme = '2-point'
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=blocks, &
        upper=[1.0_dp, 3.0_dp], bounds_keep_feasible=[.true., .true.])
    call check_status('matrix mixed with variable bounds', r, .true.)
    call check_counts('matrix and variable-bound function counts', r%constr_nfev(1:1), [0])
    call check_counts('variable-bound function counts', r%constr_nfev(3:3), [0])
    call check_counts('variable-bound Jacobian counts', r%constr_njev(3:3), [0])
    call check_counts('variable-bound Hessian counts', r%constr_nhev(3:3), [0])
    call check_values('final variable bounds', r%constr(3:4), r%x, 1.0e-14_dp)

    ! A nonlinear constraint with an analytic Hessian can coexist with a
    ! matrix block (which has an exact, identically zero Hessian).
    blocks(1) = nonlinear_constraint()
    blocks(1)%fun => circle_value
    blocks(1)%jac => circle_jac
    blocks(1)%hess => circle_hess
    blocks(1)%lower = [1.0_dp]
    blocks(1)%upper = [1.0_dp]
    blocks(2) = linear_constraint(reshape([0.0_dp, 1.0_dp], [1, 2]), &
        lower=[0.0_dp])
    r = trust_constr(circle_objective, [0.8_dp, 0.6_dp], jac=circle_gradient, &
        hess=quadratic_hess, constraint_objects=blocks, &
        gtol=1.0e-9_dp, xtol=1.0e-10_dp)
    call check_status('weighted nonlinear plus linear', r, .true.)
    call check_values('weighted nonlinear plus linear x', r%x, nonlinear_linear_x, 4.0e-5_dp)
    if (abs(r%fun - nonlinear_linear_fun) > 2.0e-6_dp) call fail('nonlinear/linear objective')
    if (max_weight < 0.5_dp) call fail('nonlinear block Hessian multiplier not passed')
    call check_counts('nonlinear/linear constant counts', r%constr_nfev(2:2), [0])
    call check_counts('nonlinear/linear constant jac counts', r%constr_njev(2:2), [0])
    call check_counts('nonlinear/linear constant hess counts', r%constr_nhev(2:2), [0])
    call check_values('nonlinear/linear final Jacobian', r%constr_jac(2, :), &
        [0.0_dp, 1.0_dp], 1.0e-14_dp)

    ! Input errors are detected before any user callback is invoked.
    invalid(1) = linear_constraint(reshape([1.0_dp], [1, 1]), lower=[0.0_dp])
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=invalid)
    call check_status('invalid matrix column count', r, .false.)

    invalid(1) = linear_constraint(reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
        lower=[0.0_dp])
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=invalid)
    call check_status('invalid matrix row count', r, .false.)

    invalid(1) = linear_constraint(reshape([nan, 1.0_dp], [1, 2]), lower=[0.0_dp])
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=invalid)
    call check_status('nonfinite matrix', r, .false.)

    invalid(1) = linear_constraint(reshape([1.0_dp, 0.0_dp], [1, 2]), lower=[0.0_dp])
    invalid(1)%fun => first_value
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=invalid)
    call check_status('callback on linear block', r, .false.)

    invalid(1) = linear_constraint(reshape([1.0_dp, 0.0_dp], [1, 2]), lower=[0.0_dp])
    invalid(1)%fd_scheme = '3-point'
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=invalid)
    call check_status('finite differences on linear block', r, .false.)

    invalid(1) = linear_constraint(reshape([1.0_dp, 0.0_dp], [1, 2]), &
        lower=[0.8_dp], keep_feasible=[.true.])
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=invalid)
    call check_status('infeasible retained linear inequality', r, .false.)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_trust_constr_linear: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_linear: PASS'

contains

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! objective
        f = (x(1) - 1.0_dp)**2 + (x(2) - 2.5_dp)**2
    end function quadratic

    subroutine quadratic_jac(x, g)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: g(:) !! objective gradient
        g = 2.0_dp * (x - [1.0_dp, 2.5_dp])
    end subroutine quadratic_jac

    subroutine quadratic_hess(x, h)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: h(:, :) !! objective Hessian
        h = 0.0_dp
        if (size(x) < 0) h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 2.0_dp
    end subroutine quadratic_hess

    subroutine first_value(x, c)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: c(:) !! first coordinate
        c(1) = x(1)
    end subroutine first_value

    subroutine first_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: a(:, :) !! first-coordinate Jacobian
        a = 0.0_dp
        a(1, 1) = 1.0_dp
        if (size(x) < 0) a = 0.0_dp
    end subroutine first_jac

    subroutine zero_hess(x, v, h)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(in) :: v(:) !! multipliers
        real(dp), intent(out) :: h(:, :) !! zero constraint Hessian
        h = 0.0_dp
        if (size(x) < 0 .or. size(v) < 0) h = 0.0_dp
    end subroutine zero_hess

    function circle_objective(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! objective
        f = (x(1) - 2.0_dp)**2 + (x(2) - 1.0_dp)**2
    end function circle_objective

    subroutine circle_gradient(x, g)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: g(:) !! objective gradient
        g = 2.0_dp * (x - [2.0_dp, 1.0_dp])
    end subroutine circle_gradient

    subroutine circle_value(x, c)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: c(:) !! squared radius
        c(1) = sum(x*x)
    end subroutine circle_value

    subroutine circle_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: a(:, :) !! circle Jacobian
        a(1, :) = 2.0_dp * x
    end subroutine circle_jac

    subroutine circle_hess(x, v, h)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(in) :: v(:) !! multiplier
        real(dp), intent(out) :: h(:, :) !! weighted circle Hessian
        max_weight = max(max_weight, abs(v(1)))
        h = 0.0_dp
        if (size(x) < 0) h = 0.0_dp
        h(1, 1) = 2.0_dp * v(1)
        h(2, 2) = 2.0_dp * v(1)
    end subroutine circle_hess

    subroutine check_solution(label, actual, expected_x, expected_fun, expected_nit, expected_nfev)
        character(len=*), intent(in) :: label !! test name
        type(multivariate_minimize_result), intent(in) :: actual !! computed result
        real(dp), intent(in) :: expected_x(:) !! reference solution
        real(dp), intent(in) :: expected_fun !! reference objective
        integer, intent(in) :: expected_nit !! reference iteration count
        integer, intent(in) :: expected_nfev !! reference objective evaluations

        call check_status(label, actual, .true.)
        call check_values(trim(label)//' x', actual%x, expected_x, 1.0e-8_dp)
        if (abs(actual%fun - expected_fun) > 1.0e-8_dp) call fail(trim(label)//' objective')
        if (actual%nit /= expected_nit) call fail(trim(label)//' iterations')
        if (actual%nfev /= expected_nfev) call fail(trim(label)//' objective evaluations')
    end subroutine check_solution

    subroutine check_status(label, actual, expected_success)
        character(len=*), intent(in) :: label !! test name
        type(multivariate_minimize_result), intent(in) :: actual !! computed result
        logical, intent(in) :: expected_success !! whether success is required

        if (actual%success .neqv. expected_success) then
            print '(4a,i0)', trim(label), ': ', trim(actual%message), ' (status) ', actual%status
            call fail(trim(label)//' success flag')
        end if
        if (.not. expected_success .and. actual%status /= -1) then
            call fail(trim(label)//' invalid-input status')
        end if
    end subroutine check_status

    subroutine check_values(label, actual, expected, tol)
        character(len=*), intent(in) :: label !! test name
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tol !! absolute tolerance

        if (size(actual) /= size(expected)) then
            call fail(trim(label)//' size')
        else if (any(abs(actual - expected) > tol)) then
            write(*, '(a,*(es24.15,1x))') trim(label)//' actual: ', actual
            call fail(trim(label)//' values')
        end if
    end subroutine check_values

    subroutine check_counts(label, actual, expected)
        character(len=*), intent(in) :: label !! test name
        integer, intent(in) :: actual(:) !! computed counts
        integer, intent(in) :: expected(:) !! reference counts

        if (size(actual) /= size(expected)) then
            call fail(trim(label)//' size')
        else if (any(actual /= expected)) then
            write(*, '(a,*(i0,1x))') trim(label)//' actual: ', actual
            call fail(trim(label)//' counts')
        end if
    end subroutine check_counts

    subroutine fail(message)
        character(len=*), intent(in) :: message !! reason for test failure
        print '(2a)', 'FAIL: ', trim(message)
        failures = failures + 1
    end subroutine fail
end program test_optimize_trust_constr_linear
