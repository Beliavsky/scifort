! SPDX-License-Identifier: MIT
! Regression tests for per-coordinate objective relative steps in trust-constr.
program test_optimize_trust_constr_objective_rel_step
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : linear_constraint, minimize, &
        multivariate_minimize_result, nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_objective_rel_step_reference
    implicit none

    type(multivariate_minimize_result) :: r
    type(nonlinear_constraint) :: objects(1)
    real(dp) :: infinity
    real(dp) :: quiet_nan
    integer :: failures
    integer :: out_of_bounds

    failures = 0
    out_of_bounds = 0
    infinity = ieee_value(0.0_dp, ieee_positive_inf)
    quiet_nan = ieee_value(0.0_dp, ieee_quiet_nan)

    ! The per-coordinate steps follow SciPy's abs(x) scaling when explicitly
    ! supplied; max(1,abs(x)) is reserved for the built-in default/fallback.
    r = trust_constr(cubic, [2.0_dp, -3.0_dp], maxiter=1, &
        objective_fd_scheme='3-point', objective_finite_diff_rel_step=[0.1_dp, 0.2_dp])
    call check_gradient('direct three point', r, three_point, 1.0e-8_dp, 5)

    r = minimize(cubic, [2.0_dp, -3.0_dp], method='trust-constr', maxiter=1, &
        objective_fd_scheme='3-point', objective_finite_diff_rel_step=[0.1_dp, 0.2_dp])
    call check_gradient('generic three point', r, three_point, 1.0e-8_dp, 5)

    r = trust_constr(cubic, [2.0_dp, -3.0_dp], maxiter=1, &
        objective_fd_scheme='2-point', objective_finite_diff_rel_step=[0.1_dp, 0.2_dp])
    call check_gradient('direct two point', r, two_point, 1.0e-8_dp, 3)

    ! With two points, the sign follows the coordinate as in SciPy. A negative
    ! relative step intentionally reverses the signed probe direction.
    r = trust_constr(cubic, [2.0_dp, -3.0_dp], maxiter=1, &
        objective_finite_diff_rel_step=[-0.1_dp, 0.2_dp])
    call check_gradient('signed two point', r, signed_two_point, 1.0e-8_dp, 3)

    objects(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[-1.0_dp], upper=[-1.0_dp])
    r = trust_constr(cubic, [2.0_dp, -3.0_dp], maxiter=1, &
        constraint_objects=objects, objective_fd_scheme='3-point', &
        objective_finite_diff_rel_step=[0.1_dp, 0.2_dp])
    call check_gradient('constraint object', r, three_point, 1.0e-8_dp, 5)
    if (.not. allocated(r%constr_nfev)) then
        call fail('constraint block diagnostics absent')
    else if (size(r%constr_nfev) /= 1) then
        call fail('constraint block diagnostics shape')
    else if (r%constr_nfev(1) /= 0) then
        call fail('linear constraint callback was evaluated')
    end if

    ! The protected lower endpoint must use a one-sided stencil, with the
    ! requested 0.2 first-coordinate absolute step and no infeasible calls.
    out_of_bounds = 0
    r = trust_constr(restricted_cubic, [2.0_dp, -3.0_dp], maxiter=1, &
        lower=[2.0_dp, -infinity], bounds_keep_feasible=[.true., .false.], &
        objective_fd_scheme='3-point', objective_finite_diff_rel_step=[0.1_dp, 0.2_dp])
    call check_gradient('retained bound', r, lower_bound, 1.0e-7_dp, 7)
    if (out_of_bounds /= 0) call fail('objective evaluated outside retained lower bound')

    ! Zero coordinates must not lead to zero-sized or unrepresentable steps.
    r = trust_constr(cubic, [0.0_dp, -3.0_dp], maxiter=1, &
        objective_fd_scheme='3-point', objective_finite_diff_rel_step=[0.1_dp, 0.2_dp])
    call check_gradient('zero coordinate fallback', r, zero_coordinate, 1.0e-6_dp, 5)

    ! Supplying an exact Jacobian always bypasses all numerical probes.
    r = trust_constr(cubic, [2.0_dp, -3.0_dp], jac=cubic_gradient, maxiter=1, &
        objective_fd_scheme='3-point', objective_finite_diff_rel_step=[0.1_dp, 0.2_dp])
    call check_gradient('analytic override', r, [12.0_dp, 54.0_dp], 1.0e-12_dp, 1)

    ! An unconstrained solve must also carry the user vector through all
    ! subsequent SQP gradient evaluations, not just the initial derivative.
    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], hess=quadratic_hessian, &
        objective_fd_scheme='3-point', &
        objective_finite_diff_rel_step=[1.0e-5_dp, 2.0e-5_dp])
    if (.not. r%success .or. maxval(abs(r%x - [1.0_dp, -2.0_dp])) > 1.0e-5_dp) then
        call fail('full unconstrained solve with custom relative steps')
    end if
    if (r%njev <= 1 .or. r%nfev <= 5) then
        call fail('custom relative steps on subsequent iterations')
    end if

    r = trust_constr(cubic, [2.0_dp, -3.0_dp], &
        objective_finite_diff_rel_step=[0.1_dp])
    call check_invalid('relative step size', r)

    r = trust_constr(cubic, [2.0_dp, -3.0_dp], &
        objective_finite_diff_rel_step=[quiet_nan, 0.2_dp])
    call check_invalid('nonfinite relative step', r)

    r = trust_constr(cubic, [2.0_dp, -3.0_dp], &
        objective_finite_diff_rel_step=[infinity, 0.2_dp], &
        constraint_objects=objects)
    call check_invalid('nonfinite relative step through blocks', r)

    r = minimize(cubic, [2.0_dp, -3.0_dp], method='bfgs', &
        objective_finite_diff_rel_step=[0.1_dp, 0.2_dp])
    call check_invalid('trust-only keyword on bfgs', r)

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_trust_constr_objective_rel_step: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_objective_rel_step: PASS'

contains

    function cubic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! nonlinear objective
        f = x(1)**3 + 2.0_dp * x(2)**3
    end function cubic

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! strictly convex quadratic
        f = (x(1) - 1.0_dp)**2 + 3.0_dp * (x(2) + 2.0_dp)**2
    end function quadratic

    subroutine quadratic_hessian(x, h)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: h(:, :) !! exact objective Hessian
        if (size(x) /= 2) error stop 'quadratic Hessian shape'
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 6.0_dp
    end subroutine quadratic_hessian

    function restricted_cubic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! strictly bounded objective
        if (x(1) < 2.0_dp) out_of_bounds = out_of_bounds + 1
        f = cubic(x)
    end function restricted_cubic

    subroutine cubic_gradient(x, g)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: g(:) !! analytic gradient
        g = [3.0_dp * x(1)**2, 6.0_dp * x(2)**2]
    end subroutine cubic_gradient

    subroutine check_gradient(name, result, expected, tolerance, expected_nfev)
        character(len=*), intent(in) :: name !! test name
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        real(dp), intent(in) :: expected(:) !! SciPy numerical gradient
        real(dp), intent(in) :: tolerance !! allowable maximum difference
        integer, intent(in) :: expected_nfev !! initial objective plus FD probes
        if (.not. allocated(result%jac)) then
            call fail(name // ': missing jac')
            return
        end if
        print '(a,a,2f15.8,a,3i6)', trim(name), ' jac=', result%jac, &
            ' nit/nfev/njev=', result%nit, result%nfev, result%njev
        if (maxval(abs(result%jac - expected)) > tolerance) then
            call fail(name // ': gradient disagrees with SciPy')
        end if
        if (result%nit /= 1 .or. result%nfev /= expected_nfev .or. result%njev /= 1) then
            call fail(name // ': evaluation accounting')
        end if
        if (result%status /= 0 .or. result%success) call fail(name // ': maxiter status')
    end subroutine check_gradient

    subroutine check_invalid(name, result)
        character(len=*), intent(in) :: name !! invalid input case
        type(multivariate_minimize_result), intent(in) :: result !! returned failure
        if (result%success .or. result%status /= -1) call fail(name)
    end subroutine check_invalid

    subroutine fail(message)
        character(len=*), intent(in) :: message !! failure description
        print '(a)', 'FAIL: ' // message
        failures = failures + 1
    end subroutine fail

end program test_optimize_trust_constr_objective_rel_step
