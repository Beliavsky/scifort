! SPDX-License-Identifier: MIT
! Numerical objective Hessian from analytic gradients, including bounded probes.
program test_optimize_trust_constr_hess_fd
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, multivariate_minimize_result, &
        nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_hess_fd_reference
    implicit none

    type(multivariate_minimize_result) :: r
    type(nonlinear_constraint) :: blocks(1)
    integer :: failures
    integer :: gradient_calls
    integer :: bound_violations
    integer :: callback_calls
    logical :: record_probes
    logical :: backward_probe_seen
    logical :: forward_probe_seen
    real(dp) :: reference_probe_x
    real(dp) :: infinity

    failures = 0
    gradient_calls = 0
    bound_violations = 0
    callback_calls = 0
    record_probes = .false.
    backward_probe_seen = .false.
    forward_probe_seen = .false.
    reference_probe_x = 0.0_dp
    infinity = ieee_value(0.0_dp, ieee_positive_inf)

    gradient_calls = 0
    r = trust_constr(rosen, [-1.2_dp, 1.0_dp], jac=rosen_grad, &
        objective_hess_fd_scheme='2-point', maxiter=200)
    call check_solution('rosen two-point', r, fd2_rosen, 5.0e-5_dp)
    call check_counts('rosen two-point', r)
    if (r%nit <= 0 .or. r%nit > 200) call fail('rosen two-point iteration bound')

    gradient_calls = 0
    r = minimize(rosen, [-1.2_dp, 1.0_dp], method='trust-constr', jac=rosen_grad, &
        objective_hess_fd_scheme='3-point', maxiter=200)
    call check_solution('rosen three-point through minimize', r, fd3_rosen, 5.0e-5_dp)
    call check_counts('rosen three-point', r)
    if (r%nit <= 0 .or. r%nit > 200) call fail('rosen three-point iteration bound')

    ! Exact weighted constraint Hessian plus objective Hessian approximated from
    ! the analytic objective gradient. This also exercises constraint objects.
    blocks(1)%fun => circle
    blocks(1)%jac => circle_jac
    blocks(1)%hess => circle_hess
    blocks(1)%lower = [1.0_dp]
    blocks(1)%upper = [1.0_dp]
    gradient_calls = 0
    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_fd_scheme='2-point', constraint_objects=blocks, maxiter=150)
    call check_solution('circle two-point', r, fd2_circle, 1.0e-5_dp)
    call check_counts('circle two-point', r)
    if (r%constr_violation > 1.0e-6_dp) call fail('circle two-point feasibility')
    if (.not. allocated(r%constr_nhev)) then
        call fail('missing constraint Hessian accounting')
    else if (r%constr_nhev(1) <= 0) then
        call fail('exact constraint Hessian never called')
    end if

    gradient_calls = 0
    r = minimize(quadratic, [0.8_dp, 0.6_dp], method='trust-constr', jac=quadratic_grad, &
        objective_hess_fd_scheme='3-point', constraint_objects=blocks, maxiter=150)
    call check_solution('circle three-point', r, fd3_circle, 1.0e-5_dp)
    call check_counts('circle three-point', r)
    if (r%constr_violation > 1.0e-6_dp) call fail('circle three-point feasibility')

    ! Inequality barrier path, including box constraints on the physical vars.
    gradient_calls = 0
    r = trust_constr(quadratic, [0.5_dp, 0.3_dp], jac=quadratic_grad, &
        lower=[0.0_dp, 0.0_dp], upper=[2.0_dp, 3.0_dp], &
        objective_hess_fd_scheme='2-point', maxiter=150)
    call check_solution('bounded two-point', r, fd2_bounds, 1.0e-3_dp)
    call check_counts('bounded two-point', r)
    if (trim(r%method) /= 'tr_interior_point') call fail('bounded barrier method')

    gradient_calls = 0
    r = trust_constr(quadratic, [0.5_dp, 0.3_dp], jac=quadratic_grad, &
        lower=[0.0_dp, 0.0_dp], upper=[2.0_dp, 3.0_dp], &
        objective_hess_fd_scheme='3-point', maxiter=150)
    call check_solution('bounded three-point', r, fd3_bounds, 1.0e-3_dp)
    call check_counts('bounded three-point', r)

    ! Two-sided finite differences switch to inward, one-sided probes at
    ! retained endpoints. An exactly fixed variable has zero feasible step.
    bound_violations = 0
    gradient_calls = 0
    r = trust_constr(quadratic, [0.1_dp, 1.0_dp], jac=domain_gradient, &
        lower=[0.1_dp, -infinity], upper=[infinity, 1.0_dp], &
        bounds_keep_feasible=[.true., .true.], &
        objective_hess_fd_scheme='3-point', maxiter=1)
    if (r%nit /= 1) call fail('one-sided init iteration')
    call check_counts('retained endpoints three-point', r)
    if (bound_violations /= 0) call fail('three-point Hessian stepped outside retained bounds')

    bound_violations = 0
    gradient_calls = 0
    r = trust_constr(quadratic, [0.1_dp, 0.0_dp], jac=domain_gradient, &
        lower=[0.1_dp, -infinity], upper=[0.1_dp, infinity], &
        objective_hess_fd_scheme='3-point', maxiter=150)
    call check_solution('fixed variable', r, [0.1_dp, -2.0_dp], 2.0e-4_dp)
    call check_counts('fixed variable', r)
    if (bound_violations /= 0) call fail('fixed variable probed outside interval')

    gradient_calls = 0
    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_fd_scheme='3-point', &
        objective_finite_diff_rel_step=[1.0e-4_dp, 2.0e-4_dp], maxiter=100)
    call check_solution('custom relative step', r, [1.0_dp, -2.0_dp], 1.0e-7_dp)
    call check_counts('custom relative step', r)

    ! SciPy's signed two-point convention follows the sign of the point,
    ! including explicit negative relative steps, unless bounds force reversal.
    call watch_probes(-1.0_dp)
    gradient_calls = 0
    r = trust_constr(quadratic, [-1.0_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_fd_scheme='2-point', maxiter=1)
    if (.not. backward_probe_seen .or. forward_probe_seen) &
        call fail('default two-point Hessian perturbation sign')
    call check_counts('default signed step', r)

    call watch_probes(1.0_dp)
    gradient_calls = 0
    r = trust_constr(quadratic, [1.0_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_fd_scheme='2-point', &
        objective_finite_diff_rel_step=[-1.0e-4_dp, 1.0e-4_dp], maxiter=1)
    if (.not. backward_probe_seen .or. forward_probe_seen) &
        call fail('explicit signed two-point Hessian perturbation')
    call check_counts('explicit signed step', r)

    call watch_probes(0.1_dp)
    gradient_calls = 0
    r = trust_constr(quadratic, [0.1_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_fd_scheme='2-point', &
        lower=[0.1_dp, -infinity], bounds_keep_feasible=[.true., .false.], &
        objective_finite_diff_rel_step=[-1.0e-4_dp, 1.0e-4_dp], maxiter=1)
    if (.not. forward_probe_seen .or. backward_probe_seen) &
        call fail('bound-reversed two-point Hessian perturbation')
    call check_counts('bound-reversed signed step', r)
    record_probes = .false.

    ! Stopping callbacks survive recursive Hessian construction.
    gradient_calls = 0
    callback_calls = 0
    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_fd_scheme='2-point', callback=stop_callback)
    if (r%status /= 3 .or. r%success .or. r%nit /= 2 .or. callback_calls /= 2) &
        call fail('callback forwarding through numerical Hessian')
    call check_counts('callback', r)

    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_fd_scheme='cs')
    call check_invalid('invalid scheme name', r)
    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], &
        objective_hess_fd_scheme='2-point')
    call check_invalid('missing analytic gradient', r)
    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        hess=quadratic_hess, objective_hess_fd_scheme='2-point')
    call check_invalid('conflicting exact Hessian', r)
    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        hessp=quadratic_hessp, objective_hess_fd_scheme='2-point')
    call check_invalid('conflicting Hessian product', r)
    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_update='sr1', objective_hess_fd_scheme='2-point')
    call check_invalid('conflicting Hessian update', r)
    r = minimize(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        method='bfgs', objective_hess_fd_scheme='2-point')
    call check_invalid('wrong minimizer', r)
    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_fd_scheme='2-point', objective_finite_diff_rel_step=[1.0e-4_dp])
    call check_invalid('invalid relative step shape', r)
    r = trust_constr(quadratic, [0.8_dp, 0.6_dp], jac=quadratic_grad, &
        objective_hess_fd_scheme='2-point', lower=[1.0_dp, 1.0_dp], &
        upper=[0.0_dp, 2.0_dp])
    call check_invalid('invalid bound interval', r)

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_hess_fd: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_hess_fd: PASS'

contains
    subroutine watch_probes(base_x)
        real(dp), intent(in) :: base_x !! expected initial first coordinate
        reference_probe_x = base_x
        backward_probe_seen = .false.
        forward_probe_seen = .false.
        record_probes = .true.
    end subroutine watch_probes

    subroutine fail(reason)
        character(len=*), intent(in) :: reason !! diagnostic
        failures = failures + 1
        print '(a)', '  FAIL: '//reason
    end subroutine fail

    subroutine check_solution(label, result, reference, tolerance)
        character(len=*), intent(in) :: label !! test name
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        real(dp), intent(in) :: reference(:) !! independent SciPy solution
        real(dp), intent(in) :: tolerance !! maximum error
        if (.not. result%success) then
            print '(a,a,i0,1x,a)', '  Status: ', label, result%status, trim(result%message)
            call fail(label//' did not converge')
        end if
        if (.not. allocated(result%x)) then
            call fail(label//' missing minimizer')
        else if (size(result%x) /= size(reference)) then
            call fail(label//' minimizer size')
        else if (maxval(abs(result%x - reference)) > tolerance) then
            print '(a,2es18.8)', '  Actual: ', result%x
            call fail(label//' incorrect minimizer')
        end if
    end subroutine check_solution

    subroutine check_counts(label, result)
        character(len=*), intent(in) :: label !! test name
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        if (gradient_calls /= result%njev) then
            print '(a,1x,a,2i7)', '  Counts:', label, gradient_calls, result%njev
            call fail(label//' gradient accounting')
        end if
        if (result%nhev /= 0) call fail(label//' numerical Hessian counted as exact Hessian')
    end subroutine check_counts

    subroutine check_invalid(label, result)
        character(len=*), intent(in) :: label !! test name
        type(multivariate_minimize_result), intent(in) :: result !! invalid result
        if (result%status /= -1 .or. result%success) call fail(label//' not rejected')
    end subroutine check_invalid

    function rosen(x) result(f)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: f !! Rosenbrock objective
        f = (1.0_dp-x(1))**2 + 100.0_dp*(x(2)-x(1)**2)**2
    end function rosen

    subroutine rosen_grad(x, g)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: g(:) !! analytic gradient
        g(1) = 2.0_dp*(x(1)-1.0_dp) - 400.0_dp*x(1)*(x(2)-x(1)**2)
        g(2) = 200.0_dp*(x(2)-x(1)**2)
        gradient_calls = gradient_calls + 1
    end subroutine rosen_grad

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: f !! strictly convex objective
        f = (x(1)-1.0_dp)**2 + 3.0_dp*(x(2)+2.0_dp)**2
    end function quadratic

    subroutine quadratic_grad(x, g)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: g(:) !! analytic gradient
        g = [2.0_dp*(x(1)-1.0_dp), 6.0_dp*(x(2)+2.0_dp)]
        gradient_calls = gradient_calls + 1
        if (record_probes) then
            if (x(1) < reference_probe_x) backward_probe_seen = .true.
            if (x(1) > reference_probe_x) forward_probe_seen = .true.
        end if
    end subroutine quadratic_grad

    subroutine domain_gradient(x, g)
        real(dp), intent(in) :: x(:) !! vector subject to explicit box bounds
        real(dp), intent(out) :: g(:) !! analytic gradient
        if (x(1) < 0.1_dp .or. x(2) > 1.0_dp) bound_violations = bound_violations + 1
        call quadratic_grad(x, g)
    end subroutine domain_gradient

    subroutine quadratic_hess(x, h)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: h(:, :) !! exact objective Hessian
        h = 0.0_dp
        h(1,1) = 2.0_dp
        h(2,2) = 6.0_dp
    end subroutine quadratic_hess

    subroutine quadratic_hessp(x, p, hp)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: p(:) !! trial direction
        real(dp), intent(out) :: hp(:) !! Hessian product
        hp = [2.0_dp*p(1), 6.0_dp*p(2)]
    end subroutine quadratic_hessp

    subroutine circle(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! scalar nonlinear constraint
        c(1) = sum(x**2)
    end subroutine circle

    subroutine circle_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: a(:, :) !! constraint Jacobian
        a(1, :) = 2.0_dp*x
    end subroutine circle_jac

    subroutine circle_hess(x, v, h)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! weighted-constraint multipliers
        real(dp), intent(out) :: h(:, :) !! exact weighted Hessian
        h = 0.0_dp
        h(1, 1) = 2.0_dp*v(1)
        h(2, 2) = 2.0_dp*v(1)
    end subroutine circle_hess

    subroutine stop_callback(x, nit, stop)
        real(dp), intent(in) :: x(:) !! current point
        integer, intent(in) :: nit !! total iteration index
        logical, intent(out) :: stop !! user termination
        callback_calls = callback_calls + 1
        stop = nit >= 2
    end subroutine stop_callback
end program test_optimize_trust_constr_hess_fd
