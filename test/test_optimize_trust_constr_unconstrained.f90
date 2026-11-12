! SPDX-License-Identifier: MIT
! Unconstrained trust-constr regressions against SciPy 1.17.0.
program test_optimize_trust_constr_unconstrained
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, multivariate_minimize_result, &
        nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_unconstrained_reference
    implicit none

    type(nonlinear_constraint) :: empty(0)
    type(multivariate_minimize_result) :: r
    integer :: failures

    failures = 0

    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], jac=quadratic_jac, hess=quadratic_hess)
    call check_result('exact quadratic', r, exact_x, exact_fun, exact_nit, &
        exact_nfev, exact_njev, exact_nhev, exact_cg_niter, 1, 1.0e-12_dp)

    r = minimize(quadratic, [-1.0_dp, 1.0_dp], method='trust-constr', &
        jac=quadratic_jac, hess=quadratic_hess)
    call check_result('generic unconstrained dispatch', r, exact_x, exact_fun, &
        exact_nit, exact_nfev, exact_njev, exact_nhev, exact_cg_niter, 1, 1.0e-12_dp)

    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], hess=quadratic_hess)
    ! At -O0/-O2, all work counts agree with SciPy. At -O3, objective
    ! cancellation at the minimizer can make a finite-difference gradient
    ! slightly noisy, causing extra rejected steps and xtol termination.
    if (r%status == 1 .and. r%nit == numeric_gradient_nit) then
        call check_result('numerical objective gradient', r, numeric_gradient_x, &
            numeric_gradient_fun, numeric_gradient_nit, numeric_gradient_nfev, &
            numeric_gradient_njev, numeric_gradient_nhev, numeric_gradient_cg_niter, &
            1, 5.0e-8_dp)
    else
        if (.not. r%success .or. r%status /= 2) call fail('numerical gradient stopping')
        if (maxval(abs(r%x - numeric_gradient_x)) > 5.0e-8_dp) then
            call fail('numerical gradient minimizer')
        end if
        if (abs(r%fun - numeric_gradient_fun) > 1.0e-12_dp) then
            call fail('numerical gradient objective')
        end if
        if (r%nfev < numeric_gradient_nfev .or. r%nit < numeric_gradient_nit .or. &
                r%njev /= numeric_gradient_njev .or. r%nhev /= numeric_gradient_nhev) then
            call fail('numerical gradient derivative counts')
        end if
        call check_empty('numerical objective gradient', r)
    end if

    r = trust_constr(quadratic, [1.0_dp, -2.0_dp], jac=quadratic_jac, hess=quadratic_hess)
    call check_result('stationary initial point', r, stationary_x, stationary_fun, &
        stationary_nit, stationary_nfev, stationary_njev, stationary_nhev, &
        stationary_cg_niter, 1, 1.0e-14_dp)

    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], jac=quadratic_jac, &
        hess=quadratic_hess, maxiter=1)
    call check_result('one-iteration limit', r, maxiter_x, maxiter_fun, &
        maxiter_nit, maxiter_nfev, maxiter_njev, maxiter_nhev, maxiter_cg_niter, &
        0, 1.0e-14_dp)

    r = trust_constr(rosenbrock, [-1.2_dp, 1.0_dp], jac=rosenbrock_jac, &
        hess=rosenbrock_hess)
    call check_result('Rosenbrock analytic Hessian', r, rosenbrock_x, &
        rosenbrock_fun, rosenbrock_nit, rosenbrock_nfev, rosenbrock_njev, &
        rosenbrock_nhev, rosenbrock_cg_niter, 1, 2.0e-8_dp)

    ! As in SciPy, the objective Hessian can be omitted and BFGS approximates it.
    ! BFGS iteration counts are intentionally not asserted (update heuristics differ).
    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], jac=quadratic_jac)
    if (.not. r%success .or. r%status /= 1) call fail('BFGS Hessian convergence')
    if (maxval(abs(r%x - [1.0_dp, -2.0_dp])) > 1.0e-7_dp) call fail('BFGS minimizer')
    if (r%nhev /= 0) call fail('BFGS should not evaluate the objective Hessian')
    call check_empty('BFGS Hessian', r)

    ! A zero-length object array is equivalent to no constraints, including
    ! when it is forwarded through generic minimize.
    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=empty)
    call check_result('empty constraint_objects', r, exact_x, exact_fun, &
        exact_nit, exact_nfev, exact_njev, exact_nhev, exact_cg_niter, 1, 1.0e-12_dp)

    r = minimize(quadratic, [-1.0_dp, 1.0_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=empty)
    call check_result('generic empty constraint_objects', r, exact_x, exact_fun, &
        exact_nit, exact_nfev, exact_njev, exact_nhev, exact_cg_niter, 1, 1.0e-12_dp)

    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=empty, lower=[0.0_dp, -3.0_dp])
    if (.not. r%success) call fail('empty constraint_objects with variable bounds')
    if (maxval(abs(r%x - [1.0_dp, -2.0_dp])) > 1.0e-3_dp) then
        call fail('empty constraint_objects with variable bounds: minimizer')
    end if
    if (.not. allocated(r%constr_nfev)) then
        call fail('bounds-only constraint counts unallocated')
    else if (size(r%constr_nfev) /= 1) then
        call fail('bounds-only constraint count size')
    else if (r%constr_nfev(1) /= 0) then
        call fail('bounds-only constraint evaluation count')
    end if

    ! With zero declared constraint outputs, a callback would be ambiguous.
    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], constraints=dummy_constraints)
    if (r%success .or. r%status /= -1) call fail('callback with zero constraint count')
    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], constraints=dummy_constraints, &
        lower=[0.0_dp, -3.0_dp])
    if (r%success .or. r%status /= -1) then
        call fail('callback with zero constraint count and variable bounds')
    end if

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_trust_constr_unconstrained: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_unconstrained: PASS'

contains

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! objective
        f = (x(1) - 1.0_dp)**2 + 3.0_dp * (x(2) + 2.0_dp)**2
    end function quadratic

    subroutine quadratic_jac(x, g)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: g(:) !! objective gradient
        g = [2.0_dp * (x(1) - 1.0_dp), 6.0_dp * (x(2) + 2.0_dp)]
    end subroutine quadratic_jac

    subroutine quadratic_hess(x, h)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: h(:, :) !! objective Hessian
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 6.0_dp
    end subroutine quadratic_hess

    function rosenbrock(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! objective
        f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
    end function rosenbrock

    subroutine rosenbrock_jac(x, g)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: g(:) !! objective gradient
        g(1) = -400.0_dp * x(1) * (x(2) - x(1)**2) - 2.0_dp * (1.0_dp - x(1))
        g(2) = 200.0_dp * (x(2) - x(1)**2)
    end subroutine rosenbrock_jac

    subroutine rosenbrock_hess(x, h)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: h(:, :) !! objective Hessian
        h(1, 1) = 1200.0_dp * x(1)**2 - 400.0_dp * x(2) + 2.0_dp
        h(1, 2) = -400.0_dp * x(1)
        h(2, 1) = h(1, 2)
        h(2, 2) = 200.0_dp
    end subroutine rosenbrock_hess

    subroutine dummy_constraints(x, c)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: c(:) !! no declared components
        if (size(c) > 0) c = x(1)
    end subroutine dummy_constraints

    subroutine check_result(name, result, expected_x, expected_f, expected_nit, &
            expected_nfev, expected_njev, expected_nhev, expected_cg, expected_status, tolerance)
        character(len=*), intent(in) :: name !! case name
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        real(dp), intent(in) :: expected_x(:) !! SciPy solution
        real(dp), intent(in) :: expected_f !! SciPy objective value
        integer, intent(in) :: expected_nit !! expected iterations
        integer, intent(in) :: expected_nfev !! expected function evaluations
        integer, intent(in) :: expected_njev !! expected gradient evaluations
        integer, intent(in) :: expected_nhev !! expected Hessian evaluations
        integer, intent(in) :: expected_cg !! expected CG iterations
        integer, intent(in) :: expected_status !! expected status
        real(dp), intent(in) :: tolerance !! absolute solution tolerance

        if (result%status /= expected_status) call fail(name // ': termination status')
        if (result%success .neqv. (expected_status == 1)) call fail(name // ': success')
        if (result%nit /= expected_nit) call fail(name // ': nit')
        if (result%nfev /= expected_nfev) call fail(name // ': nfev')
        if (result%njev /= expected_njev) call fail(name // ': njev')
        if (result%nhev /= expected_nhev) call fail(name // ': nhev')
        if (result%cg_niter /= expected_cg) call fail(name // ': cg_niter')
        if (.not. allocated(result%x)) then
            call fail(name // ': solution unallocated')
        else if (maxval(abs(result%x - expected_x)) > tolerance) then
            call fail(name // ': solution')
        end if
        if (abs(result%fun - expected_f) > max(1.0e-12_dp, tolerance)) then
            call fail(name // ': objective value')
        end if
        call check_empty(name, result)
    end subroutine check_result

    subroutine check_empty(name, result)
        character(len=*), intent(in) :: name !! case name
        type(multivariate_minimize_result), intent(in) :: result !! solver result

        if (trim(result%method) /= 'equality_constrained_sqp') call fail(name // ': method')
        if (.not. allocated(result%constr)) then
            call fail(name // ': missing constraint vector')
        else if (size(result%constr) /= 0) then
            call fail(name // ': unexpected constraint values')
        end if
        if (.not. allocated(result%constr_jac)) then
            call fail(name // ': missing constraint Jacobian')
        else if (size(result%constr_jac, 1) /= 0 .or. size(result%constr_jac, 2) /= 2) then
            call fail(name // ': incorrect empty constraint Jacobian shape')
        end if
        if (.not. allocated(result%v)) then
            call fail(name // ': missing multiplier vector')
        else if (size(result%v) /= 0) then
            call fail(name // ': unexpected multipliers')
        end if
        if (.not. allocated(result%constr_nfev)) then
            call fail(name // ': missing constraint counts')
        else if (size(result%constr_nfev) /= 0) then
            call fail(name // ': constraint function counts')
        end if
        if (.not. allocated(result%constr_njev)) then
            call fail(name // ': missing constraint Jacobian counts')
        else if (size(result%constr_njev) /= 0) then
            call fail(name // ': constraint Jacobian counts')
        end if
        if (.not. allocated(result%constr_nhev)) then
            call fail(name // ': missing constraint Hessian counts')
        else if (size(result%constr_nhev) /= 0) then
            call fail(name // ': constraint Hessian counts')
        end if
        if (result%constr_violation /= 0.0_dp .or. result%maxcv /= 0.0_dp) then
            call fail(name // ': nonzero constraint violation')
        end if
    end subroutine check_empty

    subroutine fail(message)
        character(len=*), intent(in) :: message !! failure explanation
        failures = failures + 1
        print '(a)', 'FAIL ' // message
    end subroutine fail
end program test_optimize_trust_constr_unconstrained
