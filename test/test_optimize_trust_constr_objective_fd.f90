! SPDX-License-Identifier: MIT
! Regression tests for three-point objective differentiation in trust-constr.
program test_optimize_trust_constr_objective_fd
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : linear_constraint, minimize, &
        multivariate_minimize_result, nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_objective_fd_reference
    implicit none

    type(multivariate_minimize_result) :: r
    type(nonlinear_constraint) :: blocks(1)
    real(dp) :: infinity
    integer :: failures
    integer :: n_outside
    integer :: n_upper_outside

    failures = 0
    n_outside = 0
    n_upper_outside = 0
    infinity = ieee_value(0.0_dp, ieee_positive_inf)

    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], hess=quadratic_hess, &
        objective_fd_scheme='3-point')
    call show('direct unconstrained', r)
    call good_solution('direct unconstrained', r, unconstrained_x, 1.0e-8_dp)
    call assert_counts('direct unconstrained', r, unconstrained_nit, &
        unconstrained_nfev, unconstrained_njev, unconstrained_nhev)
    call check_empty(r)

    r = minimize(quadratic, [-1.0_dp, 1.0_dp], method='trust-constr', &
        hess=quadratic_hess, objective_fd_scheme='3-point')
    call show('generic unconstrained', r)
    call good_solution('generic unconstrained', r, unconstrained_x, 1.0e-8_dp)
    call assert_counts('generic unconstrained', r, unconstrained_nit, &
        unconstrained_nfev, unconstrained_njev, unconstrained_nhev)

    ! An exact objective gradient takes priority and does not invoke FD probes.
    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], jac=quadratic_jac, &
        hess=quadratic_hess, objective_fd_scheme='3-point')
    call show('analytic gradient unchanged', r)
    call good_solution('analytic gradient unchanged', r, [1.0_dp, -2.0_dp], 1.0e-12_dp)
    call assert_counts('analytic gradient unchanged', r, 3, 3, 3, 3)

    ! Supplying eps changes the numerical objective scale; it must not be lost
    ! when the solver dispatches to equality SQP.
    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], hess=quadratic_hess, &
        objective_fd_scheme='3-point', eps=1.0e-3_dp, maxiter=1)
    call show('explicit step', r)
    if (r%nit /= 1 .or. r%nfev /= 5 .or. r%njev /= 1) then
        call fail('explicit step work counts')
    end if
    if (maxval(abs(r%jac - [-4.0_dp, 18.0_dp])) > 1.0e-9_dp) then
        call fail('explicit step central derivative')
    end if

    blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[-1.0_dp], upper=[-1.0_dp])
    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], hess=quadratic_hess, &
        constraint_objects=blocks, objective_fd_scheme='3-point')
    call show('linear block', r)
    call good_solution('linear block', r, linear_x, 1.0e-8_dp)
    call assert_counts('linear block', r, linear_nit, linear_nfev, linear_njev, linear_nhev)
    if (allocated(r%constr_nfev)) then
        if (size(r%constr_nfev) /= 1) call fail('linear block count shape')
        if (any(r%constr_nfev /= 0)) call fail('linear constraint evaluated callback')
    else
        call fail('linear block unallocated counts')
    end if

    ! Independent nonlinear constraint has an analytic weighted Hessian but a
    ! three-point Jacobian. The objective has a separate three-point gradient.
    blocks(1) = nonlinear_constraint()
    blocks(1)%fun => nonlinear_value
    blocks(1)%hess => nonlinear_hess
    blocks(1)%lower = [-1.0_dp]
    blocks(1)%upper = [-1.0_dp]
    blocks(1)%fd_scheme = '3-point'
    r = minimize(quadratic, [-1.0_dp, 1.0_dp], method='trust-constr', &
        hess=quadratic_hess, constraint_objects=blocks, objective_fd_scheme='3-point')
    call show('nonlinear block', r)
    call good_solution('nonlinear block', r, nonlinear_x, 1.0e-8_dp)
    call assert_counts('nonlinear block', r, nonlinear_nit, nonlinear_nfev, &
        nonlinear_njev, nonlinear_nhev)
    if (.not. allocated(r%constr_nfev)) then
        call fail('nonlinear block counts unallocated')
    else if (size(r%constr_nfev) /= 1) then
        call fail('nonlinear block count shape')
    else
        if (r%constr_nfev(1) <= r%njev) call fail('three-point constraint probes')
        if (r%constr_njev(1) /= 0) call fail('numerical constraint Jacobian accounting')
        if (r%constr_nhev(1) <= 0) call fail('analytic constraint Hessian accounting')
    end if

    ! Finite-difference probes may not cross a retained lower variable bound.
    n_outside = 0
    r = trust_constr(restricted_quadratic, [0.3_dp, 1.0_dp], hess=quadratic_hess, &
        lower=[0.1_dp, -infinity], bounds_keep_feasible=[.true., .false.], &
        objective_fd_scheme='3-point')
    call show('retained lower bound', r)
    call good_solution('retained lower bound', r, retained_x, 1.0e-3_dp)
    if (n_outside /= 0) call fail('objective evaluated outside strict lower bound')
    if (r%nfev < 2 * size(r%x) * r%njev) call fail('bounded objective work count')

    ! Starting exactly on a retained lower or upper endpoint forces a
    ! second-order one-sided stencil. The reported initial gradient must stay
    ! accurate without evaluating the callback beyond the endpoint.
    n_outside = 0
    r = trust_constr(restricted_quadratic, [0.1_dp, 1.0_dp], hess=quadratic_hess, &
        lower=[0.1_dp, -infinity], bounds_keep_feasible=[.true., .false.], &
        objective_fd_scheme='3-point', maxiter=1)
    call show('one-sided lower stencil', r)
    if (r%nit /= 1 .or. r%status /= 0 .or. r%nhev /= 1) then
        call fail('analytic-Hessian lower endpoint initialization')
    end if
    if (n_outside /= 0) call fail('lower one-sided stencil crossed strict bound')
    if (maxval(abs(r%jac - [-1.8_dp, 18.0_dp])) > 1.0e-7_dp) then
        call fail('incorrect one-sided lower gradient')
    end if

    n_upper_outside = 0
    r = trust_constr(restricted_upper_quadratic, [0.1_dp, 1.0_dp], hess=quadratic_hess, &
        upper=[0.1_dp, infinity], bounds_keep_feasible=[.true., .false.], &
        objective_fd_scheme='3-point', maxiter=1)
    call show('one-sided upper stencil', r)
    if (r%nit /= 1 .or. r%status /= 0 .or. r%nhev /= 1) then
        call fail('analytic-Hessian upper endpoint initialization')
    end if
    if (n_upper_outside /= 0) call fail('upper one-sided stencil crossed strict bound')
    if (maxval(abs(r%jac - [-1.8_dp, 18.0_dp])) > 1.0e-7_dp) then
        call fail('incorrect one-sided upper gradient')
    end if

    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], &
        objective_fd_scheme='cs')
    if (r%success .or. r%status /= -1) call fail('invalid numerical objective scheme')
    r = minimize(quadratic, [-1.0_dp, 1.0_dp], method='trust-constr', &
        objective_fd_scheme='unsupported')
    if (r%success .or. r%status /= -1) call fail('invalid numerical objective scheme through minimize')
    r = trust_constr(quadratic, [-1.0_dp, 1.0_dp], constraint_objects=blocks, &
        objective_fd_scheme='cs')
    if (r%success .or. r%status /= -1) call fail('invalid numerical objective scheme through objects')
    r = minimize(quadratic, [-1.0_dp, 1.0_dp], method='bfgs', &
        objective_fd_scheme='3-point')
    if (r%success .or. r%status /= -1) call fail('trust-only objective scheme on BFGS')

    if (failures > 0) then
        write(*, '(a,i0)') 'test_optimize_trust_constr_objective_fd: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_objective_fd: PASS'

contains

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! quadratic objective
        f = (x(1) - 1.0_dp)**2 + 3.0_dp * (x(2) + 2.0_dp)**2
    end function quadratic

    function restricted_quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! objective defined on the feasible domain
        if (x(1) < 0.1_dp) n_outside = n_outside + 1
        f = quadratic(x)
    end function restricted_quadratic

    function restricted_upper_quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp) :: f !! objective defined only below the strict upper bound
        if (x(1) > 0.1_dp) n_upper_outside = n_upper_outside + 1
        f = quadratic(x)
    end function restricted_upper_quadratic

    subroutine quadratic_jac(x, g)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: g(:) !! gradient
        g = [2.0_dp * (x(1) - 1.0_dp), 6.0_dp * (x(2) + 2.0_dp)]
    end subroutine quadratic_jac

    subroutine quadratic_hess(x, h)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: h(:, :) !! quadratic Hessian
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 6.0_dp
    end subroutine quadratic_hess

    subroutine nonlinear_value(x, c)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: c(:) !! constraint values
        c(1) = x(1)**2 + x(2)
    end subroutine nonlinear_value

    subroutine nonlinear_hess(x, v, h)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(in) :: v(:) !! constraint multipliers
        real(dp), intent(out) :: h(:, :) !! weighted Hessian
        if (size(x) < 2) error stop 'nonlinear Hessian shape'
        h = 0.0_dp
        h(1, 1) = 2.0_dp * v(1)
    end subroutine nonlinear_hess

    subroutine good_solution(name, result, expected, tolerance)
        character(len=*), intent(in) :: name !! case label
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        real(dp), intent(in) :: expected(:) !! expected minimizer
        real(dp), intent(in) :: tolerance !! coordinate tolerance
        if (.not. result%success) call fail(name // ': no convergence')
        if (maxval(abs(result%x - expected)) > tolerance) call fail(name // ': wrong minimizer')
    end subroutine good_solution

    subroutine assert_counts(name, result, nit, nfev, njev, nhev)
        character(len=*), intent(in) :: name !! case label
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        integer, intent(in) :: nit !! SciPy iteration count
        integer, intent(in) :: nfev !! SciPy objective count
        integer, intent(in) :: njev !! SciPy gradient count
        integer, intent(in) :: nhev !! SciPy Hessian count
        if (result%nit /= nit .or. result%nfev /= nfev .or. &
                result%njev /= njev .or. result%nhev /= nhev) then
            call fail(name // ': evaluation counts')
        end if
    end subroutine assert_counts

    subroutine check_empty(result)
        type(multivariate_minimize_result), intent(in) :: result !! unconstrained result
        if (.not. allocated(result%constr)) then
            call fail('empty constraints not allocated')
        else if (size(result%constr) /= 0) then
            call fail('unexpected unconstrained rows')
        end if
    end subroutine check_empty

    subroutine show(name, result)
        character(len=*), intent(in) :: name !! case label
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        print '(a, a, 2f15.8, a, es12.4, a, 4i6, a, i3)', &
            trim(name), ' x=', result%x, ' f=', result%fun, &
            ' nit/nfev/njev/nhev=', result%nit, result%nfev, result%njev, result%nhev, &
            ' status=', result%status
    end subroutine show

    subroutine fail(message)
        character(len=*), intent(in) :: message !! failed assertion
        print '(a)', 'FAIL: ' // message
        failures = failures + 1
    end subroutine fail

end program test_optimize_trust_constr_objective_fd
