! SPDX-License-Identifier: MIT
! Numerical Hessians of analytic nonlinear-constraint Jacobians.
program test_optimize_trust_constr_constraint_hess_fd
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : linear_constraint, minimize, &
        multivariate_minimize_result, nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_constraint_hess_fd_reference, only : &
        two_point_x, three_point_x, bounded_x
    implicit none

    type(nonlinear_constraint) :: blocks(1)
    type(multivariate_minimize_result) :: r
    type(multivariate_minimize_result) :: r2
    real(dp) :: jac_probes(2, 128)
    integer :: jac_calls
    integer :: probe_count
    integer :: bound_violations
    integer :: failures
    logical :: monitor_bounds

    failures = 0
    jac_calls = 0
    probe_count = 0
    bound_violations = 0
    monitor_bounds = .false.
    blocks(1)%fun => constraint_values
    blocks(1)%jac => constraint_jacobian
    blocks(1)%lower = [1.0_dp]
    blocks(1)%upper = [1.0_dp]

    blocks(1)%hess_fd_scheme = '2-point'
    r = trust_constr(objective, [0.8_dp, 0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, gtol=1.0e-9_dp)
    call check_solution('two-point Hessian', r, two_point_x)
    call check_counts('two-point Hessian', r)
    if (r%constr_njev(1) <= r%nit + 2) call fail('two-point Hessian did not differentiate Jacobian')

    blocks(1)%hess_fd_scheme = '3-point'
    jac_calls = 0
    r2 = minimize(objective, [0.8_dp, 0.9_dp], method='trust-constr', &
        jac=gradient, hess=objective_hessian, constraint_objects=blocks, &
        gtol=1.0e-9_dp)
    call check_solution('three-point Hessian through minimize', r2, three_point_x)
    call check_counts('three-point Hessian', r2)
    if (maxval(abs(r%x - r2%x)) > 1.0e-5_dp) call fail('two/three-point minimizers differ')

    ! Exact constraint-Hessian callbacks cannot be combined with a numerical scheme.
    blocks(1)%hess_fd_scheme = '2-point'
    blocks(1)%hess => exact_constraint_hessian
    r = trust_constr(objective, [0.8_dp, 0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=1)
    if (r%status /= -1) call fail('accepted conflicting Hessian implementations')
    nullify(blocks(1)%hess)

    ! Like SciPy NonlinearConstraint, explicit finite-difference Hessians require
    ! an analytic Jacobian, not two layers of finite differences.
    nullify(blocks(1)%jac)
    r = trust_constr(objective, [0.8_dp, 0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=1)
    if (r%status /= -1) call fail('accepted two numerical-derivative layers')
    blocks(1)%jac => constraint_jacobian

    blocks(1)%hess_fd_scheme = 'cs'
    r = trust_constr(objective, [0.8_dp, 0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=1)
    if (r%status /= -1) call fail('accepted unknown Hessian FD scheme')

    blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[1.0_dp], upper=[1.0_dp])
    blocks(1)%hess_fd_scheme = '2-point'
    r = trust_constr(objective, [0.8_dp, 0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=1)
    if (r%status /= -1) call fail('accepted numerical Hessian on linear block')

    blocks(1)%fun => constraint_values
    blocks(1)%jac => constraint_jacobian
    if (allocated(blocks(1)%matrix)) deallocate(blocks(1)%matrix)
    blocks(1)%hess_fd_scheme = '2-point'
    blocks(1)%finite_diff_rel_step = [-0.1_dp, 0.15_dp]
    jac_calls = 0
    probe_count = 0
    r = trust_constr(objective, [0.8_dp, 0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, maxiter=1)
    call check_counts('signed relative Hessian steps', r)
    if (probe_count < 3) then
        call fail('missing Hessian finite-difference probes')
    else
        if (.not. any(jac_probes(1, 1:probe_count) < 0.8_dp - 1.0e-4_dp .and. &
                abs(jac_probes(2, 1:probe_count) - 0.9_dp) < 1.0e-12_dp)) &
            call fail('negative two-point Hessian relative step ignored')
        if (.not. any(jac_probes(2, 1:probe_count) > 0.9_dp + 1.0e-4_dp .and. &
                abs(jac_probes(1, 1:probe_count) - 0.8_dp) < 1.0e-12_dp)) &
            call fail('positive two-point Hessian relative step ignored')
    end if

    ! A positive step at an active upper bound must be reflected back inside.
    blocks(1)%finite_diff_rel_step = [0.1_dp, 0.15_dp]
    monitor_bounds = .true.
    bound_violations = 0
    r = trust_constr(objective, [0.8_dp, 0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, &
        lower=[0.0_dp, -2.0_dp], upper=[0.8_dp, 2.0_dp], &
        bounds_keep_feasible=[.true., .true.], maxiter=1)
    if (r%status == -1) call fail('bounded two-point Hessian rejected')
    if (bound_violations /= 0) call fail('two-point Hessian probe violated bound')

    blocks(1)%hess_fd_scheme = '3-point'
    bound_violations = 0
    r = trust_constr(objective, [0.8_dp, 0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, &
        lower=[0.0_dp, -2.0_dp], upper=[0.8_dp, 2.0_dp], &
        bounds_keep_feasible=[.true., .true.], maxiter=1)
    if (r%status == -1) call fail('bounded three-point Hessian rejected')
    if (bound_violations /= 0) call fail('three-point Hessian probe violated bound')
    monitor_bounds = .false.

    ! The interior-point method must also handle a numerical constraint
    ! Hessian together with a BFGS objective (no exact objective Hessian).
    if (allocated(blocks(1)%upper)) deallocate(blocks(1)%upper)
    if (allocated(blocks(1)%finite_diff_rel_step)) &
        deallocate(blocks(1)%finite_diff_rel_step)
    blocks(1)%keep_feasible = [.true.]
    jac_calls = 0
    r = minimize(objective, [0.8_dp, 0.9_dp], method='trust-constr', &
        jac=gradient, constraint_objects=blocks, &
        lower=[0.0_dp, -2.0_dp], upper=[0.8_dp, 2.0_dp], &
        bounds_keep_feasible=[.true., .true.], &
        maxiter=300, gtol=1.0e-8_dp)
    if (.not. r%success) then
        print '(a,i0,1x,a)', 'BFGS barrier status=', r%status, trim(r%message)
        call fail('BFGS plus numerical constraint Hessian did not converge')
    else
        if (maxval(abs(r%x - bounded_x)) > 5.0e-3_dp) &
            call fail('BFGS barrier differs from SciPy bounded reference')
        if (r%constr_violation > 1.0e-6_dp) &
            call fail('BFGS barrier failed feasibility')
    end if
    call check_counts('BFGS interior-point Hessian', r)

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_constraint_hess_fd: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_constraint_hess_fd: PASS'

contains

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp) :: f !! quadratic objective
        f = x(1)**2 + 3.0_dp * x(2)**2 + 0.1_dp * x(1) * x(2)
    end function objective

    subroutine gradient(x, g)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: g(:) !! exact objective gradient
        g = [2.0_dp * x(1) + 0.1_dp * x(2), 6.0_dp * x(2) + 0.1_dp * x(1)]
    end subroutine gradient

    subroutine objective_hessian(x, h)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: h(:, :) !! exact objective Hessian
        if (size(x) /= 2) error stop 'unexpected objective dimension'
        h(1, :) = [2.0_dp, 0.1_dp]
        h(2, :) = [0.1_dp, 6.0_dp]
    end subroutine objective_hessian

    subroutine constraint_values(x, c)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: c(:) !! nonlinear constraint value
        c(1) = x(1)**2 + sin(x(2))
    end subroutine constraint_values

    subroutine constraint_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: a(:, :) !! analytic nonlinear constraint Jacobian
        jac_calls = jac_calls + 1
        if (probe_count < size(jac_probes, 2)) then
            probe_count = probe_count + 1
            jac_probes(:, probe_count) = x
        end if
        if (monitor_bounds) then
            if (x(1) < -1.0e-12_dp .or. x(1) > 0.8_dp + 1.0e-12_dp) &
                bound_violations = bound_violations + 1
            if (x(2) < -2.0_dp - 1.0e-12_dp .or. x(2) > 2.0_dp + 1.0e-12_dp) &
                bound_violations = bound_violations + 1
        end if
        a(1, :) = [2.0_dp * x(1), cos(x(2))]
    end subroutine constraint_jacobian

    subroutine exact_constraint_hessian(x, v, h)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(in) :: v(:) !! multiplier
        real(dp), intent(out) :: h(:, :) !! exact weighted constraint Hessian
        h = 0.0_dp
        h(1, 1) = 2.0_dp * v(1)
        h(2, 2) = -sin(x(2)) * v(1)
    end subroutine exact_constraint_hessian

    subroutine check_solution(label, result, reference_x)
        character(len=*), intent(in) :: label !! test label
        type(multivariate_minimize_result), intent(in) :: result !! optimizer result
        real(dp), intent(in) :: reference_x(:) !! SciPy 1.17 reference minimum
        if (.not. result%success) then
            print '(a,a,i0,1x,a)', trim(label), ' status=', result%status, trim(result%message)
            call fail('numerical constraint Hessian failed to converge')
            return
        end if
        if (maxval(abs(result%x - reference_x)) > 1.0e-4_dp) &
            call fail('numerical constraint Hessian differs from SciPy reference')
        if (result%constr_violation > 1.0e-7_dp) &
            call fail('numerical constraint Hessian violated equality')
    end subroutine check_solution

    subroutine check_counts(label, result)
        character(len=*), intent(in) :: label !! test label
        type(multivariate_minimize_result), intent(in) :: result !! optimizer result
        if (.not. allocated(result%constr_njev)) then
            call fail('constraint Jacobian counts unallocated')
            return
        end if
        if (result%constr_njev(1) /= jac_calls) then
            print '(a,2(1x,i0))', trim(label), result%constr_njev(1), jac_calls
            call fail('incorrect finite-difference Hessian Jacobian accounting')
        end if
        if (result%constr_nhev(1) /= 0) call fail('numeric Hessian counted as exact callback')
    end subroutine check_counts

    subroutine fail(msg)
        character(len=*), intent(in) :: msg !! assertion failure message
        failures = failures + 1
        print '(a)', 'FAIL: ' // trim(msg)
    end subroutine fail

end program test_optimize_trust_constr_constraint_hess_fd
