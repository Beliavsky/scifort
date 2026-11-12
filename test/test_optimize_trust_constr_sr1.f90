! SPDX-License-Identifier: MIT
! SciPy 1.17 trust-constr SR1 objective Hessian update regressions.
program test_optimize_trust_constr_sr1
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, multivariate_minimize_result, &
        nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_sr1_reference, only : double_well_x, &
        double_well_fd_x, equality_x, barrier_x, bounds_x
    implicit none

    type(multivariate_minimize_result) :: r
    type(multivariate_minimize_result) :: rbfgs
    type(nonlinear_constraint) :: block(1)
    integer :: failures
    integer :: hess_calls
    real(dp) :: posinf
    real(dp) :: neginf

    failures = 0
    hess_calls = 0
    posinf = ieee_value(0.0_dp, ieee_positive_inf)
    neginf = ieee_value(0.0_dp, ieee_negative_inf)

    ! SR1 permits indefinite updates in a nonconvex objective. BFGS remains
    ! the default if the caller does not select SR1.
    r = trust_constr(double_well, [0.4_dp, 1.0_dp], jac=double_well_grad, &
        objective_hess_update='SR1', gtol=1.0e-9_dp, maxiter=400)
    call check('unconstrained SR1', r, double_well_x, 2.0e-4_dp)
    if (r%nhev /= 0) call fail('SR1 objective counts')
    if (r%cg_niter <= 0) call fail('SR1 projected CG was not entered')

    r = minimize(double_well, [0.4_dp, 1.0_dp], method='trust-constr', &
        objective_hess_update='sr1', objective_fd_scheme='3-point', &
        gtol=1.0e-8_dp, maxiter=400)
    call check('numerical gradient SR1', r, double_well_fd_x, 3.0e-4_dp)
    if (r%nfev <= r%njev) call fail('SR1 gradient should be numerical')

    r = minimize(quadratic, [0.2_dp, 0.2_dp], method='trust-constr', &
        jac=quadratic_grad, constraints=circle_residual, &
        constraint_jac=circle_jac, meq=1, objective_hess_update='sr1', &
        gtol=1.0e-9_dp, maxiter=400)
    call check('nonlinear equality SR1', r, equality_x, 1.0e-4_dp)
    if (r%constr_violation > 1.0e-7_dp) call fail('SR1 equality infeasible')
    if (r%nhev /= 0) call fail('SR1 equality objective counts')

    ! The direct legacy equality interface also accepts an analytic weighted
    ! constraint Hessian when objective curvature is supplied by SR1.
    hess_calls = 0
    r = trust_constr(quadratic, [0.2_dp, 0.2_dp], jac=quadratic_grad, &
        constraints=circle_residual, constraint_jac=circle_jac, &
        constraint_hess=circle_hess, meq=1, objective_hess_update='sr1', &
        gtol=1.0e-9_dp, maxiter=400)
    call check('exact equality Hessian SR1', r, equality_x, 1.0e-4_dp)
    if (r%nhev /= 0 .or. hess_calls == 0) &
        call fail('SR1 equality exact Hessian was not called')
    if (r%constr_nhev(1) /= hess_calls) &
        call fail('SR1 equality exact Hessian counts')

    ! Exact nonlinear constraint curvature coexists with the SR1 objective.
    block(1)%fun => circle_value
    block(1)%jac => circle_jac
    block(1)%hess => circle_hess
    block(1)%lower = [neginf]
    block(1)%upper = [1.0_dp]
    hess_calls = 0
    r = trust_constr(quadratic, [0.2_dp, 0.2_dp], jac=quadratic_grad, &
        constraint_objects=block, lower=[0.0_dp, 0.0_dp], &
        upper=[posinf, posinf], objective_hess_update='sr1', &
        gtol=1.0e-9_dp, maxiter=400)
    call check('inequality and bounds SR1', r, barrier_x, 1.0e-3_dp)
    if (sum(r%x**2) > 1.0_dp + 1.0e-6_dp) call fail('SR1 inequality infeasible')
    if (trim(r%method) /= 'tr_interior_point') call fail('SR1 interior-point method')
    if (r%nhev /= 0) call fail('SR1 inequality objective Hessian count')
    if (size(r%constr_nhev) /= 2) then
        call fail('SR1 block Hessian count shape')
    else
        if (r%constr_nhev(1) /= hess_calls .or. hess_calls < 1) &
            call fail('SR1 exact constraint Hessian callback not used')
        if (r%constr_nhev(2) /= 0) call fail('SR1 bounds Hessian counts')
    end if

    ! The generic minimize dispatcher also forwards SR1 and constraint objects
    ! through the canonicalization and barrier wrappers.
    r = minimize(quadratic, [0.2_dp, 0.2_dp], method='trust-constr', &
        jac=quadratic_grad, constraint_objects=block, &
        lower=[0.0_dp, 0.0_dp], upper=[posinf, posinf], &
        objective_hess_update='sr1', gtol=1.0e-9_dp, maxiter=400)
    call check('generic constraint-object SR1', r, barrier_x, 1.0e-3_dp)
    if (r%nhev /= 0) call fail('generic SR1 objective Hessian count')

    ! Missing nonlinear curvature is updated separately using constraint
    ! BFGS. The objective SR1 update must not count as an exact Hessian.
    r = trust_constr(quadratic, [0.2_dp, 0.2_dp], jac=quadratic_grad, &
        constraints=circle_inside, constraint_jac=circle_inside_jac, &
        mineq=1, objective_hess_update='sr1', gtol=1.0e-9_dp, maxiter=400)
    call check('inequality BFGS plus objective SR1', r, barrier_x, 1.0e-3_dp)
    if (r%nhev /= 0) call fail('inequality objective SR1 Hessian count')

    ! A bounds-only barrier problem has no nonlinear constraint Hessian.
    r = trust_constr(quadratic, [0.0_dp, 0.0_dp], jac=quadratic_grad, &
        upper=[0.5_dp, 0.75_dp], objective_hess_update='sr1', &
        gtol=1.0e-9_dp, maxiter=400)
    call check('bounds-only SR1', r, bounds_x, 3.0e-3_dp)
    if (r%nhev /= 0) call fail('bounds-only SR1 exact Hessian count')

    ! The explicitly selected BFGS path remains byte-for-byte equivalent to
    ! the legacy default at a fixed early iteration budget.
    r = trust_constr(double_well, [0.4_dp, 1.0_dp], jac=double_well_grad, &
        maxiter=4)
    rbfgs = trust_constr(double_well, [0.4_dp, 1.0_dp], jac=double_well_grad, &
        objective_hess_update='BfGs', maxiter=4)
    if (maxval(abs(r%x - rbfgs%x)) > 1.0e-12_dp) call fail('BFGS default changed')
    if (r%status /= rbfgs%status) call fail('BFGS default status changed')

    ! At a fixed short iteration limit the SR1 and BFGS updates must not be
    ! interchangeable: this also checks that the selected strategy is active.
    r = trust_constr(double_well, [0.4_dp, 1.0_dp], jac=double_well_grad, &
        objective_hess_update='SR1', maxiter=4)
    if (maxval(abs(r%x - rbfgs%x)) < 1.0e-5_dp) &
        call fail('SR1 selection had no effect on the iteration')

    ! Input errors return normally, without calling the objective.
    r = trust_constr(quadratic, [0.2_dp, 0.2_dp], objective_hess_update='DFP')
    if (r%status /= -1 .or. r%success) call fail('invalid Hessian update not rejected')
    r = trust_constr(quadratic, [0.2_dp, 0.2_dp], hess=quadratic_hess, &
        objective_hess_update='SR1')
    if (r%status /= -1 .or. r%success) call fail('SR1/exact Hessian conflict not rejected')
    r = trust_constr(quadratic, [0.2_dp, 0.2_dp], hessp=quadratic_hess_product, &
        objective_hess_update='sr1')
    if (r%status /= -1 .or. r%success) call fail('SR1/objective hessp conflict not rejected')
    r = minimize(quadratic, [0.2_dp, 0.2_dp], method='bfgs', &
        objective_hess_update='SR1')
    if (r%status /= -1 .or. r%success) call fail('non-trust-constr SR1 not rejected')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_sr1: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_sr1: PASS'

contains

    function double_well(x) result(f)
        real(dp), intent(in) :: x(:) !! point
        real(dp) :: f !! objective
        f = (x(1)**2 - 1.0_dp)**2 + 0.5_dp * (x(2) - 0.3_dp)**2
    end function double_well

    subroutine double_well_grad(x, g)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: g(:) !! objective gradient
        g = [4.0_dp * x(1)**3 - 4.0_dp * x(1), x(2) - 0.3_dp]
    end subroutine double_well_grad

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! point
        real(dp) :: f !! objective
        f = sum((x - 1.2_dp)**2)
    end function quadratic

    subroutine quadratic_grad(x, g)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: g(:) !! objective gradient
        g = 2.0_dp * (x - 1.2_dp)
    end subroutine quadratic_grad

    subroutine quadratic_hess(x, h)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: h(:, :) !! objective Hessian
        integer :: i
        h = 0.0_dp
        do i = 1, size(x)
            h(i, i) = 2.0_dp
        end do
    end subroutine quadratic_hess

    subroutine quadratic_hess_product(x, p, hp)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(in) :: p(:) !! direction
        real(dp), intent(out) :: hp(:) !! Hessian product
        if (size(x) /= size(p)) error stop 'Hessian vector dimension mismatch'
        hp = 2.0_dp * p
    end subroutine quadratic_hess_product

    subroutine circle_value(x, c)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: c(:) !! circle radius squared
        c(1) = sum(x**2)
    end subroutine circle_value

    subroutine circle_residual(x, c)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: c(:) !! unit circle equality residual
        c(1) = sum(x**2) - 1.0_dp
    end subroutine circle_residual

    subroutine circle_jac(x, a)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: a(:, :) !! row Jacobian
        a(1, :) = 2.0_dp * x
    end subroutine circle_jac

    subroutine circle_inside(x, c)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: c(:) !! positive inside unit circle
        c(1) = 1.0_dp - sum(x**2)
    end subroutine circle_inside

    subroutine circle_inside_jac(x, a)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: a(:, :) !! inequality row Jacobian
        a(1, :) = -2.0_dp * x
    end subroutine circle_inside_jac

    subroutine circle_hess(x, v, h)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(in) :: v(:) !! circle multiplier
        real(dp), intent(out) :: h(:, :) !! weighted constraint Hessian
        integer :: i
        hess_calls = hess_calls + 1
        h = 0.0_dp
        do i = 1, size(x)
            h(i, i) = 2.0_dp * v(1)
        end do
    end subroutine circle_hess

    subroutine check(label, result, expected, tol)
        character(len=*), intent(in) :: label !! test label
        type(multivariate_minimize_result), intent(in) :: result !! solver output
        real(dp), intent(in) :: expected(:) !! SciPy reference optimizer
        real(dp), intent(in) :: tol !! componentwise tolerance
        print '(a,a,2f16.10,a,i0,a,i0,a,es11.3)', trim(label), ':', result%x, &
            ' status=', result%status, ' nit=', result%nit, ' optimality=', result%optimality
        if (.not. result%success) call fail(trim(label)//' convergence')
        if (maxval(abs(result%x - expected)) > tol) call fail(trim(label)//' optimizer')
    end subroutine check

    subroutine fail(message)
        character(len=*), intent(in) :: message !! failing condition
        failures = failures + 1
        print '(2a)', 'FAIL: ', trim(message)
    end subroutine fail
end program test_optimize_trust_constr_sr1
