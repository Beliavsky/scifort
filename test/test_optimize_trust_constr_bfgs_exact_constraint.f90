! SPDX-License-Identifier: MIT
! SciPy 1.17.0 reference: optimize.minimize(method='trust-constr'), BFGS
! objective Hessian plus exact nonlinear constraint Hessian.
program test_optimize_trust_constr_bfgs_exact_constraint
    use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, multivariate_minimize_result, &
        nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_bfgs_exact_constraint_reference, only : &
        equality_x, objective_fd_x, upper_x, bounds_x, mixed_x
    implicit none

    type(nonlinear_constraint) :: block2(1)
    type(nonlinear_constraint) :: block3(2)
    type(multivariate_minimize_result) :: r
    real(dp) :: posinf
    real(dp) :: neginf
    integer :: failures
    integer :: hessian_calls
    integer :: first_hessian_calls
    integer :: first_product_calls
    integer :: second_jacobian_calls

    failures = 0
    hessian_calls = 0
    first_hessian_calls = 0
    first_product_calls = 0
    second_jacobian_calls = 0
    posinf = ieee_value(0.0_dp, ieee_positive_inf)
    neginf = ieee_value(0.0_dp, ieee_negative_inf)

    ! Equalities can use an exact weighted constraint Hessian even when the
    ! objective Hessian is the default BFGS update. No objective Hessian
    ! evaluations should be reported.
    r = trust_constr(objective_2d, [0.2_dp, 0.2_dp], jac=gradient_2d, &
        constraints=circle_equality, constraint_jac=circle_jacobian, &
        constraint_hess=circle_hessian, meq=1, gtol=1.0e-10_dp, maxiter=500)
    call check_2d('direct equality', r, equality_x, 1.0e-5_dp)
    call check_hessian_counts('direct equality', r, .false.)
    if (trim(r%method) /= 'equality_constrained_sqp') call fail('direct equality method')

    ! The objective gradient may be numerical as well. Three-point objective
    ! differences must not disable the supplied exact constraint Hessian.
    hessian_calls = 0
    r = trust_constr(objective_2d, [0.2_dp, 0.2_dp], &
        constraints=circle_equality, constraint_jac=circle_jacobian, &
        constraint_hess=circle_hessian, meq=1, objective_fd_scheme='3-point', &
        gtol=1.0e-9_dp, maxiter=500)
    call check_2d('numerical objective gradient', r, objective_fd_x, 1.0e-5_dp)
    call check_hessian_counts('numerical objective gradient', r, .false.)
    if (r%nfev <= r%njev) call fail('objective finite-difference evaluations')

    ! A nonlinear inequality enters the barrier path. Its exact Hessian is
    ! still used while the objective Hessian remains BFGS.
    hessian_calls = 0
    r = minimize(objective_2d, [0.2_dp, 0.2_dp], method='trust-constr', &
        jac=gradient_2d, constraints=circle_value, constraint_jac=circle_jacobian, &
        constraint_hess=circle_hessian, constraint_upper=[1.0_dp], &
        gtol=1.0e-10_dp, maxiter=500)
    call check_2d('direct upper inequality', r, upper_x, 8.0e-4_dp)
    call check_hessian_counts('direct upper inequality', r, .false.)
    if (trim(r%method) /= 'tr_interior_point') call fail('direct inequality method')
    if (sum(r%x**2) > 1.0_dp + 1.0e-5_dp) call fail('direct inequality feasibility')

    ! A SciPy-style constraint block plus variable bounds exercises the
    ! second wrapper, with a zero-work bounds diagnostics block.
    block2(1)%fun => circle_value
    block2(1)%jac => circle_jacobian
    block2(1)%hess => circle_hessian
    block2(1)%lower = [neginf]
    block2(1)%upper = [1.0_dp]
    hessian_calls = 0
    r = trust_constr(objective_2d, [0.2_dp, 0.2_dp], jac=gradient_2d, &
        constraint_objects=block2, lower=[0.0_dp, 0.0_dp], &
        upper=[posinf, posinf], gtol=1.0e-10_dp, maxiter=500)
    call check_2d('object inequality and bounds', r, bounds_x, 1.0e-3_dp)
    call check_hessian_counts('object inequality and bounds', r, .true.)
    if (any(r%x < -1.0e-8_dp)) call fail('object variable bounds')

    ! The analytic Hessian of the first nonlinear block must coexist with a
    ! numerical Hessian of the second block (from its analytic Jacobian),
    ! while the objective itself has no supplied Hessian.
    block3(1)%fun => first_value
    block3(1)%jac => first_jacobian
    block3(1)%hess => first_hessian
    block3(1)%lower = [0.0_dp]
    block3(1)%upper = [0.0_dp]
    block3(2)%fun => second_value
    block3(2)%jac => second_jacobian
    block3(2)%lower = [0.0_dp]
    block3(2)%upper = [0.0_dp]
    first_hessian_calls = 0
    second_jacobian_calls = 0
    r = minimize(objective_3d, [1.0_dp, -0.3_dp, 1.0_dp], method='trust-constr', &
        jac=gradient_3d, constraint_objects=block3, gtol=1.0e-9_dp, maxiter=500)
    if (.not. r%success) then
        print '(a,i0,2a)', 'mixed BFGS: status ', r%status, ' ', trim(r%message)
        call fail('mixed BFGS status')
    end if
    if (maxval(abs(r%x - mixed_x)) > 6.0e-4_dp) call fail('mixed BFGS optimizer')
    if (r%constr_violation > 2.0e-6_dp) call fail('mixed BFGS feasibility')
    if (r%optimality > 1.0e-5_dp) call fail('mixed BFGS optimality')
    if (r%nhev /= 0) call fail('mixed BFGS objective Hessian count')
    if (size(r%constr_nhev) /= 2) then
        call fail('mixed BFGS Hessian count size')
    else
        if (r%constr_nhev(1) /= first_hessian_calls .or. first_hessian_calls <= 0) &
            call fail('mixed BFGS exact Hessian accounting')
        if (r%constr_nhev(2) /= 0) call fail('mixed BFGS numerical Hessian count')
        if (second_jacobian_calls <= 2 * r%constr_njev(1)) &
            call fail('mixed BFGS missing numerical curvature probes')
        if (second_jacobian_calls /= r%constr_njev(2)) &
            call fail('mixed BFGS Jacobian count')
    end if

    ! A constraint Hessian operator assembled from hessp must also combine
    ! with objective BFGS. The number of products is n per Hessian assembly.
    nullify(block3(1)%hess)
    block3(1)%hessp => first_hessian_product
    first_product_calls = 0
    r = trust_constr(objective_3d, [1.0_dp, -0.3_dp, 1.0_dp], &
        jac=gradient_3d, constraint_objects=block3, gtol=1.0e-9_dp, maxiter=500)
    if (.not. r%success .or. maxval(abs(r%x - mixed_x)) > 6.0e-4_dp) &
        call fail('constraint hessp with BFGS objective')
    if (first_product_calls /= 3 * r%constr_nhev(1) .or. first_product_calls == 0) &
        call fail('constraint hessp with BFGS objective count')
    if (r%nhev /= 0) call fail('constraint hessp supplied objective Hessian')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_trust_constr_bfgs_exact_constraint: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_bfgs_exact_constraint: PASS'

contains
    function objective_2d(x) result(f)
        real(dp), intent(in) :: x(:) !! current point
        real(dp) :: f !! objective value
        f = sum((x - 1.2_dp)**2)
    end function objective_2d

    subroutine gradient_2d(x, g)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(out) :: g(:) !! gradient
        g = 2.0_dp * (x - 1.2_dp)
    end subroutine gradient_2d

    subroutine circle_value(x, c)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(out) :: c(:) !! circle squared radius
        c(1) = sum(x**2)
    end subroutine circle_value

    subroutine circle_equality(x, c)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(out) :: c(:) !! zero residual at unit radius
        c(1) = sum(x**2) - 1.0_dp
    end subroutine circle_equality

    subroutine circle_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(out) :: a(:, :) !! row Jacobian
        a(1, :) = 2.0_dp * x
    end subroutine circle_jacobian

    subroutine circle_hessian(x, v, h)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(in) :: v(:) !! circle multiplier
        real(dp), intent(out) :: h(:, :) !! Hessian of weighted circle
        integer :: i
        hessian_calls = hessian_calls + 1
        h = 0.0_dp
        do i = 1, size(x)
            h(i, i) = 2.0_dp * v(1)
        end do
    end subroutine circle_hessian

    function objective_3d(x) result(f)
        real(dp), intent(in) :: x(:) !! point
        real(dp) :: f !! objective
        f = (x(1) - 1.3_dp)**2 + 2.0_dp * (x(2) + 0.8_dp)**2 + &
            0.7_dp * (x(3) - 1.5_dp)**2
    end function objective_3d

    subroutine gradient_3d(x, g)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: g(:) !! objective gradient
        g = [2.0_dp * (x(1) - 1.3_dp), 4.0_dp * (x(2) + 0.8_dp), &
            1.4_dp * (x(3) - 1.5_dp)]
    end subroutine gradient_3d

    subroutine first_value(x, c)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: c(:) !! first nonlinear value
        c(1) = x(1)**2 + x(2) - 1.2_dp
    end subroutine first_value

    subroutine first_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: a(:, :) !! first constraint Jacobian
        a(1, :) = [2.0_dp * x(1), 1.0_dp, 0.0_dp]
    end subroutine first_jacobian

    subroutine first_hessian(x, v, h)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(in) :: v(:) !! multiplier
        real(dp), intent(out) :: h(:, :) !! first constraint Hessian
        if (size(x) /= 3) error stop 'invalid first constraint dimension'
        first_hessian_calls = first_hessian_calls + 1
        h = 0.0_dp
        h(1, 1) = 2.0_dp * v(1)
    end subroutine first_hessian

    subroutine first_hessian_product(x, v, p, hp)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(in) :: v(:) !! multiplier
        real(dp), intent(in) :: p(:) !! direction
        real(dp), intent(out) :: hp(:) !! weighted Hessian times direction
        if (size(x) /= 3) error stop 'invalid first constraint dimension'
        first_product_calls = first_product_calls + 1
        hp = [2.0_dp * v(1) * p(1), 0.0_dp, 0.0_dp]
    end subroutine first_hessian_product

    subroutine second_value(x, c)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: c(:) !! second nonlinear value
        c(1) = x(2)**2 + x(3) - 0.7_dp
    end subroutine second_value

    subroutine second_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! point
        real(dp), intent(out) :: a(:, :) !! second constraint Jacobian
        second_jacobian_calls = second_jacobian_calls + 1
        a(1, :) = [0.0_dp, 2.0_dp * x(2), 1.0_dp]
    end subroutine second_jacobian

    subroutine check_2d(name, rr, reference, tolerance)
        character(len=*), intent(in) :: name !! test case
        type(multivariate_minimize_result), intent(in) :: rr !! optimization result
        real(dp), intent(in) :: reference(:) !! SciPy 1.17.0 optimizer
        real(dp), intent(in) :: tolerance !! optimizer coordinate tolerance
        if (.not. rr%success) then
            print '(3a,i0)', trim(name), ': ', trim(rr%message), rr%status
            call fail(trim(name)//' status')
        end if
        if (maxval(abs(rr%x - reference)) > tolerance) then
            print '(a,2f18.9)', trim(name)//' x:', rr%x
            call fail(trim(name)//' optimizer')
        end if
        if (rr%constr_violation > 1.0e-6_dp) call fail(trim(name)//' feasibility')
        if (rr%optimality > 1.0e-5_dp) call fail(trim(name)//' optimality')
        if (rr%nhev /= 0) call fail(trim(name)//' used objective Hessian')
    end subroutine check_2d

    subroutine check_hessian_counts(name, rr, has_bounds)
        character(len=*), intent(in) :: name !! test name
        type(multivariate_minimize_result), intent(in) :: rr !! solver result
        logical, intent(in) :: has_bounds !! includes zero-cost bounds block
        integer :: nblocks
        nblocks = 1
        if (has_bounds) nblocks = 2
        if (size(rr%constr_nhev) /= nblocks) then
            call fail(trim(name)//' Hessian block count')
            return
        end if
        if (rr%constr_nhev(1) /= hessian_calls .or. hessian_calls == 0) then
            call fail(trim(name)//' callback count')
        end if
        if (has_bounds) then
            if (rr%constr_nhev(2) /= 0) call fail(trim(name)//' bounds Hessian count')
        end if
    end subroutine check_hessian_counts

    subroutine fail(name)
        character(len=*), intent(in) :: name !! failure label
        failures = failures + 1
        print '(a)', 'FAIL: '//trim(name)
    end subroutine fail
end program test_optimize_trust_constr_bfgs_exact_constraint
