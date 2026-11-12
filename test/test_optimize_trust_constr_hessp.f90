! SPDX-License-Identifier: MIT
! Focused Hessian-vector product regression tests for trust-constr.
program test_optimize_trust_constr_hessp
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : linear_constraint, minimize, &
        multivariate_minimize_result, nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_hessp_reference
    implicit none

    type(multivariate_minimize_result) :: r
    type(nonlinear_constraint) :: blocks(1)
    integer :: failures
    integer :: hessp_calls
    integer :: hess_calls

    failures = 0
    hessp_calls = 0
    hess_calls = 0

    ! SciPy 1.17 reference: x=[1,-2], nit=nfev=njev=nhev=3.
    r = trust_constr(objective, [-1.0_dp, 1.0_dp], jac=gradient, hessp=hessian_product)
    call check_solution('unconstrained direct', r, unconstrained_x, 1.0e-8_dp)
    call check_counts('unconstrained direct', r, unconstrained_nit, unconstrained_nfev, &
        unconstrained_njev, unconstrained_nhev)
    call check_products('unconstrained direct', r)

    ! The generic dispatcher must forward hessp rather than silently using BFGS.
    hessp_calls = 0
    r = minimize(objective, [-1.0_dp, 1.0_dp], method='trust-constr', &
        jac=gradient, hessp=hessian_product)
    call check_solution('generic minimize', r, unconstrained_x, 1.0e-8_dp)
    call check_counts('generic minimize', r, unconstrained_nit, unconstrained_nfev, &
        unconstrained_njev, unconstrained_nhev)
    call check_products('generic minimize', r)

    ! An explicit Hessian overrides hessp, consistent with SciPy.
    hessp_calls = 0
    hess_calls = 0
    r = trust_constr(objective, [-1.0_dp, 1.0_dp], jac=gradient, &
        hess=full_hessian, hessp=hessian_product)
    call check_solution('hess override', r, unconstrained_x, 1.0e-8_dp)
    if (hessp_calls /= 0 .or. hess_calls /= r%nhev) call fail('hess overrides hessp')

    ! Linear constraints must retain their zero callback counts.
    blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[0.0_dp], upper=[0.0_dp])
    hessp_calls = 0
    r = trust_constr(objective, [-1.0_dp, 1.0_dp], jac=gradient, &
        hessp=hessian_product, constraint_objects=blocks)
    call check_solution('linear block', r, linear_x, 1.0e-7_dp)
    call check_counts('linear block', r, linear_nit, linear_nfev, linear_njev, linear_nhev)
    call check_products('linear block', r)
    if (.not. allocated(r%constr_nfev)) then
        call fail('linear block missing counters')
    else if (size(r%constr_nfev) /= 1) then
        call fail('linear block counter size')
    else if (r%constr_nfev(1) /= 0 .or. r%constr_njev(1) /= 0 .or. &
            r%constr_nhev(1) /= 0) then
        call fail('linear block should require zero constraint callbacks')
    end if

    ! Also check dispatch of both hessp and an independent constraint object
    ! through the generic minimize interface.
    hessp_calls = 0
    r = minimize(objective, [-1.0_dp, 1.0_dp], method='trust-constr', &
        jac=gradient, hessp=hessian_product, constraint_objects=blocks)
    call check_solution('generic linear block', r, linear_x, 1.0e-7_dp)
    call check_products('generic linear block', r)

    ! With finite bounds, the same adapter must work in the barrier path.
    hessp_calls = 0
    r = trust_constr(objective, [-1.0_dp, 1.0_dp], jac=gradient, &
        hessp=hessian_product, lower=[0.0_dp, -2.0_dp], upper=[2.0_dp, 0.0_dp])
    call check_solution('barrier bounds', r, bounds_x, 2.0e-3_dp)
    call check_products('barrier bounds', r)

    ! A nonlinear equality with its own exact weighted Hessian exercises
    ! the Lagrangian-Hessian addition (not only the objective Hessian).
    hessp_calls = 0
    r = trust_constr(objective, [-1.0_dp, 1.0_dp], jac=gradient, &
        hessp=hessian_product, constraints=nonlinear, constraint_jac=nonlinear_jac, &
        constraint_hess=nonlinear_hess, meq=1, mineq=0)
    call check_solution('nonlinear equality', r, nonlinear_x, 2.0e-5_dp)
    call check_counts('nonlinear equality', r, nonlinear_nit, nonlinear_nfev, &
        nonlinear_njev, nonlinear_nhev)
    call check_products('nonlinear equality', r)
    if (abs(r%x(1)**2 + r%x(2)) > 1.0e-6_dp) call fail('nonlinear equality violation')

    ! The same nonlinear constraint provided as an independent object must
    ! route its own Hessian callback while using hessp for the objective.
    if (allocated(blocks(1)%matrix)) deallocate(blocks(1)%matrix)
    blocks(1)%fun => nonlinear
    blocks(1)%jac => nonlinear_jac
    blocks(1)%hess => nonlinear_hess
    blocks(1)%lower = [0.0_dp]
    blocks(1)%upper = [0.0_dp]
    hessp_calls = 0
    r = trust_constr(objective, [-1.0_dp, 1.0_dp], jac=gradient, &
        hessp=hessian_product, constraint_objects=blocks)
    call check_solution('nonlinear block', r, nonlinear_x, 2.0e-5_dp)
    call check_products('nonlinear block', r)
    if (.not. allocated(r%constr_nhev)) then
        call fail('nonlinear block missing Hessian counts')
    else if (size(r%constr_nhev) /= 1) then
        call fail('nonlinear block Hessian counter size')
    else if (r%constr_nhev(1) <= 0) then
        call fail('nonlinear block Hessian callback not used')
    end if

    ! Approximate objective gradients also work with an exact Hessian product.
    hessp_calls = 0
    r = trust_constr(objective, [-1.0_dp, 1.0_dp], hessp=hessian_product, &
        objective_fd_scheme='3-point')
    call check_solution('numerical gradient', r, unconstrained_x, 1.0e-6_dp)
    call check_products('numerical gradient', r)
    if (r%nfev <= r%njev) call fail('numerical gradient finite differences unused')

    ! Off-diagonal and state-dependent Hessian products: the Rosenbrock Hessian
    ! changes at every iterate, unlike the earlier quadratic reference.
    hessp_calls = 0
    r = trust_constr(rosenbrock, [-1.2_dp, 1.0_dp], jac=rosenbrock_gradient, &
        hessp=rosenbrock_product)
    call check_solution('Rosenbrock product', r, rosenbrock_x, 1.0e-5_dp)
    call check_counts('Rosenbrock product', r, rosenbrock_nit, rosenbrock_nfev, &
        rosenbrock_njev, rosenbrock_nhev)
    call check_products('Rosenbrock product', r)

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_hessp: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_hessp: PASS'

contains

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: f !! quadratic objective value
        f = (x(1) - 1.0_dp)**2 + 2.0_dp * (x(2) + 2.0_dp)**2
    end function objective

    subroutine gradient(x, g)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: g(:) !! objective gradient
        g = [2.0_dp * (x(1) - 1.0_dp), 4.0_dp * (x(2) + 2.0_dp)]
    end subroutine gradient

    subroutine full_hessian(x, hh)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hh(:, :) !! explicit Hessian
        if (size(x) /= 2) error stop 'unexpected hessian size'
        hess_calls = hess_calls + 1
        hh = 0.0_dp
        hh(1, 1) = 2.0_dp
        hh(2, 2) = 4.0_dp
    end subroutine full_hessian

    subroutine hessian_product(x, p, hp)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: p(:) !! product direction
        real(dp), intent(out) :: hp(:) !! Hessian-vector product
        if (size(x) /= 2) error stop 'unexpected product size'
        hessp_calls = hessp_calls + 1
        hp = [2.0_dp * p(1), 4.0_dp * p(2)]
    end subroutine hessian_product

    subroutine nonlinear(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! nonlinear equality value
        c(1) = x(1)**2 + x(2)
    end subroutine nonlinear

    subroutine nonlinear_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: a(:, :) !! equality Jacobian
        a(1, :) = [2.0_dp * x(1), 1.0_dp]
    end subroutine nonlinear_jac

    subroutine nonlinear_hess(x, v, hh)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! constraint multiplier
        real(dp), intent(out) :: hh(:, :) !! weighted constraint Hessian
        if (size(x) /= 2) error stop 'unexpected constraint Hessian size'
        hh = 0.0_dp
        hh(1, 1) = 2.0_dp * v(1)
    end subroutine nonlinear_hess

    function rosenbrock(x) result(f)
        real(dp), intent(in) :: x(:) !! Rosenbrock parameters
        real(dp) :: f !! Rosenbrock objective
        f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
    end function rosenbrock

    subroutine rosenbrock_gradient(x, g)
        real(dp), intent(in) :: x(:) !! Rosenbrock parameters
        real(dp), intent(out) :: g(:) !! Rosenbrock gradient
        g = [-400.0_dp * x(1) * (x(2) - x(1)**2) - 2.0_dp * (1.0_dp - x(1)), &
            200.0_dp * (x(2) - x(1)**2)]
    end subroutine rosenbrock_gradient

    subroutine rosenbrock_product(x, p, hp)
        real(dp), intent(in) :: x(:) !! Rosenbrock parameters
        real(dp), intent(in) :: p(:) !! product direction
        real(dp), intent(out) :: hp(:) !! Hessian-vector product
        hessp_calls = hessp_calls + 1
        hp(1) = (1200.0_dp * x(1)**2 - 400.0_dp * x(2) + 2.0_dp) * p(1) - &
            400.0_dp * x(1) * p(2)
        hp(2) = -400.0_dp * x(1) * p(1) + 200.0_dp * p(2)
    end subroutine rosenbrock_product

    subroutine check_solution(name, result, expected, tolerance)
        character(len=*), intent(in) :: name !! regression description
        type(multivariate_minimize_result), intent(in) :: result !! returned solution
        real(dp), intent(in) :: expected(:) !! SciPy 1.17 solution
        real(dp), intent(in) :: tolerance !! absolute coordinate error limit
        print '(a,2f16.9,a,4i6,a,i0)', trim(name)//' x=', result%x, &
            ' nit/nfev/njev/nhev=', result%nit, result%nfev, result%njev, result%nhev, &
            ' status=', result%status
        if (.not. result%success) call fail(name//': solver did not converge')
        if (maxval(abs(result%x - expected)) > tolerance) call fail(name//': SciPy solution difference')
    end subroutine check_solution

    subroutine check_counts(name, result, expected_nit, expected_nfev, expected_njev, expected_nhev)
        character(len=*), intent(in) :: name !! regression description
        type(multivariate_minimize_result), intent(in) :: result !! returned solution
        integer, intent(in) :: expected_nit !! expected iteration count
        integer, intent(in) :: expected_nfev !! expected objective evaluation count
        integer, intent(in) :: expected_njev !! expected gradient evaluation count
        integer, intent(in) :: expected_nhev !! expected Hessian evaluation count
        if (result%nit /= expected_nit .or. result%nfev /= expected_nfev .or. &
                result%njev /= expected_njev .or. result%nhev /= expected_nhev) then
            call fail(name//': SciPy counts differ')
        end if
    end subroutine check_counts

    subroutine check_products(name, result)
        character(len=*), intent(in) :: name !! regression description
        type(multivariate_minimize_result), intent(in) :: result !! solver result
        if (result%nhev <= 0 .or. hessp_calls /= 2 * result%nhev) then
            call fail(name//': Hessian-vector product accounting')
        end if
    end subroutine check_products

    subroutine fail(message)
        character(len=*), intent(in) :: message !! failure reason
        print '(a)', 'FAIL: '//message
        failures = failures + 1
    end subroutine fail
end program test_optimize_trust_constr_hessp
