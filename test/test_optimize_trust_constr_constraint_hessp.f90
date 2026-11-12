! SPDX-License-Identifier: MIT
! Regression tests for multiplier-weighted nonlinear constraint Hessian products.
program test_optimize_trust_constr_constraint_hessp
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : linear_constraint, minimize, &
        multivariate_minimize_result, nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_constraint_hessp_reference
    implicit none

    type(multivariate_minimize_result) :: r
    type(multivariate_minimize_result) :: dense
    type(nonlinear_constraint) :: blocks(2)
    integer :: failures
    integer :: block_products
    integer :: full_calls
    integer :: objective_products

    failures = 0
    block_products = 0
    full_calls = 0
    objective_products = 0

    call setup_cubic_block(blocks(1))
    r = trust_constr(objective, [0.6_dp, -1.4_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks(1:1))
    call check_solution('nonconstant equality', r, equality_x, 1.0e-5_dp)
    call check_counts('nonconstant equality', r, equality_nit, equality_nfev, &
        equality_njev, equality_nhev)
    call check_block_products('nonconstant equality', r, 1)
    if (r%constr_nhev(1) /= equality_constr_nhev) call fail('equality constraint Hessian count')

    ! A full weighted Hessian must produce the same result, not just a
    ! convergent solution. The off-diagonal and x-dependent entries matter.
    blocks(1)%hess => cubic_hessian
    nullify(blocks(1)%hessp)
    full_calls = 0
    dense = trust_constr(objective, [0.6_dp, -1.4_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks(1:1))
    call check_solution('full weighted Hessian', dense, equality_x, 1.0e-5_dp)
    call check_equivalent('products versus full', r, dense, 1.0e-7_dp)
    if (full_calls /= dense%constr_nhev(1)) call fail('full Hessian callback count')

    ! When both callbacks are supplied, a full Hessian takes precedence.
    blocks(1)%hessp => cubic_product
    block_products = 0
    full_calls = 0
    dense = trust_constr(objective, [0.6_dp, -1.4_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks(1:1))
    if (block_products /= 0) call fail('full constraint Hessian should override hessp')
    if (full_calls /= dense%constr_nhev(1)) call fail('Hessian precedence count')

    ! The generic minimize dispatcher preserves the independent block pointer.
    nullify(blocks(1)%hess)
    block_products = 0
    r = minimize(objective, [0.6_dp, -1.4_dp], method='trust-constr', &
        jac=gradient, hess=objective_hessian, constraint_objects=blocks(1:1))
    call check_solution('generic dispatcher', r, equality_x, 1.0e-5_dp)
    call check_block_products('generic dispatcher', r, 1)

    ! Positive-slack barrier path and the simultaneous objective hessp adapter.
    blocks(1)%lower = [-0.1_dp]
    blocks(1)%upper = [0.6_dp]
    block_products = 0
    r = trust_constr(objective, [0.6_dp, -1.4_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks(1:1))
    call check_solution('nonlinear inequality', r, inequality_x, 3.0e-3_dp)
    call check_block_products('nonlinear inequality', r, 1)
    if (r%constr_nhev(1) /= inequality_constr_nhev) then
        call fail('inequality constraint Hessian count')
    end if
    block_products = 0
    objective_products = 0
    dense = trust_constr(objective, [0.6_dp, -1.4_dp], jac=gradient, &
        hessp=objective_product, constraint_objects=blocks(1:1))
    call check_solution('objective and constraint products', dense, inequality_x, 3.0e-3_dp)
    call check_block_products('objective and constraint products', dense, 1)
    if (objective_products /= 2 * dense%nhev) call fail('objective product assembly count')

    ! An independently supplied LinearConstraint contributes no Hessian calls.
    call setup_cubic_block(blocks(1))
    blocks(2) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[0.0_dp], upper=[0.0_dp])
    block_products = 0
    r = trust_constr(objective, [1.4_dp, -1.6_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks)
    call check_solution('linear and product blocks', r, &
        [(1.0_dp + sqrt(5.0_dp)) / 2.0_dp, -(1.0_dp + sqrt(5.0_dp)) / 2.0_dp], 1.0e-5_dp)
    if (size(r%constr_nhev) /= 2) call fail('mixed block Hessian count size')
    if (r%constr_nhev(2) /= 0 .or. r%constr_nfev(2) /= 0 .or. &
            r%constr_njev(2) /= 0) call fail('linear block should have zero evaluation counts')
    call check_block_products('linear and product blocks', r, 1)

    ! Two nonlinear blocks may supply their exact Hessians in different forms.
    call setup_first_quadratic(blocks(1))
    call setup_second_quadratic(blocks(2))
    block_products = 0
    full_calls = 0
    r = trust_constr(objective, [1.1_dp, -0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks)
    call check_solution('mixed nonlinear Hessian forms', r, [1.0_dp, -1.0_dp], 1.0e-6_dp)
    call check_block_products('mixed nonlinear Hessian forms', r, 1)
    if (r%constr_nhev(2) <= 0 .or. full_calls /= r%constr_nhev(2)) then
        call fail('second nonlinear full Hessian not used')
    end if

    ! A two-output block must receive the full multiplier vector in order.
    blocks(1) = nonlinear_constraint()
    blocks(1)%fun => two_values
    blocks(1)%jac => two_jacobian
    blocks(1)%hessp => two_product
    blocks(1)%lower = [0.0_dp, 0.0_dp]
    blocks(1)%upper = [0.0_dp, 0.0_dp]
    block_products = 0
    r = trust_constr(objective, [1.1_dp, -0.9_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks(1:1))
    call check_solution('two-output multiplier vector', r, [1.0_dp, -1.0_dp], 1.0e-6_dp)
    call check_block_products('two-output multiplier vector', r, 1)

    ! A constant-matrix block cannot request a nonzero Hessian callback.
    blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[0.0_dp], upper=[0.0_dp])
    blocks(1)%hessp => cubic_product
    r = trust_constr(objective, [0.6_dp, -1.4_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks(1:1))
    if (r%success .or. r%status >= 0) call fail('linear block hessp must be rejected')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_constraint_hessp: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_constraint_hessp: PASS'

contains

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! optimization parameters
        real(dp) :: f !! quadratic objective value
        f = (x(1) - 1.0_dp)**2 + 2.0_dp * (x(2) + 2.0_dp)**2
    end function objective

    subroutine gradient(x, g)
        real(dp), intent(in) :: x(:) !! optimization parameters
        real(dp), intent(out) :: g(:) !! gradient
        g = [2.0_dp * (x(1) - 1.0_dp), 4.0_dp * (x(2) + 2.0_dp)]
    end subroutine gradient

    subroutine objective_hessian(x, hh)
        real(dp), intent(in) :: x(:) !! optimization parameters
        real(dp), intent(out) :: hh(:, :) !! objective Hessian
        if (size(x) /= 2) error stop 'unexpected objective dimension'
        hh = 0.0_dp
        hh(1, 1) = 2.0_dp
        hh(2, 2) = 4.0_dp
    end subroutine objective_hessian

    subroutine objective_product(x, p, hp)
        real(dp), intent(in) :: x(:) !! optimization parameters
        real(dp), intent(in) :: p(:) !! product direction
        real(dp), intent(out) :: hp(:) !! objective Hessian-vector product
        if (size(x) /= 2) error stop 'unexpected objective dimension'
        objective_products = objective_products + 1
        hp = [2.0_dp * p(1), 4.0_dp * p(2)]
    end subroutine objective_product

    subroutine cubic_values(x, c)
        real(dp), intent(in) :: x(:) !! optimization parameters
        real(dp), intent(out) :: c(:) !! nonlinear constraint value
        c(1) = x(1)**3 + x(1) * x(2) + x(2)
    end subroutine cubic_values

    subroutine cubic_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! optimization parameters
        real(dp), intent(out) :: a(:, :) !! nonlinear constraint Jacobian
        a(1, :) = [3.0_dp * x(1)**2 + x(2), x(1) + 1.0_dp]
    end subroutine cubic_jacobian

    subroutine cubic_hessian(x, v, hh)
        real(dp), intent(in) :: x(:) !! optimization parameters
        real(dp), intent(in) :: v(:) !! multiplier for cubic constraint
        real(dp), intent(out) :: hh(:, :) !! weighted Hessian
        full_calls = full_calls + 1
        hh = reshape([6.0_dp * x(1), 1.0_dp, 1.0_dp, 0.0_dp], [2, 2]) * v(1)
    end subroutine cubic_hessian

    subroutine cubic_product(x, v, p, hp)
        real(dp), intent(in) :: x(:) !! optimization parameters
        real(dp), intent(in) :: v(:) !! multiplier for cubic constraint
        real(dp), intent(in) :: p(:) !! product direction
        real(dp), intent(out) :: hp(:) !! weighted Hessian-vector product
        block_products = block_products + 1
        hp(1) = v(1) * (6.0_dp * x(1) * p(1) + p(2))
        hp(2) = v(1) * p(1)
    end subroutine cubic_product

    subroutine setup_cubic_block(block)
        type(nonlinear_constraint), intent(out) :: block !! configured cubic equality block
        block%fun => cubic_values
        block%jac => cubic_jacobian
        block%hessp => cubic_product
        block%lower = [0.0_dp]
        block%upper = [0.0_dp]
    end subroutine setup_cubic_block

    subroutine first_quadratic(x, c)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: c(:) !! first constraint value
        c(1) = x(1)**2 + x(2)
    end subroutine first_quadratic

    subroutine first_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: a(:, :) !! first constraint Jacobian
        a(1, :) = [2.0_dp * x(1), 1.0_dp]
    end subroutine first_jacobian

    subroutine first_product(x, v, p, hp)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(in) :: v(:) !! first constraint multiplier
        real(dp), intent(in) :: p(:) !! product direction
        real(dp), intent(out) :: hp(:) !! weighted Hessian product
        if (size(x) /= 2) error stop 'unexpected Hessian dimension'
        block_products = block_products + 1
        hp = [2.0_dp * v(1) * p(1), 0.0_dp]
    end subroutine first_product

    subroutine setup_first_quadratic(block)
        type(nonlinear_constraint), intent(out) :: block !! first nonlinear block
        block%fun => first_quadratic
        block%jac => first_jacobian
        block%hessp => first_product
        block%lower = [0.0_dp]
        block%upper = [0.0_dp]
    end subroutine setup_first_quadratic

    subroutine second_quadratic(x, c)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: c(:) !! second constraint value
        c(1) = x(1) + x(2)**2 - 2.0_dp
    end subroutine second_quadratic

    subroutine second_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: a(:, :) !! second constraint Jacobian
        a(1, :) = [1.0_dp, 2.0_dp * x(2)]
    end subroutine second_jacobian

    subroutine second_hessian(x, v, hh)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(in) :: v(:) !! second constraint multiplier
        real(dp), intent(out) :: hh(:, :) !! weighted Hessian
        if (size(x) /= 2) error stop 'unexpected Hessian dimension'
        full_calls = full_calls + 1
        hh = 0.0_dp
        hh(2, 2) = 2.0_dp * v(1)
    end subroutine second_hessian

    subroutine setup_second_quadratic(block)
        type(nonlinear_constraint), intent(out) :: block !! second nonlinear block
        block%fun => second_quadratic
        block%jac => second_jacobian
        block%hess => second_hessian
        block%lower = [0.0_dp]
        block%upper = [0.0_dp]
    end subroutine setup_second_quadratic

    subroutine two_values(x, c)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: c(:) !! two original constraint values
        c(1) = x(1)**2 + x(2)
        c(2) = x(1) + x(2)**2 - 2.0_dp
    end subroutine two_values

    subroutine two_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(out) :: a(:, :) !! two-row Jacobian
        a(1, :) = [2.0_dp * x(1), 1.0_dp]
        a(2, :) = [1.0_dp, 2.0_dp * x(2)]
    end subroutine two_jacobian

    subroutine two_product(x, v, p, hp)
        real(dp), intent(in) :: x(:) !! parameters
        real(dp), intent(in) :: v(:) !! two distinct constraint multipliers
        real(dp), intent(in) :: p(:) !! product direction
        real(dp), intent(out) :: hp(:) !! weighted block Hessian product
        if (size(x) /= 2) error stop 'unexpected Hessian dimension'
        block_products = block_products + 1
        hp = [2.0_dp * v(1) * p(1), 2.0_dp * v(2) * p(2)]
    end subroutine two_product

    subroutine check_solution(label, result, expected, tolerance)
        character(len=*), intent(in) :: label !! test description
        type(multivariate_minimize_result), intent(in) :: result !! computed solution
        real(dp), intent(in) :: expected(:) !! expected coordinates
        real(dp), intent(in) :: tolerance !! allowed maximum absolute error
        print '(a,2f15.9,a,i0,a,i0)', trim(label)//' x=', result%x, &
            ' nit=', result%nit, ' status=', result%status
        if (.not. result%success) call fail(label//': solver did not converge')
        if (maxval(abs(result%x - expected)) > tolerance) call fail(label//': solution mismatch')
    end subroutine check_solution

    subroutine check_counts(label, result, nit, nfev, njev, nhev)
        character(len=*), intent(in) :: label !! test description
        type(multivariate_minimize_result), intent(in) :: result !! computed solution
        integer, intent(in) :: nit !! reference iteration count
        integer, intent(in) :: nfev !! reference objective function count
        integer, intent(in) :: njev !! reference gradient count
        integer, intent(in) :: nhev !! reference Hessian count
        if (result%nit /= nit .or. result%nfev /= nfev .or. &
                result%njev /= njev .or. result%nhev /= nhev) then
            call fail(label//': SciPy work-count mismatch')
        end if
    end subroutine check_counts

    subroutine check_block_products(label, result, i)
        character(len=*), intent(in) :: label !! test description
        type(multivariate_minimize_result), intent(in) :: result !! computed solution
        integer, intent(in) :: i !! index of Hessian-product constraint block
        if (.not. allocated(result%constr_nhev)) then
            call fail(label//': missing constraint counters')
            return
        end if
        if (result%constr_nhev(i) <= 0 .or. &
                block_products /= 2 * result%constr_nhev(i)) then
            call fail(label//': wrong product callback accounting')
        end if
    end subroutine check_block_products

    subroutine check_equivalent(label, first, second, tolerance)
        character(len=*), intent(in) :: label !! test description
        type(multivariate_minimize_result), intent(in) :: first !! product-based result
        type(multivariate_minimize_result), intent(in) :: second !! full-Hessian result
        real(dp), intent(in) :: tolerance !! allowed absolute difference
        if (maxval(abs(first%x - second%x)) > tolerance) call fail(label//': coordinate difference')
        if (abs(first%fun - second%fun) > tolerance) call fail(label//': objective difference')
        if (first%nit /= second%nit .or. first%nhev /= second%nhev) then
            call fail(label//': evaluation-count difference')
        end if
    end subroutine check_equivalent

    subroutine fail(message)
        character(len=*), intent(in) :: message !! failed assertion
        print '(a)', 'FAIL: '//message
        failures = failures + 1
    end subroutine fail
end program test_optimize_trust_constr_constraint_hessp
