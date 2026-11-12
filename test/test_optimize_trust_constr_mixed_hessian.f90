! SPDX-License-Identifier: MIT
! Mixed exact and numerically differentiated constraint Hessians.
! Reference minimizer obtained with SciPy 1.17.0 optimize.minimize(trust-constr).
program test_optimize_trust_constr_mixed_hessian
    use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, multivariate_minimize_result, &
        nonlinear_constraint, trust_constr
    implicit none

    type(nonlinear_constraint) :: blocks(2)
    type(multivariate_minimize_result) :: r
    type(multivariate_minimize_result) :: exact
    real(dp), parameter :: reference_x(3) = &
        [1.29147583_dp, -0.46790981_dp, 0.48106041_dp]
    real(dp) :: neginf
    integer :: failures
    integer :: block1_products
    integer :: objective_products
    integer :: block1_jac_calls
    integer :: block2_jac_calls
    integer :: block2_value_calls
    integer :: block2_hess_calls
    logical :: outside_retained_bound
    logical :: check_retained_bound

    failures = 0
    block1_products = 0
    objective_products = 0
    block1_jac_calls = 0
    block2_jac_calls = 0
    block2_value_calls = 0
    block2_hess_calls = 0
    outside_retained_bound = .false.
    check_retained_bound = .false.

    blocks(1)%fun => first_values
    blocks(1)%jac => first_jacobian
    blocks(1)%hessp => first_hessian_product
    blocks(1)%lower = [0.0_dp]
    blocks(1)%upper = [0.0_dp]
    blocks(2)%fun => second_values
    blocks(2)%jac => second_jacobian
    blocks(2)%lower = [0.0_dp]
    blocks(2)%upper = [0.0_dp]

    ! The second block has an analytic Jacobian but no Hessian. Its numerical
    ! Hessian must not trigger Jacobian work for the first, exact-Hessian block.
    r = trust_constr(objective, [1.0_dp, -0.3_dp, 1.0_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, gtol=1.0e-9_dp)
    call check_solution('mixed analytic Jacobians', r, 5.0e-5_dp)
    if (size(r%constr_nhev) /= 2) call fail('Hessian count block size')
    if (r%constr_nhev(1) <= 0 .or. r%constr_nhev(2) /= 0) then
        call fail('exact and numerical Hessian accounting')
    end if
    if (block1_products /= 3 * r%constr_nhev(1)) call fail('Hessian product calls')
    if (block1_jac_calls /= r%constr_njev(1)) call fail('first block Jacobian count')
    if (block2_jac_calls /= r%constr_njev(2)) call fail('second block Jacobian count')
    if (block2_jac_calls <= 2 * block1_jac_calls) then
        call fail('second block Hessian did not differentiate Jacobian')
    end if
    if (block2_value_calls /= r%constr_nfev(2)) call fail('second block values count')

    ! An objective Hessian assembled from products must also enable the mixed
    ! exact/numerical constraint-Hessian path through the public dispatcher.
    objective_products = 0
    r = minimize(objective, [1.0_dp, -0.3_dp, 1.0_dp], method='trust-constr', &
        jac=gradient, hessp=objective_hessian_product, constraint_objects=blocks, &
        gtol=1.0e-9_dp)
    call check_solution('objective hessp with mixed blocks', r, 5.0e-5_dp)
    if (objective_products /= 3 * r%nhev) call fail('objective product count')
    if (r%constr_nhev(2) /= 0) call fail('mixed constraint Hessian count with hessp')

    ! An explicit full Hessian on the second block must override the fallback.
    blocks(2)%hess => second_hessian
    block2_hess_calls = 0
    exact = trust_constr(objective, [1.0_dp, -0.3_dp, 1.0_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, gtol=1.0e-9_dp)
    call check_solution('fully analytic blocks', exact, 5.0e-5_dp)
    if (maxval(abs(exact%x - r%x)) > 5.0e-6_dp) then
        call fail('mixed and exact Hessian minimizers differ')
    end if
    if (block2_hess_calls /= exact%constr_nhev(2) .or. block2_hess_calls == 0) then
        call fail('explicit second block Hessian count')
    end if

    ! Two-point/three-point numerical Jacobians also work for a block whose
    ! Hessian must be approximated. No unrelated block is differentiated.
    nullify(blocks(2)%jac)
    nullify(blocks(2)%hess)
    blocks(2)%fd_scheme = '3-point'
    block1_jac_calls = 0
    block2_jac_calls = 0
    block2_value_calls = 0
    r = minimize(objective, [1.0_dp, -0.3_dp, 1.0_dp], method='trust-constr', &
        jac=gradient, hess=objective_hessian, constraint_objects=blocks, &
        gtol=1.0e-9_dp)
    call check_solution('three-point numerical Jacobian', r, 2.0e-3_dp)
    if (r%constr_nhev(2) /= 0 .or. r%constr_njev(2) /= 0) then
        call fail('numerical Jacobian/Hessian callback accounting')
    end if
    if (block2_value_calls /= r%constr_nfev(2) .or. block2_value_calls < 10) then
        call fail('numerical constraint evaluation count')
    end if
    if (block1_jac_calls /= r%constr_njev(1)) call fail('unrelated block derivative count')

    ! Retained lower bounds force one-sided outer Hessian and inner Jacobian
    ! stencils at the starting point. Neither stencil may cross that boundary.
    neginf = ieee_value(0.0_dp, ieee_negative_inf)
    check_retained_bound = .true.
    outside_retained_bound = .false.
    blocks(2)%fd_scheme = '2-point'
    blocks(2)%finite_diff_rel_step = [1.0e-5_dp, 1.0e-5_dp, 1.0e-5_dp]
    r = trust_constr(objective, [1.0_dp, -0.5_dp, 1.0_dp], jac=gradient, &
        hess=objective_hessian, constraint_objects=blocks, &
        lower=[neginf, -0.5_dp, neginf], &
        bounds_keep_feasible=[.false., .true., .false.], gtol=1.0e-8_dp)
    if (outside_retained_bound) call fail('finite-difference probe violated retained bound')
    call check_solution('one-sided bound stencils', r, 5.0e-3_dp)
    if (r%constr_nhev(2) /= 0) call fail('bounded numerical Hessian count')
    check_retained_bound = .false.

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_mixed_hessian: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_mixed_hessian: PASS'

contains

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp) :: f !! scalar objective value
        f = (x(1) - 1.3_dp)**2 + 2.0_dp * (x(2) + 0.8_dp)**2 + &
            0.7_dp * (x(3) - 1.5_dp)**2
    end function objective

    subroutine gradient(x, g)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(out) :: g(:) !! analytic objective gradient
        g = [2.0_dp * (x(1) - 1.3_dp), 4.0_dp * (x(2) + 0.8_dp), &
            1.4_dp * (x(3) - 1.5_dp)]
    end subroutine gradient

    subroutine objective_hessian(x, h)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(out) :: h(:, :) !! constant objective Hessian
        if (size(x) /= 3) error stop 'unexpected dimension'
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 4.0_dp
        h(3, 3) = 1.4_dp
    end subroutine objective_hessian

    subroutine objective_hessian_product(x, p, hp)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(in) :: p(:) !! Hessian-vector product direction
        real(dp), intent(out) :: hp(:) !! objective Hessian times direction
        if (size(x) /= 3) error stop 'unexpected dimension'
        objective_products = objective_products + 1
        hp = [2.0_dp * p(1), 4.0_dp * p(2), 1.4_dp * p(3)]
    end subroutine objective_hessian_product

    subroutine first_values(x, c)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(out) :: c(:) !! first nonlinear equality value
        c(1) = x(1)**2 + x(2) - 1.2_dp
    end subroutine first_values

    subroutine first_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(out) :: a(:, :) !! first constraint Jacobian
        block1_jac_calls = block1_jac_calls + 1
        a(1, :) = [2.0_dp * x(1), 1.0_dp, 0.0_dp]
    end subroutine first_jacobian

    subroutine first_hessian_product(x, v, p, hp)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(in) :: v(:) !! first constraint multiplier
        real(dp), intent(in) :: p(:) !! product direction
        real(dp), intent(out) :: hp(:) !! multiplier-weighted Hessian product
        if (size(x) /= 3) error stop 'unexpected dimension'
        block1_products = block1_products + 1
        hp = [2.0_dp * v(1) * p(1), 0.0_dp, 0.0_dp]
    end subroutine first_hessian_product

    subroutine second_values(x, c)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(out) :: c(:) !! second nonlinear equality value
        block2_value_calls = block2_value_calls + 1
        if (check_retained_bound .and. x(2) < -0.5_dp - 1.0e-12_dp) then
            outside_retained_bound = .true.
        end if
        c(1) = x(2)**2 + x(3) - 0.7_dp
    end subroutine second_values

    subroutine second_jacobian(x, a)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(out) :: a(:, :) !! second constraint Jacobian
        block2_jac_calls = block2_jac_calls + 1
        a(1, :) = [0.0_dp, 2.0_dp * x(2), 1.0_dp]
    end subroutine second_jacobian

    subroutine second_hessian(x, v, h)
        real(dp), intent(in) :: x(:) !! current parameter vector
        real(dp), intent(in) :: v(:) !! second constraint multiplier
        real(dp), intent(out) :: h(:, :) !! second weighted constraint Hessian
        if (size(x) /= 3) error stop 'unexpected dimension'
        block2_hess_calls = block2_hess_calls + 1
        h = 0.0_dp
        h(2, 2) = 2.0_dp * v(1)
    end subroutine second_hessian

    subroutine check_solution(name, rr, tolerance)
        character(len=*), intent(in) :: name !! case label
        type(multivariate_minimize_result), intent(in) :: rr !! solver result
        real(dp), intent(in) :: tolerance !! coordinate absolute tolerance
        if (.not. rr%success) then
            print '(3a,i0)', trim(name), ': ', trim(rr%message), rr%status
            call fail(trim(name)//' solver status')
        end if
        if (maxval(abs(rr%x - reference_x)) > tolerance) then
            write(*, '(a,3f18.9)') trim(name)//' x:', rr%x
            call fail(trim(name)//' reference optimizer')
        end if
        if (rr%constr_violation > 2.0e-5_dp) call fail(trim(name)//' feasibility')
        if (rr%optimality > 3.0e-5_dp) call fail(trim(name)//' optimality')
    end subroutine check_solution

    subroutine fail(message)
        character(len=*), intent(in) :: message !! case failure description
        print '(a)', 'FAIL: '//trim(message)
        failures = failures + 1
    end subroutine fail
end program test_optimize_trust_constr_mixed_hessian
