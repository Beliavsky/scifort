! SPDX-License-Identifier: MIT
! Regression against SciPy 1.17 NonlinearConstraint finite_diff_jac_sparsity.
program test_optimize_trust_constr_sparse_jac
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, nonlinear_constraint, &
        linear_constraint, multivariate_minimize_result, trust_constr
    use test_optimize_trust_constr_sparse_jac_reference
    implicit none

    type(nonlinear_constraint) :: blocks(1)
    type(nonlinear_constraint) :: linear(1)
    type(multivariate_minimize_result) :: dense
    type(multivariate_minimize_result) :: colored
    type(multivariate_minimize_result) :: result
    integer :: failures
    integer :: scheme
    integer :: exact_hess_calls
    real(dp) :: max_probe
    real(dp) :: infinity
    real(dp) :: x0(4)
    real(dp) :: expected(4)
    real(dp) :: exact_jac(4, 4)

    failures = 0
    exact_hess_calls = 0
    infinity = ieee_value(0.0_dp, ieee_positive_inf)
    x0 = [1.1_dp, 2.1_dp, 3.1_dp, 4.1_dp]
    blocks(1)%fun => constraint_values
    blocks(1)%lower = [1.0_dp, 13.0_dp, 9.0_dp, 16.0_dp]
    blocks(1)%upper = blocks(1)%lower
    allocate(blocks(1)%finite_diff_jac_sparsity(4, 4))
    blocks(1)%finite_diff_jac_sparsity = .false.
    blocks(1)%finite_diff_jac_sparsity(1, 1) = .true.
    blocks(1)%finite_diff_jac_sparsity(2, 2:3) = .true.
    blocks(1)%finite_diff_jac_sparsity(3, 3) = .true.
    blocks(1)%finite_diff_jac_sparsity(4, 4) = .true.

    do scheme = 1, 2
        if (scheme == 1) then
            blocks(1)%fd_scheme = '2-point'
            expected = two_colored_x
        else
            blocks(1)%fd_scheme = '3-point'
            expected = three_colored_x
        end if
        deallocate(blocks(1)%finite_diff_jac_sparsity)
        dense = minimize(objective, x0, method='trust-constr', jac=gradient, &
            hess=hessian, constraint_objects=blocks, gtol=1.0e-9_dp, maxiter=150)
        call check_result('dense reference', dense, expected)
        allocate(blocks(1)%finite_diff_jac_sparsity(4, 4))
        blocks(1)%finite_diff_jac_sparsity = .false.
        blocks(1)%finite_diff_jac_sparsity(1, 1) = .true.
        blocks(1)%finite_diff_jac_sparsity(2, 2:3) = .true.
        blocks(1)%finite_diff_jac_sparsity(3, 3) = .true.
        blocks(1)%finite_diff_jac_sparsity(4, 4) = .true.
        colored = minimize(objective, x0, method='trust-constr', jac=gradient, &
            hess=hessian, constraint_objects=blocks, gtol=1.0e-9_dp, maxiter=150)
        call check_result('colored reference', colored, expected)
        if (.not. allocated(colored%constr_nfev) .or. &
                .not. allocated(dense%constr_nfev)) then
            call fail('constraint call counts missing')
        else
            if (colored%constr_nfev(1) >= dense%constr_nfev(1)) &
                call fail('column coloring did not reduce function evaluations')
            if (scheme == 1) then
                if (colored%constr_nfev(1) /= two_colored_nfev .or. &
                        dense%constr_nfev(1) /= two_dense_nfev) &
                    call fail('two-point sparse function counts differ from SciPy')
            else
                if (colored%constr_nfev(1) /= three_colored_nfev .or. &
                        dense%constr_nfev(1) /= three_dense_nfev) &
                    call fail('three-point sparse function counts differ from SciPy')
            end if
        end if
        exact_jac = 0.0_dp
        exact_jac(1, 1) = 2.0_dp * colored%x(1)
        exact_jac(2, 2) = 2.0_dp * colored%x(2)
        exact_jac(2, 3) = 2.0_dp * colored%x(3)
        exact_jac(3, 3) = 2.0_dp * colored%x(3)
        exact_jac(4, 4) = 2.0_dp * colored%x(4)
        if (.not. allocated(colored%constr_jac)) then
            call fail('colored constraint Jacobian missing')
        else if (maxval(abs(colored%constr_jac - exact_jac)) > 1.0e-6_dp) then
            call fail('colored nonlinear Jacobian inaccurate')
        end if
    end do

    ! An exact weighted constraint Hessian can accompany a grouped numerical
    ! Jacobian. The multiplier vector still has the full block dimension.
    blocks(1)%hess => weighted_hessian
    blocks(1)%fd_scheme = '2-point'
    result = trust_constr(objective, x0, jac=gradient, hess=hessian, &
        constraint_objects=blocks, gtol=1.0e-9_dp, maxiter=150)
    call check_result('colored Jacobian and exact Hessian', result, two_colored_x)
    if (exact_hess_calls <= 0) call fail('weighted Hessian callback not evaluated')
    if (allocated(result%constr_nhev)) then
        if (result%constr_nhev(1) /= exact_hess_calls) &
            call fail('exact Hessian callback count incorrect')
    else
        call fail('exact Hessian count missing')
    end if
    nullify(blocks(1)%hess)

    ! A structurally zero column never needs a probe. It remains zero even
    ! when the corresponding objective variable is unconstrained.
    blocks(1)%fun => independent_of_fourth
    blocks(1)%finite_diff_jac_sparsity(:, 4) = .false.
    result = trust_constr(objective, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
        jac=gradient, hess=hessian, constraint_objects=blocks, maxiter=1)
    if (.not. allocated(result%constr_jac)) then
        call fail('zero-column sparse Jacobian missing')
    else if (any(result%constr_jac(:, 4) /= 0.0_dp)) then
        call fail('structurally zero column has nonzero derivative')
    end if
    if (allocated(result%constr_nfev)) then
        if (result%constr_nfev(1) /= 3) &
            call fail('zero-column probe count should be initial plus two groups')
    else
        call fail('zero-column function count missing')
    end if
    blocks(1)%fun => constraint_values
    blocks(1)%finite_diff_jac_sparsity(4, 4) = .true.

    ! With strict upper bound at the initial point, both combined probes
    ! must switch to one-sided stencils without evaluating outside the box.
    blocks(1)%fd_scheme = '3-point'
    blocks(1)%finite_diff_rel_step = [1.0e-4_dp, 1.0e-4_dp, &
        1.0e-4_dp, 1.0e-4_dp]
    max_probe = -infinity
    result = trust_constr(objective, [1.0_dp, 2.1_dp, 3.1_dp, 4.1_dp], &
        jac=gradient, hess=hessian, constraint_objects=blocks, &
        upper=[1.0_dp, infinity, infinity, infinity], &
        bounds_keep_feasible=[.true., .false., .false., .false.], maxiter=1)
    if (max_probe > 1.0_dp + 2.0e-14_dp) &
        call fail('colored three-point probe crossed retained upper bound')
    if (.not. allocated(result%constr_jac)) then
        call fail('bounded constraint Jacobian missing')
    else if (abs(result%constr_jac(1, 1) - 2.0_dp) > 1.0e-5_dp) then
        call fail('one-sided colored derivative inaccurate')
    end if

    ! An invalid pattern must be rejected before any user callback runs.
    deallocate(blocks(1)%finite_diff_jac_sparsity)
    allocate(blocks(1)%finite_diff_jac_sparsity(3, 4))
    result = trust_constr(objective, x0, jac=gradient, hess=hessian, &
        constraint_objects=blocks)
    if (result%status /= -1 .or. result%success) &
        call fail('invalid sparse Jacobian row shape accepted')
    deallocate(blocks(1)%finite_diff_jac_sparsity)
    allocate(blocks(1)%finite_diff_jac_sparsity(4, 3))
    result = trust_constr(objective, x0, jac=gradient, hess=hessian, &
        constraint_objects=blocks)
    if (result%status /= -1 .or. result%success) &
        call fail('invalid sparse Jacobian column shape accepted')

    ! Linear blocks already have exact constant Jacobians; they cannot
    ! request any of the numerical derivative options.
    linear(1) = linear_constraint(reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [1, 4]), &
        lower=[0.0_dp])
    allocate(linear(1)%finite_diff_jac_sparsity(1, 4))
    linear(1)%finite_diff_jac_sparsity = .true.
    result = trust_constr(objective, x0, jac=gradient, hess=hessian, &
        constraint_objects=linear)
    if (result%status /= -1 .or. result%success) &
        call fail('numerical sparsity accepted on linear block')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_sparse_jac: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_sparse_jac: PASS'

contains

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! current variables
        real(dp) :: f !! squared Euclidean norm
        f = sum(x**2)
    end function objective

    subroutine gradient(x, g)
        real(dp), intent(in) :: x(:) !! current variables
        real(dp), intent(out) :: g(:) !! objective gradient
        g = 2.0_dp * x
    end subroutine gradient

    subroutine hessian(x, h)
        real(dp), intent(in) :: x(:) !! current variables
        real(dp), intent(out) :: h(:, :) !! objective Hessian
        integer :: i
        h = 0.0_dp
        do i = 1, size(x)
            h(i, i) = 2.0_dp
        end do
    end subroutine hessian

    subroutine constraint_values(x, c)
        real(dp), intent(in) :: x(:) !! current variables
        real(dp), intent(out) :: c(:) !! sparse nonlinear equalities
        max_probe = max(max_probe, x(1))
        c = [x(1)**2, x(2)**2 + x(3)**2, x(3)**2, x(4)**2]
    end subroutine constraint_values

    subroutine weighted_hessian(x, v, h)
        real(dp), intent(in) :: x(:) !! current variables
        real(dp), intent(in) :: v(:) !! output multipliers
        real(dp), intent(out) :: h(:, :) !! Hessian of weighted nonlinear outputs
        if (size(x) /= 4) error stop 'unexpected Hessian input dimension'
        exact_hess_calls = exact_hess_calls + 1
        h = 0.0_dp
        h(1, 1) = 2.0_dp * v(1)
        h(2, 2) = 2.0_dp * v(2)
        h(3, 3) = 2.0_dp * (v(2) + v(3))
        h(4, 4) = 2.0_dp * v(4)
    end subroutine weighted_hessian

    subroutine independent_of_fourth(x, c)
        real(dp), intent(in) :: x(:) !! current variables
        real(dp), intent(out) :: c(:) !! values unaffected by fourth variable
        c = [x(1)**2, x(2)**2 + x(3)**2, x(3)**2, 16.0_dp]
    end subroutine independent_of_fourth

    subroutine check_result(label, r, target)
        character(len=*), intent(in) :: label !! test scenario
        type(multivariate_minimize_result), intent(in) :: r !! optimizer result
        real(dp), intent(in) :: target(:) !! SciPy 1.17 reference optimizer variables
        if (.not. r%success .or. r%status /= 1) then
            print '(a,a,a,i0,1x,a)', 'FAILED ', label, ' status=', r%status, trim(r%message)
            failures = failures + 1
        else if (maxval(abs(r%x - target)) > 2.0e-6_dp .or. &
                r%constr_violation > 1.0e-6_dp) then
            print '(a,a,1x,4(es12.4,1x))', 'FAILED solution ', label, r%x
            failures = failures + 1
        end if
    end subroutine check_result

    subroutine fail(label)
        character(len=*), intent(in) :: label !! failed assertion
        print '(a,a)', 'FAILED: ', label
        failures = failures + 1
    end subroutine fail
end program test_optimize_trust_constr_sparse_jac
