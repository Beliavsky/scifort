! SPDX-License-Identifier: MIT
! Signed objective finite differences at retained variable bounds.
program test_optimize_trust_constr_objective_signed_bounds
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : linear_constraint, minimize, &
        multivariate_minimize_result, nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_objective_signed_bounds_reference, only : &
        interior_grad, lower_first_grad, lower_second_grad, &
        positive_rel_grad
    implicit none

    type(nonlinear_constraint) :: blocks(1)
    type(multivariate_minimize_result) :: r
    integer :: failures
    integer :: observed
    integer :: outside
    real(dp) :: points(2, 50)
    real(dp) :: retained_lower(2)
    real(dp) :: retained_upper(2)
    real(dp) :: displacement(2)
    logical :: record_points

    failures = 0
    observed = 0
    outside = 0
    record_points = .false.
    blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
        lower=[-1.0_dp], upper=[-1.0_dp])

    call check_initial('interior, direct', [2.0_dp, -3.0_dp], &
        [1.0_dp, -4.0_dp], [3.0_dp, -2.0_dp], [-0.1_dp, 0.2_dp], &
        interior_grad, [1.8_dp, -3.0_dp], [2.0_dp, -3.6_dp])
    call check_initial('interior, minimize and linear block', [2.0_dp, -3.0_dp], &
        [1.0_dp, -4.0_dp], [3.0_dp, -2.0_dp], [-0.1_dp, 0.2_dp], &
        interior_grad, [1.8_dp, -3.0_dp], [2.0_dp, -3.6_dp], &
        use_blocks=.true., use_minimize=.true.)
    call check_initial('first coordinate at lower bound', [1.0_dp, -3.0_dp], &
        [1.0_dp, -4.0_dp], [3.0_dp, -2.0_dp], [-0.1_dp, 0.2_dp], &
        lower_first_grad, [1.1_dp, -3.0_dp], [1.0_dp, -3.6_dp])
    call check_initial('second coordinate at lower bound', [2.0_dp, -4.0_dp], &
        [1.0_dp, -4.0_dp], [3.0_dp, -2.0_dp], [-0.1_dp, 0.2_dp], &
        lower_second_grad, [1.8_dp, -4.0_dp], [2.0_dp, -3.2_dp])
    call check_initial('positive relative steps retain coordinate sign', &
        [2.0_dp, -3.0_dp], [1.0_dp, -4.0_dp], [3.0_dp, -2.0_dp], &
        [0.1_dp, 0.2_dp], positive_rel_grad, &
        [2.2_dp, -3.0_dp], [2.0_dp, -3.6_dp])

    ! Go beyond initial derivative caching: the equality-SQP kernel must use
    ! the requested signed stencil again after an accepted barrier step.
    retained_lower = [1.0_dp, -4.0_dp]
    retained_upper = [3.0_dp, -2.0_dp]
    record_points = .true.
    outside = 0
    observed = 0
    r = trust_constr(quadratic, [2.0_dp, -3.0_dp], &
        hess=quadratic_hessian, constraint_objects=blocks, &
        lower=retained_lower, upper=retained_upper, &
        bounds_keep_feasible=[.true., .true.], &
        objective_finite_diff_rel_step=[-0.1_dp, 0.2_dp], maxiter=2)
    record_points = .false.
    if (r%status == -1 .or. r%njev < 2 .or. r%nfev <= 3) then
        call fail('subsequent barrier/SQP objective gradient not evaluated')
    else
        if (maxval(abs(r%x - [2.0_dp, -3.0_dp])) < 1.0e-5_dp) then
            call fail('short barrier solve made no progress')
        else
            displacement = [-0.1_dp * abs(r%x(1)), -0.2_dp * abs(r%x(2))]
            ! Keep the full requested displacement when the final iterate is
            ! interior enough; otherwise the bound-aware code may reverse it.
            if (r%x(1) + displacement(1) >= retained_lower(1) .and. &
                    r%x(2) + displacement(2) >= retained_lower(2)) then
                if (maxval(abs(r%jac - [2.0_dp * (r%x(1) - 1.2_dp) + displacement(1), &
                        4.0_dp * (r%x(2) + 2.5_dp) + 2.0_dp * displacement(2)])) > 2.0e-7_dp) &
                    call fail('accepted SQP iterate lost signed gradient steps')
            else
                call fail('accepted iterate does not exercise interior signed steps')
            end if
        end if
    end if
    if (outside /= 0) call fail('subsequent objective probes crossed a retained bound')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_trust_constr_objective_signed_bounds: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_objective_signed_bounds: PASS'

contains

    subroutine check_initial(label, x0, lo, hi, rel, grad, first, second, &
            use_blocks, use_minimize)
        character(len=*), intent(in) :: label !! test description
        real(dp), intent(in) :: x0(:) !! initial parameters
        real(dp), intent(in) :: lo(:) !! retained lower bounds
        real(dp), intent(in) :: hi(:) !! retained upper bounds
        real(dp), intent(in) :: rel(:) !! signed relative steps
        real(dp), intent(in) :: grad(:) !! SciPy reference gradient
        real(dp), intent(in) :: first(:) !! first finite-difference probe
        real(dp), intent(in) :: second(:) !! second finite-difference probe
        logical, intent(in), optional :: use_blocks !! add an exact linear block
        logical, intent(in), optional :: use_minimize !! use generic dispatcher
        logical :: with_blocks
        logical :: generic

        with_blocks = .false.
        generic = .false.
        if (present(use_blocks)) with_blocks = use_blocks
        if (present(use_minimize)) generic = use_minimize
        observed = 0
        outside = 0
        retained_lower = lo
        retained_upper = hi
        record_points = .true.
        if (with_blocks) then
            if (generic) then
                r = minimize(cubic, x0, method='trust-constr', constraint_objects=blocks, &
                    lower=lo, upper=hi, bounds_keep_feasible=[.true., .true.], &
                    objective_finite_diff_rel_step=rel, maxiter=1)
            else
                r = trust_constr(cubic, x0, constraint_objects=blocks, &
                    lower=lo, upper=hi, bounds_keep_feasible=[.true., .true.], &
                    objective_finite_diff_rel_step=rel, maxiter=1)
            end if
        else
            r = trust_constr(cubic, x0, lower=lo, upper=hi, &
                bounds_keep_feasible=[.true., .true.], &
                objective_finite_diff_rel_step=rel, maxiter=1)
        end if
        record_points = .false.
        if (r%status == -1) call fail(label // ': optimizer rejected input')
        ! SciFort currently makes extra cached-state preparation calls before
        ! the first barrier stage; the first three evaluations are exactly the
        ! SciPy finite-difference probes, and the count must reflect all work.
        if (r%nfev /= observed .or. r%njev /= 1 .or. r%nit /= 1) &
            call fail(label // ': objective evaluation counts')
        if (observed < 3) then
            call fail(label // ': expected three initial finite-difference calls')
        else
            if (maxval(abs(points(:, 1) - x0)) > 1.0e-12_dp .or. &
                    maxval(abs(points(:, 2) - first)) > 1.0e-12_dp .or. &
                    maxval(abs(points(:, 3) - second)) > 1.0e-12_dp) then
                print '(a)', trim(label)
                print '(a,2f12.6)', ' actual probe 1:', points(:, 2)
                print '(a,2f12.6)', ' actual probe 2:', points(:, 3)
                call fail(label // ': signed probe directions')
            end if
        end if
        if (maxval(abs(r%jac - grad)) > 1.0e-7_dp) then
            print '(a,2f12.6)', ' gradient:', r%jac
            call fail(label // ': gradient disagrees with SciPy')
        end if
        if (outside /= 0) call fail(label // ': probe outside retained bounds')
        if (with_blocks) then
            if (size(r%constr_nfev) /= 2) call fail(label // ': constraint block layout')
            if (r%constr_nfev(1) /= 0) call fail(label // ': linear block callback accounting')
        end if
    end subroutine check_initial

    function cubic(x) result(f)
        real(dp), intent(in) :: x(:) !! decision variables
        real(dp) :: f !! cubic objective value
        call record(x)
        f = x(1)**3 + 2.0_dp * x(2)**3
    end function cubic

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! decision variables
        real(dp) :: f !! convex objective value
        call record(x)
        f = (x(1) - 1.2_dp)**2 + 2.0_dp * (x(2) + 2.5_dp)**2
    end function quadratic

    subroutine quadratic_hessian(x, h)
        real(dp), intent(in) :: x(:) !! decision variables
        real(dp), intent(out) :: h(:, :) !! exact objective Hessian
        if (size(x) /= 2) error stop 'quadratic Hessian argument size'
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 4.0_dp
    end subroutine quadratic_hessian

    subroutine record(x)
        real(dp), intent(in) :: x(:) !! objective evaluation point
        if (.not. record_points) return
        observed = observed + 1
        if (observed <= size(points, 2)) points(:, observed) = x
        if (any(x < retained_lower - 1.0e-12_dp) .or. &
                any(x > retained_upper + 1.0e-12_dp)) outside = outside + 1
    end subroutine record

    subroutine fail(label)
        character(len=*), intent(in) :: label !! failed assertion
        print '(a)', 'FAIL: ' // trim(label)
        failures = failures + 1
    end subroutine fail
end program test_optimize_trust_constr_objective_signed_bounds
