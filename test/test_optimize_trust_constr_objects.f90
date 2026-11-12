! SPDX-License-Identifier: MIT
! Regression cases against SciPy 1.17.0's independent NonlinearConstraint objects.
program test_optimize_trust_constr_objects
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, multivariate_minimize_result, &
        nonlinear_constraint, trust_constr
    use test_optimize_trust_constr_objects_reference
    implicit none

    type(nonlinear_constraint) :: c(2)
    type(nonlinear_constraint) :: bad(1)
    type(multivariate_minimize_result) :: r
    type(multivariate_minimize_result) :: r_direct
    integer :: failures
    real(dp) :: max_weight
    real(dp) :: max_clip_fun
    real(dp) :: max_clip_constr
    real(dp) :: infinity

    failures = 0
    max_weight = 0.0_dp
    infinity = ieee_value(0.0_dp, ieee_positive_inf)

    c(1)%fun => sum_value
    c(1)%jac => sum_jac
    c(1)%hess => zero_hess
    c(1)%lower = [2.0_dp]
    c(1)%upper = [2.0_dp]
    c(2)%fun => first_value
    c(2)%jac => first_jac
    c(2)%hess => zero_hess
    c(2)%lower = [0.5_dp]

    r_direct = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=c, gtol=1.0e-8_dp)
    call check_solution('two analytic blocks', r_direct, objects_analytic_x, &
        objects_analytic_fun, objects_analytic_nit, objects_analytic_nfev)
    call check_values('analytic block values', r_direct%constr, objects_analytic_constr, 1.0e-10_dp)
    call check_counts('analytic nfev', r_direct%constr_nfev, objects_analytic_constr_nfev)
    call check_counts('analytic njev', r_direct%constr_njev, objects_analytic_constr_njev)
    call check_counts('analytic nhev', r_direct%constr_nhev, objects_analytic_constr_nhev)
    call check_values('analytic first Jacobian', r_direct%constr_jac(1, :), &
        [1.0_dp, 1.0_dp], 1.0e-14_dp)
    call check_values('analytic second Jacobian', r_direct%constr_jac(2, :), &
        [1.0_dp, 0.0_dp], 1.0e-14_dp)

    ! A new constraint block alone selects trust-constr in generic minimize.
    r = minimize(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=c, gtol=1.0e-8_dp)
    call check_solution('generic minimize objects', r, objects_analytic_x, &
        objects_analytic_fun, objects_analytic_nit, objects_analytic_nfev)
    call check_counts('generic block counts', r%constr_nfev, objects_analytic_constr_nfev)

    ! SciPy records zero user-Jacobian calls for a numerically differentiated
    ! constraint, while still counting its numerical function probes.
    nullify(c(2)%jac)
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=c, gtol=1.0e-8_dp)
    call check_solution('mixed analytic/FD Jacobians', r, objects_mixed_x, &
        objects_mixed_fun, objects_mixed_nit, objects_mixed_nfev)
    call check_values('mixed block values', r%constr, objects_mixed_constr, 1.0e-10_dp)
    call check_values('mixed second Jacobian', r%constr_jac(2, :), &
        [1.0_dp, 0.0_dp], 2.0e-7_dp)
    call check_counts('mixed nfev', r%constr_nfev, objects_mixed_constr_nfev)
    call check_counts('mixed njev', r%constr_njev, objects_mixed_constr_njev)
    call check_counts('mixed nhev', r%constr_nhev, objects_mixed_constr_nhev)

    ! SciPy's three-point scheme evaluates two probes per variable, without
    ! reporting user-supplied Jacobian evaluations for the numerical block.
    c(2)%fd_scheme = '3-point'
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=c, gtol=1.0e-8_dp)
    call check_solution('mixed analytic/three-point Jacobians', r, objects_three_point_x, &
        objects_three_point_fun, objects_three_point_nit, objects_three_point_nfev)
    call check_values('three-point block values', r%constr, objects_three_point_constr, 1.0e-10_dp)
    call check_values('three-point second Jacobian', r%constr_jac(2, :), &
        [1.0_dp, 0.0_dp], 1.0e-10_dp)
    call check_counts('three-point nfev', r%constr_nfev, objects_three_point_constr_nfev)
    call check_counts('three-point njev', r%constr_njev, objects_three_point_constr_njev)
    call check_counts('three-point nhev', r%constr_nhev, objects_three_point_constr_nhev)
    c(2)%fd_scheme = '2-point'

    ! Variable bounds are reported as a separate, zero-evaluation block.
    c(2)%jac => first_jac
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=c, upper=[1.0_dp, 3.0_dp], &
        bounds_keep_feasible=[.true., .true.])
    call check_status('blocks with bounds', r, .true.)
    if (allocated(r%constr_nfev)) then
        if (size(r%constr_nfev) == 3) then
            call check_counts('bounds block nfev', r%constr_nfev(3:3), [0])
            call check_counts('bounds block njev', r%constr_njev(3:3), [0])
            call check_counts('bounds block nhev', r%constr_nhev(3:3), [0])
        else
            call fail('bounds block count size')
        end if
    else
        call fail('bounds block diagnostics missing')
    end if
    if (allocated(r%constr)) then
        if (size(r%constr) == 4) then
            call check_values('bounds are last', r%constr(3:4), r%x, 1.0e-14_dp)
        else
            call fail('bounds block values size')
        end if
    else
        call fail('bounds block values missing')
    end if

    ! A genuinely nonlinear first constraint exercises weighted Hessian
    ! assembly and multiplier mapping across separately supplied objects.
    c(1)%fun => circle_value
    c(1)%jac => circle_jac
    c(1)%hess => circle_hess
    c(1)%lower = [1.0_dp]
    c(1)%upper = [1.0_dp]
    c(2)%fun => second_value
    c(2)%jac => second_jac
    c(2)%hess => zero_hess
    c(2)%lower = [0.0_dp]
    r = trust_constr(circle_objective, [0.8_dp, 0.6_dp], &
        jac=circle_gradient, hess=quadratic_hess, constraint_objects=c, &
        gtol=1.0e-9_dp, xtol=1.0e-10_dp)
    call check_status('nonlinear weighted Hessian', r, .true.)
    call check_values('nonlinear solution', r%x, &
        [0.89442591_dp, 0.44721616_dp], 4.0e-5_dp)
    if (max_weight < 0.5_dp) call fail('nonlinear multiplier did not reach block Hessian')
    if (allocated(r%constr_nhev)) then
        if (size(r%constr_nhev) == 2) then
            if (any(r%constr_nhev <= 1)) call fail('per-block nonlinear Hessian count')
        else
            call fail('per-block nonlinear Hessian size')
        end if
    else
        call fail('per-block nonlinear Hessian diagnostics missing')
    end if

    ! A user-defined relative step applies separately to each variable. The
    ! three-point derivative of a quadratic constraint is exact to roundoff.
    nullify(c(1)%jac)
    c(1)%fd_scheme = '3-point'
    c(1)%finite_diff_rel_step = [1.0e-4_dp, 4.0e-4_dp]
    r = trust_constr(circle_objective, [0.8_dp, 0.6_dp], &
        jac=circle_gradient, hess=quadratic_hess, constraint_objects=c, &
        gtol=1.0e-9_dp, xtol=1.0e-10_dp)
    call check_status('nonlinear three-point constraint', r, .true.)
    call check_values('nonlinear three-point solution', r%x, &
        [0.89442591_dp, 0.44721616_dp], 4.0e-5_dp)
    call check_values('nonlinear three-point Jacobian', r%constr_jac(1, :), &
        2.0_dp * r%x, 1.0e-9_dp)
    if (allocated(r%constr_nfev)) then
        if (r%constr_nfev(1) <= r%constr_nfev(2)) then
            call fail('three-point nonlinear block did not count numerical probes')
        end if
    else
        call fail('three-point nonlinear diagnostics missing')
    end if

    ! Strict variable bounds must constrain numerical probes even when only
    ! the second of two heterogeneous constraint blocks lacks a Jacobian.
    c(1)%fun => first_value
    c(1)%jac => scalar_identity_jac
    c(1)%lower = [0.0_dp]
    c(1)%upper = [infinity]
    nullify(c(1)%hess)
    c(1)%fd_scheme = '2-point'
    if (allocated(c(1)%finite_diff_rel_step)) deallocate(c(1)%finite_diff_rel_step)
    c(2)%fun => clip_value
    nullify(c(2)%jac, c(2)%hess)
    c(2)%lower = [0.0_dp]
    if (allocated(c(2)%upper)) deallocate(c(2)%upper)
    max_clip_fun = -infinity
    max_clip_constr = -infinity
    r = trust_constr(clip_objective, [1.0_dp - 1.0e-10_dp], &
        constraint_objects=c, upper=[1.0_dp], bounds_keep_feasible=[.true.], &
        maxiter=300)
    call check_status('strict clipping with objects', r, .true.)
    if (max_clip_fun > 1.0_dp + 2.0e-14_dp) call fail('objective probe crossed strict bound')
    if (max_clip_constr > 1.0_dp + 2.0e-14_dp) call fail('constraint probe crossed strict bound')

    ! At the upper endpoint a centered stencil is impossible; a backward
    ! second-order stencil must stay within the retained variable bounds.
    c(2)%fd_scheme = '3-point'
    c(2)%finite_diff_rel_step = [2.0e-4_dp]
    max_clip_fun = -infinity
    max_clip_constr = -infinity
    r = trust_constr(clip_objective, [1.0_dp], &
        constraint_objects=c, upper=[1.0_dp], bounds_keep_feasible=[.true.], &
        maxiter=300)
    call check_status('one-sided three-point strict clipping', r, .true.)
    if (max_clip_fun > 1.0_dp + 2.0e-14_dp) call fail('three-point objective crossed strict bound')
    if (max_clip_constr > 1.0_dp + 2.0e-14_dp) then
        write(*, '(a,es22.14)') 'maximum constraint probe: ', max_clip_constr
        call fail('three-point constraint crossed strict bound')
    end if
    call check_values('one-sided three-point Jacobian', r%constr_jac(2, :), &
        2.0_dp * r%x, 1.0e-8_dp)

    ! Invalid block definitions are rejected before invoking callbacks.
    bad(1)%lower = [0.0_dp]
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=bad)
    call check_status('missing block procedure', r, .false.)
    bad(1)%fun => first_value
    bad(1)%upper = [-1.0_dp]
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=bad)
    call check_status('reversed block bounds', r, .false.)
    bad(1)%upper = [1.0_dp]
    bad(1)%keep_feasible = [.true., .false.]
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=bad)
    call check_status('bad block keep shape', r, .false.)
    bad(1)%keep_feasible = [.true.]
    bad(1)%lower = [1.0_dp]
    bad(1)%upper = [infinity]
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=bad)
    call check_status('retained infeasible block start', r, .false.)
    bad(1)%lower = [0.0_dp]
    bad(1)%keep_feasible = [.false.]
    r = minimize(quadratic, [0.6_dp, 1.4_dp], method='bfgs', constraint_objects=bad)
    call check_status('objects with non trust method', r, .false.)
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=bad, &
        constraints=first_value, meq=1)
    call check_status('objects and legacy constraints conflict', r, .false.)
    ! An empty constraint list is equivalent to an unconstrained problem.
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], jac=quadratic_jac, &
        hess=quadratic_hess, constraint_objects=bad(1:0))
    call check_status('empty constraint object array', r, .true.)
    call check_values('empty constraint object minimizer', r%x, [1.0_dp, 2.5_dp], 1.0e-10_dp)
    if (size(r%constr_nfev) /= 0) call fail('empty constraint object diagnostics')

    bad(1)%fd_scheme = 'central'
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=bad)
    call check_status('invalid finite-difference scheme', r, .false.)
    bad(1)%fd_scheme = '3-point'
    bad(1)%finite_diff_rel_step = [0.001_dp]
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=bad)
    call check_status('wrong relative step shape', r, .false.)
    bad(1)%finite_diff_rel_step = [0.001_dp, 0.0_dp]
    r = trust_constr(quadratic, [0.6_dp, 1.4_dp], constraint_objects=bad)
    ! SciPy accepts zero relative steps and substitutes its automatic scale.
    call check_status('zero relative step fallback', r, .true.)

    if (failures > 0) then
        write(*, '(a,i0)') 'test_optimize_trust_constr_objects: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_objects: PASS'

contains

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: f !! objective value
        f = (x(1) - 1.0_dp)**2 + (x(2) - 2.5_dp)**2
    end function quadratic

    subroutine quadratic_jac(x, g)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: g(:) !! objective gradient
        g = [2.0_dp * (x(1) - 1.0_dp), 2.0_dp * (x(2) - 2.5_dp)]
    end subroutine quadratic_jac

    subroutine quadratic_hess(x, h)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: h(:, :) !! objective Hessian
        if (size(x) < 0) h = 0.0_dp
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 2.0_dp
    end subroutine quadratic_hess

    subroutine sum_value(x, v)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: v(:) !! constraint value
        v(1) = x(1) + x(2)
    end subroutine sum_value

    subroutine sum_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: a(:, :) !! constraint Jacobian
        if (size(x) < 0) a = 0.0_dp
        a(1, :) = [1.0_dp, 1.0_dp]
    end subroutine sum_jac

    subroutine first_value(x, v)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: v(:) !! first-coordinate constraint value
        v(1) = x(1)
    end subroutine first_value

    subroutine first_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: a(:, :) !! first-coordinate constraint Jacobian
        if (size(x) < 0) a = 0.0_dp
        a(1, :) = [1.0_dp, 0.0_dp]
    end subroutine first_jac

    subroutine scalar_identity_jac(x, a)
        real(dp), intent(in) :: x(:) !! scalar parameter
        real(dp), intent(out) :: a(:, :) !! scalar identity Jacobian
        if (size(x) < 0) a = 0.0_dp
        a(1, 1) = 1.0_dp
    end subroutine scalar_identity_jac

    subroutine second_value(x, v)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: v(:) !! second-coordinate constraint value
        v(1) = x(2)
    end subroutine second_value

    subroutine second_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: a(:, :) !! second-coordinate constraint Jacobian
        if (size(x) < 0) a = 0.0_dp
        a(1, :) = [0.0_dp, 1.0_dp]
    end subroutine second_jac

    subroutine zero_hess(x, v, h)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! vector of weights
        real(dp), intent(out) :: h(:, :) !! zero weighted Hessian
        if (size(x) < 0 .or. size(v) < 0) h = 0.0_dp
        h = 0.0_dp
    end subroutine zero_hess

    function circle_objective(x) result(f)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: f !! circle objective value
        f = (x(1) - 2.0_dp)**2 + (x(2) - 1.0_dp)**2
    end function circle_objective

    subroutine circle_gradient(x, g)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: g(:) !! circle objective gradient
        g = [2.0_dp * (x(1) - 2.0_dp), 2.0_dp * (x(2) - 1.0_dp)]
    end subroutine circle_gradient

    subroutine circle_value(x, v)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: v(:) !! squared norm constraint value
        v(1) = sum(x*x)
    end subroutine circle_value

    subroutine circle_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: a(:, :) !! squared norm Jacobian
        a(1, :) = 2.0_dp * x
    end subroutine circle_jac

    subroutine circle_hess(x, v, h)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! squared norm multiplier
        real(dp), intent(out) :: h(:, :) !! weighted Hessian of squared norm
        if (size(x) < 0) h = 0.0_dp
        max_weight = max(max_weight, abs(v(1)))
        h = 0.0_dp
        h(1, 1) = 2.0_dp * v(1)
        h(2, 2) = 2.0_dp * v(1)
    end subroutine circle_hess

    function clip_objective(x) result(f)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: f !! objective tracking max evaluated coordinate
        max_clip_fun = max(max_clip_fun, x(1))
        f = (x(1) - 2.0_dp)**2
    end function clip_objective

    subroutine clip_value(x, v)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: v(:) !! constraint tracking max evaluated coordinate
        max_clip_constr = max(max_clip_constr, x(1))
        v(1) = x(1)**2
    end subroutine clip_value

    subroutine check_solution(name, actual, target_x, target_f, target_nit, target_nfev)
        character(len=*), intent(in) :: name !! identifying test case
        type(multivariate_minimize_result), intent(in) :: actual !! optimizer result
        real(dp), intent(in) :: target_x(:) !! reference minimizer
        real(dp), intent(in) :: target_f !! reference function value
        integer, intent(in) :: target_nit !! reference outer iterations
        integer, intent(in) :: target_nfev !! reference function evaluations

        call check_status(name, actual, .true.)
        call check_values(trim(name)//' x', actual%x, target_x, 1.0e-8_dp)
        if (abs(actual%fun - target_f) > 1.0e-8_dp) call fail(trim(name)//' objective mismatch')
        if (actual%nit /= target_nit) call fail(trim(name)//' outer iteration count')
        if (actual%nfev /= target_nfev) call fail(trim(name)//' objective evaluation count')
    end subroutine check_solution

    subroutine check_status(name, actual, expected_success)
        character(len=*), intent(in) :: name !! identifying test case
        type(multivariate_minimize_result), intent(in) :: actual !! optimizer result
        logical, intent(in) :: expected_success !! expected convergence state

        if (actual%success .neqv. expected_success) then
            write(*, '(a,1x,i0,2a)') trim(name), actual%status, ': ', trim(actual%message)
            call fail(trim(name)//' success flag')
        end if
        if (.not. expected_success .and. actual%status /= -1) then
            call fail(trim(name)//' invalid-input status')
        end if
    end subroutine check_status

    subroutine check_values(name, actual, expected, tol)
        character(len=*), intent(in) :: name !! identifying test case
        real(dp), intent(in) :: actual(:) !! actual vector
        real(dp), intent(in) :: expected(:) !! reference vector
        real(dp), intent(in) :: tol !! absolute error threshold

        if (size(actual) /= size(expected)) then
            call fail(trim(name)//' length mismatch')
        else if (any(abs(actual - expected) > tol)) then
            call fail(trim(name)//' value mismatch')
        end if
    end subroutine check_values

    subroutine check_counts(name, actual, expected)
        character(len=*), intent(in) :: name !! identifying test case
        integer, intent(in) :: actual(:) !! actual callback counts
        integer, intent(in) :: expected(:) !! reference callback counts

        if (size(actual) /= size(expected)) then
            call fail(trim(name)//' length mismatch')
        else if (any(actual /= expected)) then
            call fail(trim(name)//' count mismatch')
        end if
    end subroutine check_counts

    subroutine fail(message)
        character(len=*), intent(in) :: message !! reason for regression failure
        print '(a)', 'FAIL: '//trim(message)
        failures = failures + 1
    end subroutine fail

end program test_optimize_trust_constr_objects
