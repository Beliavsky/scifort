program test_optimize_minimize_constrained
    use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize, only : lbfgsb, minimize, multivariate_minimize_result, slsqp, tnc
    use test_optimize_minimize_constrained_reference
    implicit none

    real(dp), parameter :: lower_bounded(2) = [-1.0_dp, -2.0_dp]
    real(dp), parameter :: upper_bounded(2) = [1.0_dp, 0.0_dp]
    real(dp), parameter :: lower_nonnegative(2) = [0.0_dp, 0.0_dp]
    type(multivariate_minimize_result) :: r
    integer :: failures
    real(dp) :: lower_tnc(2)
    real(dp) :: neg_inf

    failures = 0
    neg_inf = ieee_value(0.0_dp, ieee_negative_inf)

    r = lbfgsb(bounded_quadratic, [0.0_dp, 0.0_dp], jac=bounded_quadratic_grad, &
        lower=lower_bounded, upper=upper_bounded)
    call check_result('l-bfgs-b analytic', r, lbfgsb_jac_x_ref, lbfgsb_jac_fun_ref, &
        lbfgsb_jac_nit_ref, lbfgsb_jac_nfev_ref, lbfgsb_jac_njev_ref, &
        lbfgsb_jac_status_ref, 2.0e-8_dp, failures)

    r = lbfgsb(bounded_quadratic, [0.0_dp, 0.0_dp], lower=lower_bounded, &
        upper=upper_bounded)
    call check_result('l-bfgs-b finite difference', r, lbfgsb_fd_x_ref, &
        lbfgsb_fd_fun_ref, lbfgsb_fd_nit_ref, lbfgsb_fd_nfev_ref, &
        lbfgsb_fd_njev_ref, lbfgsb_fd_status_ref, 5.0e-8_dp, failures)

    r = lbfgsb(bounded_quadratic, [0.0_dp, 0.0_dp], jac=bounded_quadratic_grad, &
        lower=lower_bounded, upper=upper_bounded, maxiter=1)
    call check_result('l-bfgs-b maxiter', r, lbfgsb_limit_x_ref, lbfgsb_limit_fun_ref, &
        lbfgsb_limit_nit_ref, lbfgsb_limit_nfev_ref, lbfgsb_limit_njev_ref, &
        lbfgsb_limit_status_ref, 2.0e-8_dp, failures)

    r = minimize(bounded_quadratic, [0.0_dp, 0.0_dp], jac=bounded_quadratic_grad, &
        lower=lower_bounded, upper=upper_bounded)
    call check_result('default bounded minimize', r, lbfgsb_default_x_ref, &
        lbfgsb_default_fun_ref, lbfgsb_default_nit_ref, lbfgsb_default_nfev_ref, &
        lbfgsb_default_njev_ref, lbfgsb_default_status_ref, 2.0e-8_dp, failures)

    r = tnc(bounded_quadratic, [0.0_dp, 0.0_dp], jac=bounded_quadratic_grad, &
        lower=lower_bounded, upper=upper_bounded)
    call check_tnc_result('tnc analytic', r, tnc_jac_x_ref, tnc_jac_fun_ref, &
        2.0e-7_dp, 2.0e-10_dp, failures)

    r = tnc(bounded_quadratic, [0.0_dp, 0.0_dp], lower=lower_bounded, &
        upper=upper_bounded)
    call check_tnc_result('tnc finite difference', r, tnc_fd_x_ref, tnc_fd_fun_ref, &
        2.0e-7_dp, 2.0e-10_dp, failures)

    r = minimize(bounded_quadratic, [0.0_dp, 0.0_dp], method='tnc', &
        jac=bounded_quadratic_grad, lower=lower_bounded, upper=upper_bounded)
    call check_tnc_result('tnc minimize dispatch', r, tnc_jac_x_ref, tnc_jac_fun_ref, &
        2.0e-7_dp, 2.0e-10_dp, failures)

    lower_tnc = [neg_inf, -1.5_dp]
    r = tnc(rosenbrock, [-2.0_dp, 1.0_dp], jac=rosenbrock_grad, lower=lower_tnc, &
        maxfun=200)
    call check_tnc_result('tnc Rosenbrock', r, tnc_rosen_x_ref, tnc_rosen_fun_ref, &
        5.0e-5_dp, 1.0e-8_dp, failures)

    lower_tnc = [neg_inf, 1.5_dp]
    r = tnc(rosenbrock, [-2.0_dp, 1.0_dp], jac=rosenbrock_grad, lower=lower_tnc, &
        maxfun=200)
    call check_tnc_result('tnc active bound', r, tnc_active_x_ref, tnc_active_fun_ref, &
        5.0e-5_dp, 2.0e-7_dp, failures)

    r = tnc(bounded_quadratic, [0.0_dp, 0.0_dp], jac=bounded_quadratic_grad, &
        lower=[1.0_dp, -2.0_dp], upper=[1.0_dp, 0.0_dp])
    call check_tnc_result('tnc mixed fixed bounds', r, tnc_mixed_fixed_x_ref, &
        tnc_mixed_fixed_fun_ref, 2.0e-7_dp, 2.0e-10_dp, failures)

    r = tnc(bounded_quadratic, [0.0_dp, 0.0_dp], jac=bounded_quadratic_grad, &
        lower=lower_bounded, upper=upper_bounded, maxfun=1)
    call check_int('tnc maxfun status', r%status, 3, failures)
    if (r%success) then
        write(*, '(a)') 'tnc maxfun unexpectedly succeeded'
        failures = failures + 1
    end if

    r = tnc(bounded_quadratic, [0.0_dp, 0.0_dp], jac=bounded_quadratic_grad, &
        lower=[1.0_dp, -1.0_dp], upper=[1.0_dp, -1.0_dp])
    call check_int('tnc constant-bounds status', r%status, 5, failures)
    if (r%success) then
        write(*, '(a)') 'tnc constant-bounds unexpectedly succeeded'
        failures = failures + 1
    end if

    r = tnc(bounded_quadratic, [0.0_dp, 0.0_dp], jac=bounded_quadratic_grad, &
        lower=[1.0_dp, 0.0_dp], upper=[0.0_dp, 1.0_dp])
    call check_int('tnc infeasible-bounds status', r%status, -1, failures)
    if (r%success) then
        write(*, '(a)') 'tnc infeasible-bounds unexpectedly succeeded'
        failures = failures + 1
    end if

    r = slsqp(slsqp_eq_fun, [0.5_dp, 1.5_dp], jac=slsqp_eq_grad, &
        lower=lower_nonnegative, constraints=eq_constraint, &
        constraint_jac=eq_constraint_jac, meq=1)
    call check_result('slsqp equality analytic', r, slsqp_eq_jac_x_ref, &
        slsqp_eq_jac_fun_ref, slsqp_eq_jac_nit_ref, slsqp_eq_jac_nfev_ref, &
        slsqp_eq_jac_njev_ref, slsqp_eq_jac_status_ref, 5.0e-8_dp, failures)
    call check_vector('slsqp analytic final gradient', r%jac, &
        2.0_dp * (r%x - [1.0_dp, 2.5_dp]), 1.0e-12_dp, failures)

    r = slsqp(slsqp_eq_fun, [0.5_dp, 1.5_dp], lower=lower_nonnegative, &
        constraints=eq_constraint, meq=1)
    call check_result('slsqp equality finite difference', r, slsqp_eq_fd_x_ref, &
        slsqp_eq_fd_fun_ref, slsqp_eq_fd_nit_ref, slsqp_eq_fd_nfev_ref, &
        slsqp_eq_fd_njev_ref, slsqp_eq_fd_status_ref, 5.0e-7_dp, failures)
    call check_vector('slsqp numerical final gradient', r%jac, &
        2.0_dp * (r%x - [1.0_dp, 2.5_dp]), 5.0e-7_dp, failures)

    r = slsqp(slsqp_eq_fun, [0.5_dp, 1.5_dp], jac=slsqp_eq_grad, &
        lower=lower_nonnegative, constraints=eq_constraint, meq=1)
    call check_result('slsqp numeric constraint jac', r, slsqp_eq_mixed_x_ref, &
        slsqp_eq_mixed_fun_ref, slsqp_eq_mixed_nit_ref, slsqp_eq_mixed_nfev_ref, &
        slsqp_eq_mixed_njev_ref, slsqp_eq_mixed_status_ref, 5.0e-8_dp, failures)

    r = slsqp(slsqp_eq_fun, [0.5_dp, 1.5_dp], jac=slsqp_eq_grad, &
        lower=lower_nonnegative, constraints=eq_constraint, &
        constraint_jac=eq_constraint_jac, meq=1, maxiter=1)
    call check_result('slsqp maxiter', r, slsqp_limit_x_ref, slsqp_limit_fun_ref, &
        slsqp_limit_nit_ref, slsqp_limit_nfev_ref, slsqp_limit_njev_ref, &
        slsqp_limit_status_ref, 5.0e-8_dp, failures)

    r = slsqp(slsqp_eq_fun, [0.5_dp, 1.5_dp], jac=slsqp_eq_grad, &
        lower=lower_nonnegative, constraints=eq_constraint, &
        constraint_jac=eq_constraint_jac, meq=1, maxiter=0)
    call check_int('slsqp zero maxiter status', r%status, 9, failures)
    call check_int('slsqp zero maxiter nit', r%nit, 0, failures)
    call check_vector('slsqp zero maxiter x', r%x, [0.5_dp, 1.5_dp], &
        1.0e-14_dp, failures)

    r = minimize(squared_norm, [0.2_dp, 0.2_dp], method='slsqp', jac=norm2_grad, &
        lower=lower_nonnegative, constraints=ineq_constraint, &
        constraint_jac=ineq_constraint_jac, mineq=1)
    call check_result('slsqp inequality', r, slsqp_ineq_x_ref, slsqp_ineq_fun_ref, &
        slsqp_ineq_nit_ref, slsqp_ineq_nfev_ref, slsqp_ineq_njev_ref, &
        slsqp_ineq_status_ref, 5.0e-8_dp, failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_minimize_constrained: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_minimize_constrained: PASS'

contains

    function bounded_quadratic(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = (x(1) - 2.0_dp)**2 + (x(2) + 1.0_dp)**2
    end function bounded_quadratic

    subroutine bounded_quadratic_grad(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad = [2.0_dp * (x(1) - 2.0_dp), 2.0_dp * (x(2) + 1.0_dp)]
    end subroutine bounded_quadratic_grad

    function rosenbrock(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
    end function rosenbrock

    subroutine rosenbrock_grad(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad(1) = -400.0_dp * x(1) * (x(2) - x(1)**2) - 2.0_dp * (1.0_dp - x(1))
        grad(2) = 200.0_dp * (x(2) - x(1)**2)
    end subroutine rosenbrock_grad

    function slsqp_eq_fun(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = (x(1) - 1.0_dp)**2 + (x(2) - 2.5_dp)**2
    end function slsqp_eq_fun

    subroutine slsqp_eq_grad(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad = [2.0_dp * (x(1) - 1.0_dp), 2.0_dp * (x(2) - 2.5_dp)]
    end subroutine slsqp_eq_grad

    subroutine eq_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! equality constraint values

        c(1) = x(1) + x(2) - 2.0_dp
    end subroutine eq_constraint

    subroutine eq_constraint_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: a(:, :) !! equality constraint Jacobian

        if (size(x) < 0) a = 0.0_dp
        a(1, :) = [1.0_dp, 1.0_dp]
    end subroutine eq_constraint_jac

    function squared_norm(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = dot_product(x, x)
    end function squared_norm

    subroutine norm2_grad(x, grad)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient

        grad = 2.0_dp * x
    end subroutine norm2_grad

    subroutine ineq_constraint(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! inequality values, feasible when nonnegative

        c(1) = x(1) + x(2) - 1.0_dp
    end subroutine ineq_constraint

    subroutine ineq_constraint_jac(x, a)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: a(:, :) !! inequality constraint Jacobian

        if (size(x) < 0) a = 0.0_dp
        a(1, :) = [1.0_dp, 1.0_dp]
    end subroutine ineq_constraint_jac

    subroutine check_tnc_result(name, actual, expected_x, expected_fun, xtol, ftol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(multivariate_minimize_result), intent(in) :: actual !! actual TNC result
        real(dp), intent(in) :: expected_x(:) !! SciPy reference minimizer
        real(dp), intent(in) :: expected_fun !! SciPy reference objective
        real(dp), intent(in) :: xtol !! absolute minimizer tolerance
        real(dp), intent(in) :: ftol !! absolute objective tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' x', actual%x, expected_x, xtol, nfail)
        call check_real(trim(name)//' fun', actual%fun, expected_fun, ftol, nfail)
        if (.not. actual%success) then
            write(*, '(a,1x,i0,2a)') trim(name)//' failed with status', actual%status, &
                ': ', trim(actual%message)
            nfail = nfail + 1
        end if
    end subroutine check_tnc_result

    subroutine check_result(name, actual, expected_x, expected_fun, expected_nit, &
            expected_nfev, expected_njev, expected_status, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(multivariate_minimize_result), intent(in) :: actual !! actual result
        real(dp), intent(in) :: expected_x(:) !! SciPy reference minimizer
        real(dp), intent(in) :: expected_fun !! SciPy reference objective
        integer, intent(in) :: expected_nit !! SciPy reference iteration count
        integer, intent(in) :: expected_nfev !! SciPy reference objective count
        integer, intent(in) :: expected_njev !! SciPy reference gradient count
        integer, intent(in) :: expected_status !! SciPy reference status
        real(dp), intent(in) :: tol !! absolute numerical tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' x', actual%x, expected_x, tol, nfail)
        call check_real(trim(name)//' fun', actual%fun, expected_fun, tol, nfail)
        call check_int(trim(name)//' nit', actual%nit, expected_nit, nfail)
        call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        call check_int(trim(name)//' njev', actual%njev, expected_njev, nfail)
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
    end subroutine check_result

    subroutine check_vector(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! actual vector
        real(dp), intent(in) :: expected(:) !! reference vector
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        if (size(actual) /= size(expected) .or. any(abs(actual - expected) > tol)) then
            write(*, '(a,*(1x,es24.16))') trim(name)//' actual:', actual
            write(*, '(a,*(1x,es24.16))') trim(name)//' expected:', expected
            nfail = nfail + 1
        end if
    end subroutine check_vector

    subroutine check_real(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual !! actual value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        if (abs(actual - expected) > tol) then
            write(*, '(a,2(1x,es24.16))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_real

    subroutine check_int(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual !! actual integer
        integer, intent(in) :: expected !! reference integer
        integer, intent(inout) :: nfail !! accumulated failure count

        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_int

end program test_optimize_minimize_constrained
