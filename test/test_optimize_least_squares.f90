program test_optimize_least_squares
    use scifort_kinds, only : dp
    use scifort_optimize, only : curve_fit, curve_fit_result, least_squares, &
        least_squares_result, leastsq, leastsq_result
    use test_optimize_least_squares_reference
    implicit none

    real(dp), parameter :: xdata(4) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: ydata(4) = [1.1_dp, 2.9_dp, 5.2_dp, 6.8_dp]
    real(dp), parameter :: sigma(4) = [1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp]
    type(least_squares_result) :: ls
    type(leastsq_result) :: lq
    type(curve_fit_result) :: cf
    integer :: failures

    failures = 0

    ls = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, method='lm')
    call check_ls('least_squares finite difference', ls, ls_fd_x_ref, ls_fd_cost_ref, &
        ls_fd_fun_ref, ls_fd_jac_ref, ls_fd_grad_ref, ls_fd_optimality_ref, &
        ls_fd_nfev_ref, ls_fd_njev_ref, ls_fd_status_ref, 5.0e-7_dp, failures)

    ls = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='lm', diff_scheme='3-point')
    call check_ls('least_squares three-point', ls, ls_3_x_ref, ls_3_cost_ref, &
        ls_3_fun_ref, ls_3_jac_ref, ls_3_grad_ref, ls_3_optimality_ref, &
        ls_3_nfev_ref, ls_3_njev_ref, ls_3_status_ref, 5.0e-9_dp, failures)

    ls = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, &
        method='lm', diff_scheme='cs', complex_fun=rosen_residual_complex)
    call check_ls('least_squares complex-step', ls, ls_cs_x_ref, ls_cs_cost_ref, &
        ls_cs_fun_ref, ls_cs_jac_ref, ls_cs_grad_ref, ls_cs_optimality_ref, &
        ls_cs_nfev_ref, ls_cs_njev_ref, ls_cs_status_ref, 5.0e-12_dp, failures)

    ls = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, jac=rosen_jac, method='lm')
    call check_ls('least_squares analytic', ls, ls_jac_x_ref, ls_jac_cost_ref, &
        ls_jac_fun_ref, ls_jac_jac_ref, ls_jac_grad_ref, ls_jac_optimality_ref, &
        ls_jac_nfev_ref, ls_jac_njev_ref, ls_jac_status_ref, 5.0e-10_dp, failures)

    ls = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, jac=rosen_jac, &
        method='lm', max_nfev=2)
    call check_ls('least_squares max_nfev', ls, ls_limit_x_ref, ls_limit_cost_ref, &
        ls_limit_fun_ref, ls_limit_jac_ref, ls_limit_grad_ref, ls_limit_optimality_ref, &
        ls_limit_nfev_ref, ls_limit_njev_ref, ls_limit_status_ref, 5.0e-10_dp, failures)

    ls = least_squares(rosen_residual, [-2.0_dp, 1.0_dp], 2, jac=rosen_jac, &
        method='lm', x_scale=[2.0_dp, 0.5_dp])
    call check_ls('least_squares x_scale', ls, ls_scale_x_ref, ls_scale_cost_ref, &
        ls_scale_fun_ref, ls_scale_jac_ref, ls_scale_grad_ref, ls_scale_optimality_ref, &
        ls_scale_nfev_ref, ls_scale_njev_ref, ls_scale_status_ref, 5.0e-10_dp, failures)

    lq = leastsq(rosen_residual, [-2.0_dp, 1.0_dp], 2)
    call check_lq('leastsq finite difference', lq, lq_fd_x_ref, lq_fd_cov_ref, &
        lq_fd_fvec_ref, lq_fd_nfev_ref, lq_fd_njev_ref, lq_fd_status_ref, &
        5.0e-7_dp, failures)

    lq = leastsq(rosen_residual, [-2.0_dp, 1.0_dp], 2, jac=rosen_jac)
    call check_lq('leastsq analytic', lq, lq_jac_x_ref, lq_jac_cov_ref, &
        lq_jac_fvec_ref, lq_jac_nfev_ref, lq_jac_njev_ref, lq_jac_status_ref, &
        5.0e-10_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [2.0_dp, 1.0_dp])
    call check_curve('curve_fit unweighted', cf, cf_plain_popt_ref, cf_plain_pcov_ref, &
        cf_plain_fvec_ref, cf_plain_nfev_ref, cf_plain_status_ref, 5.0e-7_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [2.0_dp, 1.0_dp], sigma=sigma)
    call check_curve('curve_fit relative sigma', cf, cf_rel_popt_ref, cf_rel_pcov_ref, &
        cf_rel_fvec_ref, cf_rel_nfev_ref, cf_rel_status_ref, 5.0e-7_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [2.0_dp, 1.0_dp], sigma=sigma, &
        absolute_sigma=.true.)
    call check_curve('curve_fit absolute sigma', cf, cf_abs_popt_ref, cf_abs_pcov_ref, &
        cf_abs_fvec_ref, cf_abs_nfev_ref, cf_abs_status_ref, 5.0e-7_dp, failures)

    cf = curve_fit(linear_model, xdata, ydata, [2.0_dp, 1.0_dp], jac=linear_model_jac, &
        sigma=sigma, absolute_sigma=.true.)
    call check_curve('curve_fit analytic jac', cf, cf_jac_popt_ref, cf_jac_pcov_ref, &
        cf_jac_fvec_ref, cf_jac_nfev_ref, cf_jac_status_ref, 5.0e-10_dp, failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_least_squares: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_least_squares: PASS'

contains

    subroutine rosen_residual(x, fvec)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: fvec(:) !! Rosenbrock residuals

        fvec = [10.0_dp * (x(2) - x(1)**2), 1.0_dp - x(1)]
    end subroutine rosen_residual

    subroutine rosen_residual_complex(x, fvec)
        complex(dp), intent(in) :: x(:) !! complex-step parameter vector
        complex(dp), intent(out) :: fvec(:) !! analytic Rosenbrock residual continuation

        fvec = [10.0_dp * (x(2) - x(1)**2), 1.0_dp - x(1)]
    end subroutine rosen_residual_complex

    subroutine rosen_jac(x, fjac)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: fjac(:, :) !! Rosenbrock residual Jacobian

        fjac(1, :) = [-20.0_dp * x(1), 10.0_dp]
        fjac(2, :) = [-1.0_dp, 0.0_dp]
    end subroutine rosen_jac

    subroutine linear_model(x, params, ymodel)
        real(dp), intent(in) :: x(:) !! independent data
        real(dp), intent(in) :: params(:) !! slope and intercept
        real(dp), intent(out) :: ymodel(:) !! model values

        ymodel = params(1) * x + params(2)
    end subroutine linear_model

    subroutine linear_model_jac(x, params, fjac)
        real(dp), intent(in) :: x(:) !! independent data
        real(dp), intent(in) :: params(:) !! slope and intercept
        real(dp), intent(out) :: fjac(:, :) !! model Jacobian

        if (size(params) < 0) fjac = 0.0_dp
        fjac(:, 1) = x
        fjac(:, 2) = 1.0_dp
    end subroutine linear_model_jac

    subroutine check_ls(name, actual, expected_x, expected_cost, expected_fun, &
            expected_jac, expected_grad, expected_optimality, expected_nfev, &
            expected_njev, expected_status, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(least_squares_result), intent(in) :: actual !! actual result
        real(dp), intent(in) :: expected_x(:), expected_cost, expected_fun(:)
        real(dp), intent(in) :: expected_jac(:, :), expected_grad(:), expected_optimality
        integer, intent(in) :: expected_nfev, expected_njev, expected_status
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' x', actual%x, expected_x, tol, nfail)
        call check_real(trim(name)//' cost', actual%cost, expected_cost, tol, nfail)
        call check_vector(trim(name)//' fun', actual%fun, expected_fun, tol, nfail)
        call check_matrix(trim(name)//' jac', actual%jac, expected_jac, tol, nfail)
        call check_vector(trim(name)//' grad', actual%grad, expected_grad, tol, nfail)
        call check_real(trim(name)//' optimality', actual%optimality, expected_optimality, &
            tol, nfail)
        call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        call check_int(trim(name)//' njev', actual%njev, expected_njev, nfail)
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
    end subroutine check_ls

    subroutine check_lq(name, actual, expected_x, expected_cov, expected_fvec, &
            expected_nfev, expected_njev, expected_status, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(leastsq_result), intent(in) :: actual !! actual result
        real(dp), intent(in) :: expected_x(:), expected_cov(:, :), expected_fvec(:)
        integer, intent(in) :: expected_nfev, expected_njev, expected_status
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' x', actual%x, expected_x, tol, nfail)
        if (.not. allocated(actual%cov_x)) then
            write(*, '(a)') trim(name)//' covariance missing'
            nfail = nfail + 1
        else
            call check_matrix(trim(name)//' cov', actual%cov_x, expected_cov, tol, nfail)
        end if
        call check_vector(trim(name)//' fvec', actual%fvec, expected_fvec, tol, nfail)
        call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        call check_int(trim(name)//' njev', actual%njev, expected_njev, nfail)
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
    end subroutine check_lq

    subroutine check_curve(name, actual, expected_popt, expected_pcov, expected_fvec, &
            expected_nfev, expected_status, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(curve_fit_result), intent(in) :: actual !! actual result
        real(dp), intent(in) :: expected_popt(:), expected_pcov(:, :), expected_fvec(:)
        integer, intent(in) :: expected_nfev, expected_status
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' popt', actual%popt, expected_popt, tol, nfail)
        call check_matrix(trim(name)//' pcov', actual%pcov, expected_pcov, tol, nfail)
        call check_vector(trim(name)//' fvec', actual%fvec, expected_fvec, tol, nfail)
        call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
    end subroutine check_curve

    subroutine check_vector(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:), expected(:), tol
        integer, intent(inout) :: nfail !! accumulated failure count

        if (size(actual) /= size(expected) .or. any(abs(actual - expected) > tol)) then
            write(*, '(a,*(1x,es24.16))') trim(name)//' actual:', actual
            write(*, '(a,*(1x,es24.16))') trim(name)//' expected:', expected
            nfail = nfail + 1
        end if
    end subroutine check_vector

    subroutine check_matrix(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:, :), expected(:, :), tol
        integer, intent(inout) :: nfail !! accumulated failure count

        if (any(shape(actual) /= shape(expected)) .or. any(abs(actual - expected) > tol)) then
            write(*, '(a)') trim(name)//' matrix mismatch'
            nfail = nfail + 1
        end if
    end subroutine check_matrix

    subroutine check_real(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual, expected, tol
        integer, intent(inout) :: nfail !! accumulated failure count

        if (abs(actual - expected) > tol) then
            write(*, '(a,2(1x,es24.16))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_real

    subroutine check_int(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual, expected
        integer, intent(inout) :: nfail !! accumulated failure count

        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_int

end program test_optimize_least_squares
