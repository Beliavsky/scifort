program test_differentiate
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_differentiate, only : derivative, derivative_array_result, derivative_result, &
        differentiate_status_invalid_input, differentiate_status_maxiter, &
        differentiate_status_nonfinite, differentiate_status_success, hessian, hessian_result, &
        jacobian, jacobian_result, complex_step_derivative, complex_step_gradient, &
        complex_step_jacobian
    use test_differentiate_reference
    implicit none

    real(dp), allocatable :: complex_grad(:)
    real(dp), allocatable :: complex_jac(:, :)
    type(derivative_result) :: d
    type(derivative_array_result) :: da
    integer :: failures
    type(hessian_result) :: h
    type(jacobian_result) :: j
    integer :: status

    failures = 0

    d = derivative(exp_function, 1.0_dp)
    call check_derivative_case('exp default', d, exp_default_df_ref, exp_default_error_ref, &
        exp_default_status_ref, exp_default_nit_ref, exp_default_nfev_ref, 2.0e-12_dp, failures)

    d = derivative(exp_function, 1.0_dp, order=4)
    call check_derivative_case('exp order4', d, exp_order4_df_ref, exp_order4_error_ref, &
        exp_order4_status_ref, exp_order4_nit_ref, exp_order4_nfev_ref, 2.0e-12_dp, failures)

    d = derivative(exp_function, 1.0_dp, order=3)
    call check_close('odd order rounded upward', d%df, exp_order4_df_ref, 2.0e-12_dp, failures)
    call check_int('odd order nit', d%nit, exp_order4_nit_ref, failures)

    d = derivative(exp_function, 1.0_dp, step_direction=1)
    call check_derivative_case('exp right', d, exp_right_df_ref, exp_right_error_ref, &
        exp_right_status_ref, exp_right_nit_ref, exp_right_nfev_ref, 3.0e-12_dp, failures)

    d = derivative(exp_function, 1.0_dp, step_direction=-1)
    call check_derivative_case('exp left', d, exp_left_df_ref, exp_left_error_ref, &
        exp_left_status_ref, exp_left_nit_ref, exp_left_nfev_ref, 3.0e-12_dp, failures)

    d = derivative(exp_function, 1.0_dp, maxiter=1)
    call check_close('exp one iteration df', d%df, exp_one_iteration_df_ref, 2.0e-12_dp, failures)
    call check_true('exp one iteration error nan', ieee_is_nan(d%error), failures)
    call check_int('exp one iteration status', d%status, exp_one_iteration_status_ref, failures)
    call check_int('exp one iteration nit', d%nit, exp_one_iteration_nit_ref, failures)
    call check_int('exp one iteration nfev', d%nfev, exp_one_iteration_nfev_ref, failures)
    call check_true('exp one iteration failure', .not. d%success, failures)

    da = derivative(exp_function, exp_array_x_ref)
    call check_array('exp array', da%df, exp_array_df_ref, 2.0e-12_dp, failures)
    call check_true('exp array success', all(da%success), failures)
    call check_true('exp array status', all(da%status == differentiate_status_success), failures)

    j = jacobian(vector_function, [0.2_dp, 0.3_dp, 0.4_dp])
    call check_matrix('vector jacobian', j%df, vector_jacobian_ref, 5.0e-12_dp, failures)
    call check_true('vector jacobian finite status', all(j%status == differentiate_status_success), failures)

    j = jacobian(rosenbrock, [0.5_dp, 0.5_dp, 0.5_dp])
    call check_array('rosen gradient', j%df(1, :), rosen_gradient_ref, 5.0e-11_dp, failures)
    call check_true('rosen gradient success', all(j%success), failures)

    h = hessian(rosenbrock, [0.5_dp, 0.5_dp, 0.5_dp])
    call check_matrix('rosen hessian', h%ddf, rosen_hessian_ref, 1.0e-8_dp, failures)
    call check_int('rosen hessian 11 nfev', h%nfev(1, 1), 121, failures)
    call check_int('rosen hessian 12 nfev', h%nfev(1, 2), 121, failures)
    call check_true('rosen hessian main success', h%success(1, 1) .and. &
        h%success(1, 2) .and. h%success(2, 1) .and. h%success(2, 2), failures)

    call check_close('complex-step derivative', &
        complex_step_derivative(complex_exp_function, 1.0_dp, status=status), &
        exp(1.0_dp), 2.0e-15_dp, failures)
    call check_int('complex-step derivative status', status, differentiate_status_success, failures)

    complex_grad = complex_step_gradient(complex_rosenbrock, &
        [0.5_dp, 0.5_dp, 0.5_dp], status=status)
    call check_array('complex-step gradient', complex_grad, [-51.0_dp, -1.0_dp, 50.0_dp], &
        2.0e-13_dp, failures)
    call check_int('complex-step gradient status', status, differentiate_status_success, failures)

    complex_jac = complex_step_jacobian(complex_vector_function, &
        [0.2_dp, 0.3_dp, 0.4_dp], status=status)
    call check_matrix('complex-step jacobian', complex_jac, vector_jacobian_ref, &
        5.0e-12_dp, failures)
    call check_int('complex-step jacobian status', status, differentiate_status_success, failures)

    d = derivative(nonfinite_function, 1.0_dp)
    call check_int('nonfinite function status', d%status, differentiate_status_nonfinite, failures)
    call check_true('nonfinite function df', ieee_is_nan(d%df), failures)

    d = derivative(exp_function, ieee_value(0.0_dp, ieee_quiet_nan))
    call check_int('nonfinite x status', d%status, differentiate_status_nonfinite, failures)

    d = derivative(exp_function, 1.0_dp, maxiter=0)
    call check_int('invalid maxiter', d%status, differentiate_status_invalid_input, failures)

    d = derivative(exp_function, 1.0_dp, step_factor=1.0_dp)
    call check_int('singular step factor', d%status, differentiate_status_invalid_input, failures)

    call check_int('maxiter status constant', differentiate_status_maxiter, -2, failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_differentiate: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_differentiate: PASS'

contains

    function exp_function(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = exp(x)
    end function exp_function

    function nonfinite_function(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = ieee_value(x, ieee_quiet_nan)
    end function nonfinite_function

    function vector_function(x) result(y)
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        allocate(y(4))
        y = [x(1), 5.0_dp * x(3), 4.0_dp * x(2)**2 - 2.0_dp * x(3), &
            x(3) * sin(x(1))]
    end function vector_function

    function rosenbrock(x) result(y)
        real(dp), intent(in) :: x(:)
        real(dp) :: y
        integer :: i
        y = 0.0_dp
        do i = 1, size(x) - 1
            y = y + 100.0_dp * (x(i + 1) - x(i)**2)**2 + (1.0_dp - x(i))**2
        end do
    end function rosenbrock

    function complex_exp_function(z) result(y)
        complex(dp), intent(in) :: z
        complex(dp) :: y
        y = exp(z)
    end function complex_exp_function

    function complex_rosenbrock(x) result(y)
        complex(dp), intent(in) :: x(:)
        complex(dp) :: y
        integer :: i
        y = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, size(x) - 1
            y = y + 100.0_dp * (x(i + 1) - x(i)**2)**2 + (1.0_dp - x(i))**2
        end do
    end function complex_rosenbrock

    function complex_vector_function(x) result(y)
        complex(dp), intent(in) :: x(:)
        complex(dp), allocatable :: y(:)
        allocate(y(4))
        y = [x(1), 5.0_dp * x(3), 4.0_dp * x(2)**2 - 2.0_dp * x(3), &
            x(3) * sin(x(1))]
    end function complex_vector_function

    subroutine check_derivative_case(name, actual, df_expected, error_expected, &
            status_expected, nit_expected, nfev_expected, tol, nfail)
        character(len=*), intent(in) :: name
        type(derivative_result), intent(in) :: actual
        real(dp), intent(in) :: df_expected
        real(dp), intent(in) :: error_expected
        integer, intent(in) :: status_expected
        integer, intent(in) :: nit_expected
        integer, intent(in) :: nfev_expected
        real(dp), intent(in) :: tol
        integer, intent(inout) :: nfail

        call check_close(trim(name)//' df', actual%df, df_expected, tol, nfail)
        call check_close(trim(name)//' error', actual%error, error_expected, tol, nfail)
        call check_int(trim(name)//' status', actual%status, status_expected, nfail)
        call check_int(trim(name)//' nit', actual%nit, nit_expected, nfail)
        call check_int(trim(name)//' nfev', actual%nfev, nfev_expected, nfail)
        call check_true(trim(name)//' success', actual%success, nfail)
    end subroutine check_derivative_case

    subroutine check_close(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tol
        integer, intent(inout) :: nfail
        if (abs(actual - expected) > tol) then
            write(*, '(a,2(1x,es24.16))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_close

    subroutine check_array(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        real(dp), intent(in) :: tol
        integer, intent(inout) :: nfail
        if (size(actual) /= size(expected)) then
            write(*, '(a)') trim(name)//' size mismatch'
            nfail = nfail + 1
        else if (maxval(abs(actual - expected)) > tol) then
            write(*, '(a,1x,es24.16)') trim(name)//' max error:', maxval(abs(actual - expected))
            nfail = nfail + 1
        end if
    end subroutine check_array

    subroutine check_matrix(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :)
        real(dp), intent(in) :: expected(:, :)
        real(dp), intent(in) :: tol
        integer, intent(inout) :: nfail
        if (any(shape(actual) /= shape(expected))) then
            write(*, '(a)') trim(name)//' shape mismatch'
            nfail = nfail + 1
        else if (maxval(abs(actual - expected)) > tol) then
            write(*, '(a,1x,es24.16)') trim(name)//' max error:', maxval(abs(actual - expected))
            nfail = nfail + 1
        end if
    end subroutine check_matrix

    subroutine check_int(name, actual, expected, nfail)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual
        integer, intent(in) :: expected
        integer, intent(inout) :: nfail
        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_int

    subroutine check_true(name, condition, nfail)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition
        integer, intent(inout) :: nfail
        if (.not. condition) then
            write(*, '(a)') trim(name)//' failed'
            nfail = nfail + 1
        end if
    end subroutine check_true

end program test_differentiate
