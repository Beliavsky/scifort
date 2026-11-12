program test_optimize_root_scalar
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize, only : newton, root_result, root_scalar, &
        root_status_convergence_error, root_status_sign_error, &
        root_status_success, root_status_value_error
    use test_optimize_root_scalar_reference
    implicit none

    integer :: failures
    type(root_result) :: r

    failures = 0

    call check_bracket_method('bisect', cubic_bisect_root_ref, &
        cubic_bisect_iterations_ref, cubic_bisect_calls_ref, failures)
    call check_bracket_method('brentq', cubic_brentq_root_ref, &
        cubic_brentq_iterations_ref, cubic_brentq_calls_ref, failures)
    call check_bracket_method('brenth', cubic_brenth_root_ref, &
        cubic_brenth_iterations_ref, cubic_brenth_calls_ref, failures)
    call check_bracket_method('ridder', cubic_ridder_root_ref, &
        cubic_ridder_iterations_ref, cubic_ridder_calls_ref, failures)
    call check_bracket_method('toms748', cubic_toms748_root_ref, &
        cubic_toms748_iterations_ref, cubic_toms748_calls_ref, failures)

    call check_cos_method('bisect', cos_bisect_root_ref, cos_bisect_iterations_ref, &
        cos_bisect_calls_ref, failures)
    call check_cos_method('brentq', cos_brentq_root_ref, cos_brentq_iterations_ref, &
        cos_brentq_calls_ref, failures)
    call check_cos_method('brenth', cos_brenth_root_ref, cos_brenth_iterations_ref, &
        cos_brenth_calls_ref, failures)
    call check_cos_method('ridder', cos_ridder_root_ref, cos_ridder_iterations_ref, &
        cos_ridder_calls_ref, failures)
    call check_cos_method('toms748', cos_toms748_root_ref, cos_toms748_iterations_ref, &
        cos_toms748_calls_ref, failures)

    r = root_scalar(cubic, method='secant', x0=0.5_dp, x1=2.0_dp)
    call check_result('secant', r, cubic_secant_root_ref, cubic_secant_iterations_ref, &
        cubic_secant_calls_ref, 2.0e-15_dp, failures)

    r = root_scalar(cubic, method='newton', x0=0.5_dp, fprime=cubic_prime)
    call check_result('newton', r, cubic_newton_root_ref, cubic_newton_iterations_ref, &
        cubic_newton_calls_ref, 2.0e-15_dp, failures)

    r = root_scalar(cubic, method='halley', x0=0.5_dp, fprime=cubic_prime, &
        fprime2=cubic_prime2)
    call check_result('halley', r, cubic_halley_root_ref, cubic_halley_iterations_ref, &
        cubic_halley_calls_ref, 2.0e-15_dp, failures)

    r = root_scalar(cubic, method='newton', x0=0.5_dp)
    call check_result('numerical newton', r, cubic_numerical_newton_root_ref, &
        cubic_numerical_newton_iterations_ref, cubic_numerical_newton_calls_ref, &
        2.0e-14_dp, failures)

    r = newton(cubic, 0.5_dp)
    call check_result('direct newton secant', r, direct_secant_root_ref, &
        direct_secant_iterations_ref, direct_secant_calls_ref, 2.0e-14_dp, failures)
    call check_true('direct newton selected secant', trim(r%method) == 'secant', failures)

    r = root_scalar(cubic, bracket=[0.0_dp, 3.0_dp])
    call check_close('automatic bracket method root', r%root, cubic_brentq_root_ref, &
        2.0e-15_dp, failures)
    call check_true('automatic bracket method', trim(r%method) == 'brentq', failures)

    r = root_scalar(cubic, x0=0.5_dp, fprime=cubic_prime, fprime2=cubic_prime2)
    call check_true('automatic halley method', trim(r%method) == 'halley', failures)
    call check_close('automatic halley root', r%root, 1.0_dp, 2.0e-15_dp, failures)

    r = root_scalar(cubic, method='brentq', bracket=[0.0_dp, 3.0_dp], maxiter=1)
    call check_close('brentq maxiter root', r%root, cubic_brentq_maxiter1_root_ref, &
        2.0e-15_dp, failures)
    call check_int('brentq maxiter iterations', r%iterations, &
        cubic_brentq_maxiter1_iterations_ref, failures)
    call check_int('brentq maxiter calls', r%function_calls, &
        cubic_brentq_maxiter1_calls_ref, failures)
    call check_int('brentq maxiter status', r%status, root_status_convergence_error, failures)
    call check_true('brentq maxiter failure', .not. r%converged, failures)

    r = root_scalar(always_positive, method='brentq', bracket=[-1.0_dp, 1.0_dp])
    call check_int('sign error status', r%status, root_status_sign_error, failures)
    call check_true('sign error failure', .not. r%converged, failures)

    r = root_scalar(nonfinite_function, method='brentq', bracket=[-1.0_dp, 1.0_dp])
    call check_int('nonfinite status', r%status, root_status_value_error, failures)

    r = root_scalar(cubic, method='unknown', x0=0.5_dp)
    call check_true('unknown method failure', .not. r%converged, failures)

    r = root_scalar(cubic, method='toms748', bracket=[0.0_dp, 3.0_dp], k=2)
    call check_close('toms748 k2 root', r%root, 1.0_dp, 2.0e-14_dp, failures)
    call check_int('success status constant', root_status_success, 0, failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_root_scalar: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_root_scalar: PASS'

contains

    subroutine check_bracket_method(method, expected_root, expected_iterations, &
            expected_calls, nfail)
        character(len=*), intent(in) :: method
        real(dp), intent(in) :: expected_root
        integer, intent(in) :: expected_iterations
        integer, intent(in) :: expected_calls
        integer, intent(inout) :: nfail
        type(root_result) :: actual

        actual = root_scalar(cubic, method=method, bracket=[0.0_dp, 3.0_dp])
        call check_result(trim(method)//' cubic', actual, expected_root, &
            expected_iterations, expected_calls, 2.0e-12_dp, nfail)
    end subroutine check_bracket_method

    subroutine check_cos_method(method, expected_root, expected_iterations, &
            expected_calls, nfail)
        character(len=*), intent(in) :: method
        real(dp), intent(in) :: expected_root
        integer, intent(in) :: expected_iterations
        integer, intent(in) :: expected_calls
        integer, intent(inout) :: nfail
        type(root_result) :: actual

        actual = root_scalar(cos_minus_x, method=method, bracket=[0.0_dp, 1.0_dp])
        call check_result(trim(method)//' cos', actual, expected_root, &
            expected_iterations, expected_calls, 2.0e-12_dp, nfail)
    end subroutine check_cos_method

    subroutine check_result(name, actual, expected_root, expected_iterations, &
            expected_calls, tol, nfail)
        character(len=*), intent(in) :: name
        type(root_result), intent(in) :: actual
        real(dp), intent(in) :: expected_root
        integer, intent(in) :: expected_iterations
        integer, intent(in) :: expected_calls
        real(dp), intent(in) :: tol
        integer, intent(inout) :: nfail

        call check_close(trim(name)//' root', actual%root, expected_root, tol, nfail)
        call check_int(trim(name)//' iterations', actual%iterations, expected_iterations, nfail)
        call check_int(trim(name)//' calls', actual%function_calls, expected_calls, nfail)
        call check_int(trim(name)//' status', actual%status, root_status_success, nfail)
        call check_true(trim(name)//' converged', actual%converged, nfail)
    end subroutine check_result

    function cubic(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = x**3 - 1.0_dp
    end function cubic

    function cubic_prime(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = 3.0_dp * x*x
    end function cubic_prime

    function cubic_prime2(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = 6.0_dp * x
    end function cubic_prime2

    function cos_minus_x(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = cos(x) - x
    end function cos_minus_x

    function always_positive(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = x*x + 1.0_dp
    end function always_positive

    function nonfinite_function(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = ieee_value(x, ieee_quiet_nan)
    end function nonfinite_function

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

end program test_optimize_root_scalar
