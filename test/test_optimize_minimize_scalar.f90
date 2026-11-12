program test_optimize_minimize_scalar
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize, only : bounded, brent, golden, minimize_result, &
        minimize_scalar, minimize_status_convergence_error, &
        minimize_status_invalid_input, minimize_status_success, &
        minimize_status_value_error
    use test_optimize_minimize_scalar_reference
    implicit none

    integer :: failures
    type(minimize_result) :: r

    failures = 0

    r = minimize_scalar(objective)
    call check_result('default brent', r, brent_default_x_ref, brent_default_fun_ref, &
        brent_default_iterations_ref, brent_default_calls_ref, 5.0e-14_dp, &
        5.0e-13_dp, failures)
    call check_true('default selects brent', trim(r%method) == 'brent', failures)

    r = brent(objective, bracket=[0.0_dp, 1.0_dp])
    call check_result('brent pair', r, brent_pair_x_ref, brent_pair_fun_ref, &
        brent_pair_iterations_ref, brent_pair_calls_ref, 5.0e-14_dp, &
        5.0e-13_dp, failures)

    r = minimize_scalar(objective, method='BrEnT', bracket=[0.0_dp, 1.0_dp, 3.0_dp])
    call check_result('brent triple', r, brent_triple_x_ref, brent_triple_fun_ref, &
        brent_triple_iterations_ref, brent_triple_calls_ref, 5.0e-14_dp, &
        5.0e-13_dp, failures)

    r = golden(objective, bracket=[0.0_dp, 1.0_dp])
    call check_result('golden pair', r, golden_pair_x_ref, golden_pair_fun_ref, &
        golden_pair_iterations_ref, golden_pair_calls_ref, 5.0e-14_dp, &
        5.0e-13_dp, failures)

    r = minimize_scalar(objective, method='golden', bracket=[0.0_dp, 1.0_dp, 3.0_dp])
    call check_result('golden triple', r, golden_triple_x_ref, golden_triple_fun_ref, &
        golden_triple_iterations_ref, golden_triple_calls_ref, 5.0e-14_dp, &
        5.0e-13_dp, failures)

    r = bounded(objective, -3.0_dp, -1.0_dp)
    call check_result('bounded default', r, bounded_default_x_ref, bounded_default_fun_ref, &
        bounded_default_iterations_ref, bounded_default_calls_ref, 5.0e-12_dp, &
        5.0e-12_dp, failures)

    r = minimize_scalar(objective, bounds=[-3.0_dp, -1.0_dp], xatol=1.0e-10_dp)
    call check_result('bounded tight', r, bounded_tight_x_ref, bounded_tight_fun_ref, &
        bounded_tight_iterations_ref, bounded_tight_calls_ref, 5.0e-10_dp, &
        5.0e-12_dp, failures)
    call check_true('bounds select bounded', trim(r%method) == 'bounded', failures)

    r = minimize_scalar(objective, method='brent', bracket=[0.0_dp, 1.0_dp], maxiter=1)
    call check_failure_result('brent maxiter', r, brent_maxiter1_x_ref, &
        brent_maxiter1_fun_ref, brent_maxiter1_iterations_ref, brent_maxiter1_calls_ref, &
        minimize_status_convergence_error, failures)

    r = minimize_scalar(objective, method='golden', bracket=[0.0_dp, 1.0_dp], maxiter=1)
    call check_failure_result('golden maxiter', r, golden_maxiter1_x_ref, &
        golden_maxiter1_fun_ref, golden_maxiter1_iterations_ref, &
        golden_maxiter1_calls_ref, minimize_status_convergence_error, failures)

    r = minimize_scalar(objective, method='bounded', bounds=[-3.0_dp, -1.0_dp], &
        maxiter=2)
    call check_failure_result('bounded maxiter', r, bounded_maxiter2_x_ref, &
        bounded_maxiter2_fun_ref, bounded_maxiter2_iterations_ref, &
        bounded_maxiter2_calls_ref, minimize_status_convergence_error, failures)

    r = minimize_scalar(nan_objective, method='brent', bracket=[0.0_dp, 1.0_dp])
    call check_int('nonfinite calls', r%function_calls, 3, failures)
    call check_int('nonfinite status', r%status, minimize_status_value_error, failures)
    call check_true('nonfinite x', ieee_is_nan(r%x), failures)
    call check_true('nonfinite fun', ieee_is_nan(r%fun), failures)
    call check_true('nonfinite failure', .not. r%success, failures)

    r = minimize_scalar(objective, method='bounded', bounds=[1.0_dp, -1.0_dp])
    call check_int('invalid bounds status', r%status, minimize_status_invalid_input, failures)
    call check_true('invalid bounds failure', .not. r%success, failures)

    r = minimize_scalar(objective, method='unknown')
    call check_int('unknown method status', r%status, minimize_status_invalid_input, failures)
    call check_true('unknown method failure', .not. r%success, failures)

    call check_int('success status constant', minimize_status_success, 0, failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_minimize_scalar: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_minimize_scalar: PASS'

contains

    function objective(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y

        y = (x - 2.0_dp) * x * (x + 2.0_dp)**2
    end function objective

    function nan_objective(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y

        y = ieee_value(x, ieee_quiet_nan)
    end function nan_objective

    subroutine check_result(name, actual, expected_x, expected_fun, expected_iterations, &
            expected_calls, x_tol, fun_tol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(minimize_result), intent(in) :: actual !! actual optimization result
        real(dp), intent(in) :: expected_x !! reference minimizer
        real(dp), intent(in) :: expected_fun !! reference objective value
        integer, intent(in) :: expected_iterations !! reference iteration count
        integer, intent(in) :: expected_calls !! reference function-call count
        real(dp), intent(in) :: x_tol !! absolute minimizer tolerance
        real(dp), intent(in) :: fun_tol !! absolute objective tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_close(trim(name)//' x', actual%x, expected_x, x_tol, nfail)
        call check_close(trim(name)//' fun', actual%fun, expected_fun, fun_tol, nfail)
        call check_int(trim(name)//' iterations', actual%iterations, expected_iterations, nfail)
        call check_int(trim(name)//' calls', actual%function_calls, expected_calls, nfail)
        call check_int(trim(name)//' status', actual%status, minimize_status_success, nfail)
        call check_true(trim(name)//' success', actual%success, nfail)
    end subroutine check_result

    subroutine check_failure_result(name, actual, expected_x, expected_fun, &
            expected_iterations, expected_calls, expected_status, nfail)
        character(len=*), intent(in) :: name !! check label
        type(minimize_result), intent(in) :: actual !! actual optimization result
        real(dp), intent(in) :: expected_x !! reference minimizer at termination
        real(dp), intent(in) :: expected_fun !! reference objective at termination
        integer, intent(in) :: expected_iterations !! reference iteration count
        integer, intent(in) :: expected_calls !! reference function-call count
        integer, intent(in) :: expected_status !! expected SciFort failure status
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_close(trim(name)//' x', actual%x, expected_x, 5.0e-13_dp, nfail)
        call check_close(trim(name)//' fun', actual%fun, expected_fun, 5.0e-13_dp, nfail)
        call check_int(trim(name)//' iterations', actual%iterations, expected_iterations, nfail)
        call check_int(trim(name)//' calls', actual%function_calls, expected_calls, nfail)
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
        call check_true(trim(name)//' failure', .not. actual%success, nfail)
    end subroutine check_failure_result

    subroutine check_close(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual !! actual value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        if (abs(actual - expected) > tol) then
            write(*, '(a,2(1x,es24.16))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_close

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

    subroutine check_true(name, condition, nfail)
        character(len=*), intent(in) :: name !! check label
        logical, intent(in) :: condition !! condition that must hold
        integer, intent(inout) :: nfail !! accumulated failure count

        if (.not. condition) then
            write(*, '(a)') trim(name)//' failed'
            nfail = nfail + 1
        end if
    end subroutine check_true

end program test_optimize_minimize_scalar
