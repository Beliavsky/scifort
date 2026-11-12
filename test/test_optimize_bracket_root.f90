! SPDX-License-Identifier: MIT
! Regression for scalar scipy.optimize.elementwise.bracket_root parity.
program test_optimize_bracket_root
    use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_root_scalar, only : bracket_root, bracket_root_result, &
        bracket_root_status_limit, root_result, root_scalar, root_status_convergence_error, &
        root_status_invalid_input, root_status_success, root_status_value_error
    use test_optimize_bracket_root_reference
    implicit none

    integer :: calls
    integer :: failures
    type(bracket_root_result) :: r
    type(root_result) :: solved

    failures = 0
    calls = 0
    r = bracket_root(cubic, 0.0_dp)
    call check_case('cubic', r, [cubic_xl_ref, cubic_xr_ref, cubic_fl_ref, cubic_fr_ref], &
        cubic_nit_ref, cubic_nfev_ref, cubic_status_ref)

    ! A bracket produced by the new helper must be consumable by root_scalar.
    solved = root_scalar(cubic, bracket=[r%xl, r%xr])
    call check_true('bracket feeds root_scalar', solved%converged)
    call check_true('cubic root matches SciPy', &
        abs(solved%root - 2.0945514815423389_dp) < 5.0e-13_dp)

    calls = 0
    r = bracket_root(right, 0.0_dp)
    call check_case('right', r, [right_xl_ref, right_xr_ref, right_fl_ref, right_fr_ref], &
        right_nit_ref, right_nfev_ref, right_status_ref)

    calls = 0
    r = bracket_root(left, 0.0_dp)
    call check_case('left', r, [left_xl_ref, left_xr_ref, left_fl_ref, left_fr_ref], &
        left_nit_ref, left_nfev_ref, left_status_ref)

    calls = 0
    r = bracket_root(right, 0.0_dp, xmax=10.0_dp)
    call check_case('bounded right', r, &
        [bounded_right_xl_ref, bounded_right_xr_ref, bounded_right_fl_ref, bounded_right_fr_ref], &
        bounded_right_nit_ref, bounded_right_nfev_ref, bounded_right_status_ref)

    calls = 0
    r = bracket_root(left, 0.0_dp, xmin=-5.0_dp)
    call check_case('bounded left', r, &
        [bounded_left_xl_ref, bounded_left_xr_ref, bounded_left_fl_ref, bounded_left_fr_ref], &
        bounded_left_nit_ref, bounded_left_nfev_ref, bounded_left_status_ref)

    calls = 0
    r = bracket_root(right, 0.0_dp, factor=3.0_dp)
    call check_case('factor 3', r, [factor3_xl_ref, factor3_xr_ref, factor3_fl_ref, factor3_fr_ref], &
        factor3_nit_ref, factor3_nfev_ref, factor3_status_ref)

    calls = 0
    r = bracket_root(initial, 0.0_dp)
    call check_case('initial bracket', r, &
        [initial_xl_ref, initial_xr_ref, initial_fl_ref, initial_fr_ref], &
        initial_nit_ref, initial_nfev_ref, initial_status_ref)

    calls = 0
    r = bracket_root(zero_left, 0.0_dp)
    call check_case('zero at left', r, &
        [zero_left_xl_ref, zero_left_xr_ref, zero_left_fl_ref, zero_left_fr_ref], &
        zero_left_nit_ref, zero_left_nfev_ref, zero_left_status_ref)

    calls = 0
    r = bracket_root(zero_right, 0.0_dp)
    call check_case('zero at right', r, &
        [zero_right_xl_ref, zero_right_xr_ref, zero_right_fl_ref, zero_right_fr_ref], &
        zero_right_nit_ref, zero_right_nfev_ref, zero_right_status_ref)

    calls = 0
    r = bracket_root(right, 0.0_dp, maxiter=1)
    call check_case('short iteration cap', r, &
        [shortcap_xl_ref, shortcap_xr_ref, shortcap_fl_ref, shortcap_fr_ref], &
        shortcap_nit_ref, shortcap_nfev_ref, shortcap_status_ref)

    calls = 0
    r = bracket_root(right, 0.2_dp, xr0=0.8_dp, xmin=0.0_dp, xmax=1.0_dp, maxiter=4)
    call check_case('bounded iteration cap', r, &
        [boundcap_xl_ref, boundcap_xr_ref, boundcap_fl_ref, boundcap_fr_ref], &
        boundcap_nit_ref, boundcap_nfev_ref, boundcap_status_ref)

    calls = 0
    r = bracket_root(two, 0.0_dp, xmin=0.0_dp)
    call check_case('left already clamped', r, &
        [leftclamped_xl_ref, leftclamped_xr_ref, leftclamped_fl_ref, leftclamped_fr_ref], &
        leftclamped_nit_ref, leftclamped_nfev_ref, leftclamped_status_ref)

    calls = 0
    r = bracket_root(right, 0.0_dp, xmin=0.0_dp)
    call check_case('only right active', r, &
        [onlyright_xl_ref, onlyright_xr_ref, onlyright_fl_ref, onlyright_fr_ref], &
        onlyright_nit_ref, onlyright_nfev_ref, onlyright_status_ref)

    calls = 0
    r = bracket_root(left, 0.0_dp, xmax=1.0_dp)
    call check_case('default right shrinks to bound', r, &
        [rightlimited_default_xl_ref, rightlimited_default_xr_ref, &
        rightlimited_default_fl_ref, rightlimited_default_fr_ref], &
        rightlimited_default_nit_ref, rightlimited_default_nfev_ref, rightlimited_default_status_ref)

    calls = 0
    r = bracket_root(right, 0.0_dp, &
        xmin=ieee_value(0.0_dp, ieee_negative_inf), &
        xmax=ieee_value(0.0_dp, ieee_positive_inf))
    call check_case('explicit infinite bounds', r, &
        [right_xl_ref, right_xr_ref, right_fl_ref, right_fr_ref], &
        right_nit_ref, right_nfev_ref, right_status_ref)

    calls = 0
    r = bracket_root(right, 0.0_dp, maxiter=0)
    call check_int('maxiter zero status', r%status, root_status_convergence_error)
    call check_int('maxiter zero evals', r%nfev, 2)
    call check_int('maxiter zero callback count', calls, 2)

    calls = 0
    r = bracket_root(right, 0.0_dp, xr0=0.0_dp)
    call check_int('invalid reversed bracket', r%status, root_status_invalid_input)
    call check_int('invalid reversed no evaluations', calls, 0)

    calls = 0
    r = bracket_root(right, 0.0_dp, factor=1.0_dp)
    call check_int('invalid growth factor', r%status, root_status_invalid_input)
    call check_int('invalid factor no evaluations', calls, 0)

    calls = 0
    r = bracket_root(right, 0.0_dp, xmax=-0.5_dp)
    call check_int('invalid bound', r%status, root_status_invalid_input)
    call check_int('invalid bound no evaluations', calls, 0)

    calls = 0
    r = bracket_root(right, 0.0_dp, maxiter=-1)
    call check_int('invalid iteration limit', r%status, root_status_invalid_input)
    call check_int('invalid limit no evaluations', calls, 0)

    calls = 0
    r = bracket_root(always_positive, 0.0_dp, xr0=1.0_dp, xmin=0.0_dp, xmax=1.0_dp)
    call check_int('exhausted both limits', r%status, bracket_root_status_limit)
    call check_int('exhausted no expansions', r%nit, 0)
    call check_int('exhausted evaluations', calls, 2)

    calls = 0
    r = bracket_root(not_finite, 0.0_dp)
    call check_int('nonfinite objective', r%status, root_status_value_error)
    call check_int('nonfinite initial evaluations', calls, 2)

    calls = 0
    r = bracket_root(right, 0.0_dp, factor=3.0_dp)
    call check_int('success status constant', r%status, root_status_success)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_bracket_root: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_bracket_root: PASS'

contains

    function right(x) result(y)
        real(dp), intent(in) :: x !! coordinate
        real(dp) :: y
        calls = calls + 1
        y = x - 8.0_dp
    end function right

    function left(x) result(y)
        real(dp), intent(in) :: x !! coordinate
        real(dp) :: y
        calls = calls + 1
        y = x + 3.0_dp
    end function left

    function cubic(x) result(y)
        real(dp), intent(in) :: x !! coordinate
        real(dp) :: y
        calls = calls + 1
        y = x**3 - 2.0_dp*x - 5.0_dp
    end function cubic

    function initial(x) result(y)
        real(dp), intent(in) :: x !! coordinate
        real(dp) :: y
        calls = calls + 1
        y = x - 0.5_dp
    end function initial

    function zero_left(x) result(y)
        real(dp), intent(in) :: x !! coordinate
        real(dp) :: y
        calls = calls + 1
        y = x
    end function zero_left

    function zero_right(x) result(y)
        real(dp), intent(in) :: x !! coordinate
        real(dp) :: y
        calls = calls + 1
        y = x - 1.0_dp
    end function zero_right

    function two(x) result(y)
        real(dp), intent(in) :: x !! coordinate
        real(dp) :: y
        calls = calls + 1
        y = x - 2.0_dp
    end function two

    function always_positive(x) result(y)
        real(dp), intent(in) :: x !! coordinate
        real(dp) :: y
        calls = calls + 1
        y = x*x + 1.0_dp
    end function always_positive

    function not_finite(x) result(y)
        real(dp), intent(in) :: x !! coordinate
        real(dp) :: y
        calls = calls + 1
        y = ieee_value(x, ieee_quiet_nan)
    end function not_finite

    subroutine check_case(label, actual, values, nit, nfev, status)
        character(len=*), intent(in) :: label !! descriptive name
        type(bracket_root_result), intent(in) :: actual !! computed bracket
        real(dp), intent(in) :: values(4) !! SciPy reference endpoints and function values
        integer, intent(in) :: nit !! expected expansion steps
        integer, intent(in) :: nfev !! expected objective evaluations
        integer, intent(in) :: status !! expected SciPy termination status
        real(dp), parameter :: tol = 2.0e-13_dp
        real(dp) :: observed(4)
        observed = [actual%xl, actual%xr, actual%fl, actual%fr]
        if (any(abs(observed - values) > tol * max(1.0_dp, abs(values)))) then
            write(*, '(a,4(1x,es15.7))') label//' observed: ', observed
            write(*, '(a,4(1x,es15.7))') label//' reference:', values
            failures = failures + 1
        end if
        call check_int(label//' status', actual%status, status)
        call check_int(label//' iterations', actual%nit, nit)
        call check_int(label//' evaluations', actual%nfev, nfev)
        call check_int(label//' observed calls', calls, nfev)
        call check_true(label//' success', actual%success .eqv. (status == root_status_success))
        if (actual%success) then
            call check_true(label//' ordered endpoints', actual%xl <= actual%xr)
            call check_true(label//' root enclosed', &
                actual%fl == 0.0_dp .or. actual%fr == 0.0_dp .or. &
                (actual%fl < 0.0_dp .and. actual%fr > 0.0_dp) .or. &
                (actual%fl > 0.0_dp .and. actual%fr < 0.0_dp))
        end if
    end subroutine check_case

    subroutine check_int(label, actual, expected)
        character(len=*), intent(in) :: label !! descriptive name
        integer, intent(in) :: actual !! observed value
        integer, intent(in) :: expected !! reference value
        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') label//' mismatch:', actual, expected
            failures = failures + 1
        end if
    end subroutine check_int

    subroutine check_true(label, condition)
        character(len=*), intent(in) :: label !! descriptive name
        logical, intent(in) :: condition !! assertion
        if (.not. condition) then
            write(*, '(a)') label//': FAIL'
            failures = failures + 1
        end if
    end subroutine check_true
end program test_optimize_bracket_root
