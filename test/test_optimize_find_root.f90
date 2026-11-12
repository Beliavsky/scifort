! SPDX-License-Identifier: MIT
! Regression tests for scalar scipy.optimize.elementwise.find_root.
program test_optimize_find_root
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize, only : bracket_root, bracket_root_result, find_root, &
        find_root_result, find_root_status_callback, root_status_invalid_input, &
        root_status_sign_error, root_status_convergence_error, &
        root_status_success, root_status_value_error
    use test_optimize_find_root_reference
    implicit none

    type(find_root_result) :: r
    type(bracket_root_result) :: br
    integer :: failures, calls, callback_calls, last_nit
    logical :: callbacks_in_progress

    failures = 0
    calls = 0
    callback_calls = 0
    last_nit = -1
    callbacks_in_progress = .true.

    r = find_root(linear, [0.0_dp, 4.0_dp])
    call check_case('linear', r, linear_x_ref, linear_f_x_ref, linear_xl_ref, linear_xr_ref, &
        linear_fl_ref, linear_fr_ref, linear_nit_ref, linear_nfev_ref, linear_status_ref)
    r = find_root(cubic, [2.0_dp, 4.0_dp])
    call check_case('cubic', r, cubic_x_ref, cubic_f_x_ref, cubic_xl_ref, cubic_xr_ref, &
        cubic_fl_ref, cubic_fr_ref, cubic_nit_ref, cubic_nfev_ref, cubic_status_ref)
    r = find_root(cosine, [0.0_dp, 1.0_dp])
    call check_case('cosine', r, cosine_x_ref, cosine_f_x_ref, cosine_xl_ref, cosine_xr_ref, &
        cosine_fl_ref, cosine_fr_ref, cosine_nit_ref, cosine_nfev_ref, cosine_status_ref)
    r = find_root(exponential, [0.0_dp, 2.0_dp])
    call check_case('exponential', r, exponential_x_ref, exponential_f_x_ref, exponential_xl_ref, exponential_xr_ref, &
        exponential_fl_ref, exponential_fr_ref, exponential_nit_ref, exponential_nfev_ref, exponential_status_ref)
    r = find_root(cubed, [-1.0_dp, 1.0_dp])
    call check_case('cubed', r, cubed_x_ref, cubed_f_x_ref, cubed_xl_ref, cubed_xr_ref, &
        cubed_fl_ref, cubed_fr_ref, cubed_nit_ref, cubed_nfev_ref, cubed_status_ref)
    r = find_root(cubic, [4.0_dp, 2.0_dp])
    call check_case('reversed', r, reversed_x_ref, reversed_f_x_ref, reversed_xl_ref, reversed_xr_ref, &
        reversed_fl_ref, reversed_fr_ref, reversed_nit_ref, reversed_nfev_ref, reversed_status_ref)
    r = find_root(left_zero, [0.0_dp, 2.0_dp])
    call check_case('left_zero', r, left_zero_x_ref, left_zero_f_x_ref, left_zero_xl_ref, left_zero_xr_ref, &
        left_zero_fl_ref, left_zero_fr_ref, left_zero_nit_ref, left_zero_nfev_ref, left_zero_status_ref)
    r = find_root(right_zero, [0.0_dp, 2.0_dp])
    call check_case('right_zero', r, right_zero_x_ref, right_zero_f_x_ref, right_zero_xl_ref, right_zero_xr_ref, &
        right_zero_fl_ref, right_zero_fr_ref, right_zero_nit_ref, right_zero_nfev_ref, right_zero_status_ref)
    r = find_root(invalid, [-1.0_dp, 1.0_dp])
    call check_case('invalid', r, invalid_x_ref, invalid_f_x_ref, invalid_xl_ref, invalid_xr_ref, &
        invalid_fl_ref, invalid_fr_ref, invalid_nit_ref, invalid_nfev_ref, invalid_status_ref)
    call check_true('invalid root NaN', ieee_is_nan(r%x) .and. ieee_is_nan(r%f_x))
    r = find_root(cosine, [0.0_dp, 1.0_dp], maxiter=0)
    call check_case('max0', r, max0_x_ref, max0_f_x_ref, max0_xl_ref, max0_xr_ref, &
        max0_fl_ref, max0_fr_ref, max0_nit_ref, max0_nfev_ref, max0_status_ref)
    r = find_root(cosine, [0.0_dp, 1.0_dp], maxiter=1)
    call check_case('max1', r, max1_x_ref, max1_f_x_ref, max1_xl_ref, max1_xr_ref, &
        max1_fl_ref, max1_fr_ref, max1_nit_ref, max1_nfev_ref, max1_status_ref)
    r = find_root(cosine, [0.0_dp, 1.0_dp], maxiter=2)
    call check_case('max2', r, max2_x_ref, max2_f_x_ref, max2_xl_ref, max2_xr_ref, &
        max2_fl_ref, max2_fr_ref, max2_nit_ref, max2_nfev_ref, max2_status_ref)
    r = find_root(cosine, [0.0_dp, 1.0_dp], xatol=0.2_dp)
    call check_case('xatol', r, xatol_x_ref, xatol_f_x_ref, xatol_xl_ref, xatol_xr_ref, &
        xatol_fl_ref, xatol_fr_ref, xatol_nit_ref, xatol_nfev_ref, xatol_status_ref)
    r = find_root(cosine, [0.0_dp, 1.0_dp], fatol=0.2_dp)
    call check_case('fatol', r, fatol_x_ref, fatol_f_x_ref, fatol_xl_ref, fatol_xr_ref, &
        fatol_fl_ref, fatol_fr_ref, fatol_nit_ref, fatol_nfev_ref, fatol_status_ref)
    r = find_root(cosine, [0.0_dp, 1.0_dp], frtol=0.2_dp)
    call check_case('frtol', r, frtol_x_ref, frtol_f_x_ref, frtol_xl_ref, frtol_xr_ref, &
        frtol_fl_ref, frtol_fr_ref, frtol_nit_ref, frtol_nfev_ref, frtol_status_ref)
    r = find_root(cosine, [0.0_dp, 1.0_dp], xatol=0.0_dp, xrtol=0.0_dp, &
        fatol=0.0_dp, frtol=0.0_dp, maxiter=100)
    call check_case('strict', r, strict_x_ref, strict_f_x_ref, strict_xl_ref, strict_xr_ref, &
        strict_fl_ref, strict_fr_ref, strict_nit_ref, strict_nfev_ref, strict_status_ref)
    r = find_root(negative, [-4.0_dp, -1.0_dp])
    call check_case('negative', r, negative_x_ref, negative_f_x_ref, negative_xl_ref, negative_xr_ref, &
        negative_fl_ref, negative_fr_ref, negative_nit_ref, negative_nfev_ref, negative_status_ref)
    r = find_root(flat, [0.0_dp, 1.0_dp])
    call check_case('flat', r, flat_x_ref, flat_f_x_ref, flat_xl_ref, flat_xr_ref, &
        flat_fl_ref, flat_fr_ref, flat_nit_ref, flat_nfev_ref, flat_status_ref)

    calls = 0
    r = find_root(linear, [0.0_dp])
    call check_int('malformed shape', r%status, root_status_invalid_input)
    call check_int('shape never calls function', calls, 0)
    r = find_root(linear, [0.0_dp, 4.0_dp], maxiter=-1)
    call check_int('negative iterations', r%status, root_status_invalid_input)
    r = find_root(linear, [0.0_dp, 4.0_dp], xrtol=-1.0_dp)
    call check_int('negative tolerance', r%status, root_status_invalid_input)
    r = find_root(linear, [0.0_dp, 4.0_dp], fatol=ieee_value(0.0_dp, ieee_quiet_nan))
    call check_int('nan tolerance', r%status, root_status_invalid_input)
    call check_int('all invalid inputs avoid evaluation', calls, 0)

    calls = 0
    r = find_root(nonfinite, [0.0_dp, 1.0_dp])
    call check_int('nonfinite status', r%status, root_status_value_error)
    call check_int('nonfinite calls', r%nfev, 2)
    call check_true('nonfinite estimate', ieee_is_nan(r%x))
    calls = 0
    r = find_root(late_nonfinite, [0.0_dp, 1.0_dp])
    call check_int('late nonfinite status', r%status, root_status_value_error)
    call check_int('late nonfinite nit', r%nit, 1)
    call check_int('late nonfinite nfev', r%nfev, 3)
    call check_true('late nonfinite estimate', ieee_is_nan(r%x))

    callback_calls = 0
    last_nit = -1
    callbacks_in_progress = .true.
    r = find_root(cosine, [0.0_dp, 1.0_dp], callback=stop_after_two)
    call check_int('callback status', r%status, find_root_status_callback)
    call check_int('callback nit', r%nit, 2)
    call check_int('callback work', r%nfev, 4)
    call check_int('callback calls', callback_calls, 3)
    call check_int('callback last nit', last_nit, 2)
    call check_true('callback sees progress', callbacks_in_progress)
    callback_calls = 0
    r = find_root(cosine, [0.0_dp, 1.0_dp], callback=stop_initial)
    call check_int('initial callback status', r%status, find_root_status_callback)
    call check_int('initial callback iterations', r%nit, 0)
    call check_int('initial callback work', r%nfev, 2)
    call check_int('initial callback count', callback_calls, 1)
    callback_calls = 0
    r = find_root(cosine, [0.0_dp, 1.0_dp], callback=never_stop)
    call check_true('nonstopping callback success', r%success)
    call check_int('nonstopping callback count', callback_calls, r%nit+1)

    br = bracket_root(cubic, xl0=0.0_dp)
    call check_true('upstream bracket found', br%success)
    if (br%success) then
        r = find_root(cubic, [br%xl, br%xr])
        call check_true('bracket feeds find_root', r%success .and. abs(r%x-2.09455148154233_dp)<1.e-10_dp)
    end if

    call check_int('root status success', root_status_success, 0)
    call check_int('root status sign', root_status_sign_error, -1)
    call check_int('root status convergence', root_status_convergence_error, -2)
    if (failures > 0) then
        write(*, '(a,i0)') 'test_optimize_find_root: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_find_root: PASS'

contains
    subroutine check_case(name, sol, x, f_x, xl, xr, fl, fr, nit, nfev, status)
        character(len=*), intent(in) :: name !! comparison name
        type(find_root_result), intent(in) :: sol !! Fortran result
        real(dp), intent(in) :: x !! SciPy root estimate
        real(dp), intent(in) :: f_x !! SciPy value at estimate
        real(dp), intent(in) :: xl !! SciPy left endpoint
        real(dp), intent(in) :: xr !! SciPy right endpoint
        real(dp), intent(in) :: fl !! SciPy left endpoint function value
        real(dp), intent(in) :: fr !! SciPy right endpoint function value
        integer, intent(in) :: nit !! expected iterations
        integer, intent(in) :: nfev !! expected function evaluations
        integer, intent(in) :: status !! expected status

        call check_int(name//' status', sol%status, status)
        call check_int(name//' nit', sol%nit, nit)
        call check_int(name//' nfev', sol%nfev, nfev)
        call check_true(name//' success consistency', sol%success .eqv. (sol%status == 0))
        if (max(abs(sol%xl-xl), abs(sol%xr-xr)) > 1.e-7_dp) then
            write(*,'(a,4es23.13)') 'FAIL '//name//' bracket', sol%xl, sol%xr, xl, xr
            failures = failures + 1
        end if
        if (max(abs(sol%fl-fl), abs(sol%fr-fr)) > 1.e-7_dp) then
            write(*,'(a,4es23.13)') 'FAIL '//name//' fbracket', sol%fl, sol%fr, fl, fr
            failures = failures + 1
        end if
        if (status == 0) then
            call check_true(name//' root in bracket', sol%x >= sol%xl .and. sol%x <= sol%xr)
        end if
        if (status /= root_status_sign_error .and. status /= root_status_value_error) then
            call check_true(name//' estimate matches SciPy', abs(sol%x-x) < 1.e-7_dp)
            call check_true(name//' objective matches SciPy', abs(sol%f_x-f_x) < 1.e-7_dp)
        end if
    end subroutine check_case

    subroutine check_int(name, actual, expected)
        character(len=*), intent(in) :: name !! comparison name
        integer, intent(in) :: actual, expected !! observed and expected
        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') 'FAIL '//name, actual, expected
            failures = failures+1
        end if
    end subroutine check_int

    subroutine check_true(name, ok)
        character(len=*), intent(in) :: name !! comparison name
        logical, intent(in) :: ok !! predicate
        if (.not. ok) then
            write(*, '(a)') 'FAIL '//name
            failures = failures+1
        end if
    end subroutine check_true

    function linear(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        calls = calls+1
        y = x-1.25_dp
    end function linear

    function cubic(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        y = x*x*x-2.0_dp*x-5.0_dp
    end function cubic

    function cosine(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        y = cos(x)-x
    end function cosine

    function exponential(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        y = exp(x)-3.0_dp
    end function exponential

    function cubed(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        y = (x-0.3_dp)**3
    end function cubed

    function left_zero(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        y = x
    end function left_zero

    function right_zero(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        y = x-2.0_dp
    end function right_zero

    function invalid(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        y = x*x+1.0_dp
    end function invalid

    function negative(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        y = x+2.3_dp
    end function negative

    function flat(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        y = 0.0_dp*x
    end function flat

    function nonfinite(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        calls = calls+1
        if (x == 0.0_dp) then
            y = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            y = x-0.3_dp
        end if
    end function nonfinite

    function late_nonfinite(x) result(y)
        real(dp), intent(in) :: x !! test abscissa
        real(dp) :: y !! value
        calls = calls+1
        if (x > 0.4_dp .and. x < 0.6_dp) then
            y = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            y = x-0.3_dp
        end if
    end function late_nonfinite

    subroutine stop_after_two(state, stop)
        type(find_root_result), intent(in) :: state !! current result
        logical, intent(out) :: stop !! request termination
        callback_calls = callback_calls+1
        last_nit = state%nit
        callbacks_in_progress = callbacks_in_progress .and. state%status == 1
        stop = state%nit >= 2
    end subroutine stop_after_two

    subroutine stop_initial(state, stop)
        type(find_root_result), intent(in) :: state !! current result
        logical, intent(out) :: stop !! request termination
        callback_calls = callback_calls+1
        stop = .true.
    end subroutine stop_initial

    subroutine never_stop(state, stop)
        type(find_root_result), intent(in) :: state !! current result
        logical, intent(out) :: stop !! request termination
        callback_calls = callback_calls+1
        stop = .false.
    end subroutine never_stop
end program test_optimize_find_root
