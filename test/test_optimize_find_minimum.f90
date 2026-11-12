! SPDX-License-Identifier: MIT
! Scalar SciPy elementwise.find_minimum / Chandrupatla regressions.
program test_optimize_find_minimum
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize, only : bracket_minimum, bracket_minimum_result, &
        find_minimum, find_minimum_result, find_minimum_status_callback, &
        minimize_status_bracket_error, &
        minimize_status_convergence_error, minimize_status_invalid_input, &
        minimize_status_success, minimize_status_value_error
    use test_optimize_find_minimum_reference
    implicit none

    type(find_minimum_result) :: r
    type(bracket_minimum_result) :: b
    integer :: failures, calls, callback_calls, last_nit
    logical :: callbacks_in_progress

    failures = 0
    calls = 0
    callback_calls = 0
    last_nit = -1
    callbacks_in_progress = .true.

    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp])
    call check_case('shifted', r, [shifted_xl_ref, shifted_xm_ref, shifted_xr_ref], &
        [shifted_fl_ref, shifted_fm_ref, shifted_fr_ref], shifted_nit_ref, shifted_nfev_ref, shifted_status_ref)
    call check_int('shifted callback count', calls, r%nfev)

    calls = 0
    r = find_minimum(shifted, [3.0_dp, 1.0_dp, 0.0_dp])
    call check_case('reversed', r, [reversed_xl_ref, reversed_xm_ref, reversed_xr_ref], &
        [reversed_fl_ref, reversed_fm_ref, reversed_fr_ref], reversed_nit_ref, reversed_nfev_ref, reversed_status_ref)

    calls = 0
    r = find_minimum(other, [-3.0_dp, -1.0_dp, 1.0_dp])
    call check_case('other', r, [other_xl_ref, other_xm_ref, other_xr_ref], &
        [other_fl_ref, other_fm_ref, other_fr_ref], other_nit_ref, other_nfev_ref, other_status_ref)

    calls = 0
    r = find_minimum(quartic, [-2.0_dp, -1.0_dp, 1.0_dp])
    call check_case('quartic', r, [quartic_xl_ref, quartic_xm_ref, quartic_xr_ref], &
        [quartic_fl_ref, quartic_fm_ref, quartic_fr_ref], quartic_nit_ref, quartic_nfev_ref, quartic_status_ref)

    calls = 0
    r = find_minimum(kink, [0.0_dp, 1.0_dp, 3.0_dp])
    call check_case('kink', r, [kink_xl_ref, kink_xm_ref, kink_xr_ref], &
        [kink_fl_ref, kink_fm_ref, kink_fr_ref], kink_nit_ref, kink_nfev_ref, kink_status_ref)

    calls = 0
    r = find_minimum(exponential, [0.0_dp, 1.0_dp, 3.0_dp])
    call check_case('exponential', r, [exponential_xl_ref, exponential_xm_ref, exponential_xr_ref], &
        [exponential_fl_ref, exponential_fm_ref, exponential_fr_ref], &
        exponential_nit_ref, exponential_nfev_ref, exponential_status_ref)

    calls = 0
    r = find_minimum(flat, [0.0_dp, 1.0_dp, 3.0_dp])
    call check_case('flat', r, [flat_xl_ref, flat_xm_ref, flat_xr_ref], &
        [flat_fl_ref, flat_fm_ref, flat_fr_ref], flat_nit_ref, flat_nfev_ref, flat_status_ref)

    calls = 0
    r = find_minimum(invalid, [0.0_dp, 1.0_dp, 2.0_dp])
    call check_int('invalid status', r%status, invalid_status_ref)
    call check_int('invalid work', r%nfev, invalid_nfev_ref)
    call check_true('invalid NaN', ieee_is_nan(r%x) .and. ieee_is_nan(r%f_x))
    call check_true('invalid endpoints', abs(r%xl-invalid_xl_ref) < 1.e-14_dp .and. &
        abs(r%xr-invalid_xr_ref) < 1.e-14_dp)

    calls = 0
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], maxiter=0)
    call check_case('max0', r, [max0_xl_ref, max0_xm_ref, max0_xr_ref], &
        [max0_fl_ref, max0_fm_ref, max0_fr_ref], max0_nit_ref, max0_nfev_ref, max0_status_ref)

    calls = 0
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], maxiter=1)
    call check_case('max1', r, [max1_xl_ref, max1_xm_ref, max1_xr_ref], &
        [max1_fl_ref, max1_fm_ref, max1_fr_ref], max1_nit_ref, max1_nfev_ref, max1_status_ref)

    calls = 0
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], maxiter=2)
    call check_case('max2', r, [max2_xl_ref, max2_xm_ref, max2_xr_ref], &
        [max2_fl_ref, max2_fm_ref, max2_fr_ref], max2_nit_ref, max2_nfev_ref, max2_status_ref)

    calls = 0
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], xatol=0.1_dp)
    call check_case('xatol', r, [xatol_xl_ref, xatol_xm_ref, xatol_xr_ref], &
        [xatol_fl_ref, xatol_fm_ref, xatol_fr_ref], xatol_nit_ref, xatol_nfev_ref, xatol_status_ref)

    calls = 0
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], fatol=0.1_dp)
    call check_case('fatol', r, [fatol_xl_ref, fatol_xm_ref, fatol_xr_ref], &
        [fatol_fl_ref, fatol_fm_ref, fatol_fr_ref], fatol_nit_ref, fatol_nfev_ref, fatol_status_ref)

    calls = 0
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], &
        xrtol=1.e-12_dp, xatol=0.0_dp, fatol=0.0_dp, frtol=0.0_dp)
    call check_case('strict', r, [strict_xl_ref, strict_xm_ref, strict_xr_ref], &
        [strict_fl_ref, strict_fm_ref, strict_fr_ref], strict_nit_ref, strict_nfev_ref, strict_status_ref)

    calls = 0
    r = find_minimum(other, [1.0_dp, -3.0_dp, -1.0_dp])
    call check_case('sort', r, [sort_xl_ref, sort_xm_ref, sort_xr_ref], &
        [sort_fl_ref, sort_fm_ref, sort_fr_ref], sort_nit_ref, sort_nfev_ref, sort_status_ref)

    calls = 0
    r = find_minimum(shifted, [0.0_dp, 1.0_dp])
    call check_int('invalid shape', r%status, minimize_status_invalid_input)
    call check_int('invalid shape no evaluations', calls, 0)
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], xatol=-1.e-8_dp)
    call check_int('invalid tolerance', r%status, minimize_status_invalid_input)
    call check_int('invalid tolerance no evaluations', calls, 0)
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], maxiter=-1)
    call check_int('invalid iterations', r%status, minimize_status_invalid_input)
    call check_int('invalid iterations no evaluations', calls, 0)
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], frtol=ieee_value(0.0_dp, ieee_quiet_nan))
    call check_int('nan tolerance', r%status, minimize_status_invalid_input)
    call check_int('nan tolerance no evaluations', calls, 0)
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 1.0_dp])
    call check_int('repeated bracket point', r%status, minimize_status_bracket_error)

    calls = 0
    r = find_minimum(nonfinite, [0.0_dp, 1.0_dp, 3.0_dp])
    call check_int('nonfinite status', r%status, minimize_status_value_error)
    call check_int('nonfinite nfev', r%nfev, 3)
    call check_true('nonfinite nan', ieee_is_nan(r%x))

    calls = 0
    r = find_minimum(nonfinite_after_start, [0.0_dp, 1.0_dp, 3.0_dp])
    call check_int('late nonfinite status', r%status, minimize_status_value_error)
    call check_int('late nonfinite iterations', r%nit, 1)
    call check_int('late nonfinite work', r%nfev, 4)
    call check_int('late nonfinite callbacks', calls, 4)
    call check_true('late nonfinite NaN minimizer', ieee_is_nan(r%x))

    callback_calls = 0
    last_nit = -1
    callbacks_in_progress = .true.
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], callback=stop_after_two)
    call check_int('callback stopped status', r%status, find_minimum_status_callback)
    call check_int('callback stopped nit', r%nit, 2)
    call check_int('callback stopped nfev', r%nfev, 5)
    call check_int('callback calls', callback_calls, 3)
    call check_true('callback sees in progress', callbacks_in_progress)
    call check_int('callback last iteration', last_nit, 2)
    call check_true('callback result remains finite', .not. ieee_is_nan(r%x))

    callback_calls = 0
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], callback=stop_initial)
    call check_int('initial callback status', r%status, find_minimum_status_callback)
    call check_int('initial callback no iteration', r%nit, 0)
    call check_int('initial callback nfev', r%nfev, 3)
    call check_int('initial callback calls', callback_calls, 1)

    callback_calls = 0
    r = find_minimum(shifted, [0.0_dp, 1.0_dp, 3.0_dp], callback=never_stop)
    call check_true('nonstopping callback success', r%success)
    call check_int('nonstopping callback calls', callback_calls, r%nit + 1)

    b = bracket_minimum(shifted, xm0=0.0_dp)
    call check_true('upstream bracket success', b%success)
    if (b%success) then
        r = find_minimum(shifted, [b%xl, b%xm, b%xr])
        call check_true('bracket feeds minimizer', r%success .and. abs(r%x-1.2_dp) < 1.e-6_dp)
    end if

    call check_int('success status', minimize_status_success, 0)
    call check_int('convergence error status', minimize_status_convergence_error, -2)

    if (failures > 0) then
        write(*, '(a,i0)') 'test_optimize_find_minimum: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_find_minimum: PASS'

contains

    subroutine check_case(name, sol, xx, ff, nit, nfev, status)
        character(len=*), intent(in) :: name !! case name
        type(find_minimum_result), intent(in) :: sol !! solver result
        real(dp), intent(in) :: xx(3) !! SciPy bracket points
        real(dp), intent(in) :: ff(3) !! SciPy bracket values
        integer, intent(in) :: nit, nfev, status !! expected work and code
        real(dp) :: err, ferr

        call check_int(name//' status', sol%status, status)
        call check_int(name//' nit', sol%nit, nit)
        call check_int(name//' nfev', sol%nfev, nfev)
        call check_int(name//' callback nfev', calls, sol%nfev)
        call check_true(name//' success consistency', sol%success .eqv. (status == 0))
        err = maxval(abs([sol%xl, sol%xm, sol%xr] - xx))
        ferr = maxval(abs([sol%fl, sol%fm, sol%fr] - ff))
        if (err > 2.e-7_dp .or. ferr > 2.e-7_dp) then
            write(*,'(a,3es22.12)') name//' actual bracket:', sol%xl, sol%xm, sol%xr
            write(*,'(a,3es22.12)') name//' reference:', xx
            failures = failures + 1
        end if
        call check_true(name//' minimum matches middle', abs(sol%x-sol%xm) < 1.e-14_dp)
        call check_true(name//' ordered bracket', sol%xl <= sol%xm .and. sol%xm <= sol%xr)
    end subroutine check_case

    subroutine check_int(name, actual, expected)
        character(len=*), intent(in) :: name !! test case label
        integer, intent(in) :: actual, expected !! value and expected value
        if (actual /= expected) then
            write(*,'(a,2(1x,i0))') 'FAIL '//name, actual, expected
            failures = failures + 1
        end if
    end subroutine check_int

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! test case label
        logical, intent(in) :: condition !! expected true
        if (.not. condition) then
            write(*,'(a)') 'FAIL '//name
            failures = failures + 1
        end if
    end subroutine check_true

    function shifted(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y !! objective value
        calls = calls + 1
        y = (x-1.2_dp)**2 + 1.0_dp
    end function shifted

    function other(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y !! objective value
        calls = calls + 1
        y = (x+1.25_dp)**2 + 0.5_dp
    end function other

    function quartic(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y !! objective value
        calls = calls + 1
        y = (x+0.75_dp)**4 + 0.2_dp*(x+0.75_dp)**2 + 0.1_dp
    end function quartic

    function kink(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y !! objective value
        calls = calls + 1
        y = abs(x-1.5_dp) + 0.2_dp
    end function kink

    function exponential(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y !! objective value
        calls = calls + 1
        y = exp(x) - 3.0_dp*x
    end function exponential

    function flat(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y !! objective value
        calls = calls + 1
        y = 2.0_dp + 0.0_dp * x
    end function flat

    function invalid(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y !! objective value
        calls = calls + 1
        y = (x-4.0_dp)**2
    end function invalid

    function nonfinite(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y !! objective value
        calls = calls + 1
        if (x > 2.0_dp) then
            y = ieee_value(x, ieee_quiet_nan)
        else
            y = x*x
        end if
    end function nonfinite

    function nonfinite_after_start(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y !! objective value
        calls = calls + 1
        if (x > 1.5_dp .and. x < 2.0_dp) then
            y = ieee_value(x, ieee_quiet_nan)
        else
            y = (x-1.2_dp)**2 + 1.0_dp
        end if
    end function nonfinite_after_start

    subroutine stop_after_two(sol, stop)
        type(find_minimum_result), intent(in) :: sol !! current iterate
        logical, intent(out) :: stop !! request stopping at iteration two
        callback_calls = callback_calls + 1
        callbacks_in_progress = callbacks_in_progress .and. (sol%status == 1)
        last_nit = sol%nit
        stop = sol%nit >= 2
    end subroutine stop_after_two

    subroutine stop_initial(sol, stop)
        type(find_minimum_result), intent(in) :: sol !! current iterate
        logical, intent(out) :: stop !! stop immediately
        callback_calls = callback_calls + 1
        stop = sol%nit == 0
    end subroutine stop_initial

    subroutine never_stop(sol, stop)
        type(find_minimum_result), intent(in) :: sol !! current iterate
        logical, intent(out) :: stop !! do not stop
        callback_calls = callback_calls + 1
        stop = sol%nit < -1
    end subroutine never_stop

end program test_optimize_find_minimum
