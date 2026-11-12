! SPDX-License-Identifier: MIT
! Scalar scipy.optimize.elementwise.bracket_minimum compatibility regressions.
program test_optimize_bracket_minimum
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize, only : bracket_minimum, bracket_minimum_result, &
        bracket_minimum_status_limit, brent, minimize_result, &
        minimize_status_convergence_error, minimize_status_invalid_input, &
        minimize_status_success, minimize_status_value_error
    use test_optimize_bracket_minimum_reference
    implicit none

    integer :: failures, calls, out_of_bounds
    type(bracket_minimum_result) :: r
    type(minimize_result) :: m

    failures = 0
    out_of_bounds = 0
    calls = 0
    r = bracket_minimum(right, 0.0_dp)
    call check_case('right', r, [right_xl_ref, right_xm_ref, right_xr_ref], &
        [right_fl_ref, right_fm_ref, right_fr_ref], &
        right_nit_ref, right_nfev_ref, right_status_ref)
    m = brent(right, bracket=[r%xl, r%xm, r%xr])
    call check_true('bracket feeds brent', m%success .and. abs(m%x - 5.0_dp) < 1.0e-7_dp)

    calls = 0
    r = bracket_minimum(left, 0.0_dp)
    call check_case('left', r, [left_xl_ref, left_xm_ref, left_xr_ref], &
        [left_fl_ref, left_fm_ref, left_fr_ref], &
        left_nit_ref, left_nfev_ref, left_status_ref)

    calls = 0
    r = bracket_minimum(initial, 0.0_dp)
    call check_case('initial', r, [initial_xl_ref, initial_xm_ref, initial_xr_ref], &
        [initial_fl_ref, initial_fm_ref, initial_fr_ref], &
        initial_nit_ref, initial_nfev_ref, initial_status_ref)

    calls = 0
    r = bracket_minimum(right, 0.0_dp, factor=3.0_dp)
    call check_case('rightfactor3', r, &
        [rightfactor3_xl_ref, rightfactor3_xm_ref, rightfactor3_xr_ref], &
        [rightfactor3_fl_ref, rightfactor3_fm_ref, rightfactor3_fr_ref], &
        rightfactor3_nit_ref, rightfactor3_nfev_ref, rightfactor3_status_ref)

    calls = 0
    r = bracket_minimum(right_bounded, 0.0_dp, xmax=6.0_dp)
    call check_case('rightbound', r, [rightbound_xl_ref, rightbound_xm_ref, rightbound_xr_ref], &
        [rightbound_fl_ref, rightbound_fm_ref, rightbound_fr_ref], &
        rightbound_nit_ref, rightbound_nfev_ref, rightbound_status_ref)
    call check_int('right bound safety', out_of_bounds, 0)

    calls = 0
    r = bracket_minimum(left, 0.0_dp, xmin=-5.0_dp)
    call check_case('leftbound', r, [leftbound_xl_ref, leftbound_xm_ref, leftbound_xr_ref], &
        [leftbound_fl_ref, leftbound_fm_ref, leftbound_fr_ref], &
        leftbound_nit_ref, leftbound_nfev_ref, leftbound_status_ref)

    calls = 0
    r = bracket_minimum(right, 0.0_dp, xl0=-1.0_dp, xr0=2.0_dp)
    call check_case('custom', r, [custom_xl_ref, custom_xm_ref, custom_xr_ref], &
        [custom_fl_ref, custom_fm_ref, custom_fr_ref], &
        custom_nit_ref, custom_nfev_ref, custom_status_ref)

    calls = 0
    r = bracket_minimum(left, 0.0_dp, xl0=-1.0_dp, xr0=1.0_dp)
    call check_case('customleft', r, [customleft_xl_ref, customleft_xm_ref, customleft_xr_ref], &
        [customleft_fl_ref, customleft_fm_ref, customleft_fr_ref], &
        customleft_nit_ref, customleft_nfev_ref, customleft_status_ref)

    calls = 0
    r = bracket_minimum(right, 0.0_dp, maxiter=0)
    call check_case('max0', r, [max0_xl_ref, max0_xm_ref, max0_xr_ref], &
        [max0_fl_ref, max0_fm_ref, max0_fr_ref], max0_nit_ref, max0_nfev_ref, max0_status_ref)

    calls = 0
    r = bracket_minimum(right, 0.0_dp, maxiter=2)
    call check_case('max2', r, [max2_xl_ref, max2_xm_ref, max2_xr_ref], &
        [max2_fl_ref, max2_fm_ref, max2_fr_ref], max2_nit_ref, max2_nfev_ref, max2_status_ref)

    calls = 0
    r = bracket_minimum(constant, 0.0_dp, maxiter=3)
    call check_case('constant', r, [constant_xl_ref, constant_xm_ref, constant_xr_ref], &
        [constant_fl_ref, constant_fm_ref, constant_fr_ref], &
        constant_nit_ref, constant_nfev_ref, constant_status_ref)

    calls = 0
    r = bracket_minimum(plateau, 0.0_dp)
    call check_case('plateau', r, [plateau_xl_ref, plateau_xm_ref, plateau_xr_ref], &
        [plateau_fl_ref, plateau_fm_ref, plateau_fr_ref], &
        plateau_nit_ref, plateau_nfev_ref, plateau_status_ref)

    calls = 0
    r = bracket_minimum(shifted, 0.0_dp)
    call check_case('shifted', r, [shifted_xl_ref, shifted_xm_ref, shifted_xr_ref], &
        [shifted_fl_ref, shifted_fm_ref, shifted_fr_ref], &
        shifted_nit_ref, shifted_nfev_ref, shifted_status_ref)

    calls = 0
    r = bracket_minimum(initial, 0.0_dp, xmin=-0.1_dp, xmax=0.1_dp)
    call check_case('boundboth', r, [boundboth_xl_ref, boundboth_xm_ref, boundboth_xr_ref], &
        [boundboth_fl_ref, boundboth_fm_ref, boundboth_fr_ref], &
        boundboth_nit_ref, boundboth_nfev_ref, boundboth_status_ref)

    calls = 0
    r = bracket_minimum(monotone, 0.0_dp, xmax=1.0_dp)
    call check_case('rightlimit', r, [rightlimit_xl_ref, rightlimit_xm_ref, rightlimit_xr_ref], &
        [rightlimit_fl_ref, rightlimit_fm_ref, rightlimit_fr_ref], &
        rightlimit_nit_ref, rightlimit_nfev_ref, rightlimit_status_ref)
    call check_int('limit status constant', r%status, bracket_minimum_status_limit)

    calls = 0
    r = bracket_minimum(monotone_left, 0.0_dp, xmin=-1.0_dp)
    call check_case('leftlimit', r, [leftlimit_xl_ref, leftlimit_xm_ref, leftlimit_xr_ref], &
        [leftlimit_fl_ref, leftlimit_fm_ref, leftlimit_fr_ref], &
        leftlimit_nit_ref, leftlimit_nfev_ref, leftlimit_status_ref)

    calls = 0
    r = bracket_minimum(monotone, 0.0_dp, xl0=-0.2_dp, xr0=1.0_dp, xmax=1.0_dp)
    call check_case('limitstart', r, [limitstart_xl_ref, limitstart_xm_ref, limitstart_xr_ref], &
        [limitstart_fl_ref, limitstart_fm_ref, limitstart_fr_ref], &
        limitstart_nit_ref, limitstart_nfev_ref, limitstart_status_ref)

    calls = 0
    r = bracket_minimum(right, 0.0_dp, factor=1.0_dp)
    call check_int('invalid growth', r%status, minimize_status_invalid_input)
    call check_int('invalid growth callbacks', calls, 0)

    r = bracket_minimum(right, 0.0_dp, factor=ieee_value(0.0_dp, ieee_positive_inf))
    call check_int('invalid infinite factor', r%status, minimize_status_invalid_input)

    r = bracket_minimum(right, 0.0_dp, maxiter=-1)
    call check_int('invalid maxiter', r%status, minimize_status_invalid_input)
    call check_int('invalid no callbacks', calls, 0)

    r = bracket_minimum(right, 0.0_dp, xl0=0.2_dp)
    call check_int('invalid left initial', r%status, minimize_status_invalid_input)
    call check_int('invalid left no callbacks', calls, 0)

    r = bracket_minimum(right, 0.0_dp, xr0=0.0_dp)
    call check_int('invalid right initial', r%status, minimize_status_invalid_input)

    r = bracket_minimum(right, 0.0_dp, xmax=0.0_dp)
    call check_int('invalid interval to bound', r%status, minimize_status_invalid_input)

    r = bracket_minimum(right, 0.0_dp, xmin=2.0_dp)
    call check_int('invalid left bound', r%status, minimize_status_invalid_input)

    r = bracket_minimum(right, 0.0_dp, xmin=ieee_value(0.0_dp, ieee_quiet_nan))
    call check_int('nan bound', r%status, minimize_status_invalid_input)

    calls = 0
    r = bracket_minimum(not_finite, 0.0_dp)
    call check_int('nan objective status', r%status, minimize_status_value_error)
    call check_int('nan objective work', r%nfev, 3)
    call check_int('nan objective callbacks', calls, 3)

    calls = 0
    r = bracket_minimum(nan_later, 0.0_dp)
    call check_int('later nan objective status', r%status, minimize_status_value_error)
    call check_int('later nan objective callbacks', calls, r%nfev)
    call check_true('later nan objective progressed', r%nit > 0)

    calls = 0
    r = bracket_minimum(right, 0.0_dp, &
        xmin=ieee_value(0.0_dp, ieee_positive_inf))
    call check_int('wrong infinity bound', r%status, minimize_status_invalid_input)
    call check_int('wrong infinity no callbacks', calls, 0)

    call check_int('success status constant', minimize_status_success, 0)
    call check_int('maxiter status constant', minimize_status_convergence_error, -2)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_bracket_minimum: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_bracket_minimum: PASS'

contains

    function right(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        y = (x - 5.0_dp)**2
    end function right

    function right_bounded(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        if (x > 6.0_dp) out_of_bounds = out_of_bounds + 1
        y = right(x)
    end function right_bounded

    function left(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        y = (x + 3.0_dp)**2
    end function left

    function initial(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        y = x*x
    end function initial

    function constant(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        y = 0.0_dp*x
    end function constant

    function monotone(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        y = -x
    end function monotone

    function monotone_left(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        y = x
    end function monotone_left

    function plateau(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        if (x >= 1.0_dp) then
            y = 0.0_dp
        else
            y = (x - 1.0_dp)**2
        end if
    end function plateau

    function shifted(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        y = (x - 1.0_dp)**2 + 2.0_dp
    end function shifted

    function not_finite(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        y = ieee_value(x, ieee_quiet_nan)
    end function not_finite

    function nan_later(x) result(y)
        real(dp), intent(in) :: x !! scalar abscissa
        real(dp) :: y
        calls = calls + 1
        if (x > 1.0_dp) then
            y = ieee_value(x, ieee_quiet_nan)
        else
            y = (x - 5.0_dp)**2
        end if
    end function nan_later

    subroutine check_case(label, actual, ref_x, ref_f, ref_nit, ref_nfev, ref_status)
        character(len=*), intent(in) :: label !! reference test name
        type(bracket_minimum_result), intent(in) :: actual !! observed bracket result
        real(dp), intent(in) :: ref_x(3) !! SciPy bracket endpoints and middle
        real(dp), intent(in) :: ref_f(3) !! SciPy function values at bracket points
        integer, intent(in) :: ref_nit !! SciPy expansion count
        integer, intent(in) :: ref_nfev !! SciPy objective evaluation count
        integer, intent(in) :: ref_status !! SciPy termination status
        real(dp), parameter :: tol = 2.0e-12_dp

        call check_int(trim(label)//' status', actual%status, ref_status)
        call check_int(trim(label)//' iterations', actual%nit, ref_nit)
        call check_int(trim(label)//' evaluations', actual%nfev, ref_nfev)
        call check_int(trim(label)//' callbacks', calls, ref_nfev)
        call check_vector(trim(label)//' abscissae', [actual%xl, actual%xm, actual%xr], ref_x, tol)
        call check_vector(trim(label)//' ordinates', [actual%fl, actual%fm, actual%fr], ref_f, tol)
        call check_true(trim(label)//' success', actual%success .eqv. (ref_status == 0))
        call check_true(trim(label)//' order', actual%xl < actual%xm .and. actual%xm < actual%xr)
        if (actual%success) then
            call check_true(trim(label)//' valid bracket', &
                (actual%fl >= actual%fm .and. actual%fr > actual%fm) .or. &
                (actual%fl > actual%fm .and. actual%fr >= actual%fm))
        end if
    end subroutine check_case

    subroutine check_vector(label, actual, expected, tol)
        character(len=*), intent(in) :: label !! assertion label
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! SciPy reference values
        real(dp), intent(in) :: tol !! relative/absolute tolerance
        if (any(abs(actual - expected) > tol * max(1.0_dp, abs(expected)))) then
            write(*, '(a,3(1x,es23.15))') trim(label)//' observed:', actual
            write(*, '(a,3(1x,es23.15))') trim(label)//' expected:', expected
            failures = failures + 1
        end if
    end subroutine check_vector

    subroutine check_int(label, actual, expected)
        character(len=*), intent(in) :: label !! assertion label
        integer, intent(in) :: actual !! computed integer
        integer, intent(in) :: expected !! reference integer
        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') trim(label)//' observed, expected:', actual, expected
            failures = failures + 1
        end if
    end subroutine check_int

    subroutine check_true(label, condition)
        character(len=*), intent(in) :: label !! assertion label
        logical, intent(in) :: condition !! assertion to verify
        if (.not. condition) then
            write(*, '(a)') trim(label)//': FAIL'
            failures = failures + 1
        end if
    end subroutine check_true
end program test_optimize_bracket_minimum
