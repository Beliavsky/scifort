! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Validate scalar and logical firwin overloads and odd-edge multiband designs.
program test_signal_firwin_multiband
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_design, only : firwin, design_status_success, &
        design_status_invalid_input
    implicit none

    real(dp), allocatable :: taps(:), ref(:), width_ref(:)
    integer :: st
    real(dp), parameter :: tol = 3.0e-13_dp

    ! Existing array/string and omitted-argument calls must still resolve.
    ref = firwin(19, [0.28_dp], pass_zero='lowpass', window='hamming')
    taps = firwin(19, 0.28_dp, status=st)
    call check(st == design_status_success, 'scalar lowpass status')
    call check(size(taps) == 19, 'scalar cutoff length')
    call check(maxval(abs(taps-ref)) < tol, 'scalar vs vector lowpass')

    taps = firwin(19, 0.28_dp, pass_zero=.true., status=st)
    call check(st == design_status_success, 'scalar logical lowpass')
    call check(maxval(abs(taps-ref)) < tol, 'scalar logical lowpass equivalence')

    ref = firwin(19, [0.28_dp], pass_zero='highpass')
    taps = firwin(19, 0.28_dp, pass_zero=.false., status=st)
    call check(st == design_status_success, 'scalar logical highpass')
    call check(maxval(abs(taps-ref)) < tol, 'scalar logical highpass equivalence')

    ! Three transitions, with pass_zero true, stop at Nyquist; with false,
    ! pass at Nyquist (and require odd numtaps, SciPy Type-I/II restriction).
    ref = firwin(18, [0.13_dp, 0.39_dp, 0.66_dp], pass_zero='bandstop', status=st)
    call check(st == design_status_success, 'three-edge bandstop, even taps')
    taps = firwin(18, [0.13_dp, 0.39_dp, 0.66_dp], pass_zero=.true.)
    call check(maxval(abs(taps-ref)) < tol, 'three-edge logical bandstop')
    call check(maxval(abs(taps(1:9)-taps(18:10:-1))) < tol, 'even symmetric multiband')

    taps = firwin(18, [0.13_dp, 0.39_dp, 0.66_dp], pass_zero=.false., status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, &
        'even taps reject Nyquist passband on odd edge count')
    taps = firwin(19, [0.13_dp, 0.39_dp, 0.66_dp], pass_zero='bandpass', status=st)
    call check(st == design_status_success, 'odd-edge bandpass at Nyquist')
    ref = firwin(19, [0.13_dp, 0.39_dp, 0.66_dp], pass_zero=.false.)
    call check(maxval(abs(taps-ref)) < tol, 'three-edge logical bandpass')

    ! Five transition frequencies: either starting passband choice works.
    taps = firwin(24, [0.08_dp,0.16_dp,0.28_dp,0.57_dp,0.72_dp], &
        pass_zero=.true., status=st)
    call check(st == design_status_success, 'five-edge bandstop')
    call check(abs(sum(taps)-1.0_dp) < tol, 'DC unity after scaling')
    taps = firwin(25, [0.08_dp,0.16_dp,0.28_dp,0.57_dp,0.72_dp], &
        pass_zero=.false., status=st)
    call check(st == design_status_success, 'five-edge bandpass')
    call check(maxval(abs(taps-taps(25:1:-1))) < tol, 'odd multiband symmetry')

    ! With four cuts, a bandstop passes at Nyquist and requires odd taps.
    taps = firwin(24, [0.15_dp,0.25_dp,0.45_dp,0.70_dp], pass_zero=.true., status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, &
        'four-cut bandstop Nyquist restriction')
    taps = firwin(24, [0.15_dp,0.25_dp,0.45_dp,0.70_dp], pass_zero=.false., status=st)
    call check(st == design_status_success, 'four-cut bandpass even taps')

    ! width overrides the named window and may be combined with multiband
    ! and logical pass_zero; both scale choices should pass through unchanged.
    taps = firwin(21, [0.14_dp,0.3_dp,0.5_dp], pass_zero=.false., &
        fs=2.0_dp, width=0.04_dp, window='window-ignored', scale=.false., status=st)
    call check(st == design_status_success, 'width plus logical multiband')
    width_ref = firwin(21, [0.14_dp,0.3_dp,0.5_dp], pass_zero='bandpass', &
        fs=2.0_dp, width=0.04_dp, scale=.false.)
    call check(maxval(abs(taps-width_ref)) < tol, 'logical plus width tap equivalence')

    taps = firwin(15, 0.3_dp, pass_zero='bandstop', status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, &
        'scalar bandstop requires at least two cutoffs')
    taps = firwin(15, [0.2_dp,0.5_dp,0.3_dp], pass_zero=.true., status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, &
        'reject unordered multiband frequencies')
    taps = firwin(15, [0.2_dp,0.5_dp,0.5_dp], pass_zero=.true., status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, &
        'reject duplicate multiband frequencies')
    taps = firwin(15, [0.2_dp,0.5_dp,0.9_dp], pass_zero='invalid', status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, &
        'reject invalid pass_zero keyword')
    taps = firwin(15, ieee_value(0.0_dp,ieee_quiet_nan), status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, &
        'reject invalid scalar cutoff')

    print '(a)', 'PASS: scalar, logical and odd-edge multiband FIR designs'

contains

    subroutine check(condition, message)
        logical, intent(in) :: condition !! assertion to check
        character(*), intent(in) :: message !! diagnostic description
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine check

end program test_signal_firwin_multiband
