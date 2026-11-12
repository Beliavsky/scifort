! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Native regressions for transition-width Kaiser FIR and shaped firwin2.
program test_signal_fir_width_windows
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_design, only : firwin, kaiser_atten, kaiser_beta, &
        design_status_success, design_status_invalid_input
    use scifort_signal_fir_advanced, only : firwin2, fir_advanced_status_success, &
        fir_advanced_status_invalid_input
    use scifort_signal_spectral, only : get_window
    implicit none

    real(dp), parameter :: tol = 5.0e-13_dp
    real(dp), allocatable :: taps(:), expected(:), weights(:)
    integer :: st

    ! width uses the physical transition width, normalized by fs/2 before
    ! applying the Kaiser attenuation and beta formulas.
    taps = firwin(21, [3.5_dp], fs=40.0_dp, width=2.2_dp, &
        window='nonexistent-window-ignored-when-width-given', status=st)
    call check(st == design_status_success, 'firwin transition-width status')
    expected = firwin(21, [3.5_dp], fs=40.0_dp, window='kaiser', &
        window_parameter=kaiser_beta(kaiser_atten(21, 2.2_dp/20.0_dp)))
    call check(size(taps) == 21, 'transition-width tap count')
    call check(maxval(abs(taps-expected)) < tol, 'firwin width uses Kaiser beta')

    taps = firwin(18, [0.22_dp, 0.62_dp], pass_zero='bandpass', width=0.05_dp, &
        scale=.false., status=st)
    call check(st == design_status_success, 'even-tap bandpass with width')
    call check(maxval(abs(taps(1:9)-taps(18:10:-1))) < tol, 'bandpass symmetry')

    taps = firwin(13, [0.4_dp], width=0.0_dp, status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, 'reject zero width')
    taps = firwin(13, [0.4_dp], width=-0.03_dp, status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, 'reject negative width')
    taps = firwin(13, [0.4_dp], width=ieee_value(0.0_dp,ieee_quiet_nan), status=st)
    call check(st == design_status_invalid_input .and. size(taps) == 0, 'reject NaN width')
    taps = firwin(13, [0.4_dp], window='boxcar', &
        window_parameter=3.0_dp, status=st)
    call check(st == design_status_invalid_input, 'do not ignore unexpected window parameter')

    taps = firwin2(15, [0.0_dp,0.4_dp,1.0_dp], &
        [1.0_dp,1.0_dp,0.0_dp], window='kaiser', window_parameter=3.4_dp, status=st)
    call check(st == fir_advanced_status_success, 'firwin2 Kaiser window')
    weights = get_window('kaiser', 15, fftbins=.false., parameter=3.4_dp)
    expected = firwin2(15, [0.0_dp,0.4_dp,1.0_dp], [1.0_dp,1.0_dp,0.0_dp], &
        window='none', status=st)
    call check(maxval(abs(taps-expected*weights)) < tol, 'Kaiser window multiplies FIR design')

    taps = firwin2(15, [0.0_dp,0.4_dp,1.0_dp], [1.0_dp,1.0_dp,0.0_dp], &
        window_values=weights, status=st)
    call check(st == fir_advanced_status_success, 'explicit symmetric window')
    call check(maxval(abs(taps-expected*weights)) < tol, 'explicit window applied')

    taps = firwin2(13, [0.0_dp,0.2_dp,1.0_dp], [1.0_dp,0.5_dp,0.0_dp], &
        window='general_gaussian', window_parameter=1.7_dp, &
        window_parameter2=4.0_dp, status=st)
    call check(st == fir_advanced_status_success, 'general Gaussian two-parameter window')

    taps = firwin2(11, [0.0_dp,0.25_dp,1.0_dp], [1.0_dp,1.0_dp,0.0_dp], &
        window='general_cosine', window_coefficients=[0.4_dp,0.3_dp,0.2_dp,0.1_dp], &
        status=st)
    call check(st == fir_advanced_status_success, 'general cosine window')

    taps = firwin2(11, [0.0_dp,0.25_dp,1.0_dp], [1.0_dp,1.0_dp,0.0_dp], &
        window_values=weights, status=st)
    call check(st == fir_advanced_status_invalid_input .and. size(taps) == 0, &
        'explicit window length mismatch')
    taps = firwin2(11, [0.0_dp,0.25_dp,1.0_dp], [1.0_dp,1.0_dp,0.0_dp], &
        window='hann', window_parameter=2.0_dp, status=st)
    call check(st == fir_advanced_status_invalid_input, 'extraneous parameter rejected')
    taps = firwin2(11, [0.0_dp,0.25_dp,1.0_dp], [1.0_dp,1.0_dp,0.0_dp], &
        window='hann', window_values=weights(:11), status=st)
    call check(st == fir_advanced_status_invalid_input, 'conflicting named and explicit windows')
    taps = firwin2(11, [0.0_dp,0.25_dp,1.0_dp], [1.0_dp,1.0_dp,0.0_dp], &
        window='none', window_parameter=2.0_dp, status=st)
    call check(st == fir_advanced_status_invalid_input, 'no parameter for window none')
    taps = firwin2(11, [0.0_dp,0.25_dp,1.0_dp], [1.0_dp,1.0_dp,0.0_dp], &
        window='gaussian', window_parameter=3.0_dp, window_parameter2=1.2_dp, status=st)
    call check(st == fir_advanced_status_invalid_input, 'unsupported second parameter')
    print '(a)', 'PASS: firwin transition-width and firwin2 parameter-window native tests'

contains
    subroutine check(condition, message)
        logical, intent(in) :: condition !! evaluated assertion
        character(*), intent(in) :: message !! diagnostic label
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine check
end program test_signal_fir_width_windows
