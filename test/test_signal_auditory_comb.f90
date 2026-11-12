! SPDX-License-Identifier: MIT
! SciPy-compatible gammatone and IIR comb filter design regressions.
program test_signal_auditory_comb
    use scifort_signal, only : dp, butter_result, gammatone, iircomb, &
        design_status_success, design_status_invalid_input
    implicit none
    type(butter_result) :: filt
    real(dp), parameter :: expected_fir(8) = [ &
        0.0_dp, 2.2219671862468276e-7_dp, 1.649421011505555e-6_dp, &
        4.992982267588435e-6_dp, 1.0199396946176995e-5_dp, &
        1.631257696277016e-5_dp, 2.1464894019747388e-5_dp, &
        2.2994726253098723e-5_dp]
    real(dp), parameter :: eps = 1.0e-10_dp

    filt = gammatone(440.0_dp, 'FIR', numtaps=16, fs=16000.0_dp)
    call assert_status(filt, design_status_success)
    if (size(filt%b) /= 16 .or. size(filt%a) /= 1) error stop 'FIR shapes'
    if (maxval(abs(filt%b(1:8)-expected_fir)) > eps) error stop 'FIR example'
    if (filt%a(1) /= 1.0_dp) error stop 'FIR denominator'

    filt = gammatone(440.0_dp, 'iir', fs=16000.0_dp)
    call assert_status(filt, design_status_success)
    if (size(filt%b) /= 5 .or. size(filt%a) /= 9) error stop 'IIR shapes'
    if (filt%a(1) /= 1.0_dp) error stop 'IIR denominator'
    if (maxval(abs(filt%a)) < 1.0_dp) error stop 'IIR coefficients'

    filt = gammatone(300.0_dp, 'fir', fs=1000.0_dp)
    call assert_status(filt, design_status_success)
    if (size(filt%b) /= 15) error stop 'low sampling-rate default length'
    filt = gammatone(1000.0_dp, 'fir', fs=10000.0_dp, order=1, numtaps=3)
    call assert_status(filt, design_status_success)
    if (filt%b(1) <= 0.0_dp) error stop 'first order response'

    filt = iircomb(20.0_dp,30.0_dp,fs=200.0_dp)
    call assert_status(filt, design_status_success)
    if (size(filt%b) /= 11 .or. size(filt%a) /= 11) error stop 'comb order'
    if (abs(sum(filt%b)) > eps) error stop 'default notch DC rejection'
    filt = iircomb(20.0_dp,30.0_dp,fs=200.0_dp,pass_zero=.true.)
    call assert_status(filt, design_status_success)
    if (abs(sum(filt%b)) < 0.01_dp) error stop 'notch DC pass'
    filt = iircomb(25.0_dp,30.0_dp,'PEAK',fs=200.0_dp)
    call assert_status(filt, design_status_success)
    if (abs(sum(filt%b)) > eps) error stop 'peak DC rejection'
    filt = iircomb(25.0_dp,30.0_dp,'peak',fs=200.0_dp,pass_zero=.true.)
    call assert_status(filt, design_status_success)
    if (abs(sum(filt%b)) < 0.01_dp) error stop 'peak DC pass'

    filt = iircomb(23.0_dp,30.0_dp,fs=200.0_dp)
    call assert_status(filt, design_status_invalid_input)
    filt = iircomb(-1.0_dp,30.0_dp)
    call assert_status(filt, design_status_invalid_input)
    filt = iircomb(0.25_dp,0.0_dp)
    call assert_status(filt, design_status_invalid_input)
    filt = iircomb(0.25_dp,30.0_dp,ftype='wrong')
    call assert_status(filt, design_status_invalid_input)
    filt = gammatone(-3.0_dp,'fir')
    call assert_status(filt, design_status_invalid_input)
    filt = gammatone(0.3_dp,'not_a_type')
    call assert_status(filt, design_status_invalid_input)
    filt = gammatone(0.3_dp,'fir',order=0)
    call assert_status(filt, design_status_invalid_input)
    filt = gammatone(0.3_dp,'fir',numtaps=0)
    call assert_status(filt, design_status_invalid_input)
    filt = gammatone(0.5_dp,'iir',fs=1.0_dp)
    call assert_status(filt, design_status_invalid_input)
    print '(a)', 'test_signal_auditory_comb passed'
contains
    subroutine assert_status(result, wanted)
        type(butter_result), intent(in) :: result !! returned coefficients/status
        integer, intent(in) :: wanted !! expected design status
        if (result%status /= wanted) then
            print '(a,2(i0,1x))', 'Unexpected status: ', result%status, wanted
            error stop 'filter status mismatch'
        end if
    end subroutine assert_status
end program test_signal_auditory_comb
