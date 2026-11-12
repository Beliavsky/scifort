! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_design
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_signal, only : dp, firwin, freqz, sosfreqz, butter, butter_result, &
        frequency_response, design_status_success, design_status_invalid_input
    implicit none
    type(butter_result) :: design
    type(frequency_response) :: response, from_sos
    real(dp), allocatable :: taps(:)
    real(dp), parameter :: tol = 4.0e-12_dp
    complex(dp) :: cb(2), ca(2), q
    real(dp) :: frequencies(7)
    integer :: status, k, degree, mode
    character(8) :: names(4) = ['lowpass ', 'highpass', 'bandpass', 'bandstop']

    taps = firwin(7,[0.35_dp],window='hamming',status=status)
    if (status /= design_status_success) error stop 'FIR lowpass status'
    call check_real(taps,[-0.001536969_dp,0.046201278_dp,0.252776226_dp, &
        0.405118930_dp,0.252776226_dp,0.046201278_dp,-0.001536969_dp], 3.0e-8_dp)
    if (abs(sum(taps)-1.0_dp)>tol) error stop 'FIR unity DC'
    taps = firwin(7,[0.35_dp],pass_zero='highpass',status=status)
    if (status /= design_status_success) error stop 'FIR highpass status'
    if (abs(sum(taps(2:7:2))-sum(taps(1:7:2))-1.0_dp)>tol) &
        error stop 'FIR unity Nyquist'
    taps = firwin(9,[0.2_dp,0.6_dp],pass_zero='bandpass',status=status)
    if (status /= design_status_success) error stop 'FIR bandpass status'
    if (maxval(abs(taps-taps(9:1:-1)))>tol) error stop 'FIR band symmetry'
    taps = firwin(9,[0.2_dp,0.6_dp],pass_zero='bandstop',status=status)
    if (status /= design_status_success) error stop 'FIR bandstop status'
    if (abs(sum(taps)-1.0_dp)>tol) error stop 'FIR bandstop gain'
    taps = firwin(6,[0.4_dp],pass_zero='highpass',status=status)
    if (status /= design_status_invalid_input .or. size(taps)/=0) &
        error stop 'FIR even Nyquist passband'
    taps = firwin(7,[0.5_dp,0.3_dp],pass_zero='bandpass',status=status)
    if (status /= design_status_invalid_input) error stop 'FIR descending cutoffs'
    taps = firwin(7,[0.4_dp],window='unrecognized',status=status)
    if (status /= design_status_invalid_input) error stop 'FIR unknown window'

    design = butter(3,0.25_dp)
    if (design%status /= design_status_success) error stop 'Butterworth status'
    call check_real(design%b,[0.03168934384971104_dp,0.09506803154913313_dp, &
        0.09506803154913313_dp,0.03168934384971104_dp],tol)
    call check_real(design%a,[1.0_dp,-1.459029062228061_dp,0.9103690002900687_dp, &
        -0.1978251872643194_dp],tol)
    if (size(design%sos,1)/=2) error stop 'Butterworth section count'
    frequencies = [0.0_dp,0.06_dp,0.1_dp,0.17_dp,0.25_dp,0.38_dp,0.49_dp]
    do mode=1,4
        do degree=1,8
            if (mode<=2) then
                design=butter(degree,0.37_dp,btype=trim(names(mode)))
            else
                design=butter(degree,[0.18_dp,0.62_dp],btype=trim(names(mode)))
            end if
            if (design%status/=design_status_success) error stop 'Butterworth general design'
            if (.not. all(ieee_is_finite(design%a))) error stop 'Butterworth nonfinite a'
            if (.not. all(ieee_is_finite(design%b))) error stop 'Butterworth nonfinite b'
            response=freqz(design%b,design%a,w=frequencies,fs=1.0_dp)
            from_sos=sosfreqz(design%sos,w=frequencies,fs=1.0_dp)
            if (response%status/=design_status_success .or. &
                from_sos%status/=design_status_success) error stop 'frequency response status'
            if (maxval(abs(response%response-from_sos%response))>3.0e-8_dp) &
                error stop 'frequency response mismatch'
        end do
    end do
    design=butter(2,0.4_dp,btype='highpass')
    call check_real(design%b,[0.3913357725017685_dp,-0.7826715450035371_dp, &
        0.3913357725017685_dp],tol)
    call check_real(design%a,[1.0_dp,-0.3695273773512414_dp,0.1958157126558331_dp],tol)
    design=butter(2,[0.2_dp,0.6_dp],btype='bandpass')
    call check_real(design%b,[0.2065720838261479_dp,0.0_dp,-0.4131441676522958_dp, &
        0.0_dp,0.2065720838261479_dp],tol)
    design=butter(2,[0.2_dp,0.6_dp],btype='bandstop')
    call check_real(design%b,[0.3913357725017686_dp,-0.5979078563279168_dp, &
        1.011052023980213_dp,-0.5979078563279168_dp,0.3913357725017687_dp],tol)
    design=butter(0,0.4_dp)
    if (design%status/=design_status_invalid_input) error stop 'Butterworth zero order'
    design=butter(3,[0.4_dp,0.2_dp],btype='bandpass')
    if (design%status/=design_status_invalid_input) error stop 'Butterworth inverted band'
    design=butter(3,1.0_dp)
    if (design%status/=design_status_invalid_input) error stop 'Butterworth Nyquist edge'

    cb = [cmplx(0.4_dp,0.25_dp,dp),cmplx(0.2_dp,-0.1_dp,dp)]
    ca = [cmplx(1.0_dp,0.0_dp,dp),cmplx(-0.32_dp,0.1_dp,dp)]
    response=freqz(cb,ca,w=frequencies,fs=1.0_dp)
    if (response%status /= design_status_success) error stop 'complex freqz status'
    do k=1,size(frequencies)
        q=cmplx(cos(-2.0_dp*acos(-1.0_dp)*frequencies(k)), &
            sin(-2.0_dp*acos(-1.0_dp)*frequencies(k)),dp)
        if (abs(response%response(k)-(cb(1)+cb(2)*q)/(ca(1)+ca(2)*q))>tol) &
            error stop 'complex freqz numerical value'
    end do
    response=freqz([1.0_dp],worN=11,include_nyquist=.true.)
    if (response%status/=design_status_success) error stop 'freqz optional a'
    if (abs(response%frequencies(11)-acos(-1.0_dp))>tol) &
        error stop 'freqz Nyquist endpoint'
    if (any(abs(response%response-cmplx(1.0_dp,0.0_dp,dp))>tol)) &
        error stop 'constant response'
    response=sosfreqz(reshape([1.0_dp,0.0_dp,0.0_dp,2.0_dp,0.0_dp,0.0_dp],[1,6]))
    if (response%status/=design_status_invalid_input) error stop 'sos a0 validation'

    print '(a)', 'test_signal_design passed'
contains
    subroutine check_real(v,w,eps)
        real(dp), intent(in) :: v(:) !! actual values
        real(dp), intent(in) :: w(:) !! known SciPy reference values
        real(dp), intent(in) :: eps !! maximum permitted difference
        if (size(v)/=size(w)) error stop 'reference shape mismatch'
        if (any(abs(v-w)>eps)) error stop 'SciPy reference mismatch'
    end subroutine check_real
end program test_signal_design
