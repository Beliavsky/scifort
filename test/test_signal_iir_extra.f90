! SPDX-License-Identifier: MIT
! Chebyshev-II and scalar order-selection/automatic IIR design regression.
program test_signal_iir_extra
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_signal, only : dp, cheby2, iirfilter, iirdesign, buttord, &
        cheb1ord, cheb2ord, iirnotch, iirpeak, butter_result, iir_order_result, &
        frequency_response, sosfreqz, freqz, design_status_success, &
        design_status_invalid_input
    implicit none
    type(butter_result) :: design, selected
    type(iir_order_result) :: spec
    type(frequency_response) :: ba_response, sections_response
    integer :: n, kind, i
    character(8), parameter :: modes(4) = &
        ['lowpass ', 'highpass', 'bandpass', 'bandstop']
    real(dp), parameter :: tol = 3.0e-10_dp
    real(dp) :: freq(9)

    design = cheby2(3,40.0_dp,0.4_dp)
    if (design%status /= design_status_success) error stop 'cheby2 lowpass status'
    call compare(design%b, [0.022849864710835058_dp,0.014905707366380064_dp, &
        0.014905707366380061_dp,0.022849864710835055_dp], tol)
    call compare(design%a, [1.0_dp,-2.0452520774397329_dp, &
        1.498781707830726_dp,-0.3780184862365627_dp], tol)

    design = cheby2(5,50.0_dp,0.3_dp,btype='highpass')
    if (design%status /= design_status_success) error stop 'cheby2 highpass status'
    call compare(design%b, [0.090697175573627026_dp,-0.35463290186388008_dp, &
        0.63312804153832047_dp,-0.63312804153832047_dp, &
        0.35463290186388008_dp,-0.090697175573627012_dp],tol)

    design = cheby2(4,35.0_dp,[0.2_dp,0.65_dp],btype='bandpass')
    if (design%status /= design_status_success) error stop 'cheby2 bandpass status'
    call compare(design%b, [0.054603487878345756_dp,-0.047365971256972145_dp, &
        -0.042163267254786606_dp,-0.0055514749561408503_dp, &
        0.099119780439311023_dp,-0.005551474956140964_dp, &
        -0.042163267254786557_dp,-0.04736597125697211_dp, &
        0.054603487878345749_dp],5.0e-10_dp)

    do kind = 1, 4
        do n = 1, 8
            if (kind <= 2) then
                design = cheby2(n, 45.0_dp, 0.37_dp, btype=trim(modes(kind)))
            else
                design = cheby2(n, 45.0_dp, [0.19_dp,0.68_dp], &
                    btype=trim(modes(kind)))
            end if
            if (design%status /= design_status_success) error stop 'cheby2 modes'
            if (.not. all(ieee_is_finite(design%b))) error stop 'cheby2 b finite'
            if (.not. all(ieee_is_finite(design%a))) error stop 'cheby2 a finite'
            if (.not. all(ieee_is_finite(design%sos))) error stop 'cheby2 sos finite'
            freq = [(real(i-1,dp)*0.061_dp,i=1,9)]
            ba_response = freqz(design%b,design%a,w=freq)
            sections_response = sosfreqz(design%sos,w=freq)
            if (ba_response%status /= design_status_success) error stop 'freqz status'
            if (sections_response%status /= design_status_success) error stop 'sosfreqz status'
            if (maxval(abs(ba_response%response-sections_response%response)) > 2.0e-6_dp) &
                error stop 'TF and SOS mismatch'
        end do
    end do
    design = cheby2(0,40.0_dp,0.3_dp)
    if (design%status /= design_status_invalid_input) error stop 'cheby2 invalid order'
    design = cheby2(3,-5.0_dp,0.3_dp)
    if (design%status /= design_status_invalid_input) error stop 'cheby2 negative attenuation'
    design = cheby2(3,40.0_dp,[0.7_dp,0.4_dp],btype='bandpass')
    if (design%status /= design_status_invalid_input) error stop 'cheby2 reversed edges'
    design = cheby2(3,40.0_dp,1.0_dp)
    if (design%status /= design_status_invalid_input) error stop 'cheby2 Nyquist edge'

    selected = iirfilter(4,0.27_dp,ftype='butter')
    if (selected%status /= design_status_success) error stop 'iirfilter butter'
    design = iirfilter(4,0.27_dp,ftype='cheby1',rp=1.0_dp)
    if (design%status /= design_status_success) error stop 'iirfilter cheby1'
    design = iirfilter(4,[0.27_dp,0.59_dp],btype='bandstop',ftype='cheby2',rs=40.0_dp)
    if (design%status /= design_status_success) error stop 'iirfilter cheby2'
    design = iirfilter(4,0.3_dp,ftype='cheby2')
    if (design%status /= design_status_invalid_input) error stop 'missing rs'
    design = iirfilter(4,0.3_dp,ftype='ellip',rs=40.0_dp)
    if (design%status /= design_status_invalid_input) error stop 'unsupported ellip'

    spec = buttord(0.2_dp,0.34_dp,1.0_dp,42.0_dp)
    call check_order(spec,10,0.21298874786440833_dp)
    spec = cheb1ord(0.2_dp,0.34_dp,1.0_dp,42.0_dp)
    call check_order(spec,6,0.2_dp)
    spec = cheb2ord(0.2_dp,0.34_dp,1.0_dp,42.0_dp)
    call check_order(spec,6,0.3025907645473806_dp)
    spec = buttord(0.7_dp,0.55_dp,1.0_dp,30.0_dp)
    if (spec%status /= design_status_success) error stop 'highpass buttord'
    spec = cheb2ord(0.2_dp,0.3_dp,6.0_dp,2.0_dp)
    if (spec%status /= design_status_invalid_input) error stop 'stop less than pass'
    spec = cheb1ord(0.2_dp,0.2_dp,1.0_dp,40.0_dp)
    if (spec%status /= design_status_invalid_input) error stop 'coincident edges'
    selected = iirdesign(0.2_dp,0.34_dp,1.0_dp,42.0_dp,ftype='cheby2')
    if (selected%status /= design_status_success) error stop 'iirdesign cheby2'
    selected = iirdesign(0.7_dp,0.52_dp,2.0_dp,35.0_dp,ftype='cheby1')
    if (selected%status /= design_status_success) error stop 'iirdesign cheby1 highpass'
    selected = iirdesign(0.2_dp,0.34_dp,1.0_dp,42.0_dp,ftype='ellip')
    if (selected%status /= design_status_success) error stop 'iirdesign elliptic family'

    design = iirnotch(60.0_dp,30.0_dp,fs=200.0_dp)
    if (design%status /= design_status_success) error stop 'notch status'
    sections_response = sosfreqz(design%sos,w=[60.0_dp],fs=200.0_dp)
    if (abs(sections_response%response(1)) > 5.0e-12_dp) error stop 'notch center'
    design = iirpeak(60.0_dp,30.0_dp,fs=200.0_dp)
    if (design%status /= design_status_success) error stop 'peak status'
    sections_response = sosfreqz(design%sos,w=[60.0_dp],fs=200.0_dp)
    if (abs(sections_response%response(1)-cmplx(1.0_dp,0.0_dp,dp)) > 5.0e-12_dp) &
        error stop 'peak center'
    design = iirnotch(10.0_dp,0.0_dp,fs=200.0_dp)
    if (design%status /= design_status_invalid_input) error stop 'invalid notch Q'

    print '(a)', 'test_signal_iir_extra passed'
contains
    subroutine compare(actual, reference, atol)
        real(dp), intent(in) :: actual(:) !! computed coefficients
        real(dp), intent(in) :: reference(:) !! SciPy reference
        real(dp), intent(in) :: atol !! absolute tolerance
        if (size(actual) /= size(reference)) error stop 'reference length'
        if (maxval(abs(actual-reference)) > atol) error stop 'coefficient difference'
    end subroutine compare

    subroutine check_order(out, n, cutoff)
        type(iir_order_result), intent(in) :: out !! selected order and critical frequency
        integer, intent(in) :: n !! SciPy reference order
        real(dp), intent(in) :: cutoff !! SciPy reference critical frequency
        if (out%status /= design_status_success) error stop 'order status'
        if (out%order /= n) error stop 'order mismatch'
        if (abs(out%cutoff-cutoff) > 1.0e-11_dp) error stop 'natural frequency'
    end subroutine check_order
end program test_signal_iir_extra
