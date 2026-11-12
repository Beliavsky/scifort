! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Fixed SciPy 1.17.0 references for elliptic IIR design and order selection.
program test_signal_elliptic
    use scifort_signal, only : dp, ellip, ellipord, iirfilter, iirdesign, sosfreqz, &
        butter_result, iir_order_result, frequency_response, design_status_success, design_status_invalid_input
    implicit none
    type(butter_result) :: filt, expected
    type(iir_order_result) :: selection
    type(frequency_response) :: response

    filt = ellip(3,1.0_dp,40.0_dp,0.3_dp)
    call assert_coeff(filt, &
        [0.05548085891202_dp,0.09185345473382_dp,0.09185345473382_dp,0.05548085891202_dp], &
        [1.0_dp,-1.56894139395241_dp,1.26274246588815_dp,-0.39913244464406_dp])
    if (.not. allocated(filt%sos)) error stop 'missing digital SOS'
    if (size(filt%sos,1) /= 2) error stop 'odd-order SOS count'

    filt = ellip(4,1.0_dp,40.0_dp,0.3_dp,btype='highpass')
    call assert_coeff(filt, &
        [0.26680730348480_dp,-0.94818241520540_dp,1.37070681961636_dp, &
        -0.94818241520540_dp,0.26680730348480_dp], &
        [1.0_dp,-1.28814857598772_dp,1.31763796625209_dp, &
        -0.44624167493511_dp,0.21241190218729_dp])

    filt = ellip(3,1.0_dp,45.0_dp,[0.2_dp,0.55_dp],btype='bandpass')
    call assert_coeff(filt, &
        [0.0674427434972181_dp,-0.0237278188230936_dp,-0.138811757837852_dp,0.0_dp, &
        0.138811757837852_dp,0.0237278188230936_dp,-0.0674427434972181_dp], &
        [1.0_dp,-1.91701304573160_dp,2.59849851528030_dp,-2.41134601491775_dp, &
        1.93061891920096_dp,-0.92954001137792_dp,0.34151163407316_dp])

    filt = ellip(3,1.0_dp,45.0_dp,[0.2_dp,0.55_dp],btype='bandstop')
    call assert_coeff(filt, &
        [0.28421640756081_dp,-0.74841990393189_dp,1.47907674246370_dp, &
        -1.66166920674474_dp,1.47907674246370_dp,-0.74841990393189_dp, &
        0.28421640756081_dp], &
        [1.0_dp,-1.64517834388501_dp,1.65686988290951_dp,-1.32581698382817_dp, &
        0.92713480130227_dp,-0.18751368689533_dp,-0.05741838416275_dp])

    filt = ellip(4,2.0_dp,45.0_dp,10.0_dp,analog=.true.)
    call assert_coeff(filt, &
        [0.00562341325190348_dp,0.0_dp,9.13720896123334_dp,0.0_dp,2098.19380407035_dp], &
        [1.0_dp,7.08643779715383_dp,131.714830526545_dp,561.715319397038_dp, &
        2641.46949881324_dp])
    if (allocated(filt%sos)) error stop 'analog SOS must be absent'

    filt = ellip(3,1.0_dp,40.0_dp,[20.0_dp,50.0_dp],btype='bandstop',analog=.true.)
    call assert_coeff(filt, &
        [1.0_dp,0.0_dp,3118.28940894583_dp,0.0_dp,3118289.40894583_dp, &
        0.0_dp,1000000000.0_dp], &
        [1.0_dp,70.8455913508225_dp,4672.15330773301_dp,192971.614792078_dp, &
        4672153.30773301_dp,70845591.3508225_dp,1000000000.0_dp])

    selection = ellipord(0.2_dp,0.3_dp,1.0_dp,40.0_dp)
    call assert_order(selection,4,[0.2_dp])
    selection = ellipord([0.2_dp,0.5_dp],[0.1_dp,0.65_dp],1.0_dp,40.0_dp)
    call assert_order(selection,4,[0.2_dp,0.5_dp])
    selection = ellipord([20.0_dp,50.0_dp],[14.0_dp,60.0_dp],3.0_dp,40.0_dp,analog=.true.)
    call assert_order(selection,4,[20.0_dp,50.0_dp])
    selection = ellipord([20.0_dp,50.0_dp],[30.0_dp,40.0_dp],1.0_dp,40.0_dp,analog=.true.)
    if (selection%status /= design_status_success .or. selection%order /= 3) error stop 'bandstop order'

    expected = ellip(4,1.0_dp,40.0_dp,0.2_dp)
    filt = iirfilter(4,0.2_dp,ftype='ellip',rp=1.0_dp,rs=40.0_dp)
    call assert_coeff(filt,expected%b,expected%a)
    filt = ellip(4,1.0_dp,40.0_dp,10.0_dp,fs=100.0_dp)
    call assert_coeff(filt,expected%b,expected%a)
    response = sosfreqz(filt%sos,w=[0.0_dp],fs=100.0_dp)
    if (response%status /= design_status_success) error stop 'digital SOS response'
    if (abs(abs(response%response(1))-10.0_dp**(-1.0_dp/20.0_dp)) > 2.0e-12_dp) &
        error stop 'even elliptic DC gain'
    filt = iirdesign(0.2_dp,0.3_dp,1.0_dp,40.0_dp,ftype='elliptic')
    call assert_coeff(filt,expected%b,expected%a)
    response = sosfreqz(filt%sos,w=[0.2_dp,0.3_dp],fs=2.0_dp)
    if (response%status /= design_status_success) error stop 'minimum-order SOS response'
    if (-20.0_dp*log10(abs(response%response(1))) > 1.001_dp) error stop 'passband specification'
    if (-20.0_dp*log10(abs(response%response(2))) < 39.999_dp) error stop 'stopband specification'
    filt = iirdesign([20.0_dp,50.0_dp],[14.0_dp,60.0_dp], &
        3.0_dp,40.0_dp,ftype='ellip',analog=.true.)
    if (filt%status /= design_status_success) error stop 'analog iirdesign'
    filt = iirfilter(3,[0.2_dp,0.4_dp],btype='bandpass',ftype='cauer', &
        rp=1.0_dp,rs=45.0_dp)
    if (filt%status /= design_status_success) error stop 'iirfilter aliases'

    filt = ellip(3,0.0_dp,40.0_dp,0.2_dp)
    if (filt%status /= design_status_invalid_input) error stop 'zero ripple'
    filt = ellip(3,3.0_dp,2.0_dp,0.2_dp)
    if (filt%status /= design_status_invalid_input) error stop 'inverted attenuation'
    filt = ellip(21,1.0_dp,40.0_dp,0.2_dp)
    if (filt%status /= design_status_invalid_input) error stop 'excessive order'
    filt = ellip(4,1.0_dp,40.0_dp,[0.2_dp],btype='bandpass')
    if (filt%status /= design_status_invalid_input) error stop 'missing edge'
    filt = ellip(4,1.0_dp,40.0_dp,0.2_dp,analog=.true.,fs=2.0_dp)
    if (filt%status /= design_status_invalid_input) error stop 'analog fs'
    filt = iirfilter(4,0.2_dp,ftype='ellip',rp=1.0_dp)
    if (filt%status /= design_status_invalid_input) error stop 'missing rs'
    selection = ellipord(0.2_dp,0.3_dp,1.0_dp,40.0_dp,analog=.true.,fs=2.0_dp)
    if (selection%status /= design_status_invalid_input) error stop 'invalid order fs'
    selection = ellipord(0.2_dp,0.3_dp,4.0_dp,3.0_dp)
    if (selection%status /= design_status_invalid_input) error stop 'invalid order loss'

    print '(a)', 'test_signal_elliptic passed'
contains
    subroutine assert_coeff(output, b, a)
        type(butter_result), intent(in) :: output !! filter under test
        real(dp), intent(in) :: b(:) !! reference numerator coefficients
        real(dp), intent(in) :: a(:) !! reference denominator coefficients
        real(dp), parameter :: tolerance = 3.0e-10_dp
        if (output%status /= design_status_success) error stop 'filter design status'
        if (size(output%b) /= size(b) .or. size(output%a) /= size(a)) error stop 'coefficient shapes'
        if (maxval(abs(output%b-b)/(1.0_dp+abs(b))) > tolerance) error stop 'numerator coefficients'
        if (maxval(abs(output%a-a)/(1.0_dp+abs(a))) > tolerance) error stop 'denominator coefficients'
    end subroutine assert_coeff

    subroutine assert_order(output, n, cutoff)
        type(iir_order_result), intent(in) :: output !! calculated order/critical edges
        integer, intent(in) :: n !! expected order
        real(dp), intent(in) :: cutoff(:) !! reference critical frequency(s)
        if (output%status /= design_status_success) error stop 'order selection status'
        if (output%order /= n) error stop 'selected filter order'
        if (.not. allocated(output%cutoffs)) error stop 'missing cutoffs'
        if (size(output%cutoffs) /= size(cutoff)) error stop 'cutoffs length'
        if (maxval(abs(output%cutoffs-cutoff)) > 1.0e-8_dp) error stop 'cutoff values'
    end subroutine assert_order
end program test_signal_elliptic
