! SPDX-License-Identifier: MIT
! Analog and digital order selection for scalar and two-edge IIR designs.
program test_signal_iir_order_bands
    use scifort_signal, only : dp, iir_order_result, butter_result, &
        frequency_response, buttord, cheb1ord, cheb2ord, iirdesign, &
        freqs, sosfreqz, design_status_success, design_status_invalid_input
    implicit none
    type(iir_order_result) :: choice
    type(butter_result) :: filt
    type(frequency_response) :: hpass, hstop

    ! The canonical analog Butterworth bandpass example from SciPy.
    choice = buttord([20.0_dp,50.0_dp], [14.0_dp,60.0_dp], &
        3.0_dp,40.0_dp,analog=.true.)
    call assert_order(choice,13,2)
    call assert_edges(choice%cutoffs,[19.99843436_dp,50.00391441_dp],2.0e-5_dp)

    choice = cheb1ord([10.0_dp,70.0_dp], [20.0_dp,50.0_dp], &
        3.0_dp,40.0_dp,analog=.true.)
    call assert_order(choice,5,2)
    call assert_edges(choice%cutoffs,[14.28571415_dp,69.99999362_dp],2.0e-5_dp)

    choice = cheb2ord([0.1_dp,0.6_dp], [0.2_dp,0.5_dp], &
        3.0_dp,60.0_dp)
    call assert_order(choice,7,2)
    call assert_edges(choice%cutoffs,[0.19776722_dp,0.50381287_dp],5.0e-5_dp)

    choice = buttord([0.2_dp,0.5_dp], [0.1_dp,0.6_dp], &
        3.0_dp,40.0_dp)
    call assert_order(choice,9,2)
    call assert_edges(choice%cutoffs,[0.19997485_dp,0.50004279_dp],2.0e-8_dp)

    choice = cheb1ord([0.2_dp,0.5_dp], [0.1_dp,0.6_dp], &
        2.0_dp,37.0_dp)
    call assert_order(choice,5,2)
    call assert_edges(choice%cutoffs,[0.2_dp,0.5_dp],1.0e-12_dp)

    choice = cheb2ord([0.2_dp,0.5_dp], [0.1_dp,0.6_dp], &
        2.0_dp,37.0_dp)
    call assert_order(choice,5,2)
    call assert_edges(choice%cutoffs,[0.1532788_dp,0.58802976_dp],2.0e-8_dp)

    choice = buttord(10.0_dp,20.0_dp,2.0_dp,40.0_dp,analog=.true.)
    call assert_order(choice,8,1)
    call assert_edges(choice%cutoffs,[10.340885147072871_dp],1.0e-10_dp)
    choice = cheb1ord(30.0_dp,10.0_dp,1.0_dp,40.0_dp,analog=.true.)
    call assert_order(choice,4,1)
    call assert_edges(choice%cutoffs,[30.0_dp],1.0e-10_dp)
    choice = cheb2ord(30.0_dp,10.0_dp,1.0_dp,40.0_dp,analog=.true.)
    call assert_order(choice,4,1)
    call assert_edges(choice%cutoffs,[12.828465344953258_dp],1.0e-10_dp)

    ! Spec-based analog bandstop design, checked at all four original edges.
    filt = iirdesign([10.0_dp,70.0_dp], [20.0_dp,50.0_dp], &
        3.0_dp,40.0_dp,ftype='cheby1',analog=.true.)
    if (filt%status /= design_status_success) error stop 'analog iirdesign bandstop'
    hpass = freqs(filt%b,filt%a,w=[10.0_dp,70.0_dp])
    hstop = freqs(filt%b,filt%a,w=[20.0_dp,50.0_dp])
    call check_specs(hpass,hstop,3.0_dp,40.0_dp)

    ! SOS avoids numerical loss in moderately high-order digital BA forms.
    filt = iirdesign([0.1_dp,0.6_dp], [0.2_dp,0.5_dp], &
        3.0_dp,60.0_dp,ftype='cheby2')
    if (filt%status /= design_status_success) error stop 'digital iirdesign bandstop'
    if (.not. allocated(filt%sos)) error stop 'missing digital SOS'
    hpass = sosfreqz(filt%sos,w=[0.1_dp,0.6_dp],fs=2.0_dp)
    hstop = sosfreqz(filt%sos,w=[0.2_dp,0.5_dp],fs=2.0_dp)
    call check_specs(hpass,hstop,3.0_dp,60.0_dp)

    filt = iirdesign(10.0_dp,20.0_dp,2.0_dp,40.0_dp, &
        ftype='butter',analog=.true.)
    if (filt%status /= design_status_success) error stop 'scalar analog iirdesign'

    choice = buttord([0.2_dp,0.6_dp], [0.1_dp],1.0_dp,40.0_dp)
    if (choice%status /= design_status_invalid_input) error stop 'unequal edge counts'
    choice = cheb1ord([0.2_dp,0.6_dp,0.8_dp],[0.1_dp,0.7_dp,0.9_dp],1.0_dp,40.0_dp)
    if (choice%status /= design_status_invalid_input) error stop 'three edges'
    choice = cheb2ord([0.2_dp,0.5_dp],[0.3_dp,0.7_dp],1.0_dp,40.0_dp)
    if (choice%status /= design_status_invalid_input) error stop 'interleaved edges'
    choice = buttord([20.0_dp,50.0_dp], [14.0_dp,60.0_dp], &
        3.0_dp,40.0_dp,analog=.true.,fs=100.0_dp)
    if (choice%status /= design_status_invalid_input) error stop 'analog fs invalid'
    choice = buttord(0.2_dp,0.4_dp,1.0_dp,40.0_dp,analog=.true.,fs=2.0_dp)
    if (choice%status /= design_status_invalid_input) error stop 'scalar analog fs invalid'
    choice = cheb2ord([0.1_dp,0.8_dp],[0.2_dp,0.6_dp],1.0_dp,40.0_dp,fs=1.0_dp)
    if (choice%status /= design_status_invalid_input) error stop 'digital Nyquist invalid'
    filt = iirdesign([20.0_dp,50.0_dp],[14.0_dp,60.0_dp], &
        3.0_dp,40.0_dp,analog=.true.,fs=100.0_dp)
    if (filt%status /= design_status_invalid_input) error stop 'analog iirdesign fs invalid'

    print '(a)', 'test_signal_iir_order_bands passed'
contains
    subroutine assert_order(selection, expected, edge_count)
        type(iir_order_result), intent(in) :: selection !! computed order and critical edges
        integer, intent(in) :: expected !! reference prototype order
        integer, intent(in) :: edge_count !! required number of edges
        if (selection%status /= design_status_success) error stop 'order status'
        if (selection%order /= expected) error stop 'order mismatch'
        if (.not. allocated(selection%cutoffs)) error stop 'cutoffs unallocated'
        if (size(selection%cutoffs) /= edge_count) error stop 'cutoff count'
        if (selection%cutoff /= selection%cutoffs(1)) error stop 'legacy cutoff mismatch'
    end subroutine assert_order

    subroutine assert_edges(actual, expected, tolerance)
        real(dp), intent(in) :: actual(:) !! computed cutoff frequencies
        real(dp), intent(in) :: expected(:) !! SciPy reference values
        real(dp), intent(in) :: tolerance !! absolute tolerance
        if (size(actual) /= size(expected)) error stop 'edge count'
        if (maxval(abs(actual-expected)) > tolerance) error stop 'edge mismatch'
    end subroutine assert_edges

    subroutine check_specs(pass_response, stop_response, pass_loss, stop_loss)
        type(frequency_response), intent(in) :: pass_response !! responses at passband boundaries
        type(frequency_response), intent(in) :: stop_response !! responses at stopband boundaries
        real(dp), intent(in) :: pass_loss !! allowed passband attenuation in dB
        real(dp), intent(in) :: stop_loss !! required stopband attenuation in dB
        real(dp) :: losses(2), attenuations(2)
        if (pass_response%status /= design_status_success) error stop 'pass response'
        if (stop_response%status /= design_status_success) error stop 'stop response'
        losses = -20.0_dp*log10(abs(pass_response%response))
        attenuations = -20.0_dp*log10(abs(stop_response%response))
        if (maxval(losses) > pass_loss+0.2_dp) error stop 'passband attenuation exceeded'
        if (minval(attenuations) < stop_loss-0.2_dp) error stop 'stopband attenuation insufficient'
    end subroutine check_specs
end program test_signal_iir_order_bands
