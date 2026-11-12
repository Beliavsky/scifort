! SPDX-License-Identifier: MIT
! Native ShortTimeFFT support-grid and coordinate regression checks.
program test_signal_short_time_fft_helpers
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result, &
        short_fft_success, short_fft_invalid_input, short_fft_noninvertible
    implicit none
    type(short_time_fft) :: sft
    type(short_time_fft_result) :: frame
    integer :: status, first, last, bounds(2), n, i
    real(dp), allocatable :: f(:), t(:)
    real(dp) :: extent(4)
    real(dp) :: x(32)

    if (sft%invertible()) error stop 'default plan is invertible'
    if (sft%frequency_count() /= 0) error stop 'uninitialized frequency count'
    call sft%initialize_window('hann',8,8,12.0_dp)
    if (sft%status /= short_fft_noninvertible) error stop 'noninvertible plan not detected'
    call sft%initialize_window('hann',8,3,12.0_dp,mfft=8,fft_mode='twosided')
    if (.not. sft%invertible()) error stop 'twosided plan'
    if (abs(sft%delta_f()-1.5_dp) > 1.e-13_dp) error stop 'frequency step'
    if (abs(sft%delta_t()-0.25_dp) > 1.e-13_dp) error stop 'time step'
    if (sft%frequency_count() /= 8) error stop 'twosided frequency count'
    f = sft%frequency_points()
    if (maxval(abs(f-[0.0_dp,1.5_dp,3.0_dp,4.5_dp, &
        -6.0_dp,-4.5_dp,-3.0_dp,-1.5_dp])) > 1.e-13_dp) &
        error stop 'negative Nyquist frequency'
    call sft%extent(32,extent,status)
    if (status /= short_fft_invalid_input) error stop 'twosided extent accepted'
    if (sft%nearest_k_p(-7) /= -9) error stop 'negative floor'
    if (sft%nearest_k_p(-7,left=.false.) /= -6) error stop 'negative ceiling'
    if (sft%nearest_k_p(6,left=.false.) /= 6) error stop 'exact hop'
    if (sft%k_min() > 0 .or. sft%k_max(32) < 32) error stop 'signal coverage'
    bounds = sft%lower_border_end()
    if (bounds(2) < 0) error stop 'lower border'
    bounds = sft%upper_border_begin(32)
    if (bounds(2) > sft%p_max(32)) error stop 'upper border'
    call sft%p_range(32,first,last,status)
    if (status /= short_fft_success) error stop 'p_range'
    if (last-first /= sft%p_num(32)) error stop 'p_num'
    call sft%p_range(32,first,last,status,p1=sft%p_min)
    if (status /= short_fft_invalid_input) error stop 'invalid p_range'
    t = sft%time_points(32,status,k_offset=5)
    if (status /= short_fft_success) error stop 'time_points'
    x = [(sin(real(i,dp)*0.14_dp),i=1,32)]
    call sft%stft(x,frame,k_offset=5)
    if (frame%status /= short_fft_success) error stop 'STFT with offset'
    if (maxval(abs(t-frame%times)) > 2.e-13_dp) error stop 'offset timestamps'
    if (maxval(abs(f-frame%frequencies)) > 2.e-13_dp) error stop 'frequency metadata'
    call sft%initialize_window('boxcar',7,2,10.0_dp,mfft=11,fft_mode='centered')
    if (sft%status /= short_fft_success) error stop 'centered plan'
    f = sft%frequency_points()
    if (abs(f(1)+5.0_dp*10.0_dp/11.0_dp) > 1.e-13_dp) &
        error stop 'odd centered start'
    if (abs(f(11)-5.0_dp*10.0_dp/11.0_dp) > 1.e-13_dp) &
        error stop 'odd centered end'
    call sft%extent(32,extent,status,axes_seq='ft',center_bins=.true.)
    if (status /= short_fft_success) error stop 'centered extent'
    if (extent(1) >= 0.0_dp .or. extent(2) <= 0.0_dp) error stop 'extent spectrum'
    call sft%extent(32,extent,status,axes_seq='bad')
    if (status /= short_fft_invalid_input) error stop 'invalid extent axes'
    print *, 'test_signal_short_time_fft_helpers: all checks passed'
end program test_signal_short_time_fft_helpers
