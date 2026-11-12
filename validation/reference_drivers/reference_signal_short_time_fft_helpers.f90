! SPDX-License-Identifier: MIT
! Exercise ShortTimeFFT coordinate and boundary helper methods for SciPy comparison.
program reference_signal_short_time_fft_helpers
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result, &
        short_fft_success, short_fft_invalid_input
    implicit none
    type(short_time_fft) :: sft
    type(short_time_fft_result) :: frame
    character(20) :: name, mode
    real(dp), allocatable :: freq(:), times(:)
    real(dp), allocatable :: x(:)
    real(dp) :: bounds(4), fs
    integer :: uin, uout, ncases, ic, nw, hop, mfft, n, off
    integer :: first, last, status, extent_status, b0(2), b1(2)
    integer :: kleft, kright, l0, l1, i

    open(newunit=uin,file='stft_helper_cases.txt',status='old',action='read')
    open(newunit=uout,file='stft_helper_results.txt',status='replace',action='write')
    read(uin,*) ncases
    do ic=1,ncases
        read(uin,*) name,nw,hop,mfft,n,fs,mode,off
        call sft%initialize_window(trim(name),nw,hop,fs,mfft=mfft,fft_mode=trim(mode), &
            scale_to='psd')
        if (.not. sft%invertible()) error stop 'helpers not invertible'
        if (sft%status /= short_fft_success) error stop 'bad helper test plan'
        if (abs(sft%delta_t()-real(hop,dp)/fs) > 1.e-14_dp) error stop 'helper delta_t'
        if (abs(sft%delta_f()-fs/real(mfft,dp)) > 1.e-14_dp) error stop 'helper delta_f'
        b0 = sft%lower_border_end()
        b1 = sft%upper_border_begin(n)
        kleft = sft%nearest_k_p(-7)
        kright = sft%nearest_k_p(-7,left=.false.)
        call sft%p_range(n,first,last,status)
        if (status /= short_fft_success) error stop 'bad p_range'
        if (first /= sft%p_min .or. last /= sft%p_max(n)) error stop 'p_range default'
        call sft%p_range(n,l0,l1,status,p0=first-1)
        if (status /= short_fft_invalid_input) error stop 'p_range invalid accepted'
        if (last-first > 2) then
            call sft%p_range(n,l0,l1,status,p0=first+1,p1=last-1)
            if (status /= short_fft_success) error stop 'p_range restricted'
            if (l0 /= first+1 .or. l1 /= last-1) error stop 'p_range values'
        end if
        freq = sft%frequency_points()
        times = sft%time_points(n,status,k_offset=off)
        if (status /= short_fft_success) error stop 'time_points invalid'
        if (size(freq) /= sft%frequency_count()) error stop 'frequency_count'
        if (size(times) /= sft%p_num(n)) error stop 'time count'
        allocate(x(n))
        x = [(sin(real(i,dp)*0.13_dp),i=1,n)]
        call sft%stft(x,frame,k_offset=off)
        if (frame%status /= short_fft_success) error stop 'STFT metadata'
        if (maxval(abs(frame%times-times)) > 1.e-13_dp) error stop 'times mismatch'
        if (maxval(abs(frame%frequencies-freq)) > 1.e-13_dp) error stop 'freq mismatch'
        call sft%extent(n,bounds,extent_status,center_bins=.true.)
        if (mode == 'twosided') then
            if (extent_status /= short_fft_invalid_input) error stop 'twosided extent accepted'
        else
            if (extent_status /= short_fft_success) error stop 'extent invalid'
        end if
        write(uout,*) ic,sft%p_min,sft%p_max(n),sft%k_min(),sft%k_max(n), &
            b0,b1,kleft,kright,sft%frequency_count(),first,last,extent_status
        write(uout,'(6(es26.17,1x))') sft%delta_t(),sft%delta_f(),bounds
        write(uout,*) size(freq),size(times)
        write(uout,'(*(es26.17,1x))') freq
        write(uout,'(*(es26.17,1x))') times
        deallocate(x)
    end do
    close(uin)
    close(uout)
    print '(a,i0)', 'STFT helper reference cases: ',ncases
end program reference_signal_short_time_fft_helpers
