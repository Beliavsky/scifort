! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! A single-case reference driver for independent multichannel SciPy checks.
! Header: mode nchannels n nperseg noverlap nfft complex_flag pad_flag fs
! Then: boundary scaling detrend window
! Forward input is n-by-nchannels samples, channel by channel.
! Inverse input is n frames of nfreq-by-nchannels bins, frame by frame.
! All complex values are provided as real/imaginary pairs.
program reference_signal_stream_multi
    use scifort_signal, only : dp, stft_stream_multi, stft_stream_multi_init, &
        stft_stream_multi_push, stft_stream_multi_finish, stft_stream_multi_frequencies, &
        spectrogram_stream_multi_init, spectrogram_stream_multi_push, &
        spectrogram_stream_multi_finish, istft_stream_multi, istft_stream_multi_init, &
        istft_stream_multi_push, istft_stream_multi_finish, &
        spectral_status_success, istft_status_success, istft_status_noninvertible
    implicit none
    type(stft_stream_multi) :: forward
    type(istft_stream_multi) :: inverse
    integer :: mode, nch, n, nw, no, nf, complex_flag, pad_flag, st, i, j, k, first, last
    integer :: nfreq, nframes, steps
    real(dp) :: fs, re, im
    character(32) :: edge, scale, trend, win
    real(dp), allocatable :: xr(:,:), frequency(:), time(:), part_times(:), power(:,:,:)
    real(dp), allocatable :: part_real(:,:)
    complex(dp), allocatable :: xc(:,:), frames(:,:,:), part(:,:,:), inverse_frames(:,:,:)
    complex(dp), allocatable :: reconstructed(:,:), part_complex(:,:)

    read(*,*) mode,nch,n,nw,no,nf,complex_flag,pad_flag,fs
    read(*,*) edge,scale,trend,win
    allocate(time(0))
    if (mode==0 .or. mode==1) then
        allocate(xr(n,nch),xc(n,nch))
        do j=1,nch
            do i=1,n
                read(*,*) re,im
                xr(i,j)=re
                xc(i,j)=cmplx(re,im,dp)
            end do
        end do
        if (mode==0) then
            call stft_stream_multi_init(forward,nch,nw,st,fs=fs,window=trim(win), &
                noverlap=no,nfft=nf,detrend=trim(trend),boundary=trim(edge), &
                padded=pad_flag==1,scaling=trim(scale),input_complex=complex_flag==1)
        else
            call spectrogram_stream_multi_init(forward,nch,nw,st,fs=fs,window=trim(win), &
                noverlap=no,nfft=nf,detrend=trim(trend),scaling=trim(scale), &
                input_complex=complex_flag==1)
        end if
        if (st/=spectral_status_success) error stop 'multichannel forward init'
        call stft_stream_multi_frequencies(forward,frequency,st)
        if (st/=spectral_status_success) error stop 'multichannel frequency query'
        allocate(frames(size(frequency),nch,0))
        first=1
        steps=0
        do while(first<=n)
            steps=steps+1
            last=min(n,first+mod(steps*13+nw*7+no,17))
            if (mode==0) then
                if (complex_flag==1) then
                    call stft_stream_multi_push(forward,xc(first:last,:),part,st,part_times)
                else
                    call stft_stream_multi_push(forward,xr(first:last,:),part,st,part_times)
                end if
            else
                if (complex_flag==1) then
                    call spectrogram_stream_multi_push(forward,xc(first:last,:),power,st,part_times)
                else
                    call spectrogram_stream_multi_push(forward,xr(first:last,:),power,st,part_times)
                end if
                part=cmplx(power,0.0_dp,dp)
            end if
            if (st/=spectral_status_success) error stop 'multichannel forward push'
            call append_frames(frames,part)
            time=[time,part_times]
            first=last+1
        end do
        if (mode==0) then
            call stft_stream_multi_finish(forward,part,st,part_times)
        else
            call spectrogram_stream_multi_finish(forward,power,st,part_times)
            part=cmplx(power,0.0_dp,dp)
        end if
        if (st/=spectral_status_success) error stop 'multichannel forward finish'
        call append_frames(frames,part)
        time=[time,part_times]
        nfreq=size(frames,1)
        nframes=size(frames,3)
        print '(i0,1x,i0,1x,i0)',nfreq,nframes,nch
        do i=1,nfreq
            print '(es26.17)',frequency(i)
        end do
        do k=1,size(time)
            print '(es26.17)',time(k)
        end do
        do k=1,nframes
            do j=1,nch
                do i=1,nfreq
                    print '(2(es26.17,1x))',real(frames(i,j,k),dp),aimag(frames(i,j,k))
                end do
            end do
        end do
    else if (mode==2 .or. mode==3) then
        nfreq=nf
        nframes=n
        allocate(inverse_frames(nfreq,nch,nframes))
        do k=1,nframes
            do j=1,nch
                do i=1,nfreq
                    read(*,*) re,im
                    inverse_frames(i,j,k)=cmplx(re,im,dp)
                end do
            end do
        end do
        call istft_stream_multi_init(inverse,nch,nfreq,st,fs=fs,window=trim(win), &
            nperseg=nw,noverlap=no,nfft=complex_flag,boundary=pad_flag==1, &
            scaling=trim(scale),input_onesided=mode==2)
        if (st/=istft_status_success) error stop 'multichannel inverse init'
        allocate(reconstructed(0,nch))
        do k=1,nframes
            if (mode==2) then
                call istft_stream_multi_push(inverse,inverse_frames(:,:,k),part_real,st,part_times)
                part_complex=cmplx(part_real,0.0_dp,dp)
            else
                call istft_stream_multi_push(inverse,inverse_frames(:,:,k),part_complex,st,part_times)
            end if
            if (st/=istft_status_success .and. st/=istft_status_noninvertible) &
                error stop 'multichannel inverse push'
            call append_samples(reconstructed,part_complex)
            time=[time,part_times]
        end do
        if (mode==2) then
            call istft_stream_multi_finish(inverse,part_real,st,part_times)
            part_complex=cmplx(part_real,0.0_dp,dp)
        else
            call istft_stream_multi_finish(inverse,part_complex,st,part_times)
        end if
        if (st/=istft_status_success .and. st/=istft_status_noninvertible) &
            error stop 'multichannel inverse finish'
        call append_samples(reconstructed,part_complex)
        time=[time,part_times]
        print '(i0,1x,i0,1x,i0)',size(reconstructed,1),size(reconstructed,2),0
        do k=1,size(time)
            print '(es26.17)',time(k)
        end do
        do j=1,nch
            do i=1,size(reconstructed,1)
                print '(2(es26.17,1x))',real(reconstructed(i,j),dp),aimag(reconstructed(i,j))
            end do
        end do
    else
        error stop 'invalid driver mode'
    end if
contains

    subroutine append_frames(all_frames,new_frames)
        complex(dp), allocatable, intent(inout) :: all_frames(:,:,:) !! saved frames
        complex(dp), intent(in) :: new_frames(:,:,:) !! newly completed frames
        complex(dp), allocatable :: work(:,:,:)
        integer :: nold, nadd
        nold=size(all_frames,3)
        nadd=size(new_frames,3)
        allocate(work(size(all_frames,1),size(all_frames,2),nold+nadd))
        work(:,:,1:nold)=all_frames
        work(:,:,nold+1:nold+nadd)=new_frames
        call move_alloc(work,all_frames)
    end subroutine append_frames

    subroutine append_samples(all_samples,new_samples)
        complex(dp), allocatable, intent(inout) :: all_samples(:,:) !! saved samples
        complex(dp), intent(in) :: new_samples(:,:) !! newly completed samples
        complex(dp), allocatable :: work(:,:)
        integer :: nold, nadd
        nold=size(all_samples,1)
        nadd=size(new_samples,1)
        allocate(work(nold+nadd,size(all_samples,2)))
        work(1:nold,:)=all_samples
        work(nold+1:nold+nadd,:)=new_samples
        call move_alloc(work,all_samples)
    end subroutine append_samples
end program reference_signal_stream_multi
