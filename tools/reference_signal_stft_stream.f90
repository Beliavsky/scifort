! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Streaming forward time-frequency reference protocol, one case per process.
program reference_signal_stft_stream
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_stft_stream, only : stft_stream, stft_stream_init, &
        stft_stream_push, stft_stream_finish, stft_stream_frequencies, &
        spectrogram_stream_init, spectrogram_stream_push, spectrogram_stream_finish
    implicit none
    type(stft_stream) :: state
    integer :: mode, n, nw, no, nf, complex_flag, pad_flag, st, i, j, k
    integer :: b, steps
    real(dp) :: fs, re, im
    character(32) :: edge, scale, trend, win
    real(dp), allocatable :: xr(:), freq(:), time(:), times(:), part_times(:)
    real(dp), allocatable :: all_power(:,:), part_power(:,:)
    complex(dp), allocatable :: xc(:), all_spectra(:,:), part_spectra(:,:)

    read(*,*) mode,n,nw,no,nf,complex_flag,pad_flag,fs
    read(*,*) edge,scale,trend,win
    allocate(xr(n),xc(n),time(0))
    do i=1,n
        read(*,*) re,im
        xr(i)=re
        xc(i)=cmplx(re,im,dp)
    end do
    if (mode==0 .and. win=='default') then
        call stft_stream_init(state,nw,st,fs=fs,noverlap=no,nfft=nf,detrend=trim(trend), &
            boundary=trim(edge),padded=pad_flag==1,scaling=trim(scale), &
            input_complex=complex_flag==1)
    else if (mode==1 .and. win=='default') then
        call spectrogram_stream_init(state,nw,st,fs=fs,noverlap=no,nfft=nf, &
            detrend=trim(trend),scaling=trim(scale),input_complex=complex_flag==1)
    else if (mode==0) then
        call stft_stream_init(state,nw,st,fs=fs,window=trim(win),noverlap=no, &
            nfft=nf,detrend=trim(trend),boundary=trim(edge),padded=pad_flag==1, &
            scaling=trim(scale),input_complex=complex_flag==1)
    else
        call spectrogram_stream_init(state,nw,st,fs=fs,window=trim(win),noverlap=no, &
            nfft=nf,detrend=trim(trend),scaling=trim(scale),input_complex=complex_flag==1)
    end if
    if (st/=spectral_status_success) error stop 'stream init'
    call stft_stream_frequencies(state,freq,st)
    if (st/=spectral_status_success) error stop 'frequency query'
    allocate(all_spectra(size(freq),0),all_power(size(freq),0))
    b=1
    steps=0
    do while(b<=n)
        steps=steps+1
        k=min(n,b+mod(steps*13+7*nw+no,17))
        if (mode==0) then
            if (complex_flag==1) then
                call stft_stream_push(state,xc(b:k),part_spectra,st,part_times)
            else
                call stft_stream_push(state,xr(b:k),part_spectra,st,part_times)
            end if
            if (st/=spectral_status_success) error stop 'STFT push'
            all_spectra=reshape([reshape(all_spectra,[size(all_spectra)]), &
                reshape(part_spectra,[size(part_spectra)])], &
                [size(freq),size(all_spectra,2)+size(part_spectra,2)])
        else
            if (complex_flag==1) then
                call spectrogram_stream_push(state,xc(b:k),part_power,st,part_times)
            else
                call spectrogram_stream_push(state,xr(b:k),part_power,st,part_times)
            end if
            if (st/=spectral_status_success) error stop 'spectrogram push'
            all_power=reshape([reshape(all_power,[size(all_power)]), &
                reshape(part_power,[size(part_power)])], &
                [size(freq),size(all_power,2)+size(part_power,2)])
        end if
        time=[time,part_times]
        b=k+1
    end do
    if (mode==0) then
        call stft_stream_finish(state,part_spectra,st,part_times)
        if (st/=spectral_status_success) error stop 'STFT finish'
        all_spectra=reshape([reshape(all_spectra,[size(all_spectra)]), &
            reshape(part_spectra,[size(part_spectra)])], &
            [size(freq),size(all_spectra,2)+size(part_spectra,2)])
    else
        call spectrogram_stream_finish(state,part_power,st,part_times)
        if (st/=spectral_status_success) error stop 'spectrogram finish'
        all_power=reshape([reshape(all_power,[size(all_power)]), &
            reshape(part_power,[size(part_power)])], &
            [size(freq),size(all_power,2)+size(part_power,2)])
    end if
    time=[time,part_times]
    if (mode==0) then
        print '(i0,1x,i0)',size(freq),size(all_spectra,2)
    else
        print '(i0,1x,i0)',size(freq),size(all_power,2)
    end if
    do i=1,size(freq)
        print '(es26.17)',freq(i)
    end do
    do j=1,size(time)
        print '(es26.17)',time(j)
    end do
    do j=1,size(time)
        do i=1,size(freq)
            if (mode==0) then
                print '(2(es26.17,1x))',real(all_spectra(i,j),dp),aimag(all_spectra(i,j))
            else
                print '(2(es26.17,1x))',all_power(i,j),0.0_dp
            end if
        end do
    end do
end program reference_signal_stft_stream
