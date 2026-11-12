! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_stft_stream
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, stft, spectrogram, stft_result, spectrogram_result, &
        stft_stream, stft_stream_init, stft_stream_push, stft_stream_finish, &
        stft_stream_frequencies, get_window, spectrogram_stream_init, spectrogram_stream_push, &
        spectrogram_stream_finish, spectral_status_success, spectral_status_invalid_input, &
        istft_stream, istft_stream_init, istft_stream_push, istft_stream_finish
    implicit none
    integer, parameter :: window_lengths(*) = [3,8,9,16]
    integer :: nw, no, nf, n, b, s, d, t, mode, j, cases, st
    integer :: nframes
    character(8) :: scaling, trend, edge
    logical :: pad
    real(dp), allocatable :: real_values(:), power(:,:), times(:), frequency(:)
    complex(dp), allocatable :: complex_values(:), spectra(:,:)
    type(stft_stream) :: state

    cases = 0
    do nw=3,9
        do no=0,nw-1,max(1,nw/2)
            do b=0,1
                edge='none'
                if (b==1) edge='zeros'
                do s=0,1
                    scaling='spectrum'
                    if (s==1) scaling='psd'
                    do d=0,2
                        trend='none'
                        if (d==1) trend='constant'
                        if (d==2) trend='linear'
                        do t=0,1
                            pad=t==1
                            do mode=0,1
                                nf=nw+mod(nw+mode,3)
                                n=3*nw+2*no+1
                                call check_stft(n,nw,no,nf,trim(edge),pad,trim(scaling), &
                                    trim(trend),mode==1)
                                cases=cases+1
                            end do
                        end do
                    end do
                end do
            end do
            do s=0,1
                scaling='density'
                if (s==1) scaling='spectrum'
                do d=0,2
                    trend='none'
                    if (d==1) trend='constant'
                    if (d==2) trend='linear'
                    do mode=0,1
                        nf=nw+mod(nw+mode,3)
                        n=3*nw+2*no+1
                        call check_spectrogram(n,nw,no,nf,trim(scaling),trim(trend), &
                            mode==1,window='boxcar')
                        cases=cases+1
                    end do
                end do
            end do
        end do
    end do

    ! Default windows, explicit Tukey names, and custom weights must agree with batch.
    do j=1,size(window_lengths)
        nw=window_lengths(j)
        do s=0,1
            scaling='density'
            if (s==1) scaling='spectrum'
            do mode=0,1
                call check_spectrogram(4*nw+5,nw,nw/4,nw+3,trim(scaling),'constant', &
                    mode==1)
                call check_spectrogram(4*nw+5,nw,nw/4,nw+3,trim(scaling),'constant', &
                    mode==1,window='tukey')
                call check_spectrogram(4*nw+5,nw,nw/4,nw+3,trim(scaling),'constant', &
                    mode==1,window_values=get_window('tukey',nw,parameter=0.25_dp))
                cases=cases+3
            end do
        end do
    end do

    call stft_stream_init(state,8,st,window='boxcar',noverlap=4)
    if (st/=spectral_status_success) error stop 'valid initialization'
    call stft_stream_init(state,8,st,noverlap=8)
    if (st/=spectral_status_invalid_input) error stop 'invalid overlap not rejected'
    call stft_stream_init(state,0,st)
    if (st/=spectral_status_invalid_input) error stop 'zero length accepted'
    call stft_stream_init(state,8,st,nfft=7)
    if (st/=spectral_status_invalid_input) error stop 'short transform accepted'
    call stft_stream_init(state,8,st,fs=0.0_dp)
    if (st/=spectral_status_invalid_input) error stop 'zero sampling rate accepted'
    call stft_stream_init(state,8,st,window='hann',window_values=[(1.0_dp,j=1,8)])
    if (st/=spectral_status_invalid_input) error stop 'two windows accepted'
    call stft_stream_init(state,8,st,scaling='density')
    if (st/=spectral_status_invalid_input) error stop 'bad STFT scaling accepted'
    call stft_stream_init(state,8,st,detrend='bad')
    if (st/=spectral_status_invalid_input) error stop 'bad detrend accepted'
    call stft_stream_init(state,8,st,window_values=[(0.0_dp,j=1,8)])
    if (st/=spectral_status_invalid_input) error stop 'zero-normalization accepted'
    call stft_stream_init(state,8,st,window_values=[(1.0_dp,j=1,7)])
    if (st/=spectral_status_invalid_input) error stop 'bad window length accepted'
    call stft_stream_frequencies(state,frequency,st)
    if (st/=spectral_status_success .or. size(frequency)/=5) error stop 'failed init mutated stream'
    call stft_stream_push(state,[ieee_value(1.0_dp,ieee_quiet_nan)],spectra,st)
    if (st/=spectral_status_invalid_input .or. size(spectra)/=0) error stop 'invalid chunk accepted'
    call stft_stream_push(state,[(1.0_dp,j=1,7)],spectra,st,times)
    if (st/=spectral_status_success .or. size(times)/=1) error stop 'rejected chunk advanced stream'
    call stft_stream_push(state,[(cmplx(1.0_dp,0.0_dp,dp),j=1,2)],spectra,st)
    if (st/=spectral_status_invalid_input) error stop 'complex chunk in real stream'
    call stft_stream_finish(state,spectra,st)
    if (st/=spectral_status_success) error stop 'finish after rejected chunks'
    call stft_stream_finish(state,spectra,st)
    if (st/=spectral_status_invalid_input .or. size(spectra)/=0) error stop 'repeat finish accepted'
    call stft_stream_push(state,[1.0_dp],spectra,st)
    if (st/=spectral_status_invalid_input) error stop 'push after finish accepted'

    call spectrogram_stream_init(state,4,st,window='boxcar',noverlap=2,input_complex=.true.)
    if (st/=spectral_status_success) error stop 'complex spectrogram init'
    call stft_stream_frequencies(state,frequency,st)
    if (st/=spectral_status_success .or. size(frequency)/=4) error stop 'two-sided frequency count'
    if (abs(frequency(4)+0.25_dp)>1e-12_dp) error stop 'negative frequency bin'
    call spectrogram_stream_push(state,[1.0_dp,2.0_dp],power,st)
    if (st/=spectral_status_invalid_input) error stop 'real chunk in complex spectrogram'
    call spectrogram_stream_push(state,[(cmplx(1.0_dp,0.2_dp,dp),j=1,4)],power,st)
    if (st/=spectral_status_success .or. size(power,2)/=1) error stop 'valid complex spectrogram'
    call spectrogram_stream_finish(state,power,st)
    if (st/=spectral_status_success) error stop 'finish complex spectrogram'

    ! Custom windows must work without naming a window family.
    call stft_stream_init(state,5,st,window_values=[1.0_dp,2.0_dp,3.0_dp,2.0_dp,1.0_dp], &
        noverlap=2,nfft=8,boundary='none',padded=.false.)
    if (st/=spectral_status_success) error stop 'custom window initialization'
    call stft_stream_push(state,[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],spectra,st)
    if (st/=spectral_status_success .or. size(spectra,2)/=1) error stop 'custom window transform'
    if (abs(spectra(1,1)-cmplx(3.0_dp,0.0_dp,dp))>1e-12_dp) error stop 'custom DC bin'

    call check_round_trip(.false.)
    call check_round_trip(.true.)
    call check_real_twosided()
    call check_short_and_empty()
    call check_spectrogram_default()

    ! Pushing one sample at a time has constant state size, even after many frames.
    call stft_stream_init(state,8,st,window='boxcar',noverlap=4, &
        boundary='zeros',padded=.true.)
    if (st/=spectral_status_success) error stop 'stress init'
    nframes=0
    do j=1,20000
        call stft_stream_push(state,[sin(0.01_dp*real(j,dp))],spectra,st)
        if (st/=spectral_status_success) error stop 'stress push'
        nframes=nframes+size(spectra,2)
    end do
    call stft_stream_finish(state,spectra,st)
    nframes=nframes+size(spectra,2)
    if (st/=spectral_status_success .or. nframes/=5001) error stop 'stress frame count'

    print '(a,i0,a)', 'test_signal_stft_stream passed: ',cases,' batch/stream comparisons'
contains
    subroutine check_round_trip(complex_input)
        logical, intent(in) :: complex_input
        type(stft_stream) :: encoder
        type(istft_stream) :: decoder
        real(dp) :: xr(120)
        complex(dp) :: xc(120)
        complex(dp), allocatable :: values(:,:), restored_c(:), part_c(:)
        real(dp), allocatable :: restored_r(:), part_r(:)
        integer :: j, k, st, chunk_end, bin_count, offset

        do j=1,120
            xr(j)=sin(0.14_dp*j)+0.7_dp*cos(0.37_dp*j)
            xc(j)=cmplx(xr(j),sin(0.43_dp*j),dp)
        end do
        bin_count=9
        if (complex_input) bin_count=16
        call stft_stream_init(encoder,16,st,nfft=16,noverlap=12, &
            window='hann',boundary='zeros',input_complex=complex_input)
        if (st/=spectral_status_success) error stop 'roundtrip forward initialization'
        call istft_stream_init(decoder,bin_count,st,nperseg=16,nfft=16, &
            noverlap=12,window='hann',boundary=.true.,input_onesided=.not.complex_input)
        if (st/=spectral_status_success) error stop 'roundtrip inverse initialization'
        allocate(restored_r(0),restored_c(0))
        do j=1,120,11
            chunk_end=min(120,j+10)
            if (complex_input) then
                call stft_stream_push(encoder,xc(j:chunk_end),values,st)
            else
                call stft_stream_push(encoder,xr(j:chunk_end),values,st)
            end if
            if (st/=spectral_status_success) error stop 'roundtrip forward push'
            do k=1,size(values,2)
                if (complex_input) then
                    call istft_stream_push(decoder,values(:,k),part_c,st)
                    restored_c=[restored_c,part_c]
                else
                    call istft_stream_push(decoder,values(:,k),part_r,st)
                    restored_r=[restored_r,part_r]
                end if
                if (st/=spectral_status_success) error stop 'roundtrip inverse push'
            end do
        end do
        call stft_stream_finish(encoder,values,st)
        if (st/=spectral_status_success) error stop 'roundtrip forward finish'
        do k=1,size(values,2)
            if (complex_input) then
                call istft_stream_push(decoder,values(:,k),part_c,st)
                restored_c=[restored_c,part_c]
            else
                call istft_stream_push(decoder,values(:,k),part_r,st)
                restored_r=[restored_r,part_r]
            end if
            if (st/=spectral_status_success) error stop 'roundtrip inverse final frames'
        end do
        if (complex_input) then
            call istft_stream_finish(decoder,part_c,st)
            restored_c=[restored_c,part_c]
            if (st/=spectral_status_success .or. size(restored_c)/=120) error stop 'complex roundtrip length'
            if (maxval(abs(restored_c-xc))>2.0e-12_dp) error stop 'complex roundtrip values'
        else
            call istft_stream_finish(decoder,part_r,st)
            restored_r=[restored_r,part_r]
            if (st/=spectral_status_success .or. size(restored_r)/=120) error stop 'real roundtrip length'
            if (maxval(abs(restored_r-xr))>2.0e-12_dp) error stop 'real roundtrip values'
        end if
    end subroutine check_round_trip

    subroutine check_real_twosided()
        type(stft_stream) :: state
        real(dp) :: x(31)
        real(dp), allocatable :: frequencies(:)
        complex(dp), allocatable :: values(:,:), tail(:,:)
        integer :: i, st

        x=[(sin(0.07_dp*real(i,dp)),i=1,31)]
        call stft_stream_init(state,7,st,window='boxcar',nfft=9,noverlap=4, &
            boundary='none',padded=.false.,return_onesided=.false.)
        if (st/=spectral_status_success) error stop 'two-sided real init'
        call stft_stream_frequencies(state,frequencies,st)
        if (size(frequencies)/=9) error stop 'real two-sided bin count'
        if (abs(frequencies(6)+4.0_dp/9.0_dp)>1e-12_dp) error stop 'real negative frequency'
        call stft_stream_push(state,x,values,st)
        call stft_stream_finish(state,tail,st)
        if (st/=spectral_status_success .or. size(values,2)+size(tail,2)/=9) &
            error stop 'real two-sided frames'
    end subroutine check_real_twosided

    subroutine check_short_and_empty()
        type(stft_stream) :: state
        type(stft_result) :: batch
        complex(dp), allocatable :: part(:,:), tail(:,:)
        real(dp), allocatable :: times(:)
        integer :: st

        call stft_stream_init(state,1,st,window='boxcar',noverlap=0,nfft=1)
        if (st/=spectral_status_success) error stop 'single-sample stream init'
        call stft_stream_finish(state,part,st,times)
        if (st/=spectral_status_success .or. size(part)/=0 .or. size(times)/=0) &
            error stop 'empty stream finish'
        call stft_stream_init(state,1,st,window='boxcar',noverlap=0,nfft=1)
        if (st/=spectral_status_success) error stop 'empty stream reinitialization'
        call stft_stream_push(state,[2.0_dp,3.0_dp,5.0_dp],part,st,times)
        if (st/=spectral_status_success) error stop 'single-sample stream push'
        batch=stft([2.0_dp,3.0_dp,5.0_dp],window='boxcar',nperseg=1,noverlap=0,nfft=1)
        if (batch%status/=spectral_status_success) error stop 'single-sample batch STFT'
        if (any(shape(batch%values)/=shape(part))) error stop 'single-sample STFT shape'
        if (maxval(abs(part-batch%values))>1e-12_dp) error stop 'single-sample STFT'
        call stft_stream_finish(state,tail,st)
        if (st/=spectral_status_success .or. size(tail)/=0) error stop 'single-sample tail'
        call stft_stream_init(state,7,st,window='boxcar',boundary='none',padded=.false.)
        call stft_stream_push(state,[1.0_dp,2.0_dp],part,st)
        if (st/=spectral_status_success .or. size(part)/=0) error stop 'short push'
        call stft_stream_finish(state,part,st)
        if (st/=spectral_status_success .or. size(part,2)/=1) error stop 'short final frame'
    end subroutine check_short_and_empty

    subroutine check_spectrogram_default()
        type(stft_stream) :: state
        type(spectrogram_result) :: batch
        real(dp) :: x(32)
        real(dp), allocatable :: part(:,:), tail(:,:), custom(:)
        integer :: i, st

        x=[(sin(0.19_dp*real(i,dp)),i=1,32)]
        batch=spectrogram(x,nperseg=8,noverlap=2,nfft=8)
        if (batch%status/=spectral_status_success) error stop 'Tukey batch'
        call spectrogram_stream_init(state,8,st,noverlap=2,nfft=8)
        if (st/=spectral_status_success) error stop 'default Tukey initialization'
        call spectrogram_stream_push(state,x,part,st)
        call spectrogram_stream_finish(state,tail,st)
        if (st/=spectral_status_success) error stop 'Tukey finalization'
        if (size(tail,2)/=0) error stop 'spectrogram unexpected tail'
        if (any(shape(part)/=shape(batch%power))) error stop 'Tukey shapes'
        if (maxval(abs(part-batch%power))>1e-12_dp) error stop 'Tukey batch/stream mismatch'
        custom=get_window('tukey',8,parameter=0.25_dp)
        call spectrogram_stream_init(state,8,st,noverlap=2,nfft=8,window_values=custom)
        if (st/=spectral_status_success) error stop 'parameterized Tukey init'
        call spectrogram_stream_push(state,x,part,st)
        if (st/=spectral_status_success .or. size(part,2)/=5) &
            error stop 'parameterized Tukey frames'
    end subroutine check_spectrogram_default

    subroutine check_stft(n,nw,no,nf,edge,pad,scaling,trend,complex_input)
        integer, intent(in) :: n,nw,no,nf
        character(*), intent(in) :: edge,scaling,trend
        logical, intent(in) :: pad,complex_input
        real(dp) :: x(n)
        complex(dp) :: z(n)
        type(stft_result) :: reference
        type(stft_stream) :: state
        complex(dp), allocatable :: all(:,:), part(:,:)
        real(dp), allocatable :: all_t(:),chunk_t(:),freq(:)
        integer :: j,start,stop,step,st,bin_count

        do j=1,n
            x(j)=sin(0.17_dp*j)+0.1_dp*cos(0.8_dp*j)+0.01_dp*j
            z(j)=cmplx(x(j),cos(0.27_dp*j)-0.04_dp*j,dp)
        end do
        if (complex_input) then
            reference=stft(z,fs=20.0_dp,window='boxcar',nperseg=nw,noverlap=no, &
                nfft=nf,detrend=trend,boundary=edge,padded=pad,scaling=scaling)
        else
            reference=stft(x,fs=20.0_dp,window='boxcar',nperseg=nw,noverlap=no, &
                nfft=nf,detrend=trend,boundary=edge,padded=pad,scaling=scaling)
        end if
        if (reference%status/=spectral_status_success) error stop 'batch stft failed'
        call stft_stream_init(state,nw,st,fs=20.0_dp,window='boxcar',noverlap=no, &
            nfft=nf,detrend=trend,boundary=edge,padded=pad,scaling=scaling, &
            input_complex=complex_input)
        if (st/=spectral_status_success) error stop 'stream init failed'
        call stft_stream_frequencies(state,freq,st)
        if (st/=spectral_status_success) error stop 'frequency query failed'
        if (size(freq)/=size(reference%frequencies)) error stop 'frequency count'
        if (maxval(abs(freq-reference%frequencies))>2e-12_dp) error stop 'frequency mismatch'
        bin_count=size(freq)
        allocate(all(bin_count,0),all_t(0))
        start=1
        do while(start<=n)
            step=mod(start+nw+no,11)+1
            stop=min(n,start+step-1)
            if (complex_input) then
                call stft_stream_push(state,z(start:stop),part,st,chunk_t)
            else
                call stft_stream_push(state,x(start:stop),part,st,chunk_t)
            end if
            if (st/=spectral_status_success) error stop 'STFT stream push'
            all=reshape([reshape(all,[size(all)]),reshape(part,[size(part)])], &
                [bin_count,size(all,2)+size(part,2)])
            all_t=[all_t,chunk_t]
            start=stop+1
        end do
        call stft_stream_finish(state,part,st,chunk_t)
        if (st/=spectral_status_success) error stop 'STFT stream finish'
        all=reshape([reshape(all,[size(all)]),reshape(part,[size(part)])], &
            [bin_count,size(all,2)+size(part,2)])
        all_t=[all_t,chunk_t]
        if (any(shape(all)/=shape(reference%values))) error stop 'STFT output shape'
        if (size(all_t)/=size(reference%times)) error stop 'STFT times shape'
        if (size(all)>0) then
            if (maxval(abs(all-reference%values))>2e-11_dp) error stop 'STFT values mismatch'
        end if
        if (size(all_t)>0) then
            if (maxval(abs(all_t-reference%times))>1e-12_dp) then
                print *, 'TIMES DIFF',n,nw,no,nf,edge,pad,scaling,trend,complex_input
                print *, 'got=',all_t
                print *, 'want=',reference%times
                error stop 'STFT times mismatch'
            end if
        end if
    end subroutine check_stft

    subroutine check_spectrogram(n,nw,no,nf,scaling,trend,complex_input,window,window_values)
        implicit none
        integer, intent(in) :: n !! positive input sample count
        integer, intent(in) :: nw !! positive analysis window length
        integer, intent(in) :: no !! overlap samples, 0 <= no < nw
        integer, intent(in) :: nf !! FFT length, at least nw
        character(*), intent(in) :: scaling !! density or spectrum normalization
        character(*), intent(in) :: trend !! none, constant, or linear detrending
        logical, intent(in) :: complex_input !! true checks complex samples and two-sided output
        character(*), intent(in), optional :: window !! periodic window; absent uses batch default
        real(dp), intent(in), optional :: window_values(:) !! custom weights, size nw
        real(dp) :: x(n)
        complex(dp) :: z(n)
        type(spectrogram_result) :: reference
        type(stft_stream) :: state
        real(dp), allocatable :: all(:,:), part(:,:), all_t(:), chunk_t(:),freq(:)
        integer :: j,start,stop,step,st,bin_count

        do j=1,n
            x(j)=sin(0.29_dp*j)+0.03_dp*j+cos(0.16_dp*j)
            z(j)=cmplx(x(j),cos(0.61_dp*j)-0.05_dp*j,dp)
        end do
        if (complex_input) then
            reference=spectrogram(z,fs=37.0_dp,window=window,nperseg=nw, &
                noverlap=no,nfft=nf,detrend=trend,scaling=scaling,window_values=window_values)
        else
            reference=spectrogram(x,fs=37.0_dp,window=window,nperseg=nw, &
                noverlap=no,nfft=nf,detrend=trend,scaling=scaling,window_values=window_values)
        end if
        if (reference%status/=spectral_status_success) error stop 'batch spectrogram failed'
        call spectrogram_stream_init(state,nw,st,fs=37.0_dp,window=window, &
            noverlap=no,nfft=nf,detrend=trend,scaling=scaling,window_values=window_values, &
            input_complex=complex_input)
        if (st/=spectral_status_success) error stop 'spectrogram stream init'
        call stft_stream_frequencies(state,freq,st)
        if (st/=spectral_status_success) error stop 'frequency query'
        if (size(freq)/=size(reference%frequencies)) error stop 'frequency count'
        if (maxval(abs(freq-reference%frequencies))>2e-12_dp) error stop 'spectrogram frequencies'
        bin_count=size(freq)
        allocate(all(bin_count,0),all_t(0))
        start=1
        do while(start<=n)
            step=mod(start+no,13)+1
            stop=min(n,start+step-1)
            if (complex_input) then
                call spectrogram_stream_push(state,z(start:stop),part,st,chunk_t)
            else
                call spectrogram_stream_push(state,x(start:stop),part,st,chunk_t)
            end if
            if (st/=spectral_status_success) error stop 'spectrogram push failed'
            all=reshape([reshape(all,[size(all)]),reshape(part,[size(part)])], &
                [bin_count,size(all,2)+size(part,2)])
            all_t=[all_t,chunk_t]
            start=stop+1
        end do
        call spectrogram_stream_finish(state,part,st,chunk_t)
        if (st/=spectral_status_success) error stop 'spectrogram finish'
        all=reshape([reshape(all,[size(all)]),reshape(part,[size(part)])], &
            [bin_count,size(all,2)+size(part,2)])
        all_t=[all_t,chunk_t]
        if (any(shape(all)/=shape(reference%power))) error stop 'spectrogram shape'
        if (size(all)>0) then
            if (maxval(abs(all-reference%power))>2e-11_dp) error stop 'spectrogram power mismatch'
        end if
        if (size(all_t)/=size(reference%times)) error stop 'spectrogram times shape'
        if (size(all_t)>0) then
            if (maxval(abs(all_t-reference%times))>1e-12_dp) error stop 'spectrogram times mismatch'
        end if
    end subroutine check_spectrogram
end program test_signal_stft_stream
