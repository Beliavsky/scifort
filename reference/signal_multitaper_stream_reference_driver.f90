! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program signal_multitaper_stream_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_stream, &
        multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(multitaper_stream) :: stream
    type(multitaper_spectral_matrix_spectrogram_result) :: new
    complex(dp), allocatable :: x(:,:), power(:,:,:,:)
    real(dp), allocatable :: times(:)
    integer :: trials, trial, n, seg, overlap, kmax, nfft, kind, center, nframes, &
        nfreq, status, sample, last, nt, step, done, bin, frame
    real(dp) :: nw, fs, ar, ai, br, bi
    character(12) :: method
    read(*,*) trials
    do trial = 1, trials
        read(*,*) n,seg,overlap,kmax,nfft,kind,center,nw,fs,method
        allocate(x(n,2))
        do sample = 1, n
            read(*,*) ar,ai,br,bi
            x(sample,1) = cmplx(ar,ai,dp)
            x(sample,2) = cmplx(br,bi,dp)
        end do
        step = seg-overlap
        nframes = 1+(n-seg)/step
        nfreq = nfft
        if (kind == 0) nfreq = nfft/2+1
        allocate(power(nfreq,nframes,2,2),times(nframes))
        call stream%initialize(nw,seg,2,status,noverlap=overlap,kmax=kmax,nfft=nfft, &
            fs=fs,detrend=center/=0,weighting=trim(method),complex_input=kind==1)
        if (status /= 0) error stop 'init failed in reference driver'
        done = 0
        sample = 1
        do while (sample <= n)
            last = min(n,sample+mod(7*sample+trial,13))
            if (kind == 0) then
                new = stream%push(real(x(sample:last,:),dp))
            else
                new = stream%push(x(sample:last,:))
            end if
            if (new%status /= 0) error stop 'push failed in reference driver'
            nt = size(new%times)
            if (nt > 0) then
                power(:,done+1:done+nt,:,:) = new%power
                times(done+1:done+nt) = new%times
            end if
            done = done+nt
            sample = last+1
        end do
        if (done /= nframes) error stop 'unexpected frame count'
        print '(4(i0,1x))', status,nfreq,nframes,2
        do bin = 1, nfreq
            print '(es26.17)', new%frequencies(bin)
        end do
        do frame = 1, nframes
            print '(es26.17)', times(frame)
        end do
        do frame = 1, nframes
            do bin = 1, nfreq
                print '(4(es26.17,1x))', real(power(bin,frame,1,2)),aimag(power(bin,frame,1,2)), &
                    real(power(bin,frame,1,1)),real(power(bin,frame,2,2))
            end do
        end do
        deallocate(x,power,times)
    end do
end program signal_multitaper_stream_reference_driver
