! SPDX-License-Identifier: MIT
! Text protocol driver for independent time-frequency partial-coherence checks.
program reference_signal_tf_partial_coherence
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result, &
        multitaper_spectral_matrix_spectrogram
    use scifort_signal_partial_coherence, only : partial_coherence_tf_result, &
        partial_coherence_tf_matrix, multitaper_tf_partial_coherence
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectrum
    type(partial_coherence_tf_result) :: pc
    real(dp), allocatable :: rx(:,:)
    complex(dp), allocatable :: cx(:,:)
    real(dp) :: re, im, ridge, nw, fs
    integer :: mode, nf, nt, nc, n, nperseg, noverlap, kmax, nfft
    integer :: is_real, axis, i, j, f, t, ierr, valid
    character(20) :: weighting, trend

    do
        read(*,*,iostat=ierr) mode
        if (ierr /= 0) exit
        if (mode == 0) then
            read(*,*) nf, nt, nc, ridge
            if (allocated(spectrum%frequencies)) deallocate(spectrum%frequencies)
            if (allocated(spectrum%times)) deallocate(spectrum%times)
            if (allocated(spectrum%power)) deallocate(spectrum%power)
            allocate(spectrum%frequencies(nf), spectrum%times(nt), spectrum%power(nf,nt,nc,nc))
            spectrum%frequencies = [(real(f-1,dp),f=1,nf)]
            spectrum%times = [(real(t-1,dp),t=1,nt)]
            spectrum%status = 0
            do t = 1, nt
                do f = 1, nf
                    do j = 1, nc
                        do i = 1, nc
                            read(*,*) re, im
                            spectrum%power(f,t,i,j) = cmplx(re,im,dp)
                        end do
                    end do
                end do
            end do
            pc = partial_coherence_tf_matrix(spectrum,ridge=ridge)
        else
            read(*,*) n,nc,nperseg,noverlap,nw,kmax,nfft,fs,ridge,is_real,axis,weighting,trend
            allocate(cx(n,nc),rx(n,nc))
            do j = 1,nc
                do i = 1,n
                    read(*,*) re,im
                    cx(i,j) = cmplx(re,im,dp)
                    rx(i,j) = re
                end do
            end do
            if (is_real /= 0) then
                if (axis == 1) then
                    spectrum = multitaper_spectral_matrix_spectrogram(rx,nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft, &
                        weighting=trim(weighting),detrend_mode=trim(trend))
                    pc = multitaper_tf_partial_coherence(rx,nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft, &
                        weighting=trim(weighting),detrend_mode=trim(trend),ridge=ridge)
                else
                    spectrum = multitaper_spectral_matrix_spectrogram(transpose(rx),nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft,axis=2, &
                        weighting=trim(weighting),detrend_mode=trim(trend))
                    pc = multitaper_tf_partial_coherence(transpose(rx),nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft,axis=2, &
                        weighting=trim(weighting),detrend_mode=trim(trend),ridge=ridge)
                end if
            else
                if (axis == 1) then
                    spectrum = multitaper_spectral_matrix_spectrogram(cx,nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft, &
                        weighting=trim(weighting),detrend_mode=trim(trend))
                    pc = multitaper_tf_partial_coherence(cx,nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft, &
                        weighting=trim(weighting),detrend_mode=trim(trend),ridge=ridge)
                else
                    spectrum = multitaper_spectral_matrix_spectrogram(transpose(cx),nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft,axis=2, &
                        weighting=trim(weighting),detrend_mode=trim(trend))
                    pc = multitaper_tf_partial_coherence(transpose(cx),nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft,axis=2, &
                        weighting=trim(weighting),detrend_mode=trim(trend),ridge=ridge)
                end if
            end if
            deallocate(cx,rx)
        end if
        nf = size(pc%frequencies)
        nt = size(pc%times)
        nc = 0
        if (allocated(pc%coherence)) nc = size(pc%coherence,3)
        write(*,'(4(i0,1x))') pc%status,nf,nt,nc
        do t=1,nt
            do f=1,nf
                valid=0
                if (pc%valid(f,t)) valid=1
                write(*,'(2(es26.17,1x),i0)') pc%frequencies(f),pc%times(t),valid
                do j=1,nc
                    do i=1,nc
                        write(*,'(5(es26.17,1x))') pc%coherence(f,t,i,j), &
                            real(pc%precision(f,t,i,j),dp),aimag(pc%precision(f,t,i,j)), &
                            real(spectrum%power(f,t,i,j),dp),aimag(spectrum%power(f,t,i,j))
                    end do
                end do
            end do
        end do
    end do
end program reference_signal_tf_partial_coherence
