! SPDX-License-Identifier: MIT
! Plain-text driver for independent partial coherency/shrinkage comparisons.
program reference_signal_partial_coherency
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix_spectrogram
    use scifort_signal_partial_coherence, only : partial_coherence_result, &
        partial_coherence_tf_result, partial_coherence_matrix, partial_coherence_tf_matrix, &
        multitaper_tf_partial_coherence
    implicit none
    type(multitaper_matrix_result) :: static_spectrum
    type(multitaper_spectral_matrix_spectrogram_result) :: spectrum
    type(partial_coherence_result) :: pc1
    type(partial_coherence_tf_result) :: pc
    real(dp), allocatable :: rx(:,:)
    complex(dp), allocatable :: cx(:,:)
    real(dp) :: re, im, ridge, shrinkage, nw, fs, fvalue, tvalue, cv
    complex(dp) :: pv, zv
    integer :: mode, nf, nt, nc, n, nperseg, noverlap, kmax, nfft
    integer :: axis, i, j, f, t, ierr, valid, status
    character(20) :: weighting, trend

    do
        read(*,*,iostat=ierr) mode
        if (ierr /= 0) exit
        if (mode == 0 .or. mode == 1) then
            read(*,*) nf,nt,nc,ridge,shrinkage
            if (allocated(spectrum%frequencies)) deallocate(spectrum%frequencies)
            if (allocated(spectrum%times)) deallocate(spectrum%times)
            if (allocated(spectrum%power)) deallocate(spectrum%power)
            allocate(spectrum%frequencies(nf),spectrum%times(nt),spectrum%power(nf,nt,nc,nc))
            spectrum%frequencies = [(real(f-1,dp),f=1,nf)]
            spectrum%times = [(real(t-1,dp),t=1,nt)]
            spectrum%status = 0
            do t = 1, nt
                do f = 1, nf
                    do j = 1, nc
                        do i = 1, nc
                            read(*,*) re,im
                            spectrum%power(f,t,i,j) = cmplx(re,im,dp)
                        end do
                    end do
                end do
            end do
            if (mode == 0) then
                if (allocated(static_spectrum%frequencies)) deallocate(static_spectrum%frequencies)
                if (allocated(static_spectrum%power)) deallocate(static_spectrum%power)
                allocate(static_spectrum%frequencies(nf),static_spectrum%power(nf,nc,nc))
                static_spectrum%status = 0
                static_spectrum%frequencies = spectrum%frequencies
                static_spectrum%power = spectrum%power(:,1,:,:)
                pc1 = partial_coherence_matrix(static_spectrum,ridge=ridge, &
                    shrinkage=shrinkage,return_coherency=.true.)
            else
                pc = partial_coherence_tf_matrix(spectrum,ridge=ridge, &
                    shrinkage=shrinkage,return_coherency=.true.)
            end if
        else
            read(*,*) n,nc,nperseg,noverlap,nw,kmax,nfft,fs,ridge, &
                shrinkage,axis,weighting,trend
            allocate(cx(n,nc),rx(n,nc))
            do j = 1,nc
                do i = 1,n
                    read(*,*) re,im
                    cx(i,j) = cmplx(re,im,dp)
                    rx(i,j) = re
                end do
            end do
            if (mode == 2) then
                if (axis == 1) then
                    spectrum = multitaper_spectral_matrix_spectrogram(rx,nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft, &
                        weighting=trim(weighting),detrend_mode=trim(trend))
                    pc = multitaper_tf_partial_coherence(rx,nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft, &
                        weighting=trim(weighting),detrend_mode=trim(trend), &
                        ridge=ridge,shrinkage=shrinkage,return_coherency=.true.)
                else
                    spectrum = multitaper_spectral_matrix_spectrogram(transpose(rx),nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft,axis=2, &
                        weighting=trim(weighting),detrend_mode=trim(trend))
                    pc = multitaper_tf_partial_coherence(transpose(rx),nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft,axis=2, &
                        weighting=trim(weighting),detrend_mode=trim(trend), &
                        ridge=ridge,shrinkage=shrinkage,return_coherency=.true.)
                end if
            else
                if (axis == 1) then
                    spectrum = multitaper_spectral_matrix_spectrogram(cx,nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft, &
                        weighting=trim(weighting),detrend_mode=trim(trend))
                    pc = multitaper_tf_partial_coherence(cx,nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft, &
                        weighting=trim(weighting),detrend_mode=trim(trend), &
                        ridge=ridge,shrinkage=shrinkage,return_coherency=.true.)
                else
                    spectrum = multitaper_spectral_matrix_spectrogram(transpose(cx),nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft,axis=2, &
                        weighting=trim(weighting),detrend_mode=trim(trend))
                    pc = multitaper_tf_partial_coherence(transpose(cx),nw,nperseg, &
                        noverlap=noverlap,kmax=kmax,fs=fs,nfft=nfft,axis=2, &
                        weighting=trim(weighting),detrend_mode=trim(trend), &
                        ridge=ridge,shrinkage=shrinkage,return_coherency=.true.)
                end if
            end if
            deallocate(cx,rx)
        end if
        if (mode == 0) then
            status = pc1%status
            nf = size(pc1%frequencies)
            nt = 1
            nc = size(pc1%coherence,2)
        else
            status = pc%status
            nf = size(pc%frequencies)
            nt = size(pc%times)
            nc = size(pc%coherence,3)
        end if
        write(*,'(4(i0,1x))') status,nf,nt,nc
        do t = 1, nt
            do f = 1, nf
                if (mode == 0) then
                    fvalue = pc1%frequencies(f)
                    tvalue = 0.0_dp
                    valid = 0
                    if (pc1%valid(f)) valid = 1
                else
                    fvalue = pc%frequencies(f)
                    tvalue = pc%times(t)
                    valid = 0
                    if (pc%valid(f,t)) valid = 1
                end if
                write(*,'(2(es26.17,1x),i0)') fvalue,tvalue,valid
                do j = 1,nc
                    do i = 1,nc
                        if (mode == 0) then
                            cv = pc1%coherence(f,i,j)
                            zv = pc1%coherency(f,i,j)
                            pv = pc1%precision(f,i,j)
                        else
                            cv = pc%coherence(f,t,i,j)
                            zv = pc%coherency(f,t,i,j)
                            pv = pc%precision(f,t,i,j)
                        end if
                        write(*,'(7(es26.17,1x))') cv,real(zv,dp),aimag(zv), &
                            real(pv,dp),aimag(pv),real(spectrum%power(f,t,i,j),dp), &
                            aimag(spectrum%power(f,t,i,j))
                    end do
                end do
            end do
        end do
    end do
end program reference_signal_partial_coherency
