! SPDX-License-Identifier: MIT
! Machine-readable driver for independent SciPy/NumPy comparisons.
program reference_signal_partial_coherence
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, multitaper_spectral_matrix
    use scifort_signal_partial_coherence, only : partial_coherence_result, partial_coherence_matrix
    implicit none
    type(multitaper_matrix_result) :: spectral
    type(partial_coherence_result) :: pc
    complex(dp), allocatable :: cx(:,:)
    real(dp), allocatable :: rx(:,:)
    real(dp) :: re, im, ridge, nw, fs
    integer :: mode, nf, nc, n, kmax, nfft, is_real, i, j, f, ierr, valid
    character(20) :: weighting, trend

    do
        read(*,*,iostat=ierr) mode
        if (ierr /= 0) exit
        if (mode == 0) then
            read(*,*) nf, nc, ridge
            if (allocated(spectral%frequencies)) deallocate(spectral%frequencies)
            if (allocated(spectral%power)) deallocate(spectral%power)
            allocate(spectral%frequencies(nf), spectral%power(nf,nc,nc))
            spectral%frequencies = [(real(f-1,dp),f=1,nf)]
            spectral%status = 0
            do f=1,nf
                do j=1,nc
                    do i=1,nc
                        read(*,*) re, im
                        spectral%power(f,i,j) = cmplx(re,im,dp)
                    end do
                end do
            end do
        else
            read(*,*) n, nc, nw, kmax, nfft, fs, ridge, is_real, weighting, trend
            allocate(cx(n,nc),rx(n,nc))
            do j=1,nc
                do i=1,n
                    read(*,*) re, im
                    cx(i,j) = cmplx(re,im,dp)
                    rx(i,j) = re
                end do
            end do
            if (is_real /= 0) then
                spectral = multitaper_spectral_matrix(rx,nw,kmax=kmax,fs=fs, &
                    nfft=nfft,weighting=trim(weighting),detrend_mode=trim(trend))
            else
                spectral = multitaper_spectral_matrix(cx,nw,kmax=kmax,fs=fs, &
                    nfft=nfft,weighting=trim(weighting),detrend_mode=trim(trend))
            end if
            deallocate(cx,rx)
        end if
        pc = partial_coherence_matrix(spectral, ridge=ridge)
        nf=0
        nc=0
        if (allocated(pc%frequencies)) nf = size(pc%frequencies)
        if (allocated(pc%coherence)) nc = size(pc%coherence,2)
        write(*,'(3(i0,1x))') pc%status,nf,nc
        do f=1,nf
            valid=0
            if (pc%valid(f)) valid=1
            write(*,'(es26.17,1x,i0)') pc%frequencies(f),valid
            do j=1,nc
                do i=1,nc
                    write(*,'(5(es26.17,1x))') pc%coherence(f,i,j), &
                        real(pc%precision(f,i,j),dp),aimag(pc%precision(f,i,j)), &
                        real(spectral%power(f,i,j),dp),aimag(spectral%power(f,i,j))
                end do
            end do
        end do
    end do
end program reference_signal_partial_coherence
