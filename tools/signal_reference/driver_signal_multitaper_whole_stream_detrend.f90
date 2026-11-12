! SPDX-License-Identifier: MIT
! Reference driver for independent SciPy comparisons of whole-record and
! streaming linear multitaper detrending. Input samples are channel-major.
program driver_signal_multitaper_whole_stream_detrend
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_stream, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram, multitaper_psd, multitaper_psd_ci, &
        multitaper_matrix_result, multitaper_confidence_result, &
        multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_spectral, only : spectral_result
    implicit none
    type(multitaper_stream) :: stream
    type(multitaper_matrix_result) :: full
    type(multitaper_spectral_matrix_spectrogram_result) :: tf
    type(multitaper_confidence_result) :: ci
    type(spectral_result) :: psd
    integer :: ncase, icase, n, nc, seg, ov, nt, nf, iscomplex, mode, i, j, k, l, c, d, status, nr
    real(dp) :: nw, rate, r, im
    complex(dp), allocatable :: x(:,:)
    real(dp), allocatable :: y(:,:)
    character(8), parameter :: modes(2) = [character(8) :: 'equal', 'eigen']

    read(*,*) ncase
    do icase=1,ncase
        read(*,*) n,nc,seg,ov,nt,nf,iscomplex,mode,nw,rate
        allocate(x(n,nc),y(n,nc))
        do c=1,nc
            do i=1,n
                read(*,*) r,im
                x(i,c)=cmplx(r,im,dp)
            end do
        end do
        if (iscomplex == 1) then
            full=multitaper_spectral_matrix(x,nw,kmax=nt,nfft=nf,fs=rate, &
                weighting=trim(modes(mode)),detrend_mode='linear')
            call stream%initialize(nw,seg,nc,status,noverlap=ov,kmax=nt,nfft=seg+3, &
                fs=rate,weighting=trim(modes(mode)),complex_input=.true.,detrend_mode='linear')
            tf=stream%push(x)
            psd=multitaper_psd(x(:,1),nw,kmax=nt,nfft=nf,fs=rate, &
                weighting=trim(modes(mode)),detrend_mode='linear')
            ci=multitaper_psd_ci(x(:,1),nw,kmax=nt,nfft=nf,fs=rate, &
                weighting=trim(modes(mode)),detrend_mode='linear')
        else
            y=real(x,dp)
            full=multitaper_spectral_matrix(y,nw,kmax=nt,nfft=nf,fs=rate, &
                weighting=trim(modes(mode)),detrend_mode='linear')
            call stream%initialize(nw,seg,nc,status,noverlap=ov,kmax=nt,nfft=seg+3, &
                fs=rate,weighting=trim(modes(mode)),detrend_mode='linear')
            tf=stream%push(y)
            psd=multitaper_psd(y(:,1),nw,kmax=nt,nfft=nf,fs=rate, &
                weighting=trim(modes(mode)),detrend_mode='linear')
            ci=multitaper_psd_ci(y(:,1),nw,kmax=nt,nfft=nf,fs=rate, &
                weighting=trim(modes(mode)),detrend_mode='linear')
        end if
        nr=size(full%frequencies)
        write(*,'(7(i0,1x))') full%status,status,tf%status,psd%status,ci%status,nr,size(tf%times)
        do i=1,nr
            write(*,'(3(es26.17,1x))') full%frequencies(i),psd%power(i),ci%power(i)
        end do
        do d=1,nc
            do c=1,nc
                do i=1,nr
                    write(*,'(2(es26.17,1x))') real(full%power(i,c,d),dp),aimag(full%power(i,c,d))
                end do
            end do
        end do
        do j=1,size(tf%times)
            do d=1,nc
                do c=1,nc
                    do i=1,size(tf%frequencies)
                        write(*,'(2(es26.17,1x))') real(tf%power(i,j,c,d),dp), &
                            aimag(tf%power(i,j,c,d))
                    end do
                end do
            end do
        end do
        deallocate(x,y)
    end do
end program driver_signal_multitaper_whole_stream_detrend
