! SPDX-License-Identifier: MIT
! Reference driver for time-frequency linear detrending, including cross channels.
program driver_signal_multitaper_linear_detrend
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram, &
        multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: out
    integer :: ncase, icase, n, nc, seg, ov, nk, nf, iscomplex, i, j, k, frame, ch, ch2
    real(dp) :: nw, fs, re, im
    complex(dp), allocatable :: data(:,:)
    real(dp), allocatable :: real_data(:,:)

    read(*,*) ncase
    do icase=1,ncase
        read(*,*) n,nc,seg,ov,nk,nf,iscomplex,nw,fs
        allocate(data(n,nc))
        do ch=1,nc
            do i=1,n
                read(*,*) re,im
                data(i,ch) = cmplx(re,im,dp)
            end do
        end do
        if (iscomplex == 1) then
            out = multitaper_spectral_matrix_spectrogram(data,nw,seg, &
                noverlap=ov,kmax=nk,nfft=nf,fs=fs,weighting='equal',detrend_mode='linear')
        else
            real_data = real(data,dp)
            out = multitaper_spectral_matrix_spectrogram(real_data,nw,seg, &
                noverlap=ov,kmax=nk,nfft=nf,fs=fs,weighting='equal',detrend_mode='linear')
        end if
        write(*,'(4(i0,1x))') out%status,size(out%power,1),size(out%power,2),nc
        if (out%status == 0) then
            do frame=1,size(out%power,2)
                do ch2=1,nc
                    do ch=1,nc
                        do k=1,size(out%power,1)
                            write(*,'(2(es26.17,1x))') real(out%power(k,frame,ch,ch2),dp), &
                                aimag(out%power(k,frame,ch,ch2))
                        end do
                    end do
                end do
            end do
        end if
        deallocate(data)
        if (allocated(real_data)) deallocate(real_data)
    end do
end program driver_signal_multitaper_linear_detrend
