! SPDX-License-Identifier: MIT
! Serialization driver for independent NumPy/SciPy exact-band tests.
program fixed_taper_spectrogram_driver
    use scifort_kinds, only : dp
    use scifort_signal_fixed_taper_spectrogram, only : fixed_taper_spectral_matrix_spectrogram
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: out
    complex(dp), allocatable :: data(:,:)
    real(dp), allocatable :: raw(:), tapers(:,:), weights(:), freqs(:), bands(:,:)
    integer :: ncases, case_idx, ns, nc, length, nk, nf, nb, noverlap, norm, avg, complex_flag
    integer :: i, j, k, f, w

    read(*,*) ncases
    do case_idx = 1, ncases
        read(*,*) ns,nc,length,nk,nf,nb,noverlap,norm,avg,complex_flag
        allocate(raw(2*ns*nc),data(ns,nc),tapers(length,nk),weights(nk), &
            freqs(nf),bands(2,nb))
        read(*,*) raw
        do j = 1, nc
            do i = 1, ns
                k = 2*((j-1)*ns+i-1)+1
                data(i,j) = cmplx(raw(k),raw(k+1),dp)
            end do
        end do
        read(*,*) tapers
        read(*,*) weights
        if (nf > 0) read(*,*) freqs
        if (nb > 0) then
            do f = 1, nb
                read(*,*) bands(:,f)
            end do
        end if
        if (complex_flag == 1) then
            out = fixed_taper_spectral_matrix_spectrogram(data,tapers,freqs, &
                bands=bands,taper_weights=weights,noverlap=noverlap, &
                normalized=norm==1,average=avg==1)
        else
            out = fixed_taper_spectral_matrix_spectrogram(real(data,dp),tapers,freqs, &
                bands=bands,taper_weights=weights,noverlap=noverlap, &
                normalized=norm==1,average=avg==1)
        end if
        if (allocated(out%power)) then
            write(*,'(i0,1x,i0,1x,es25.16,1x)',advance='no') &
                out%status,size(out%power,2),out%effective_tapers
            do w = 1, size(out%power,2)
                do f = 1, size(out%power,1)
                    do j = 1, nc
                        do i = 1, nc
                            write(*,'(es25.16,1x,es25.16,1x)',advance='no') &
                                real(out%power(f,w,i,j),dp),aimag(out%power(f,w,i,j))
                        end do
                    end do
                end do
            end do
            print '(a)', ''
        else
            print '(i0,1x,i0,1x,es25.16)', out%status,0,0.0_dp
        end if
        deallocate(raw,data,tapers,weights,freqs,bands)
    end do
end program fixed_taper_spectrogram_driver
