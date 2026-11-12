! SPDX-License-Identifier: MIT
! Numerical oracle driver for verify_band_coherence.py; stdin -> stdout.
program band_coherence_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_band_coherence, only : integrate_spectral_bands, &
        integrate_tf_spectral_bands, band_multiple_coherence_matrix, &
        band_multiple_coherence_tf_matrix, band_partial_coherence_matrix, &
        band_partial_coherence_tf_matrix
    use scifort_signal_multiple_coherence, only : multiple_coherence_result, multiple_coherence_tf_result
    use scifort_signal_partial_coherence, only : partial_coherence_result, partial_coherence_tf_result
    implicit none
    type(multitaper_matrix_result) :: s, sb
    type(multitaper_spectral_matrix_spectrogram_result) :: st, stb
    type(multiple_coherence_result) :: r
    type(multiple_coherence_tf_result) :: rt
    type(partial_coherence_result) :: p
    type(partial_coherence_tf_result) :: pt
    real(dp), allocatable :: bands(:,:)
    real(dp) :: re, im, ridge, shrinkage
    integer :: ncase, icase, nf, nc, nb, nt, avg, k, i, j, t, ib

    read(*,*) ncase
    do icase = 1, ncase
        read(*,*) nf,nc,nb,nt,ridge,shrinkage,avg
        allocate(s%frequencies(nf),s%power(nf,nc,nc),bands(2,nb))
        allocate(st%frequencies(nf),st%times(nt),st%power(nf,nt,nc,nc))
        s%status = 0
        st%status = 0
        read(*,*) s%frequencies
        st%frequencies = s%frequencies
        do ib = 1, nb
            read(*,*) bands(:,ib)
        end do
        read(*,*) st%times
        do k = 1, nf
            do j = 1, nc
                do i = 1, nc
                    read(*,*) re,im
                    s%power(k,i,j) = cmplx(re,im,dp)
                end do
            end do
        end do
        do t = 1, nt
            do k = 1, nf
                do j = 1, nc
                    do i = 1, nc
                        read(*,*) re,im
                        st%power(k,t,i,j) = cmplx(re,im,dp)
                    end do
                end do
            end do
        end do
        sb = integrate_spectral_bands(s,bands,average=avg == 1)
        stb = integrate_tf_spectral_bands(st,bands,average=avg == 1)
        r = band_multiple_coherence_matrix(s,bands,average=avg == 1, &
            ridge=ridge,shrinkage=shrinkage)
        rt = band_multiple_coherence_tf_matrix(st,bands,average=avg == 1, &
            ridge=ridge,shrinkage=shrinkage)
        p = band_partial_coherence_matrix(s,bands,average=avg == 1, &
            ridge=ridge,shrinkage=shrinkage,return_coherency=.true.)
        pt = band_partial_coherence_tf_matrix(st,bands,average=avg == 1, &
            ridge=ridge,shrinkage=shrinkage,return_coherency=.true.)
        if (sb%status /= 0 .or. stb%status /= 0 .or. r%status /= 0 .or. &
            rt%status /= 0 .or. p%status /= 0 .or. pt%status /= 0) error stop 'driver invalid'
        print '(a,i0)', 'CASE ',icase
        do ib = 1, nb
            do j = 1, nc
                do i = 1, nc
                    call emit_complex(sb%power(ib,i,j))
                end do
            end do
        end do
        do t = 1, nt
            do ib = 1, nb
                do j = 1, nc
                    do i = 1, nc
                        call emit_complex(stb%power(ib,t,i,j))
                    end do
                end do
            end do
        end do
        do ib = 1, nb
            do i = 1, nc
                print '(3(es26.17,1x))',r%coherence(ib,i), &
                    r%conditional_power(ib,i),r%explained_power(ib,i)
            end do
            do j = 1, nc
                do i = 1, nc
                    print '(es26.17,1x)',p%coherence(ib,i,j)
                    call emit_complex(p%coherency(ib,i,j))
                end do
            end do
        end do
        do t = 1, nt
            do ib = 1, nb
                do i = 1, nc
                    print '(3(es26.17,1x))',rt%coherence(ib,t,i), &
                        rt%conditional_power(ib,t,i),rt%explained_power(ib,t,i)
                end do
                do j = 1, nc
                    do i = 1, nc
                        print '(es26.17,1x)',pt%coherence(ib,t,i,j)
                        call emit_complex(pt%coherency(ib,t,i,j))
                    end do
                end do
            end do
        end do
        deallocate(s%frequencies,s%power,bands,st%frequencies,st%times,st%power)
    end do
contains
    subroutine emit_complex(z)
        complex(dp), intent(in) :: z !! complex value to print as two real components
        print '(2(es26.17,1x))', real(z,dp),aimag(z)
    end subroutine emit_complex
end program band_coherence_driver
