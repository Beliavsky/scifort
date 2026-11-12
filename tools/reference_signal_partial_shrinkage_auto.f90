! SPDX-License-Identifier: MIT
! Reads complex Hermitian spectral matrices and emits regularized precision,
! complex partial coherency and automatic shrinkage for NumPy comparison.
program reference_signal_partial_shrinkage_auto
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_partial_coherence, only : partial_coherence_tf_matrix, partial_coherence_tf_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: s
    type(partial_coherence_tf_result) :: r
    integer :: cases, k, nf, nt, nc, i, j, f, t
    real(dp) :: effective, ridge, re, im
    read(*,*) cases
    write(*,'(i0)') cases
    do k = 1, cases
        read(*,*) nf,nt,nc,effective,ridge
        allocate(s%power(nf,nt,nc,nc),s%frequencies(nf),s%times(nt))
        s%frequencies = [(real(f-1,dp),f=1,nf)]
        s%times = [(real(t-1,dp),t=1,nt)]
        s%status = 0
        do t = 1,nt
            do f = 1,nf
                do j = 1,nc
                    do i = 1,nc
                        read(*,*) re,im
                        s%power(f,t,i,j) = cmplx(re,im,dp)
                    end do
                end do
            end do
        end do
        r = partial_coherence_tf_matrix(s, auto_shrinkage=.true., &
            effective_tapers=effective, ridge=ridge, return_coherency=.true.)
        write(*,'(i0,1x,i0,1x,i0,1x,i0)') r%status,nf,nt,nc
        do t = 1,nt
            do f = 1,nf
                write(*,'(i0,1x,es27.18e3)') merge(1,0,r%valid(f,t)),r%shrinkage_used(f,t)
                if (.not. r%valid(f,t)) cycle
                do j = 1,nc
                    do i = 1,nc
                        write(*,'(5(es27.18e3,1x))') r%coherence(f,t,i,j), &
                            real(r%precision(f,t,i,j),dp), aimag(r%precision(f,t,i,j)), &
                            real(r%coherency(f,t,i,j),dp), aimag(r%coherency(f,t,i,j))
                    end do
                end do
            end do
        end do
        deallocate(s%power,s%frequencies,s%times)
    end do
end program reference_signal_partial_shrinkage_auto
