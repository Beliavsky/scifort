! SPDX-License-Identifier: MIT
! Text protocol for NumPy/SciPy independent delete-group reference.
program canonical_block_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_jackknife, only : canonical_block_jackknife, canonical_window_jackknife_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: s
    type(canonical_window_jackknife_result) :: out
    real(dp), allocatable :: bands(:,:)
    integer, allocatable :: a(:), b(:), c(:)
    integer :: ncases, icase, nf, nt, nc, na, nb, ncond, nband, bl, f, t, i, k
    real(dp) :: ridge, shrinkage, tol

    read(*,*) ncases
    do icase = 1, ncases
        read(*,*) nf,nt,nc,na,nb,ncond,nband,bl,ridge,shrinkage,tol
        allocate(a(na),b(nb),c(ncond),bands(2,nband))
        read(*,*) a
        read(*,*) b
        if (ncond > 0) read(*,*) c
        if (nband > 0) read(*,*) bands
        allocate(s%frequencies(nf),s%times(nt),s%power(nf,nt,nc,nc))
        read(*,*) s%frequencies
        read(*,*) s%times
        do f = 1, nf
            do t = 1, nt
                do i = 1, nc
                    read(*,*) s%power(f,t,i,:)
                end do
            end do
        end do
        s%status = 0
        if (nband > 0) then
            out = canonical_block_jackknife(s,a,b,bl,conditioned_on=c, &
                ridge=ridge,shrinkage=shrinkage,bands=bands,gap_tolerance=tol)
        else
            out = canonical_block_jackknife(s,a,b,bl,conditioned_on=c, &
                ridge=ridge,shrinkage=shrinkage,gap_tolerance=tol)
        end if
        if (out%status /= 0) error stop 'unexpected block jackknife failure'
        write(*,'(5(i0,1x))') size(out%frequencies),size(out%squared_modes,2), &
            out%n_windows,out%n_blocks,out%block_length
        do f = 1, size(out%frequencies)
            write(*,'(i0,1x,2(es26.17e3,1x),i0)') f,out%coherence(f), &
                out%coherence_stderr(f),merge(1,0,out%coherence_valid(f))
            do k = 1, size(out%squared_modes,2)
                write(*,'(2(i0,1x),6(es26.17e3,1x),i0)') f,k, &
                    out%squared_modes(f,k),out%stderr(f,k), &
                    out%stderr_fisher(f,k),out%lower95(f,k), &
                    out%upper95(f,k),out%bias(f,k),merge(1,0,out%valid(f,k))
            end do
        end do
        deallocate(a,b,c,bands,s%frequencies,s%times,s%power)
    end do
end program canonical_block_driver
