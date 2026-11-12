! SPDX-License-Identifier: MIT
! Text protocol driver for independent NumPy reference on window jackknife.
program canonical_jackknife_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_jackknife, only : canonical_window_jackknife, &
        canonical_window_jackknife_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: s
    type(canonical_window_jackknife_result) :: answer
    real(dp), allocatable :: bands(:,:)
    integer, allocatable :: a(:),b(:),c(:)
    integer :: ncases, icase, nf, nt, nc, na, nb, ncond, nband, f, t, i, j, k
    real(dp) :: ridge, shrinkage, tol

    read(*,*) ncases
    do icase = 1, ncases
        read(*,*) nf,nt,nc,na,nb,ncond,nband,ridge,shrinkage,tol
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
            answer = canonical_window_jackknife(s,a,b,conditioned_on=c, &
                ridge=ridge,shrinkage=shrinkage,bands=bands,gap_tolerance=tol)
        else
            answer = canonical_window_jackknife(s,a,b,conditioned_on=c, &
                ridge=ridge,shrinkage=shrinkage,gap_tolerance=tol)
        end if
        if (answer%status /= 0) error stop 'unexpected case failure'
        write(*,'(3(i0,1x))') size(answer%frequencies),size(answer%squared_modes,2),answer%n_windows
        do f = 1, size(answer%frequencies)
            write(*,'(i0,1x,2(es26.17e3,1x),i0)') f,answer%coherence(f), &
                answer%coherence_stderr(f),merge(1,0,answer%coherence_valid(f))
            do k = 1, size(answer%squared_modes,2)
                write(*,'(2(i0,1x),6(es26.17e3,1x),i0)') f,k, &
                    answer%squared_modes(f,k),answer%stderr(f,k), &
                    answer%stderr_fisher(f,k),answer%lower95(f,k), &
                    answer%upper95(f,k),answer%bias(f,k),merge(1,0,answer%valid(f,k))
            end do
        end do
        deallocate(a,b,c,bands,s%frequencies,s%times,s%power)
    end do
end program canonical_jackknife_driver
