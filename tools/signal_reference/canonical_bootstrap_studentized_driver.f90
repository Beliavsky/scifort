! SPDX-License-Identifier: MIT
! Text protocol for independent NumPy bootstrap-t with block jackknife SEs.
program canonical_bootstrap_studentized_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(canonical_moving_block_bootstrap_result) :: out
    real(dp), allocatable :: bands(:,:)
    integer, allocatable :: a(:), b(:), c(:)
    integer :: ncases, icase, nf, nt, nc, na, nb, ncond, nband, bl, nr, seed
    integer :: circular, f, t, i, k
    real(dp) :: ridge, shrinkage, tol

    read(*,*) ncases
    do icase = 1, ncases
        read(*,*) nf,nt,nc,na,nb,ncond,nband,bl,nr,seed,circular,ridge,shrinkage,tol
        allocate(a(na),b(nb),c(ncond),bands(2,nband))
        read(*,*) a
        read(*,*) b
        if (ncond > 0) read(*,*) c
        if (nband > 0) read(*,*) bands
        allocate(spectra%frequencies(nf),spectra%times(nt),spectra%power(nf,nt,nc,nc))
        read(*,*) spectra%frequencies
        read(*,*) spectra%times
        do f = 1,nf
            do t = 1,nt
                do i = 1,nc
                    read(*,*) spectra%power(f,t,i,:)
                end do
            end do
        end do
        spectra%status = 0
        if (nband > 0) then
            out = canonical_moving_block_bootstrap(spectra,a,b,bl,n_resamples=nr, &
                conditioned_on=c,ridge=ridge,shrinkage=shrinkage,bands=bands, &
                gap_tolerance=tol,seed=seed,circular=circular /= 0,interval_method='studentized')
        else
            out = canonical_moving_block_bootstrap(spectra,a,b,bl,n_resamples=nr, &
                conditioned_on=c,ridge=ridge,shrinkage=shrinkage, &
                gap_tolerance=tol,seed=seed,circular=circular /= 0,interval_method='studentized')
        end if
        if (out%status == 1) error stop 'unexpected invalid input in studentized reference'
        write(*,'(7(i0,1x))') out%status,size(out%frequencies),size(out%squared_modes,2), &
            out%n_windows,out%block_length,out%n_resamples,out%blocks_per_resample
        do f = 1,size(out%frequencies)
            write(*,'(i0,1x,6(es26.17e3,1x),i0)') f,out%coherence(f), &
                out%coherence_stderr(f),out%coherence_bias(f), &
                out%coherence_studentizing_stderr(f),out%coherence_lower95(f), &
                out%coherence_upper95(f),merge(1,0,out%coherence_valid(f))
            do k = 1,size(out%squared_modes,2)
                write(*,'(2(i0,1x),6(es26.17e3,1x),i0)') f,k, &
                    out%squared_modes(f,k),out%stderr(f,k),out%bias(f,k), &
                    out%studentizing_stderr(f,k),out%lower95(f,k),out%upper95(f,k), &
                    merge(1,0,out%valid(f,k))
            end do
        end do
        deallocate(a,b,c,bands,spectra%frequencies,spectra%times,spectra%power)
    end do
end program canonical_bootstrap_studentized_driver
