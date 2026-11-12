! SPDX-License-Identifier: MIT
! Text protocol for independent NumPy/SciPy moving-block bootstrap reference.
program canonical_bootstrap_intervals_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: s
    type(canonical_moving_block_bootstrap_result) :: out
    real(dp), allocatable :: bands(:,:)
    integer, allocatable :: a(:), b(:), c(:)
    integer :: ncases, icase, nf, nt, nc, na, nb, ncond, nband, bl, nr, seed, wrap, f, t, i, k
    real(dp) :: ridge, shrinkage, tol
    character(len=10) :: method

    read(*,*) ncases
    do icase = 1, ncases
        read(*,*) nf,nt,nc,na,nb,ncond,nband,bl,nr,seed,wrap,ridge,shrinkage,tol,method
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
            out = canonical_moving_block_bootstrap(s,a,b,bl,n_resamples=nr, &
                conditioned_on=c,ridge=ridge,shrinkage=shrinkage,bands=bands, &
                gap_tolerance=tol,seed=seed,circular=wrap /= 0,interval_method=method)
        else
            out = canonical_moving_block_bootstrap(s,a,b,bl,n_resamples=nr, &
                conditioned_on=c,ridge=ridge,shrinkage=shrinkage, &
                gap_tolerance=tol,seed=seed,circular=wrap /= 0,interval_method=method)
        end if
        if (out%status == 1) error stop 'unexpected invalid canonical bootstrap input'
        write(*,'(8(i0,1x))') out%status, size(out%frequencies),size(out%squared_modes,2), &
            out%n_windows,out%block_length,out%n_resamples,out%blocks_per_resample, &
            merge(1,0,out%circular)
        do f = 1, size(out%frequencies)
            write(*,'(i0,1x,5(es26.17e3,1x),i0)') f,out%coherence(f), &
                out%coherence_stderr(f),out%coherence_bias(f), &
                out%coherence_lower95(f),out%coherence_upper95(f), &
                merge(1,0,out%coherence_valid(f))
            do k = 1, size(out%squared_modes,2)
                write(*,'(2(i0,1x),5(es26.17e3,1x),i0)') f,k, &
                    out%squared_modes(f,k),out%stderr(f,k),out%bias(f,k), &
                    out%lower95(f,k),out%upper95(f,k),merge(1,0,out%valid(f,k))
            end do
        end do
        deallocate(a,b,c,bands,s%frequencies,s%times,s%power)
    end do
end program canonical_bootstrap_intervals_driver
