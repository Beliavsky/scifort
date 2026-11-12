! SPDX-License-Identifier: MIT
! Fixed-seed coverage diagnostic; not a claim of calibrated bootstrap intervals.
program canonical_influence_coverage_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spec
    type(canonical_moving_block_bootstrap_result) :: spectral, influence, independent
    integer :: count, trial, nw, nc, seed, t, i, j, offset
    real(dp) :: raw(50)

    read(*,*) count
    do trial = 1, count
        read(*,*) nw,nc,seed
        allocate(spec%frequencies(1),spec%times(nw),spec%power(1,nw,nc,nc))
        spec%frequencies = [0.0_dp]
        spec%times = [(real(t,dp),t=1,nw)]
        spec%status = 0
        do t = 1, nw
            read(*,*) raw(1:2*nc*nc)
            do j = 1, nc
                do i = 1, nc
                    offset = 2*((j-1)*nc+i-1)+1
                    spec%power(1,t,i,j) = cmplx(raw(offset),raw(offset+1),dp)
                end do
            end do
        end do
        spectral = canonical_moving_block_bootstrap(spec,[1],[2], &
            n_resamples=99,seed=seed,circular=.true.)
        influence = canonical_moving_block_bootstrap(spec,[1],[2], &
            n_resamples=99,seed=seed,circular=.true.,selection_method='influence')
        independent = canonical_moving_block_bootstrap(spec,[1],[2],block_length=1, &
            n_resamples=99,seed=seed,circular=.true.)
        if (spectral%status /= 0 .or. influence%status /= 0 .or. independent%status /= 0) &
            error stop 'invalid AR spectral bootstrap fit'
        print '(2(i0,1x),7(es25.16,1x))', spectral%block_length,influence%block_length, &
            spectral%coherence(1),spectral%coherence_lower95(1),spectral%coherence_upper95(1), &
            influence%coherence_lower95(1),influence%coherence_upper95(1), &
            independent%coherence_lower95(1),independent%coherence_upper95(1)
        deallocate(spec%frequencies,spec%times,spec%power)
    end do
end program canonical_influence_coverage_driver
