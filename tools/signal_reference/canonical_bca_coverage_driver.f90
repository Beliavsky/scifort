! SPDX-License-Identifier: MIT
! Exploratory VAR coverage diagnostics for alternative pointwise intervals.
program canonical_bca_coverage_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spec
    type(canonical_moving_block_bootstrap_result) :: percentile, basic, bc, bca, studentized
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
        percentile = canonical_moving_block_bootstrap(spec,[1],[2],block_length=4, &
            n_resamples=199,seed=seed,circular=.true.)
        basic = canonical_moving_block_bootstrap(spec,[1],[2],block_length=4, &
            n_resamples=199,seed=seed,circular=.true.,interval_method='basic')
        bc = canonical_moving_block_bootstrap(spec,[1],[2],block_length=4, &
            n_resamples=199,seed=seed,circular=.true.,interval_method='bc')
        bca = canonical_moving_block_bootstrap(spec,[1],[2],block_length=4, &
            n_resamples=199,seed=seed,circular=.true.,interval_method='bca')
        studentized = canonical_moving_block_bootstrap(spec,[1],[2],block_length=4, &
            n_resamples=199,seed=seed,circular=.true.,interval_method='studentized')
        if (percentile%status /= 0 .or. basic%status /= 0) &
            error stop 'unexpected invalid percentile/basic bootstrap'
        if (percentile%coherence(1) /= basic%coherence(1) .or. &
            percentile%coherence(1) /= bc%coherence(1)) &
            error stop 'interval methods changed full-sample estimate'
        write(*,'(10(es25.16,1x),3(i0,1x))') &
            percentile%coherence_lower95(1),percentile%coherence_upper95(1), &
            basic%coherence_lower95(1),basic%coherence_upper95(1), &
            bc%coherence_lower95(1),bc%coherence_upper95(1), &
            bca%coherence_lower95(1),bca%coherence_upper95(1), &
            studentized%coherence_lower95(1),studentized%coherence_upper95(1), &
            merge(1,0,bc%coherence_valid(1)), &
            merge(1,0,bca%coherence_valid(1)), &
            merge(1,0,studentized%coherence_valid(1))
        deallocate(spec%frequencies,spec%times,spec%power)
    end do
end program canonical_bca_coverage_driver
