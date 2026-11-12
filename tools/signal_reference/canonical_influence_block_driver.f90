! SPDX-License-Identifier: MIT
! Text protocol for independent NumPy testing of statistic-specific block pilot.
program canonical_influence_block_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_block_selection, only : canonical_block_length_result, &
        select_canonical_influence_block_length
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spec
    type(canonical_block_length_result) :: fit
    integer :: ncases, caseid, nw, nf, nc, mode, frequency, cond, band, &
        lag, cap, f, t, i, j, offset
    real(dp) :: ridge, shrinkage, raw(50)

    read(*,*) ncases
    do caseid = 1, ncases
        read(*,*) nw,nf,nc,mode,frequency,cond,band,lag,cap,ridge,shrinkage
        allocate(spec%times(nw),spec%frequencies(nf),spec%power(nf,nw,nc,nc))
        read(*,*) spec%frequencies
        read(*,*) spec%times
        do t = 1, nw
            do f = 1, nf
                read(*,*) raw(1:2*nc*nc)
                do j = 1, nc
                    do i = 1, nc
                        offset = 2*((j-1)*nc+i-1)+1
                        spec%power(f,t,i,j) = cmplx(raw(offset),raw(offset+1),dp)
                    end do
                end do
            end do
        end do
        spec%status = 0
        if (cond == 0 .and. band == 0) then
            fit = select_canonical_influence_block_length(spec,[1,2],[3,4],ridge=ridge, &
                shrinkage=shrinkage,mode_index=mode,frequency_index=frequency, &
                max_lag=lag,max_block_length=cap)
        else if (cond == 1 .and. band == 0) then
            fit = select_canonical_influence_block_length(spec,[1,2],[3,4],conditioned_on=[5], &
                ridge=ridge,shrinkage=shrinkage,mode_index=mode,frequency_index=frequency, &
                max_lag=lag,max_block_length=cap)
        else if (cond == 0 .and. band == 1) then
            fit = select_canonical_influence_block_length(spec,[1,2],[3,4],ridge=ridge, &
                shrinkage=shrinkage,bands=reshape([0.05_dp,0.42_dp],[2,1]), &
                average=.true.,mode_index=mode,frequency_index=frequency, &
                max_lag=lag,max_block_length=cap)
        else
            fit = select_canonical_influence_block_length(spec,[1,2],[3,4],conditioned_on=[5], &
                ridge=ridge,shrinkage=shrinkage,bands=reshape([0.05_dp,0.42_dp],[2,1]), &
                average=.true.,mode_index=mode,frequency_index=frequency, &
                max_lag=lag,max_block_length=cap)
        end if
        if (fit%status /= 0) error stop 'reference influence selection failed'
        print '(3(i0,1x),2(es26.17,1x))', fit%status,fit%block_length,fit%lags_used, &
            fit%full_statistic,fit%first_lag_correlation
        print '(*(es26.17,1x))', fit%influence
        print '(*(es26.17,1x))', fit%autocorrelations
        deallocate(spec%times,spec%frequencies,spec%power)
    end do
end program canonical_influence_block_driver
