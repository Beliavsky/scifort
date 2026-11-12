! SPDX-License-Identifier: MIT
! Exchange spectral-window matrices with independent NumPy validation.
program canonical_auto_block_driver
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, iostat_end
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_block_selection, only : canonical_block_length_result, &
        select_canonical_block_length
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(canonical_block_length_result) :: choice
    type(canonical_moving_block_bootstrap_result) :: auto, independent
    integer :: iteration, nf, nw, nc, mode, seed, r, f, w, i, j, ios
    real(dp) :: re, im

    read(input_unit,*,iostat=ios) iteration
    if (ios /= 0) error stop 'missing case count'
    do r = 1, iteration
        read(input_unit,*,iostat=ios) mode, nf, nw, nc, seed
        if (ios /= 0) error stop 'missing case header'
        allocate(spectra%frequencies(nf),spectra%times(nw),spectra%power(nf,nw,nc,nc))
        spectra%frequencies = [(real(f-1,dp)/real(nf,dp),f=1,nf)]
        spectra%times = [(real(w,dp),w=1,nw)]
        spectra%status = 0
        do f = 1, nf
            do w = 1, nw
                do i = 1, nc
                    do j = 1, nc
                        read(input_unit,*,iostat=ios) re, im
                        if (ios /= 0) error stop 'missing spectral matrix component'
                        spectra%power(f,w,i,j) = cmplx(re,im,dp)
                    end do
                end do
            end do
        end do
        choice = select_canonical_block_length(spectra)
        if (choice%status /= 0) error stop 'selection unsuccessful'
        if (mode == 0) then
            write(output_unit,'(i0,1x,i0,1x,i0,1x,i0,1x,es24.16,1x,es24.16,1x,es24.16)') &
                choice%block_length, choice%lags_used, choice%max_lag, &
                choice%maximum_block_length, choice%first_lag_correlation, &
                choice%long_run_variance_ratio, choice%weighted_lag_ratio
        else
            auto = canonical_moving_block_bootstrap(spectra,[1],[2],n_resamples=99,seed=seed,circular=.true.)
            independent = canonical_moving_block_bootstrap(spectra,[1],[2],block_length=1, &
                n_resamples=99,seed=seed,circular=.true.)
            if (auto%status /= 0 .or. independent%status /= 0) error stop 'bootstrap unsuccessful'
            write(output_unit,'(i0,1x,5(es24.16,1x),i0,1x,2(es24.16,1x),i0)') &
                choice%block_length, auto%coherence(1), &
                auto%coherence_lower95(1),auto%coherence_upper95(1), &
                auto%coherence_stderr(1), independent%coherence_stderr(1), &
                merge(1,0,auto%coherence_valid(1)), &
                independent%coherence_lower95(1),independent%coherence_upper95(1), &
                merge(1,0,independent%coherence_valid(1))
        end if
        deallocate(spectra%frequencies,spectra%times,spectra%power)
    end do
end program canonical_auto_block_driver
