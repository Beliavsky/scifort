! SPDX-License-Identifier: MIT
! Native deterministic fixed-DPSS-taper band estimator followed by all
! five existing conditional canonical moving-block bootstrap intervals.
! Input is a simulated real time series and a fixed taper bank, not a
! precomputed window spectral matrix. No FFT-grid band approximation.
program canonical_varp_multitaper_coverage_driver
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_fixed_taper_spectrogram, only : fixed_taper_spectral_matrix_spectrogram
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap
    use scifort_signal_partial_coherence, only : partial_coherence_invalid_input
    implicit none

    character(len=12), parameter :: methods(5) = [character(len=12) :: &
        'percentile','basic','bc','bca','studentized']
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(canonical_moving_block_bootstrap_result) :: result(5)
    real(dp), allocatable :: data(:,:), tapers(:,:), weights(:), frequencies(:), bands(:,:), row(:)
    real(dp) :: value
    integer :: count, trial, nw, length, nk, nc, nf, nb, block_length, nr, seed
    integer :: f, t, j, m, method, k, nm

    read(*,*) count
    if (count < 1) error stop 'invalid number of datasets'
    do trial = 1, count
        read(*,*) nw,nc,length,nk,nf,nb,block_length,nr,seed
        if (nw < 8 .or. nc /= 5 .or. length < 3 .or. nk < 1 .or. &
            nf+nb < 1 .or. block_length < 1 .or. block_length >= nw .or. nr < 20) &
            error stop 'invalid input dimensions'
        allocate(data(nw*length,nc),tapers(length,nk),weights(nk), &
            frequencies(nf),bands(2,nb))
        do j = 1, nc
            read(*,*) data(:,j)
        end do
        do j = 1, nk
            read(*,*) tapers(:,j)
        end do
        read(*,*) weights
        if (nf > 0) read(*,*) frequencies
        do f = 1, nb
            read(*,*) bands(:,f)
        end do
        spectra = fixed_taper_spectral_matrix_spectrogram(data,tapers,frequencies, &
            bands=bands,taper_weights=weights,normalized=.true.,average=.true.)
        if (spectra%status /= spectral_status_success) &
            error stop 'unexpected invalid fixed-taper spectral matrix'
        do method = 1, size(methods)
            result(method) = canonical_moving_block_bootstrap(spectra,[1,2],[3,4], &
                block_length=block_length,n_resamples=nr,seed=seed,circular=.true., &
                conditioned_on=[5],interval_method=trim(methods(method)))
            if (result(method)%status == partial_coherence_invalid_input) &
                error stop 'unexpected invalid canonical bootstrap input'
            if (.not. allocated(result(method)%squared_modes) .or. &
                .not. allocated(result(method)%coherence)) &
                error stop 'canonical outputs missing'
            if (method > 1) then
                if (any(result(method)%squared_modes /= result(1)%squared_modes) .or. &
                    any(result(method)%coherence /= result(1)%coherence)) &
                    error stop 'interval method altered canonical estimates'
                if (any(result(method)%stderr /= result(1)%stderr)) &
                    error stop 'interval method altered resamples'
            end if
        end do
        nm = size(result(1)%squared_modes,2)
        allocate(row((nf+nb)*(nm+1)*(1+3*size(methods))))
        k = 0
        do f = 1, nf+nb
            do m = 1, nm+1
                if (m <= nm) then
                    value = result(1)%squared_modes(f,m)
                else
                    value = result(1)%coherence(f)
                end if
                k = k+1
                row(k) = value
                do method = 1, size(methods)
                    if (m <= nm) then
                        row(k+1) = result(method)%lower95(f,m)
                        row(k+2) = result(method)%upper95(f,m)
                        row(k+3) = merge(1.0_dp,0.0_dp,result(method)%valid(f,m))
                    else
                        row(k+1) = result(method)%coherence_lower95(f)
                        row(k+2) = result(method)%coherence_upper95(f)
                        row(k+3) = merge(1.0_dp,0.0_dp,result(method)%coherence_valid(f))
                    end if
                    k = k+3
                end do
            end do
        end do
        if (k /= size(row)) error stop 'serialization size mismatch'
        write(*,'(*(es25.16,1x))') row
        deallocate(data,tapers,weights,frequencies,bands,row)
    end do
end program canonical_varp_multitaper_coverage_driver
