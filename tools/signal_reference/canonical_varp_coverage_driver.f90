! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Compare all five canonical moving-block bootstrap interval methods for
! multiple frequencies, two canonical roots and a conditioned channel.
! This is a numerical-study driver, not a public Fortran library API.
program canonical_varp_coverage_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap
    use scifort_signal_partial_coherence, only : partial_coherence_invalid_input
    implicit none

    character(len=12), parameter :: methods(5) = [character(len=12) :: &
        'percentile', 'basic', 'bc', 'bca', 'studentized']
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(canonical_moving_block_bootstrap_result) :: result(5)
    real(dp), allocatable :: raw(:), row(:)
    real(dp) :: value
    integer :: count, trial, nw, nc, nf, length, nr, seed, f, t, i, j, m, method, k, nm

    read(*,*) count
    if (count < 1) error stop 'invalid number of datasets'
    do trial = 1, count
        read(*,*) nw, nc, nf, length, nr, seed
        if (nw < 8 .or. nc /= 5 .or. nf < 1 .or. length < 1 .or. &
            length >= nw .or. nr < 20) error stop 'invalid dataset dimensions'
        allocate(raw(2*nc*nc), spectra%frequencies(nf), spectra%times(nw), &
            spectra%power(nf,nw,nc,nc))
        read(*,*) spectra%frequencies
        spectra%times = [(real(t,dp), t=1,nw)]
        spectra%status = 0
        do t = 1, nw
            do f = 1, nf
                read(*,*) raw
                do j = 1, nc
                    do i = 1, nc
                        k = 2*((j-1)*nc+i-1)+1
                        spectra%power(f,t,i,j) = cmplx(raw(k),raw(k+1),kind=dp)
                    end do
                end do
            end do
        end do

        do method = 1, size(methods)
            result(method) = canonical_moving_block_bootstrap(spectra,[1,2],[3,4], &
                block_length=length, n_resamples=nr, seed=seed, circular=.true., &
                conditioned_on=[5], interval_method=trim(methods(method)))
            if (result(method)%status == partial_coherence_invalid_input) &
                error stop 'unexpected invalid canonical input'
            if (.not. allocated(result(method)%squared_modes) .or. &
                .not. allocated(result(method)%coherence)) &
                error stop 'canonical outputs not allocated'
            if (method > 1) then
                if (any(result(method)%squared_modes /= result(1)%squared_modes) .or. &
                    any(result(method)%coherence /= result(1)%coherence)) &
                    error stop 'interval method changed point estimate'
                if (any(result(method)%stderr /= result(1)%stderr)) &
                    error stop 'interval method changed resamples'
            end if
        end do
        nm = size(result(1)%squared_modes,2)
        allocate(row(nf*(nm+1)*(1+3*size(methods))))
        k = 0
        do f = 1, nf
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
                        k = k+1
                        row(k) = result(method)%lower95(f,m)
                        k = k+1
                        row(k) = result(method)%upper95(f,m)
                        k = k+1
                        row(k) = merge(1.0_dp,0.0_dp,result(method)%valid(f,m))
                    else
                        k = k+1
                        row(k) = result(method)%coherence_lower95(f)
                        k = k+1
                        row(k) = result(method)%coherence_upper95(f)
                        k = k+1
                        row(k) = merge(1.0_dp,0.0_dp,result(method)%coherence_valid(f))
                    end if
                end do
            end do
        end do
        if (k /= size(row)) error stop 'internal output size mismatch'
        write(*,'(*(es25.16,1x))') row
        deallocate(raw,row,spectra%frequencies,spectra%times,spectra%power)
    end do
end program canonical_varp_coverage_driver
