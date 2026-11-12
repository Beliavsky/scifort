! SPDX-License-Identifier: MIT
program test_signal_canonical_influence_block
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_block_selection, only : canonical_block_length_result, &
        select_canonical_block_length, select_canonical_influence_block_length
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap, multitaper_canonical_moving_block_bootstrap
    implicit none
    integer, parameter :: nw = 96
    type(multitaper_spectral_matrix_spectrogram_result) :: spec, bad
    type(canonical_block_length_result) :: spectral, influence, target, chosen, invalid
    type(canonical_moving_block_bootstrap_result) :: automatic, explicit, direct
    real(dp) :: pi, latent, rho, exact, loo, x(288,3)
    complex(dp) :: cross, z(288,3)
    integer :: t, f, i, j

    pi = acos(-1.0_dp)
    allocate(spec%frequencies(3),spec%times(nw),spec%power(3,nw,3,3))
    spec%status = 0
    spec%frequencies = [0.0_dp,0.25_dp,0.50_dp]
    spec%times = [(real(t,dp),t=1,nw)]
    spec%power = (0.0_dp,0.0_dp)
    latent = 0.0_dp
    do t = 1, nw
        latent = 0.97_dp*latent+sin(0.05_dp*real(t,dp))
        rho = 0.37_dp+0.08_dp*sin(0.5_dp*pi*real(t,dp))
        cross = cmplx(rho,0.03_dp,dp)
        do f = 1, 3
            spec%power(f,t,1,1) = (1.0_dp,0.0_dp)
            spec%power(f,t,2,2) = (1.0_dp,0.0_dp)
            spec%power(f,t,3,3) = cmplx(250.0_dp+10.0_dp*latent,0.0_dp,dp)
            spec%power(f,t,1,2) = cross
            spec%power(f,t,2,1) = conjg(cross)
        end do
    end do
    spectral = select_canonical_block_length(spec)
    influence = select_canonical_influence_block_length(spec,[1],[2])
    call require(spectral%status == 0 .and. influence%status == 0,'both selectors succeed')
    call require(spectral%block_length > 1 .and. influence%block_length == 1, &
        'unrelated slowly varying channel does not control influence length')
    call require(abs(influence%first_lag_correlation) < 0.15_dp,'target influence lag one nearly zero')
    call require(size(influence%influence) == nw .and. influence%mode_index == 0 .and. &
        influence%frequency_index == 0,'diagnostics and target metadata')
    call require(abs(influence%full_statistic-(0.37_dp**2+0.03_dp**2)) < 1.0e-12_dp, &
        'full-sample group coherence')
    do t = 1, nw
        rho = 0.37_dp+0.08_dp*sin(0.5_dp*pi*real(t,dp))
        loo = (real(nw,dp)*0.37_dp-rho)/real(nw-1,dp)
        exact = real(nw-1,dp)*(0.37_dp**2-loo**2)
        call require(abs(influence%influence(t)-exact) < 2.0e-11_dp, &
            'delete-one influence matches analytic Hermitian coherence')
    end do
    target = select_canonical_influence_block_length(spec,[1],[2],mode_index=1,frequency_index=2)
    call require(target%status == 0 .and. target%block_length == influence%block_length .and. &
        all(abs(target%influence-influence%influence) < 1.0e-11_dp), &
        'singleton canonical root agrees with group coherence')
    chosen = select_canonical_influence_block_length(spec,[1],[2],max_lag=1,max_block_length=2)
    call require(chosen%status == 0 .and. chosen%max_lag == 1 .and. &
        chosen%maximum_block_length == 2,'pilot controls forwarded')
    chosen = select_canonical_influence_block_length(spec,[1],[2], &
        bands=reshape([0.02_dp,0.48_dp],[2,1]),frequency_index=1)
    call require(chosen%status == 0 .and. chosen%block_length == influence%block_length, &
        'band-integrated influence selection')

    automatic = canonical_moving_block_bootstrap(spec,[1],[2],n_resamples=31,seed=321, &
        selection_method='influence')
    explicit = canonical_moving_block_bootstrap(spec,[1],[2],influence%block_length, &
        n_resamples=31,seed=321)
    call require(automatic%status == 0 .and. explicit%status == 0, &
        'automatic influence and fixed bootstrap fits')
    call require(automatic%automatic_block_length .and. automatic%selection_method == 'influence' &
        .and. automatic%selection_lags == influence%lags_used .and. &
        automatic%block_length == influence%block_length, 'influence selector metadata')
    call require(all(automatic%selection_influence == influence%influence) .and. &
        all(automatic%coherence_stderr == explicit%coherence_stderr) .and. &
        all(automatic%coherence_lower95 == explicit%coherence_lower95), &
        'selected and explicitly fixed block lengths give identical seeded results')
    automatic = canonical_moving_block_bootstrap(spec,[1],[2],n_resamples=25,seed=321, &
        selection_method='influence',selection_mode=1,selection_frequency=2)
    call require(automatic%status == 0 .and. automatic%selection_mode == 1 .and. &
        automatic%selection_frequency == 2,'target-mode and frequency API')

    invalid = select_canonical_influence_block_length(spec,[1],[2],mode_index=2)
    call require(invalid%status == 1,'invalid mode rejected')
    invalid = select_canonical_influence_block_length(spec,[1],[2],frequency_index=4)
    call require(invalid%status == 1,'invalid frequency rejected')
    bad = spec
    bad%power(1,1,1,1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    invalid = select_canonical_influence_block_length(bad,[1],[2])
    call require(invalid%status == 1,'nonfinite spectral matrices rejected')
    bad = spec
    bad%times(2) = bad%times(1)
    invalid = select_canonical_influence_block_length(bad,[1],[2])
    call require(invalid%status == 1,'nonincreasing times rejected')
    automatic = canonical_moving_block_bootstrap(spec,[1],[2],n_resamples=20, &
        selection_method='influenc')
    call require(automatic%status == 1,'unknown method rejected')
    automatic = canonical_moving_block_bootstrap(spec,[1],[2],block_length=2,n_resamples=20, &
        selection_method='influence')
    call require(automatic%status == 1,'automatic influence method invalid with manual length')
    automatic = canonical_moving_block_bootstrap(spec,[1],[2],n_resamples=20, &
        selection_method='spectral',selection_mode=1)
    call require(automatic%status == 1,'target mode invalid for spectral selection')

    do i = 1, 288
        do j = 1, 3
            x(i,j) = sin(0.075_dp*i*j)+0.24_dp*cos(0.038_dp*i+0.17_dp*j)
            z(i,j) = cmplx(x(i,j),0.15_dp*sin(0.055_dp*i*j),dp)
        end do
    end do
    direct = multitaper_canonical_moving_block_bootstrap(x,3.0_dp,24,[1],[2], &
        n_resamples=20,seed=3123,kmax=4,weighting='equal',ridge=0.05_dp, &
        selection_method='influence')
    call require(direct%status == 0 .and. direct%automatic_block_length .and. &
        direct%selection_method == 'influence','real-signal influence API')
    direct = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,24,[1],[2], &
        n_resamples=20,seed=3123,kmax=4,weighting='adaptive',ridge=0.05_dp, &
        bands=reshape([0.01_dp,0.20_dp],[2,1]),selection_method='influence')
    call require(direct%status == 0 .and. direct%automatic_block_length .and. &
        size(direct%frequencies) == 1,'complex adaptive band influence API')
    print '(a)', 'PASS influence-targeted canonical block selection native regressions'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! tested proposition
        character(*), intent(in) :: message !! failure reason

        if (.not. ok) then
            print '(a)', 'FAIL '//message
            error stop 1
        end if
    end subroutine require
end program test_signal_canonical_influence_block
