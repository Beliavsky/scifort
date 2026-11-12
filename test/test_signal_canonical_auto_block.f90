! SPDX-License-Identifier: MIT
program test_signal_canonical_auto_block
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_block_selection, only : canonical_block_length_result, &
        select_canonical_block_length
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap, multitaper_canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spec, copy_spec
    type(canonical_block_length_result) :: selected, second, capped
    type(canonical_moving_block_bootstrap_result) :: automatic, explicit, direct
    integer, parameter :: nw = 80
    real(dp) :: latent, x(192,3)
    complex(dp) :: z(192,3), cross
    integer :: w, f, j, i

    allocate(spec%power(3,nw,3,3),spec%times(nw),spec%frequencies(3))
    spec%status = 0
    spec%frequencies = [0.0_dp,0.25_dp,0.50_dp]
    spec%times = [(real(w,dp),w=1,nw)]
    spec%power = (0.0_dp,0.0_dp)
    latent = 0.0_dp
    do w = 1, nw
        latent = 0.93_dp*latent + sin(0.17_dp*real(w,dp))
        do f = 1, 3
            spec%power(f,w,1,1) = cmplx(3.0_dp+0.10_dp*latent,0.0_dp,dp)
            spec%power(f,w,2,2) = cmplx(2.5_dp+0.16_dp*latent,0.0_dp,dp)
            spec%power(f,w,3,3) = cmplx(2.0_dp+0.12_dp*latent,0.0_dp,dp)
            cross = cmplx(0.30_dp+0.035_dp*latent,0.12_dp+0.013_dp*latent,dp)
            spec%power(f,w,1,2) = cross
            spec%power(f,w,2,1) = conjg(cross)
            spec%power(f,w,1,3) = (0.20_dp,0.0_dp)
            spec%power(f,w,3,1) = (0.20_dp,0.0_dp)
            spec%power(f,w,2,3) = (0.08_dp,0.0_dp)
            spec%power(f,w,3,2) = (0.08_dp,0.0_dp)
        end do
    end do
    selected = select_canonical_block_length(spec)
    call require(selected%status == 0,'selection succeeds')
    call require(selected%n_windows == nw .and. selected%max_lag > 0, 'pilot metadata')
    call require(selected%lags_used > 0 .and. selected%first_lag_correlation > 0.5_dp, &
        'serial dependence detected by spectral pilot')
    call require(selected%block_length > 1 .and. selected%block_length <= nw/3, &
        'dependent data receive a longer bounded block')
    call require(all(abs(selected%autocorrelations) <= 1.01_dp), &
        'normalized short-lag spectral correlations')

    copy_spec = spec
    copy_spec%power = copy_spec%power*1.0e50_dp
    second = select_canonical_block_length(copy_spec)
    call require(second%status == 0 .and. second%block_length == selected%block_length, &
        'common spectral rescaling preserves block selection')
    copy_spec%power = spec%power*1.0e-50_dp
    second = select_canonical_block_length(copy_spec)
    call require(second%status == 0 .and. second%block_length == selected%block_length, &
        'small spectral magnitudes preserve block selection')
    capped = select_canonical_block_length(spec,max_block_length=2)
    call require(capped%block_length == 2 .and. capped%maximum_block_length == 2, &
        'caller cap is obeyed')
    capped = select_canonical_block_length(spec,max_lag=1)
    call require(capped%max_lag == 1 .and. capped%lags_used <= 1, &
        'pilot lag cap is obeyed')

    automatic = canonical_moving_block_bootstrap(spec,[1],[2],n_resamples=49,seed=321)
    explicit = canonical_moving_block_bootstrap(spec,[1],[2],selected%block_length, &
        n_resamples=49,seed=321)
    call require(automatic%status == 0 .and. explicit%status == 0, &
        'automatic and fixed resampling fit')
    call require(automatic%automatic_block_length .and. .not. explicit%automatic_block_length, &
        'automatic selection metadata')
    call require(automatic%selection_lags == selected%lags_used .and. &
        automatic%block_length == selected%block_length, 'selector feeds resampler')
    call require(all(automatic%squared_modes == explicit%squared_modes) .and. &
        all(automatic%coherence_stderr == explicit%coherence_stderr) .and. &
        all(automatic%lower95 == explicit%lower95), &
        'selected and explicitly fixed block length produce identical seeded results')

    copy_spec = spec
    copy_spec%power = cmplx(2.0_dp,0.0_dp,dp)
    do i = 1, 3
        do j = 1, 3
            if (i /= j) copy_spec%power(:,:,i,j) = (0.0_dp,0.0_dp)
        end do
    end do
    second = select_canonical_block_length(copy_spec)
    call require(second%status == 0 .and. second%block_length == 1 .and. &
        second%lags_used == 0, 'constant window spectra use length one')

    copy_spec = spec
    copy_spec%times(3) = copy_spec%times(2)
    second = select_canonical_block_length(copy_spec)
    call require(second%status == 1,'invalid window times rejected')
    second = select_canonical_block_length(spec,max_lag=nw)
    call require(second%status == 1,'invalid pilot lag cap rejected')
    second = select_canonical_block_length(spec,max_block_length=nw)
    call require(second%status == 1,'invalid block cap rejected')
    copy_spec = spec
    copy_spec%power(1,3,1,1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    second = select_canonical_block_length(copy_spec)
    call require(second%status == 1,'nonfinite spectra rejected')

    do i = 1, 192
        do j = 1, 3
            x(i,j) = sin(0.08_dp*i*j)+0.3_dp*cos(0.11_dp*i+0.17_dp*j)
            z(i,j) = cmplx(x(i,j),0.12_dp*sin(0.19_dp*i*j),dp)
        end do
    end do
    direct = multitaper_canonical_moving_block_bootstrap(x,3.0_dp,16,[1],[2], &
        n_resamples=20,seed=234,ridge=0.02_dp,kmax=4,weighting='equal')
    call require(direct%status == 0 .and. direct%automatic_block_length .and. &
        direct%block_length >= 1,'direct real automatically selected bootstrap')
    direct = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,16,[1],[2], &
        n_resamples=20,seed=235,ridge=0.02_dp,kmax=4,weighting='adaptive', &
        circular=.true.,bands=reshape([0.01_dp,0.12_dp],[2,1]))
    call require(direct%status == 0 .and. direct%automatic_block_length .and. &
        direct%block_length >= 1,'direct complex adaptive band automatically selected bootstrap')
    print '(a)', 'PASS automatic canonical bootstrap block selection and native regressions'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! checked proposition
        character(*), intent(in) :: message !! test failure description
        if (.not. ok) then
            print '(a)', 'FAIL '//message
            error stop 1
        end if
    end subroutine require
end program test_signal_canonical_auto_block
