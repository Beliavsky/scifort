! SPDX-License-Identifier: MIT
program test_signal_canonical_interval_methods
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap, multitaper_canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra, zero_spectra, extreme
    type(canonical_moving_block_bootstrap_result) :: original, percentile, basic, bc, invalid, direct
    real(dp) :: x(384,3), expected_low, expected_high
    real(dp), parameter :: correlations(6) = &
        [0.83_dp,-0.74_dp,0.42_dp,-0.39_dp,0.28_dp,-0.40_dp]
    complex(dp) :: factor(3,3), z(384,3)
    integer :: f, t, i, j

    allocate(spectra%frequencies(3),spectra%times(18),spectra%power(3,18,3,3))
    spectra%frequencies = [0.0_dp,0.25_dp,0.5_dp]
    spectra%times = [(real(t,dp),t=1,18)]
    spectra%status = 0
    do f = 1, 3
        do t = 1, 18
            factor = (0.0_dp,0.0_dp)
            do i = 1, 3
                factor(i,i) = cmplx(1.0_dp+0.09_dp*i+0.07_dp*t,0.0_dp,dp)
                do j = 1, i-1
                    factor(i,j) = cmplx(0.31_dp*cos(0.21_dp*i*j+0.29_dp*t), &
                        0.22_dp*sin(0.15_dp*j*f+0.12_dp*t),dp)
                end do
            end do
            spectra%power(f,t,:,:) = matmul(factor,transpose(conjg(factor)))*(0.7_dp+0.08_dp*f)
        end do
    end do

    original = canonical_moving_block_bootstrap(spectra,[1],[2],3,n_resamples=199, &
        conditioned_on=[3],ridge=0.03_dp,shrinkage=0.1_dp,seed=2468)
    percentile = canonical_moving_block_bootstrap(spectra,[1],[2],3,n_resamples=199, &
        conditioned_on=[3],ridge=0.03_dp,shrinkage=0.1_dp,seed=2468, &
        interval_method='percentile')
    basic = canonical_moving_block_bootstrap(spectra,[1],[2],3,n_resamples=199, &
        conditioned_on=[3],ridge=0.03_dp,shrinkage=0.1_dp,seed=2468,interval_method='basic')
    bc = canonical_moving_block_bootstrap(spectra,[1],[2],3,n_resamples=199, &
        conditioned_on=[3],ridge=0.03_dp,shrinkage=0.1_dp,seed=2468,interval_method='bc')
    call require(original%status == 0 .and. percentile%status == 0 .and. &
        basic%status == 0,'percentile and basic successful')
    call require(original%interval_method == 'percentile' .and. &
        basic%interval_method == 'basic' .and. bc%interval_method == 'bc', &
        'explicit interval metadata')
    call require(all(original%lower95 == percentile%lower95) .and. &
        all(original%coherence_upper95 == percentile%coherence_upper95), &
        'explicit percentile precisely reproduces default')
    call require(all(original%stderr == basic%stderr) .and. &
        all(original%stderr == bc%stderr) .and. &
        all(original%coherence_bias == basic%coherence_bias), &
        'interval method does not change bootstrap samples or moments')
    call require(all(basic%lower95 >= 0.0_dp) .and. all(basic%upper95 <= 1.0_dp) &
        .and. all(basic%valid),'basic confidence bounds in the parameter space')
    do f = 1, size(basic%frequencies)
        expected_low = max(0.0_dp,2.0_dp*basic%coherence(f)- &
            percentile%coherence_upper95(f))
        expected_high = min(1.0_dp,2.0_dp*basic%coherence(f)- &
            percentile%coherence_lower95(f))
        call require(abs(basic%coherence_lower95(f)-expected_low) < 1.0e-14_dp &
            .and. abs(basic%coherence_upper95(f)-expected_high) < 1.0e-14_dp, &
            'basic interval is reflected percentile interval')
    end do
    call require(count(bc%coherence_valid) > 0 .and. count(bc%valid) > 0, &
        'BC intervals identifiable in regular cases')
    call require(all(pack(bc%lower95 <= bc%upper95,bc%valid)), &
        'BC intervals have ordered endpoints')
    call require(maxval(abs(basic%coherence_lower95-percentile%coherence_lower95)) &
        > 1.0e-7_dp,'basic bounds generally differ from percentile bounds')
    invalid = canonical_moving_block_bootstrap(spectra,[1],[2],3,interval_method='unknown')
    call require(invalid%status == 1,'reject unsupported interval option')
    invalid = canonical_moving_block_bootstrap(spectra,[1],[2],3,interval_method='normal')
    call require(invalid%status == 1,'reject unsupported normal/studentized option')

    zero_spectra = spectra
    zero_spectra%power = (0.0_dp,0.0_dp)
    do t = 1, 18
        do f = 1, 3
            do i = 1, 3
                zero_spectra%power(f,t,i,i) = (1.0_dp,0.0_dp)
            end do
        end do
    end do
    invalid = canonical_moving_block_bootstrap(zero_spectra,[1],[2],3, &
        n_resamples=30,interval_method='bc')
    call require(invalid%status == 0 .and. .not. any(invalid%valid) &
        .and. all(ieee_is_nan(invalid%lower95)), &
        'boundary canonical roots remain invalid for BC intervals')

    ! All bootstrap replicates have nonzero association but the original
    ! mean cross-spectrum cancels, so the BC bias-offset rank is zero.
    allocate(extreme%frequencies(1),extreme%times(6),extreme%power(1,6,2,2))
    extreme%frequencies = [0.0_dp]
    extreme%times = [(real(t,dp),t=1,6)]
    extreme%power = (0.0_dp,0.0_dp)
    extreme%status = 0
    do t = 1, 6
        extreme%power(1,t,1,1) = (1.0_dp,0.0_dp)
        extreme%power(1,t,2,2) = (1.0_dp,0.0_dp)
        extreme%power(1,t,1,2) = cmplx(correlations(t),0.0_dp,dp)
        extreme%power(1,t,2,1) = cmplx(correlations(t),0.0_dp,dp)
    end do
    invalid = canonical_moving_block_bootstrap(extreme,[1],[2],2, &
        n_resamples=30,seed=3,interval_method='bc')
    call require(invalid%status == 2 .and. .not. any(invalid%coherence_valid), &
        'unidentifiable BC bias correction returns invalid aggregate interval')
    call require(all(ieee_is_nan(invalid%coherence_lower95)) .and. &
        all(ieee_is_nan(invalid%coherence_upper95)), &
        'undefined BC endpoints are NaNs, not invented tail probabilities')
    call require(all(invalid%coherence_stderr > 0.0_dp), &
        'bootstrap moment estimates remain available when BC fails')

    do i = 1, 384
        do j = 1, 3
            x(i,j) = sin(0.09_dp*i*j)+0.31_dp*cos(0.06_dp*i+0.2_dp*j)
            z(i,j) = cmplx(x(i,j),0.16_dp*sin(0.04_dp*i*j),dp)
        end do
    end do
    direct = multitaper_canonical_moving_block_bootstrap(x,3.0_dp,24,[1],[2],3, &
        n_resamples=50,kmax=4,weighting='eigen',ridge=0.08_dp, &
        interval_method='basic')
    call require(direct%status == 0 .and. direct%interval_method == 'basic', &
        'real direct multitaper method forwarding')
    direct = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,24,[1],[2],3, &
        n_resamples=50,kmax=4,weighting='adaptive',ridge=0.08_dp, &
        bands=reshape([0.03_dp,0.18_dp],[2,1]), &
        interval_method='bc')
    call require(direct%interval_method == 'bc' .and. direct%n_windows == 16 &
        .and. allocated(direct%coherence_lower95), &
        'complex adaptive band method forwarding')
    print '(a)', 'PASS canonical bootstrap percentile, basic and BC interval methods'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! condition to verify
        character(*), intent(in) :: message !! assertion description
        if (.not. ok) then
            print '(a)', 'FAIL '//message
            error stop 1
        end if
    end subroutine require
end program test_signal_canonical_interval_methods
