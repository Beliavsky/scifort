! SPDX-License-Identifier: MIT
program test_signal_canonical_moving_bootstrap
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap, multitaper_canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra, bad
    type(canonical_moving_block_bootstrap_result) :: a, b, c, d, band, direct
    real(dp) :: x(256,5), bands(2,2), r1, r2
    complex(dp) :: factor(5,5), z(256,5)
    integer, allocatable :: random_seed_state(:)
    integer :: f, t, i, j, nseed

    allocate(spectra%power(4,15,5,5),spectra%frequencies(4),spectra%times(15))
    spectra%frequencies = [0.0_dp,0.25_dp,0.5_dp,0.75_dp]
    spectra%times = [(real(t,dp),t=1,15)]
    spectra%status = 0
    do f = 1, 4
        do t = 1, 15
            factor = (0.0_dp,0.0_dp)
            do i = 1, 5
                factor(i,i) = cmplx(1.0_dp+0.16_dp*i+0.035_dp*t,0.0_dp,dp)
                do j = 1, i-1
                    factor(i,j) = cmplx(0.24_dp*cos(0.17_dp*i*j+0.11_dp*t), &
                        0.19_dp*sin(0.21_dp*j*f+0.09_dp*t),dp)
                end do
            end do
            spectra%power(f,t,:,:) = matmul(factor,transpose(conjg(factor)))*(0.7_dp+0.08_dp*f)
        end do
    end do

    a = canonical_moving_block_bootstrap(spectra,[1,2],[3,4],4,n_resamples=83, &
        conditioned_on=[5],ridge=0.04_dp,shrinkage=0.08_dp,seed=1729)
    call require(a%status == 0,'ordinary MBB success')
    call require(a%n_windows == 15 .and. a%n_resamples == 83 .and. &
        a%block_length == 4 .and. a%blocks_per_resample == 4,'non-divisible block metadata')
    call require(a%seed == 1729 .and. .not. a%circular,'seed and ordinary flag')
    call require(all(a%coherence_valid),'aggregate valid across ordinary MBB replicates')
    call require(all(a%valid),'isolated interior modes')
    call require(all(a%stderr >= 0.0_dp) .and. all(a%coherence_stderr >= 0.0_dp), &
        'nonnegative standard errors')
    call require(all(a%lower95 >= 0.0_dp) .and. all(a%upper95 <= 1.0_dp), &
        'bounded bootstrap percentile intervals')
    call require(all(a%lower95 <= a%upper95),'ordered percentile intervals')
    b = canonical_moving_block_bootstrap(spectra,[1,2],[3,4],4,n_resamples=83, &
        conditioned_on=[5],ridge=0.04_dp,shrinkage=0.08_dp,seed=1729)
    call require(all(a%stderr == b%stderr) .and. all(a%lower95 == b%lower95) .and. &
        all(a%coherence_bias == b%coherence_bias),'same seed reproduces every statistic')
    c = canonical_moving_block_bootstrap(spectra,[1,2],[3,4],4,n_resamples=83, &
        conditioned_on=[5],ridge=0.04_dp,shrinkage=0.08_dp,seed=1730)
    call require(maxval(abs(a%stderr-c%stderr)) > 1.0e-10_dp,'new seed changes draws')
    d = canonical_moving_block_bootstrap(spectra,[1,2],[3,4],4,n_resamples=83, &
        conditioned_on=[5],ridge=0.04_dp,shrinkage=0.08_dp,seed=1729,circular=.true.)
    call require(d%status == 0 .and. d%circular,'circular moving blocks')
    call require(maxval(abs(a%stderr-d%stderr)) > 1.0e-10_dp,'circular sampling changes draws')
    call require(all(a%squared_modes == d%squared_modes),'full sample unaffected by block geometry')

    bands(:,1) = [0.04_dp,0.37_dp]
    bands(:,2) = [0.42_dp,0.71_dp]
    band = canonical_moving_block_bootstrap(spectra,[1,2],[3,4],3,n_resamples=60, &
        conditioned_on=[5],ridge=0.04_dp,shrinkage=0.08_dp,seed=1201,bands=bands)
    call require(band%status == 0 .and. size(band%frequencies) == 2,'band integration MBB')
    d = canonical_moving_block_bootstrap(spectra,[1,2],[3,4],3,n_resamples=60, &
        conditioned_on=[5],ridge=0.04_dp,shrinkage=0.08_dp,seed=1201, &
        bands=bands,average=.true.)
    call require(maxval(abs(band%squared_modes-d%squared_modes)) < 1.0e-12_dp, &
        'band averaging and integration give equivalent coherence')
    call require(maxval(abs(band%stderr-d%stderr)) < 1.0e-12_dp, &
        'band averaging and integration give equivalent bootstrap variability')

    call random_seed(size=nseed)
    allocate(random_seed_state(nseed))
    random_seed_state = 13579
    call random_seed(put=random_seed_state)
    call random_number(r1)
    d = canonical_moving_block_bootstrap(spectra,[1,2],[3,4],1,n_resamples=20,seed=1, &
        ridge=0.04_dp)
    call random_seed(put=random_seed_state)
    call random_number(r2)
    call require(r1 == r2,'bootstrap does not modify intrinsic RNG state')

    direct = canonical_moving_block_bootstrap(spectra,[1],[2],0,n_resamples=20)
    call require(direct%status == 1,'reject zero block length')
    direct = canonical_moving_block_bootstrap(spectra,[1],[2],15,n_resamples=20)
    call require(direct%status == 1,'reject full-window nonidentifiable block')
    direct = canonical_moving_block_bootstrap(spectra,[1],[2],16,n_resamples=20)
    call require(direct%status == 1,'reject too-long block')
    direct = canonical_moving_block_bootstrap(spectra,[1],[2],2,n_resamples=19)
    call require(direct%status == 1,'reject too few bootstrap replicates')
    direct = canonical_moving_block_bootstrap(spectra,[1],[2],2,n_resamples=20,seed=0)
    call require(direct%status == 1,'reject invalid RNG seed')
    direct = canonical_moving_block_bootstrap(spectra,[1],[1],2,n_resamples=20)
    call require(direct%status == 1,'reject overlapping channel groups')
    bad = spectra
    bad%times(4) = bad%times(3)
    direct = canonical_moving_block_bootstrap(bad,[1],[2],2,n_resamples=20)
    call require(direct%status == 1,'reject repeated window times')

    bad = spectra
    bad%power = (0.0_dp,0.0_dp)
    do t = 1, 15
        do f = 1, 4
            do i = 1, 5
                bad%power(f,t,i,i) = (1.0_dp,0.0_dp)
            end do
        end do
    end do
    direct = canonical_moving_block_bootstrap(bad,[1],[2],3,n_resamples=40)
    call require(direct%status == 0 .and. .not. any(direct%valid), &
        'zero modes produce no rankwise inference')
    call require(all(ieee_is_nan(direct%stderr)) .and. all(direct%coherence_valid), &
        'zero modes retain valid aggregate but not rankwise interval')
    bad%power = (0.0_dp,0.0_dp)
    direct = canonical_moving_block_bootstrap(bad,[1],[2],3,n_resamples=40)
    call require(direct%status == 2 .and. .not. any(direct%coherence_valid), &
        'singular spectra invalidate bootstrap inference')

    do i = 1, 256
        do j = 1, 5
            x(i,j) = sin(0.095_dp*i*j)+0.23_dp*cos(0.16_dp*i+0.09_dp*j) &
                +0.04_dp*sin(0.17_dp*i*i+0.17_dp*j)
            z(i,j) = cmplx(x(i,j),0.16_dp*cos(0.08_dp*i*j+0.12_dp*j),dp)
        end do
    end do
    direct = multitaper_canonical_moving_block_bootstrap(x,3.0_dp,16,[1,2],[3,4],3, &
        n_resamples=30,conditioned_on=[5],kmax=5,weighting='eigen', &
        ridge=0.10_dp,shrinkage=0.14_dp,seed=2345)
    call require(direct%status == 0 .and. direct%n_windows == 16, &
        'direct real DPSS moving-block bootstrap')
    direct = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,16,[1,2],[3,4],3, &
        n_resamples=30,conditioned_on=[5],kmax=5,weighting='adaptive', &
        ridge=0.10_dp,shrinkage=0.14_dp,seed=2345,circular=.true., &
        bands=reshape([0.04_dp,0.18_dp,0.24_dp,0.42_dp],[2,2]))
    call require(direct%status == 0 .and. direct%n_windows == 16 .and. &
        size(direct%frequencies) == 2,'direct complex adaptive circular band bootstrap')
    call require(all(direct%coherence_valid),'complex adaptive circular aggregate fit')
    print '(a)', 'PASS canonical moving-block bootstrap, boundaries, seed, bands and real/complex signals'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! condition to check
        character(*), intent(in) :: message !! assertion description
        if (.not. ok) then
            print '(a)', 'FAIL '//message
            error stop 1
        end if
    end subroutine require
end program test_signal_canonical_moving_bootstrap
