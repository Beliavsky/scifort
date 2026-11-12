! SPDX-License-Identifier: MIT
program test_signal_canonical_studentized
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap, multitaper_canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: s, constant_s
    type(canonical_moving_block_bootstrap_result) :: t, again, default_result, short, bad, real_result, complex_result
    real(dp) :: x(384,3)
    complex(dp) :: z(384,3), matrix(3,3)
    integer :: f, j, k, w

    allocate(s%frequencies(3),s%times(15),s%power(3,15,3,3))
    s%frequencies = [0.1_dp,0.3_dp,0.6_dp]
    s%times = [(real(w,dp),w=1,15)]
    s%status = 0
    do f = 1, 3
        do w = 1, 15
            matrix = (0.0_dp,0.0_dp)
            do j = 1, 3
                matrix(j,j) = cmplx(1.0_dp+0.03_dp*j+0.09_dp*w,0.0_dp,dp)
                do k = 1,j-1
                    matrix(j,k) = cmplx(0.4_dp+0.06_dp*cos(0.21_dp*w*j), &
                        0.19_dp*sin(0.27_dp*w*f+k),dp)
                end do
            end do
            s%power(f,w,:,:) = matmul(matrix,transpose(conjg(matrix)))
        end do
    end do
    t = canonical_moving_block_bootstrap(s,[1],[2],3,n_resamples=59, &
        conditioned_on=[3],ridge=0.02_dp,shrinkage=0.08_dp,seed=778, &
        interval_method='studentized')
    call check(t%status == 0,'studentized basic fit succeeds')
    call check(t%interval_method == 'studentized','method metadata')
    call check(allocated(t%studentizing_stderr) .and. &
        allocated(t%coherence_studentizing_stderr),'se allocated only on demand')
    call check(all(t%coherence_valid) .and. all(t%valid),'regular studentized intervals exist')
    call check(all(t%studentizing_stderr > 0.0_dp) .and. &
        all(t%coherence_studentizing_stderr > 0.0_dp),'positive block standard errors')
    call check(all(t%lower95 >= 0.0_dp) .and. all(t%upper95 <= 1.0_dp) .and. &
        all(t%coherence_lower95 <= t%coherence_upper95),'bounded ordered intervals')
    again = canonical_moving_block_bootstrap(s,[1],[2],3,n_resamples=59, &
        conditioned_on=[3],ridge=0.02_dp,shrinkage=0.08_dp,seed=778, &
        interval_method='studentized')
    call check(all(again%lower95 == t%lower95) .and. &
        all(again%coherence_upper95 == t%coherence_upper95), &
        'identical seeds reproduce studentized bounds exactly')
    default_result = canonical_moving_block_bootstrap(s,[1],[2],3,n_resamples=59, &
        conditioned_on=[3],ridge=0.02_dp,shrinkage=0.08_dp,seed=778)
    call check(.not. allocated(default_result%studentizing_stderr) .and. &
        .not. allocated(default_result%coherence_studentizing_stderr), &
        'default percentile avoids studentizing allocations')
    call check(all(default_result%stderr == t%stderr) .and. &
        all(default_result%bias == t%bias) .and. &
        all(default_result%coherence_bias == t%coherence_bias), &
        'studentization does not change bootstrap replicates or moments')
    call check(maxval(abs(default_result%coherence_lower95-t%coherence_lower95)) > 1.0e-8_dp, &
        'studentized endpoints differ from percentile endpoints')
    ! Five blocks, the last containing only three of fifteen windows: a
    ! length-four test exercises the weighted unequal-block pseudovalues.
    short = canonical_moving_block_bootstrap(s,[1],[2],4,n_resamples=35, &
        seed=17,interval_method='studentized',circular=.true.,ridge=0.03_dp)
    call check(short%status == 0 .and. all(short%coherence_valid), &
        'circular resampling with a shorter final block')
    short = canonical_moving_block_bootstrap(s,[1],[2],5,n_resamples=30, &
        seed=17,interval_method='studentized')
    call check(short%status == 2 .and. .not. any(short%coherence_valid), &
        'fewer than four blocks do not create untrustworthy t intervals')
    call check(all(ieee_is_nan(short%coherence_studentizing_stderr)) .and. &
        all(ieee_is_finite(short%coherence_stderr)), &
        'invalid studentization retains bootstrap standard deviations')
    constant_s = s
    do w = 2,15
        constant_s%power(:,w,:,:) = constant_s%power(:,1,:,:)
    end do
    short = canonical_moving_block_bootstrap(constant_s,[1],[2],3, &
        n_resamples=30,interval_method='studentized')
    call check(.not. any(short%coherence_valid) .and. &
        all(ieee_is_nan(short%coherence_lower95)), &
        'constant windows cannot provide a studentizing denominator')
    bad = canonical_moving_block_bootstrap(s,[1],[2],3,interval_method='bootstrap-t')
    call check(bad%status == 1,'reject unknown interval method')
    do j = 1, 384
        do k = 1, 3
            x(j,k) = sin(0.058_dp*j*k)+0.24_dp*cos(0.073_dp*j+0.8_dp*k)
            z(j,k) = cmplx(x(j,k),0.17_dp*sin(0.06_dp*j*k),dp)
        end do
    end do
    real_result = multitaper_canonical_moving_block_bootstrap(x,3.0_dp,24,[1],[2],3, &
        n_resamples=30,kmax=4,weighting='eigen',seed=19,ridge=0.04_dp, &
        interval_method='studentized')
    call check(real_result%status /= 1 .and. real_result%interval_method == 'studentized' .and. &
        allocated(real_result%studentizing_stderr), 'real multitaper forwarding')
    complex_result = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,24,[1],[2],3, &
        n_resamples=30,kmax=4,weighting='adaptive',seed=29,ridge=0.08_dp, &
        bands=reshape([0.03_dp,0.18_dp],[2,1]),interval_method='studentized')
    call check(complex_result%status /= 1 .and. complex_result%interval_method == 'studentized' &
        .and. size(complex_result%frequencies) == 1, 'complex adaptive band forwarding')
    print '(a)', 'PASS canonical studentized moving-block bootstrap and edge cases'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! asserted Boolean
        character(*), intent(in) :: label !! diagnostic label
        if (.not. ok) then
            print '(a)', 'FAIL '//label
            error stop 1
        end if
    end subroutine check
end program test_signal_canonical_studentized
