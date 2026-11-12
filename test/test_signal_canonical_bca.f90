! SPDX-License-Identifier: MIT
program test_signal_canonical_bca
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_bootstrap, only : canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap, multitaper_canonical_moving_block_bootstrap
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra, constant_spectra
    type(canonical_moving_block_bootstrap_result) :: bca, previous, again, short, bad, direct
    real(dp) :: x(384,3)
    complex(dp) :: z(384,3), factor(3,3)
    integer :: f, w, i, j

    allocate(spectra%frequencies(3),spectra%times(18),spectra%power(3,18,3,3))
    spectra%frequencies = [0.0_dp,0.25_dp,0.50_dp]
    spectra%times = [(real(w,dp),w=1,18)]
    spectra%status = 0
    do f = 1, 3
        do w = 1, 18
            factor = (0.0_dp,0.0_dp)
            do i = 1, 3
                factor(i,i) = cmplx(1.0_dp+0.04_dp*i+0.08_dp*w,0.0_dp,dp)
                do j = 1, i-1
                    factor(i,j) = cmplx(0.36_dp+0.08_dp*cos(0.16_dp*w*i), &
                        0.16_dp*sin(0.19_dp*f*w+j),dp)
                end do
            end do
            spectra%power(f,w,:,:) = matmul(factor,transpose(conjg(factor)))
        end do
    end do
    bca = canonical_moving_block_bootstrap(spectra,[1],[2],3,n_resamples=99, &
        conditioned_on=[3],ridge=0.03_dp,shrinkage=0.08_dp,seed=141,interval_method='bca')
    call check(bca%status == 0 .and. bca%interval_method == 'bca','BCa option accepted')
    call check(allocated(bca%bca_acceleration) .and. &
        allocated(bca%coherence_bca_acceleration),'BCa acceleration allocated')
    call check(all(ieee_is_finite(bca%bca_acceleration)) .and. &
        all(ieee_is_finite(bca%coherence_bca_acceleration)),'finite nondegenerate acceleration')
    call check(all(bca%coherence_valid) .and. all(bca%valid),'BCa valid in regular case')
    call check(all(bca%lower95 >= 0.0_dp) .and. all(bca%upper95 <= 1.0_dp), &
        'BCa intervals respect parameter space')
    again = canonical_moving_block_bootstrap(spectra,[1],[2],3,n_resamples=99, &
        conditioned_on=[3],ridge=0.03_dp,shrinkage=0.08_dp,seed=141,interval_method='bca')
    call check(all(bca%lower95 == again%lower95) .and. &
        all(bca%coherence_bca_acceleration == again%coherence_bca_acceleration), &
        'seeded BCa intervals reproducible')
    previous = canonical_moving_block_bootstrap(spectra,[1],[2],3,n_resamples=99, &
        conditioned_on=[3],ridge=0.03_dp,shrinkage=0.08_dp,seed=141)
    call check(all(bca%stderr == previous%stderr) .and. &
        all(bca%coherence_bias == previous%coherence_bias), &
        'BCa does not change bootstrap draws or summary moments')
    call check(.not. allocated(previous%bca_acceleration) .and. &
        .not. allocated(previous%coherence_bca_acceleration), &
        'default avoids BCa allocations')
    short = canonical_moving_block_bootstrap(spectra,[1],[2],4,n_resamples=50, &
        seed=191,ridge=0.04_dp,interval_method='bca',circular=.true.)
    call check(short%status == 0 .and. all(short%coherence_valid), &
        'weighted incomplete final deletion block handled')
    short = canonical_moving_block_bootstrap(spectra,[1],[2],6,n_resamples=35, &
        seed=191,interval_method='bca')
    call check(short%status == 2 .and. .not. any(short%coherence_valid) .and. &
        all(ieee_is_nan(short%coherence_bca_acceleration)), &
        'fewer than four deletion blocks invalidates BCa acceleration')
    call check(all(ieee_is_finite(short%coherence_stderr)) .and. &
        all(ieee_is_nan(short%coherence_lower95)), &
        'invalid BCa retains valid bootstrap moments')
    constant_spectra = spectra
    do w = 2,18
        constant_spectra%power(:,w,:,:) = constant_spectra%power(:,1,:,:)
    end do
    short = canonical_moving_block_bootstrap(constant_spectra,[1],[2],3, &
        n_resamples=40,interval_method='bca')
    call check(.not. any(short%coherence_valid) .and. &
        all(ieee_is_nan(short%coherence_bca_acceleration)), &
        'constant windows cannot identify BCa acceleration')
    bad = canonical_moving_block_bootstrap(spectra,[1],[2],3,interval_method='bcaa')
    call check(bad%status == 1,'invalid interval method rejected')
    do i = 1,384
        do j = 1,3
            x(i,j) = sin(0.045_dp*i*j)+0.33_dp*cos(0.06_dp*i+0.21_dp*j)
            z(i,j) = cmplx(x(i,j),0.13_dp*sin(0.04_dp*i*j),dp)
        end do
    end do
    direct = multitaper_canonical_moving_block_bootstrap(x,3.0_dp,24,[1],[2],3, &
        n_resamples=30,kmax=4,weighting='eigen',ridge=0.08_dp,interval_method='bca')
    call check(direct%status /= 1 .and. allocated(direct%bca_acceleration), &
        'real multitaper BCa forwarding')
    direct = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,24,[1],[2],3, &
        n_resamples=30,kmax=4,weighting='adaptive',ridge=0.08_dp, &
        bands=reshape([0.03_dp,0.18_dp],[2,1]),interval_method='bca')
    call check(direct%status /= 1 .and. size(direct%frequencies) == 1 .and. &
        allocated(direct%coherence_bca_acceleration), &
        'complex adaptive band BCa forwarding')
    print '(a)', 'PASS canonical block-aware BCa bootstrap and edge cases'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! asserted condition
        character(*), intent(in) :: label !! diagnostic label
        if (.not. ok) then
            print '(a)', 'FAIL '//label
            error stop 1
        end if
    end subroutine check
end program test_signal_canonical_bca
