! SPDX-License-Identifier: MIT
program test_signal_canonical_bca_public
    use scifort_signal, only : dp, canonical_moving_block_bootstrap_result, &
        multitaper_canonical_moving_block_bootstrap
    implicit none
    real(dp) :: x(384,3)
    complex(dp) :: z(384,3)
    type(canonical_moving_block_bootstrap_result) :: a, b
    integer :: i, j

    do i = 1,384
        do j = 1,3
            x(i,j) = sin(0.057_dp*i*j)+0.27_dp*cos(0.019_dp*i+0.33_dp*j)
            z(i,j) = cmplx(x(i,j),0.19_dp*sin(0.041_dp*i*j),dp)
        end do
    end do
    a = multitaper_canonical_moving_block_bootstrap(x,3.0_dp,24,[1],[2], &
        block_length=3,n_resamples=30,kmax=4,weighting='eigen', &
        ridge=0.05_dp,interval_method='bca')
    if (a%status == 1 .or. .not. allocated(a%bca_acceleration) .or. &
        a%interval_method /= 'bca') error stop 'real public BCa interface'
    b = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,24,[1],[2], &
        block_length=3,n_resamples=30,kmax=4,weighting='adaptive', &
        ridge=0.06_dp,interval_method='bca', &
        bands=reshape([0.03_dp,0.15_dp],[2,1]))
    if (b%status == 1 .or. .not. allocated(b%coherence_bca_acceleration) .or. &
        size(b%frequencies) /= 1) error stop 'complex adaptive public BCa band interface'
    print '(a)', 'PASS public canonical BCa bootstrap interfaces'
end program test_signal_canonical_bca_public
