! SPDX-License-Identifier: MIT
program test_signal_canonical_interval_methods_public
    use scifort_signal, only : dp, canonical_moving_block_bootstrap_result, &
        multitaper_canonical_moving_block_bootstrap
    implicit none
    real(dp) :: x(384,3)
    complex(dp) :: z(384,3)
    type(canonical_moving_block_bootstrap_result) :: a, b
    integer :: i, j

    do i = 1, 384
        do j = 1, 3
            x(i,j) = sin(0.068_dp*i*j)+0.21_dp*cos(0.025_dp*i*j)
            z(i,j) = cmplx(x(i,j),0.11_dp*sin(0.05_dp*i*j),dp)
        end do
    end do
    a = multitaper_canonical_moving_block_bootstrap(x,3.0_dp,24,[1],[2], &
        block_length=3,n_resamples=39,kmax=4,weighting='eigen', &
        ridge=0.06_dp,interval_method='basic')
    if (a%status /= 0 .or. a%interval_method /= 'basic' .or. &
        .not. all(a%coherence_valid)) error stop 'public basic bootstrap interface'
    b = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,24,[1],[2], &
        block_length=3,n_resamples=39,kmax=4,weighting='adaptive', &
        ridge=0.06_dp,interval_method='bc')
    if (b%interval_method /= 'bc' .or. .not. allocated(b%frequencies)) &
        error stop 'public BC bootstrap interface'
    print '(a)', 'PASS public canonical bootstrap interval methods'
end program test_signal_canonical_interval_methods_public
