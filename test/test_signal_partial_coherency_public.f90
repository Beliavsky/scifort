! SPDX-License-Identifier: MIT
! The new optional arguments must be visible through the top-level module.
program test_signal_partial_coherency_public
    use scifort_signal, only : dp, partial_coherence_matrix, &
        multitaper_partial_coherence, partial_coherence_tf_matrix, &
        multitaper_tf_partial_coherence, partial_coherence_result, &
        partial_coherence_tf_result, multitaper_matrix_result, &
        multitaper_spectral_matrix, multitaper_spectral_matrix_spectrogram_result, &
        partial_coherence_success
    implicit none
    real(dp) :: x(72,3)
    complex(dp) :: z(72,3)
    type(partial_coherence_result) :: p, q
    type(partial_coherence_tf_result) :: tf
    type(multitaper_matrix_result) :: s
    type(multitaper_spectral_matrix_spectrogram_result) :: ts
    integer :: i, j

    do j = 1, 3
        do i = 1, 72
            x(i,j) = sin(0.19_dp*real(i*j,dp))+0.05_dp*cos(0.13_dp*real(i+j,dp))
            z(i,j) = cmplx(x(i,j),0.10_dp*cos(0.07_dp*real(i*j,dp)),dp)
        end do
    end do
    p = multitaper_partial_coherence(x,3.0_dp,kmax=6, &
        shrinkage=0.15_dp,return_coherency=.true.)
    if (p%status /= partial_coherence_success) error stop 'real public partial coherency'
    s = multitaper_spectral_matrix(x,3.0_dp,kmax=6)
    q = partial_coherence_matrix(s,shrinkage=0.15_dp,return_coherency=.true.)
    if (maxval(abs(p%coherency-q%coherency)) > 1.0e-12_dp) &
        error stop 'public static adapter discrepancy'
    tf = multitaper_tf_partial_coherence(z,2.5_dp,24,noverlap=12,kmax=5, &
        shrinkage=0.08_dp,return_coherency=.true.)
    if (tf%status /= partial_coherence_success) error stop 'complex public TF partial coherency'
    allocate(ts%frequencies(1),ts%times(1),ts%power(1,1,3,3))
    ts%status = 0
    ts%frequencies = [0.0_dp]
    ts%times = [0.0_dp]
    ts%power = (0.0_dp,0.0_dp)
    do i = 1,3
        ts%power(1,1,i,i) = (1.0_dp,0.0_dp)
    end do
    tf = partial_coherence_tf_matrix(ts,shrinkage=0.5_dp,return_coherency=.true.)
    if (tf%status /= partial_coherence_success) error stop 'public TF matrix interface'
    if (maxval(abs(tf%coherency(1,1,:,:)-cmplx(reshape([ &
        1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp], &
        [3,3]),0.0_dp,dp))) > 1.0e-12_dp) error stop 'public diagonal coherence'
    print '(a)', 'PASS public complex partial coherency and shrinkage'
end program test_signal_partial_coherency_public
