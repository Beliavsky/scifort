! SPDX-License-Identifier: MIT
program test_signal_tf_partial_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result, &
        multitaper_spectral_matrix_spectrogram
    use scifort_signal_partial_coherence, only : partial_coherence_tf_result, &
        partial_coherence_tf_matrix, multitaper_tf_partial_coherence, &
        partial_coherence_success, partial_coherence_singular, partial_coherence_invalid_input
    implicit none

    type(multitaper_spectral_matrix_spectrogram_result) :: matrix
    type(partial_coherence_tf_result) :: pc, loaded, scaled, direct, from_data
    real(dp) :: data(97,3), t
    complex(dp) :: cdata(97,3)
    integer :: i, j, f, k

    allocate(matrix%frequencies(3), matrix%times(2), matrix%power(3,2,3,3))
    matrix%status = 0
    matrix%frequencies = [0.0_dp, 0.3_dp, 0.6_dp]
    matrix%times = [1.0_dp, 2.0_dp]
    matrix%power = (0.0_dp, 0.0_dp)
    do f = 1, 3
        do k = 1, 2
            do i = 1, 3
                matrix%power(f,k,i,i) = (1.0_dp, 0.0_dp)
            end do
        end do
    end do

    ! Markov chain: conditioning on channel 2 removes channels 1-3 coherence.
    matrix%power(1,1,1,2) = (0.4_dp,0.0_dp)
    matrix%power(1,1,2,1) = (0.4_dp,0.0_dp)
    matrix%power(1,1,2,3) = (0.4_dp,0.0_dp)
    matrix%power(1,1,3,2) = (0.4_dp,0.0_dp)
    matrix%power(1,1,1,3) = (0.16_dp,0.0_dp)
    matrix%power(1,1,3,1) = (0.16_dp,0.0_dp)

    ! One rank-deficient frequency/time bin.
    matrix%power(3,1,:,:) = (1.0_dp,0.0_dp)
    ! The second time bin includes complex off-diagonal elements.
    matrix%power(1,2,1,2) = (0.0_dp,0.25_dp)
    matrix%power(1,2,2,1) = (0.0_dp,-0.25_dp)

    pc = partial_coherence_tf_matrix(matrix)
    if (pc%status /= partial_coherence_singular) error stop 'missing singular-bin status'
    if (.not. all(pc%valid(:2,:))) error stop 'valid time-frequency bins lost'
    if (pc%valid(3,1)) error stop 'singular bin reported valid'
    if (.not. pc%valid(3,2)) error stop 'valid bin after singular not processed'
    if (.not. ieee_is_nan(pc%coherence(3,1,1,2))) error stop 'missing NaN'
    if (abs(pc%coherence(1,1,1,3)) > 1.0e-25_dp) error stop 'conditional independence'
    if (maxval(abs(pc%coherence(1,2,:,:)-transpose(pc%coherence(1,2,:,:)))) > 1.0e-13_dp) &
        error stop 'Hermitian partial-coherence symmetry'
    if (any(shape(pc%coherence) /= [3,2,3,3])) error stop 'wrong shape'
    if (maxval(abs(pc%times-matrix%times)) > 0.0_dp) error stop 'time grid mismatch'

    loaded = partial_coherence_tf_matrix(matrix, ridge=0.05_dp)
    if (loaded%status /= partial_coherence_success) error stop 'loading failed'
    if (.not. all(loaded%valid)) error stop 'loaded bin not invertible'
    if (any(loaded%coherence < 0.0_dp) .or. any(loaded%coherence > 1.0_dp)) &
        error stop 'coherence outside unit interval'
    do i = 1, 3
        if (maxval(abs(loaded%coherence(:,:,i,i)-1.0_dp)) > 1.0e-14_dp) &
            error stop 'nonunit diagonal'
    end do
    matrix%power = matrix%power*1.0e6_dp
    scaled = partial_coherence_tf_matrix(matrix,ridge=0.05_dp)
    if (maxval(abs(scaled%coherence-loaded%coherence)) > 1.0e-12_dp) &
        error stop 'not invariant under channel scaling'
    matrix%power = matrix%power/1.0e6_dp

    ! Reject invalid inputs with correctly initialized empty outputs.
    pc = partial_coherence_tf_matrix(matrix,ridge=-0.5_dp)
    if (pc%status /= partial_coherence_invalid_input) error stop 'negative ridge accepted'
    if (size(pc%coherence) /= 0) error stop 'invalid result has coefficients'
    matrix%power(1,1,1,2) = (0.0_dp,2.0_dp)
    pc = partial_coherence_tf_matrix(matrix)
    if (pc%status /= partial_coherence_invalid_input) error stop 'non-Hermitian accepted'
    matrix%power(1,1,1,2) = (0.4_dp,0.0_dp)
    matrix%times(1) = ieee_value(0.0_dp,ieee_quiet_nan)
    pc = partial_coherence_tf_matrix(matrix)
    if (pc%status /= partial_coherence_invalid_input) error stop 'nonfinite timestamps accepted'
    matrix%times(1) = 1.0_dp

    do j = 1, 3
        do i = 1, 97
            t = real(i,dp)
            data(i,j) = sin((0.12_dp+0.05_dp*j)*t) + 0.19_dp*cos(0.27_dp*t*j)
            cdata(i,j) = cmplx(data(i,j),0.5_dp*sin(0.14_dp*j*t),dp)
        end do
    end do

    from_data = multitaper_tf_partial_coherence(data,3.0_dp,32, &
        noverlap=16,kmax=5,nfft=40,fs=100.0_dp,ridge=0.04_dp)
    if (from_data%status /= partial_coherence_success) error stop 'real signals failed'
    matrix = multitaper_spectral_matrix_spectrogram(data,3.0_dp,32, &
        noverlap=16,kmax=5,nfft=40,fs=100.0_dp)
    direct = partial_coherence_tf_matrix(matrix,ridge=0.04_dp)
    if (maxval(abs(from_data%coherence-direct%coherence)) > 1.0e-14_dp) &
        error stop 'real convenience/direct discrepancy'
    if (size(from_data%times) /= 5) error stop 'number of windows incorrect'
    if (size(from_data%frequencies) /= 21) error stop 'real spectral grid incorrect'
    from_data = multitaper_tf_partial_coherence(transpose(data),3.0_dp,32, &
        noverlap=16,kmax=5,nfft=40,fs=100.0_dp,axis=2,ridge=0.04_dp)
    if (maxval(abs(from_data%coherence-direct%coherence)) > 1.0e-14_dp) &
        error stop 'sample-axis mismatch'
    from_data = multitaper_tf_partial_coherence(cdata,3.0_dp,32, &
        noverlap=16,kmax=6,nfft=40,ridge=0.06_dp,weighting='eigen')
    if (from_data%status /= partial_coherence_success) error stop 'complex signals failed'
    matrix = multitaper_spectral_matrix_spectrogram(cdata,3.0_dp,32, &
        noverlap=16,kmax=6,nfft=40,weighting='eigen')
    direct = partial_coherence_tf_matrix(matrix,ridge=0.06_dp)
    if (maxval(abs(from_data%coherence-direct%coherence)) > 1.0e-14_dp) &
        error stop 'complex convenience/direct discrepancy'
    if (size(from_data%frequencies) /= 40) error stop 'complex Fourier grid incorrect'

    from_data = multitaper_tf_partial_coherence(data(1:97:2,1:3:2),2.5_dp,24, &
        noverlap=12,kmax=4,nfft=32,ridge=0.02_dp)
    if (from_data%status /= partial_coherence_success) error stop 'strided input failed'
    if (size(from_data%times) /= 3) error stop 'strided time grid incorrect'
    from_data = multitaper_tf_partial_coherence(data,3.0_dp,200,ridge=0.1_dp)
    if (from_data%status /= partial_coherence_invalid_input) error stop 'invalid window accepted'
    from_data = multitaper_tf_partial_coherence(data,3.0_dp,32,axis=3,ridge=0.1_dp)
    if (from_data%status /= partial_coherence_invalid_input) error stop 'invalid axis accepted'

    print '(a)', 'PASS multitaper time-frequency partial coherence native tests'
end program test_signal_tf_partial_coherence
