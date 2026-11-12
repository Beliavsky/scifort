! SPDX-License-Identifier: MIT
! Regression tests for the Gaussian-Wishart automatic spectral shrinkage.
program test_signal_partial_shrinkage_auto
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram
    use scifort_signal_partial_coherence, only : partial_coherence_matrix, &
        partial_coherence_tf_matrix, multitaper_partial_coherence, &
        multitaper_tf_partial_coherence, partial_coherence_result, &
        partial_coherence_tf_result, partial_coherence_success, &
        partial_coherence_invalid_input
    implicit none
    type(multitaper_matrix_result) :: s, raw
    type(multitaper_spectral_matrix_spectrogram_result) :: tf_s, raw_tf
    type(partial_coherence_result) :: a, b, c
    type(partial_coherence_tf_result) :: tf, tf_signal
    real(dp) :: x(64,3), t
    complex(dp) :: z(64,3)
    integer :: i, j, f

    allocate(s%frequencies(3), s%power(3,3,3))
    s%status = 0
    s%frequencies = [0.0_dp, 1.0_dp, 2.0_dp]
    s%power = cmplx(0.0_dp,0.0_dp,dp)
    do f = 1, 3
        do i = 1, 3
            s%power(f,i,i) = 1.0_dp + real(i,dp)/3.0_dp
        end do
        s%power(f,1,2) = cmplx(0.3_dp,0.15_dp,dp)
        s%power(f,2,1) = conjg(s%power(f,1,2))
        s%power(f,2,3) = cmplx(-0.1_dp,0.08_dp,dp)
        s%power(f,3,2) = conjg(s%power(f,2,3))
    end do

    a = partial_coherence_matrix(s, effective_tapers=4.0_dp, &
        auto_shrinkage=.true., return_coherency=.true.)
    if (a%status /= partial_coherence_success) error stop 'automatic matrix invalid'
    if (.not. all(a%valid)) error stop 'automatic matrix marked singular'
    if (any(a%shrinkage_used < 0.0_dp) .or. any(a%shrinkage_used > 1.0_dp)) &
        error stop 'shrinkage outside [0,1]'
    if (maxval(abs(a%shrinkage_used - a%shrinkage_used(1))) > 0.0_dp) &
        error stop 'identical bins should share alpha'
    if (maxval(abs(abs(a%coherency)**2 - a%coherence)) > 1.0e-14_dp) &
        error stop 'coherence vs coherency mismatch'
    b = partial_coherence_matrix(s, shrinkage=a%shrinkage_used(1), return_coherency=.true.)
    if (maxval(abs(a%precision-b%precision)) > 1.0e-13_dp) &
        error stop 'automatic coefficient differs from manual shrinkage'
    if (any(abs(b%shrinkage_used-a%shrinkage_used)>1.0e-15_dp)) &
        error stop 'manual shrinkage report wrong'
    c = partial_coherence_matrix(s)
    if (maxval(abs(c%shrinkage_used)) /= 0.0_dp) error stop 'default shrinkage changed'
    if (size(c%coherency) /= 0) error stop 'coherency allocated by default'

    a = partial_coherence_matrix(s, auto_shrinkage=.true., effective_tapers=1.0_dp)
    if (a%status /= partial_coherence_success) error stop 'one taper should diagonally shrink'
    if (maxval(abs(a%shrinkage_used-1.0_dp)) /= 0.0_dp) error stop 'one taper alpha not 1'
    if (maxval(abs(a%coherence(:,1,2))) > 1.0e-15_dp) error stop 'one taper offdiag nonzero'

    s%power(2,:,:) = cmplx(0.0_dp,0.0_dp,dp)
    a = partial_coherence_matrix(s, auto_shrinkage=.true., effective_tapers=4.0_dp)
    if (a%valid(2)) error stop 'zero spectral cell incorrectly valid'
    if (.not. ieee_is_nan(a%shrinkage_used(2))) error stop 'singular alpha not NaN'
    s%power(2,:,:) = s%power(1,:,:)
    a = partial_coherence_matrix(s, auto_shrinkage=.true.)
    if (a%status /= partial_coherence_invalid_input) error stop 'missing effective count accepted'
    a = partial_coherence_matrix(s, auto_shrinkage=.true., effective_tapers=0.0_dp)
    if (a%status /= partial_coherence_invalid_input) error stop 'zero effective count accepted'
    a = partial_coherence_matrix(s, auto_shrinkage=.true., effective_tapers=4.0_dp, shrinkage=0.2_dp)
    if (a%status /= partial_coherence_invalid_input) error stop 'automatic/manual conflict accepted'
    a = partial_coherence_matrix(s, effective_tapers=4.0_dp)
    if (a%status /= partial_coherence_invalid_input) error stop 'unused effective count accepted'

    allocate(tf_s%frequencies(3), tf_s%times(2), tf_s%power(3,2,3,3))
    tf_s%frequencies = s%frequencies
    tf_s%times = [0.25_dp, 0.75_dp]
    tf_s%status = 0
    tf_s%power(:,1,:,:) = s%power
    tf_s%power(:,2,:,:) = 1.5_dp*s%power
    tf = partial_coherence_tf_matrix(tf_s, auto_shrinkage=.true., effective_tapers=4.0_dp, &
        return_coherency=.true.)
    if (tf%status /= partial_coherence_success) error stop 'TF automatic invalid'
    if (any(shape(tf%shrinkage_used) /= [3,2])) error stop 'TF alpha shape'
    if (maxval(abs(tf%shrinkage_used(:,1)-tf%shrinkage_used(:,2))) > 1.0e-14_dp) &
        error stop 'TF scalar rescaling changes alpha'
    if (maxval(abs(tf%coherence(:,1,:,:)-tf%coherence(:,2,:,:))) > 1.0e-13_dp) &
        error stop 'TF scalar rescaling changes coherence'
    tf = partial_coherence_tf_matrix(tf_s, auto_shrinkage=.true., effective_tapers=2.0_dp, &
        shrinkage=0.2_dp)
    if (tf%status /= partial_coherence_invalid_input) error stop 'TF conflicting shrinkage'

    do j = 1, 3
        do i = 1, 64
            t = real(i,dp)
            x(i,j) = sin(0.07_dp*t+0.3_dp*j)+0.4_dp*cos(0.26_dp*t-0.5_dp*j)
        end do
    end do
    z = cmplx(x,0.25_dp*x,dp)
    raw = multitaper_spectral_matrix(x, 2.5_dp, kmax=4, weighting='equal')
    a = multitaper_partial_coherence(x, 2.5_dp, kmax=4, weighting='equal', &
        auto_shrinkage=.true., ridge=0.03_dp)
    b = partial_coherence_matrix(raw, auto_shrinkage=.true., effective_tapers=4.0_dp, ridge=0.03_dp)
    if (a%status /= partial_coherence_success) error stop 'real high-level automatic failed'
    if (maxval(abs(a%shrinkage_used-b%shrinkage_used)) > 1.0e-14_dp) &
        error stop 'real high-level shrinkage does not match matrix'
    if (maxval(abs(a%coherence-b%coherence)) > 1.0e-13_dp) &
        error stop 'real high-level coherence does not match matrix'
    raw = multitaper_spectral_matrix(transpose(z), 2.5_dp, kmax=4, axis=2, weighting='equal')
    a = multitaper_partial_coherence(transpose(z), 2.5_dp, kmax=4, axis=2, &
        auto_shrinkage=.true., ridge=0.03_dp, return_coherency=.true.)
    b = partial_coherence_matrix(raw, auto_shrinkage=.true., effective_tapers=4.0_dp, &
        ridge=0.03_dp, return_coherency=.true.)
    if (a%status /= partial_coherence_success) error stop 'complex high-level automatic failed'
    if (maxval(abs(a%coherency-b%coherency)) > 1.0e-13_dp) &
        error stop 'complex high-level coherency does not match matrix'
    a = multitaper_partial_coherence(x, 2.5_dp, kmax=4, weighting='eigen', auto_shrinkage=.true.)
    if (a%status /= partial_coherence_success) error stop 'weighted auto shrinkage failed'
    a = multitaper_partial_coherence(z, 2.5_dp, kmax=4, shrinkage=0.2_dp, auto_shrinkage=.true.)
    if (a%status /= partial_coherence_invalid_input) error stop 'conflicting high-level args not rejected'

    raw_tf = multitaper_spectral_matrix_spectrogram(z, 2.5_dp, 32, noverlap=16, &
        kmax=4, weighting='equal')
    tf_signal = multitaper_tf_partial_coherence(z, 2.5_dp, 32, noverlap=16, kmax=4, &
        auto_shrinkage=.true., ridge=0.03_dp, return_coherency=.true.)
    tf = partial_coherence_tf_matrix(raw_tf, auto_shrinkage=.true., effective_tapers=4.0_dp, &
        ridge=0.03_dp, return_coherency=.true.)
    if (tf_signal%status /= partial_coherence_success) error stop 'TF high-level automatic failed'
    if (maxval(abs(tf_signal%shrinkage_used-tf%shrinkage_used)) > 1.0e-14_dp) &
        error stop 'TF high-level shrinkage does not match matrix'
    if (maxval(abs(tf_signal%coherency-tf%coherency)) > 1.0e-13_dp) &
        error stop 'TF high-level coherency does not match matrix'
    tf_signal = multitaper_tf_partial_coherence(x, 2.5_dp, 32, noverlap=16, &
        kmax=4, auto_shrinkage=.true., ridge=0.03_dp)
    if (tf_signal%status /= partial_coherence_success) error stop 'real TF automatic failed'
    tf_signal = multitaper_tf_partial_coherence(x, 2.5_dp, 32, noverlap=16, &
        kmax=4, weighting='adaptive', auto_shrinkage=.true.)
    if (tf_signal%status /= partial_coherence_invalid_input) error stop 'adaptive TF auto accepted'

    print *, 'PASS automatic diagonal spectral shrinkage native tests'
end program test_signal_partial_shrinkage_auto
