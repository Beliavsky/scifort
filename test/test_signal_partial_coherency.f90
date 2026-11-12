! SPDX-License-Identifier: MIT
! Complex conditional spectra and diagonal-target shrinkage regression tests.
program test_signal_partial_coherency
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_partial_coherence, only : partial_coherence_result, &
        partial_coherence_tf_result, partial_coherence_matrix, &
        partial_coherence_tf_matrix, multitaper_partial_coherence, &
        multitaper_tf_partial_coherence, partial_coherence_success, &
        partial_coherence_singular, partial_coherence_invalid_input
    implicit none
    type(multitaper_matrix_result) :: s, scaled_s
    type(multitaper_spectral_matrix_spectrogram_result) :: ts
    type(partial_coherence_result) :: p, q, r, default
    type(partial_coherence_tf_result) :: tf, tf_default, by_signal
    real(dp) :: x(96,3), d(3), t, nan
    complex(dp) :: z(96,3)
    integer :: f, i, j

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    allocate(s%frequencies(3), s%power(3,3,3))
    s%status = 0
    s%frequencies = [0.0_dp, 0.1_dp, 0.2_dp]
    s%power = cmplx(0.0_dp,0.0_dp,dp)
    do f = 1, 3
        do i = 1, 3
            s%power(f,i,i) = cmplx(1.0_dp+0.1_dp*i,0.0_dp,dp)
        end do
        s%power(f,1,2) = cmplx(0.30_dp,0.12_dp,dp)
        s%power(f,2,1) = conjg(s%power(f,1,2))
        s%power(f,2,3) = cmplx(-0.16_dp,0.08_dp,dp)
        s%power(f,3,2) = conjg(s%power(f,2,3))
    end do

    default = partial_coherence_matrix(s)
    if (size(default%coherency) /= 0) error stop 'unrequested coherency allocated'
    p = partial_coherence_matrix(s, shrinkage=0.2_dp, return_coherency=.true.)
    if (p%status /= partial_coherence_success) error stop 'complex coherency invalid'
    if (any(shape(p%coherency) /= [3,3,3])) error stop 'coherency shape'
    do f = 1, 3
        if (maxval(abs(p%coherency(f,:,:)-transpose(conjg(p%coherency(f,:,:))))) > 1.0e-14_dp) &
            error stop 'coherency must be Hermitian'
        if (maxval(abs(abs(p%coherency(f,:,:))**2-p%coherence(f,:,:))) > 1.0e-14_dp) &
            error stop 'complex magnitude inconsistent with coherence'
        do i = 1, 3
            if (abs(p%coherency(f,i,i)-1.0_dp) > 0.0_dp) error stop 'complex diagonal not one'
        end do
    end do
    if (abs(aimag(p%coherency(1,1,2))) < 1.0e-3_dp) error stop 'complex phase lost'

    ! Positive rescaling of each channel independently commutes with diagonal
    ! shrinkage (unlike trace-ridge regularization).
    d = [0.17_dp, 3.5_dp, 12.0_dp]
    scaled_s = s
    do j = 1, 3
        do i = 1, 3
            scaled_s%power(:,i,j) = s%power(:,i,j)*d(i)*d(j)
        end do
    end do
    q = partial_coherence_matrix(scaled_s, shrinkage=0.2_dp, return_coherency=.true.)
    if (q%status /= partial_coherence_success) error stop 'scaled channels singular'
    if (maxval(abs(q%coherency-p%coherency)) > 3.0e-14_dp) &
        error stop 'coherency changes under independent positive channel scaling'

    r = partial_coherence_matrix(s, shrinkage=1.0_dp, return_coherency=.true.)
    if (r%status /= partial_coherence_success) error stop 'full diagonal shrinkage failed'
    if (maxval(abs(r%coherency(:,1,2))) > 1.0e-14_dp) error stop 'full diagonal shrinkage not zero'
    if (maxval(abs(r%coherence(:,1,2))) > 1.0e-14_dp) error stop 'full diagonal coherence not zero'

    ! Rank-one input must be singular without shrinkage but recover with it.
    s%power(3,:,:) = cmplx(1.0_dp,0.0_dp,dp)
    p = partial_coherence_matrix(s, return_coherency=.true.)
    if (p%status /= partial_coherence_singular) error stop 'rank deficient not detected'
    if (p%valid(3)) error stop 'rank deficient bin marked valid'
    if (.not. ieee_is_nan(real(p%coherency(3,1,2),dp))) error stop 'singular NaN missing'
    q = partial_coherence_matrix(s, shrinkage=0.1_dp, return_coherency=.true.)
    if (q%status /= partial_coherence_success) error stop 'shrinkage did not rescue singular input'
    if (.not. all(q%valid)) error stop 'shrinkage failed at one bin'

    ! Invalid shrinkage must reject with empty arrays; default behavior is kept.
    p = partial_coherence_matrix(s, shrinkage=-0.01_dp, return_coherency=.true.)
    if (p%status /= partial_coherence_invalid_input) error stop 'negative shrinkage accepted'
    if (size(p%coherency) /= 0) error stop 'invalid result has coherency'
    p = partial_coherence_matrix(s, shrinkage=1.01_dp)
    if (p%status /= partial_coherence_invalid_input) error stop 'shrinkage greater than one accepted'
    p = partial_coherence_matrix(s, shrinkage=nan)
    if (p%status /= partial_coherence_invalid_input) error stop 'NaN shrinkage accepted'

    allocate(ts%frequencies(3),ts%times(2),ts%power(3,2,3,3))
    ts%status = 0
    ts%frequencies = s%frequencies
    ts%times = [0.5_dp, 1.0_dp]
    ts%power(:,1,:,:) = s%power
    ts%power(:,2,:,:) = scaled_s%power
    tf = partial_coherence_tf_matrix(ts, shrinkage=0.15_dp, return_coherency=.true.)
    if (tf%status /= partial_coherence_success) error stop 'time-frequency shrinkage failed'
    if (any(shape(tf%coherency) /= [3,2,3,3])) error stop 'time-frequency coherency shape'
    do j = 1, 2
        if (j == 1) then
            q = partial_coherence_matrix(s, shrinkage=0.15_dp, return_coherency=.true.)
        else
            q = partial_coherence_matrix(scaled_s, shrinkage=0.15_dp, return_coherency=.true.)
        end if
        if (maxval(abs(q%coherency-tf%coherency(:,j,:,:))) > 1.0e-13_dp) &
            error stop 'time-frequency and static coherency disagree'
    end do
    tf_default = partial_coherence_tf_matrix(ts)
    if (size(tf_default%coherency) /= 0) error stop 'unrequested TF coherency allocated'
    ts%power(3,1,:,:) = cmplx(1.0_dp,0.0_dp,dp)
    tf = partial_coherence_tf_matrix(ts, return_coherency=.true.)
    if (tf%status /= partial_coherence_singular) error stop 'TF singular not flagged'
    if (.not. ieee_is_nan(real(tf%coherency(3,1,1,2),dp))) error stop 'TF NaN missing'
    if (.not. tf%valid(3,2)) error stop 'valid neighboring TF cell lost'

    ! Both signal convenience overloads and the noncontiguous sample axis.
    do j = 1, 3
        do i = 1, 96
            t = real(i,dp)
            x(i,j) = sin((0.13_dp+0.035_dp*j)*t)+0.2_dp*cos((0.07_dp*j)*t)
            z(i,j) = cmplx(x(i,j),0.3_dp*cos(t*(0.11_dp+0.04_dp*j)),dp)
        end do
    end do
    p = multitaper_partial_coherence(x, 3.0_dp, kmax=5, &
        shrinkage=0.1_dp, return_coherency=.true.)
    if (p%status /= partial_coherence_success) error stop 'real signal static failed'
    if (size(p%coherency) == 0) error stop 'static signal coherency missing'
    p = multitaper_partial_coherence(z, 3.0_dp, kmax=5, &
        shrinkage=0.1_dp, return_coherency=.true.)
    if (p%status /= partial_coherence_success) error stop 'complex signal static failed'
    by_signal = multitaper_tf_partial_coherence(x,3.0_dp,32,noverlap=16, &
        kmax=5,shrinkage=0.15_dp,return_coherency=.true.)
    if (by_signal%status /= partial_coherence_success) error stop 'real signal TF failed'
    if (maxval(abs(abs(by_signal%coherency)**2-by_signal%coherence)) > 1.0e-13_dp) &
        error stop 'signal TF complex power mismatch'
    tf = multitaper_tf_partial_coherence(transpose(z),3.0_dp,32, &
        noverlap=16,kmax=5,axis=2,shrinkage=0.2_dp,return_coherency=.true.)
    if (tf%status /= partial_coherence_success) error stop 'complex signal TF axis=2 failed'
    tf = multitaper_tf_partial_coherence(z(1:95:2,1:3:2),2.5_dp,24, &
        noverlap=12,kmax=4,shrinkage=0.1_dp,return_coherency=.true.)
    if (tf%status /= partial_coherence_success) error stop 'strided complex signal TF failed'
    tf = multitaper_tf_partial_coherence(x,3.0_dp,32,shrinkage=2.0_dp)
    if (tf%status /= partial_coherence_invalid_input) error stop 'TF invalid shrinkage accepted'

    print '(a)', 'PASS complex partial coherency and diagonal spectral shrinkage'
end program test_signal_partial_coherency
