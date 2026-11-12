! SPDX-License-Identifier: MIT
program test_signal_multiple_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram
    use scifort_signal_multiple_coherence, only : multiple_coherence_result, &
        multiple_coherence_tf_result, multiple_coherence_matrix, multiple_coherence_tf_matrix, &
        multitaper_multiple_coherence, multitaper_tf_multiple_coherence
    use scifort_signal_partial_coherence, only : partial_coherence_success, &
        partial_coherence_invalid_input, partial_coherence_singular
    implicit none

    type(multitaper_matrix_result) :: s, spect, sc
    type(multitaper_spectral_matrix_spectrogram_result) :: tf, sp_tf
    type(multiple_coherence_result) :: r, rr, rc, real_direct, complex_direct
    type(multiple_coherence_tf_result) :: rt, direct_tf
    real(dp) :: x(96,3), expected, alpha, lambda, theta
    complex(dp) :: cx(96,3)
    integer :: i, j, t

    ! Exact two-channel coherence with a phase-bearing cross spectrum.
    allocate(s%frequencies(2),s%power(2,2,2))
    s%status = 0
    s%frequencies = [0.25_dp,0.5_dp]
    s%power = (0.0_dp,0.0_dp)
    s%power(1,1,1) = (1.0_dp,0.0_dp)
    s%power(1,2,2) = (2.0_dp,0.0_dp)
    s%power(1,1,2) = (0.6_dp,0.4_dp)
    s%power(1,2,1) = conjg(s%power(1,1,2))
    r = multiple_coherence_matrix(s)
    call require(r%status == partial_coherence_singular,'zero bin reports singular')
    call require(all(r%valid .eqv. [.true.,.false.]),'valid mask')
    expected = (0.6_dp**2+0.4_dp**2)/2.0_dp
    call near(r%coherence(1,1),expected,2.0e-13_dp,'two-channel squared multiple coherence')
    call near(r%coherence(1,2),expected,2.0e-13_dp,'two-channel symmetry')
    call near(r%conditional_power(1,1),1.0_dp-expected,2.0e-13_dp,'first residual PSD')
    call near(r%conditional_power(1,2),2.0_dp*(1.0_dp-expected),2.0e-13_dp, &
        'second residual PSD')
    call near(r%explained_power(1,1),expected,2.0e-13_dp,'explained PSD')
    call require(ieee_is_nan(r%coherence(2,1)), 'singular bin returns NaN')

    alpha = 0.25_dp
    lambda = 0.10_dp
    rr = multiple_coherence_matrix(s, ridge=lambda, shrinkage=alpha)
    expected = (0.6_dp**2+0.4_dp**2)*(1.0_dp-alpha)**2 / &
        ((1.0_dp+1.5_dp*lambda)*(2.0_dp+1.5_dp*lambda))
    call near(rr%coherence(1,1),expected,2.0e-13_dp,'ridge plus shrinkage')
    call near(rr%conditional_power(1,1), &
        (1.0_dp+1.5_dp*lambda)*(1.0_dp-expected),2.0e-13_dp,'regularized residual PSD')
    call near(rr%explained_power(1,1)+rr%conditional_power(1,1), &
        1.0_dp+1.5_dp*lambda,2.0e-13_dp,'power decomposition')
    call require(.not. rr%valid(2),'zero spectral matrix remains singular with trace ridge')


    ! A nonzero rank-one cross-spectral matrix needs ridge or shrinkage.
    allocate(sc%frequencies(1),sc%power(1,2,2))
    sc%status = 0
    sc%frequencies = [0.3_dp]
    sc%power = (1.0_dp,0.0_dp)
    rc = multiple_coherence_matrix(sc)
    call require(rc%status == partial_coherence_singular,'nonzero rank deficient matrix')
    rc = multiple_coherence_matrix(sc,ridge=0.03_dp)
    call require(rc%status == partial_coherence_success,'ridge rescues nonzero rank deficiency')
    expected = 1.0_dp/(1.0_dp+0.03_dp)**2
    call near(rc%coherence(1,1),expected,3.0e-13_dp,'ridge on rank-one matrix')
    rc = multiple_coherence_matrix(sc,shrinkage=1.0_dp)
    call require(rc%status == partial_coherence_success,'full shrinkage rescues rank deficiency')
    call near(rc%coherence(1,1),0.0_dp,3.0e-13_dp,'full diagonal is independent')
    ! Malformed Hermitian data must not produce plausible multiple coherence.
    sc%power(1,1,2) = cmplx(0.4_dp,0.2_dp,dp)
    rc = multiple_coherence_matrix(sc)
    call require(rc%status == partial_coherence_invalid_input,'reject non-Hermitian PSD')

    ! Full diagonal shrinkage destroys all interchannel dependencies.
    rr = multiple_coherence_matrix(s,shrinkage=1.0_dp)
    call near(maxval(abs(rr%coherence(1,:))),0.0_dp,5.0e-14_dp,'full shrinkage')
    rr = multiple_coherence_matrix(s,shrinkage=-0.001_dp)
    call require(rr%status == partial_coherence_invalid_input,'invalid shrinkage')
    call require(size(rr%coherence) == 0,'invalid output allocation')
    rr = multiple_coherence_matrix(s,cv_shrinkage=.true.)
    call require(rr%status == partial_coherence_invalid_input,'missing CV metadata')
    rr = multiple_coherence_matrix(s,ridge=-0.1_dp)
    call require(rr%status == partial_coherence_invalid_input,'invalid ridge')
    rr = multiple_coherence_matrix(s,auto_shrinkage=.true.)
    call require(rr%status == partial_coherence_invalid_input,'missing taper count')
    rr = multiple_coherence_matrix(s,auto_shrinkage=.true.,effective_tapers=6.0_dp)
    call require(rr%status == partial_coherence_singular,'auto shrinkage with zero bin')
    call require(rr%shrinkage_used(1) >= 0.0_dp .and. &
        rr%shrinkage_used(1) <= 1.0_dp, 'automatic coefficient bounds')

    s%cv_shrinkage = [0.45_dp,0.95_dp]
    rr = multiple_coherence_matrix(s,cv_shrinkage=.true.)
    rc = multiple_coherence_matrix(s,shrinkage=s%cv_shrinkage(1))
    call near(rr%coherence(1,1),rc%coherence(1,1),3.0e-13_dp,'CV coefficient applied')
    call near(rr%shrinkage_used(1),0.45_dp,0.0_dp,'CV coefficient retained')
    rr = multiple_coherence_matrix(s,cv_shrinkage=.true.,shrinkage=0.4_dp)
    call require(rr%status == partial_coherence_invalid_input,'CV and manual conflict')
    s%cv_shrinkage(1) = -0.01_dp
    rr = multiple_coherence_matrix(s,cv_shrinkage=.true.)
    call require(rr%status == partial_coherence_invalid_input,'reject invalid CV coefficients')
    s%cv_shrinkage(1) = 0.45_dp

    ! Time-frequency matrices must give the same answer as isolated frames.
    allocate(tf%frequencies(2),tf%times(3),tf%power(2,3,2,2),tf%cv_shrinkage(2,3))
    tf%status = 0
    tf%frequencies = s%frequencies
    tf%times = [0.5_dp,1.0_dp,1.5_dp]
    do t = 1, 3
        tf%power(:,t,:,:) = s%power*real(t,dp)
        tf%cv_shrinkage(:,t) = min(1.0_dp,0.15_dp*real(t,dp))
    end do
    rt = multiple_coherence_tf_matrix(tf,ridge=0.02_dp,cv_shrinkage=.true.)
    call require(rt%status == partial_coherence_singular,'time-frequency singular bin')
    do t = 1, 3
        sc%status = 0
        sc%frequencies = tf%frequencies
        sc%power = tf%power(:,t,:,:)
        sc%cv_shrinkage = tf%cv_shrinkage(:,t)
        rc = multiple_coherence_matrix(sc,ridge=0.02_dp,cv_shrinkage=.true.)
        call near(maxval(abs(rt%coherence(1,t,:)-rc%coherence(1,:))), &
            0.0_dp,2.0e-13_dp,'time-frequency matrix vs frame')
        call near(maxval(abs(rt%conditional_power(1,t,:)-rc%conditional_power(1,:))), &
            0.0_dp,2.0e-13_dp,'time-frequency conditional PSD')
    end do
    tf%cv_shrinkage(1,1) = 1.01_dp
    rt = multiple_coherence_tf_matrix(tf,cv_shrinkage=.true.)
    call require(rt%status == partial_coherence_invalid_input,'time-frequency invalid metadata')

    ! Real/complex end-to-end multitaper, comparing precomputed matrix API.
    do i = 1, 96
        theta = real(i,dp)
        x(i,1) = sin(0.13_dp*theta)+0.3_dp*cos(0.29_dp*theta)
        x(i,2) = 0.4_dp*x(i,1)+sin(0.31_dp*theta)
        x(i,3) = 0.6_dp*x(i,2)+0.4_dp*cos(0.41_dp*theta)
        do j = 1, 3
            cx(i,j) = cmplx(x(i,j),0.2_dp*sin(real(j,dp)*0.073_dp*theta),dp)
        end do
    end do
    spect = multitaper_spectral_matrix(x,3.5_dp,kmax=6,nfft=128,weighting='eigen')
    real_direct = multitaper_multiple_coherence(x,3.5_dp,kmax=6,nfft=128, &
        weighting='eigen',ridge=0.01_dp,auto_shrinkage=.true.)
    rc = multiple_coherence_matrix(spect,ridge=0.01_dp,auto_shrinkage=.true., &
        effective_tapers=spect%effective_tapers)
    call require(real_direct%status == partial_coherence_success,'static real direct status')
    call near(maxval(abs(real_direct%coherence-rc%coherence)), &
        0.0_dp,2.0e-13_dp,'real static direct vs spectral')
    call near(maxval(abs(real_direct%conditional_power-rc%conditional_power)), &
        0.0_dp,2.0e-12_dp,'real static conditional PSD')

    rc = multitaper_multiple_coherence(cx,3.5_dp,kmax=6,nfft=128, &
        weighting='adaptive',auto_shrinkage=.true.)
    call require(rc%status == partial_coherence_invalid_input, &
        'reject Gaussian-Wishart automatic shrinkage for adaptive weights')
    spect = multitaper_spectral_matrix(cx,3.5_dp,kmax=6,nfft=128, &
        weighting='adaptive',cv_shrinkage=.true.)
    complex_direct = multitaper_multiple_coherence(cx,3.5_dp,kmax=6,nfft=128, &
        weighting='adaptive',ridge=0.02_dp,cv_shrinkage=.true.)
    rc = multiple_coherence_matrix(spect,ridge=0.02_dp,cv_shrinkage=.true.)
    call require(complex_direct%status == partial_coherence_success,'static adaptive CV')
    call near(maxval(abs(complex_direct%coherence-rc%coherence)), &
        0.0_dp,2.0e-13_dp,'complex static adaptive CV match')
    call require(all(complex_direct%coherence >= 0.0_dp) .and. &
        all(complex_direct%coherence <= 1.0_dp),'multiple coherence is bounded')

    sp_tf = multitaper_spectral_matrix_spectrogram(cx,3.5_dp,48,noverlap=24, &
        kmax=5,nfft=64,weighting='eigen')
    direct_tf = multitaper_tf_multiple_coherence(cx,3.5_dp,48,noverlap=24, &
        kmax=5,nfft=64,weighting='eigen',ridge=0.01_dp,auto_shrinkage=.true.)
    rt = multiple_coherence_tf_matrix(sp_tf,ridge=0.01_dp, &
        auto_shrinkage=.true.,effective_tapers=sp_tf%effective_tapers)
    call require(direct_tf%status == partial_coherence_success,'time-frequency complex direct')
    call near(maxval(abs(direct_tf%coherence-rt%coherence)), &
        0.0_dp,2.0e-13_dp,'complex time-frequency direct vs spectral')
    call near(maxval(abs(direct_tf%conditional_power-rt%conditional_power)), &
        0.0_dp,2.0e-12_dp,'complex time-frequency PSD')

    sp_tf = multitaper_spectral_matrix_spectrogram(x,3.5_dp,48,noverlap=24, &
        kmax=5,nfft=64,weighting='adaptive',cv_shrinkage=.true.)
    direct_tf = multitaper_tf_multiple_coherence(x,3.5_dp,48,noverlap=24, &
        kmax=5,nfft=64,weighting='adaptive',ridge=0.05_dp,cv_shrinkage=.true.)
    rt = multiple_coherence_tf_matrix(sp_tf,ridge=0.05_dp,cv_shrinkage=.true.)
    call require(direct_tf%status == partial_coherence_success,'time-frequency real CV')
    call near(maxval(abs(direct_tf%coherence-rt%coherence)), &
        0.0_dp,2.0e-13_dp,'real time-frequency adaptive CV')
    call near(maxval(abs(direct_tf%conditional_power-rt%conditional_power)), &
        0.0_dp,2.0e-12_dp,'real time-frequency adaptive CV residual PSD')

    print *, 'PASS multiple coherence native regressions'

contains
    subroutine require(condition,description)
        logical, intent(in) :: condition !! test predicate
        character(*), intent(in) :: description !! failure label
        if (.not. condition) then
            print *, 'FAIL: ',description
            error stop 1
        end if
    end subroutine require
    subroutine near(actual,reference,tol,description)
        real(dp), intent(in) :: actual !! measured scalar
        real(dp), intent(in) :: reference !! expected scalar
        real(dp), intent(in) :: tol !! allowed absolute error
        character(*), intent(in) :: description !! assertion label
        call require(abs(actual-reference) <= tol,description)
    end subroutine near
end program test_signal_multiple_coherence
