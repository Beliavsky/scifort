! SPDX-License-Identifier: MIT
! Covariance-metric canonical subspace tracking and conditional spectral metrics.
program test_signal_canonical_covmetric_tracking
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, &
        canonical_coherence_tf_result, canonical_coherence_matrix, canonical_coherence_tf_matrix, &
        band_canonical_coherence_matrix, band_canonical_coherence_tf_matrix
    use scifort_signal_canonical_subspace_tracking, only : track_canonical_subspaces, &
        canonical_subspace_tracks_result, canonical_subspace_tracks_tf_result, &
        canonical_subspace_success, canonical_subspace_invalid_input
    implicit none
    type(canonical_coherence_result) :: base, scaled, invalid, result, band
    type(canonical_coherence_tf_result) :: tf, tf_result
    type(canonical_subspace_tracks_result) :: tbase, tscaled, tlegacy
    type(canonical_subspace_tracks_tf_result) :: ttf
    type(multitaper_matrix_result) :: spectrum
    type(multitaper_spectral_matrix_spectrogram_result) :: spec_tf
    complex(dp) :: d(3), phase, input(5,5), cs(5,5), reg(5,5), a(2,2), b(2,2)
    real(dp) :: theta, bands(2,2), mean_power
    integer :: f, i, j, k

    ! Two nearly tied nonzero modes, with complex rotation and one direction
    ! changing its joint A/B subspace between bins.
    allocate(base%frequencies(3),base%valid(3),base%squared_modes(3,2), &
        base%coefficients_a(3,3,2),base%coefficients_b(3,3,2), &
        base%covariance_a(3,3,3),base%covariance_b(3,3,3))
    base%frequencies = [0.0_dp,1.0_dp,2.0_dp]
    base%valid = .true.
    base%squared_modes(:,1) = 0.750000001_dp
    base%squared_modes(:,2) = 0.75_dp
    base%coefficients_a = (0.0_dp,0.0_dp)
    base%coefficients_b = (0.0_dp,0.0_dp)
    base%covariance_a = (0.0_dp,0.0_dp)
    base%covariance_b = (0.0_dp,0.0_dp)
    do f = 1, 3
        do i = 1, 3
            base%covariance_a(f,i,i) = (1.0_dp,0.0_dp)
            base%covariance_b(f,i,i) = (1.0_dp,0.0_dp)
        end do
        base%coefficients_a(f,2,2) = (1.0_dp,0.0_dp)
        base%coefficients_b(f,2,2) = (1.0_dp,0.0_dp)
    end do
    theta = 0.32_dp
    phase = cmplx(cos(0.6_dp),sin(0.6_dp),dp)
    base%coefficients_a(1,1,1) = (1.0_dp,0.0_dp)
    base%coefficients_b(1,1,1) = (1.0_dp,0.0_dp)
    do f = 2, 3
        base%coefficients_a(f,1,1) = cmplx(cos(theta),0.0_dp,dp)*phase
        base%coefficients_a(f,3,1) = cmplx(sin(theta),0.0_dp,dp)*phase
        base%coefficients_b(f,1,1) = cmplx(cos(theta),0.0_dp,dp)*phase
        base%coefficients_b(f,3,1) = cmplx(sin(theta),0.0_dp,dp)*phase
    end do
    tbase = track_canonical_subspaces(base,min_similarity=0.8_dp,covariance_metric=.true.)
    call require(tbase%status == canonical_subspace_success,'identity metric accepted')
    call require(tbase%covariance_metric_used,'metric diagnostic true')
    call require(all(tbase%matched(2:,:)),'continuity with metric')
    call require(abs(tbase%principal_cosines(2,2)-cos(theta)) < 1.0e-12_dp, &
        'metric principal angle of moving direction')
    tlegacy = track_canonical_subspaces(base,min_similarity=0.8_dp)
    call require(tlegacy%status == canonical_subspace_success,'legacy default')
    call require(.not. tlegacy%covariance_metric_used,'legacy diagnostic false')
    call require(maxval(abs(tlegacy%basis_a-tbase%basis_a)) < 1.0e-12_dp, &
        'identity metric reproduces legacy aligned bases')

    ! Rescale each original input channel by a different unit, with a complex
    ! phase on the third channel. Contravariant coefficient transformation and
    ! covariant conditional spectra leave the invariant comparisons unchanged.
    scaled = base
    d(1) = (100.0_dp,0.0_dp)
    d(2) = (0.01_dp,0.0_dp)
    d(3) = (0.0_dp,30.0_dp)
    do f = 1, 3
        do i = 1, 3
            scaled%coefficients_a(f,i,:) = base%coefficients_a(f,i,:)/conjg(d(i))
            scaled%coefficients_b(f,i,:) = base%coefficients_b(f,i,:)/conjg(d(i))
            do j = 1, 3
                scaled%covariance_a(f,i,j) = d(i)*base%covariance_a(f,i,j)*conjg(d(j))
                scaled%covariance_b(f,i,j) = d(i)*base%covariance_b(f,i,j)*conjg(d(j))
            end do
        end do
    end do
    tscaled = track_canonical_subspaces(scaled,min_similarity=0.8_dp,covariance_metric=.true.)
    call require(tscaled%status == canonical_subspace_success,'scaled metric accepted')
    call require(all(tscaled%track_id == tbase%track_id),'scaling leaves block labels')
    call require(all(tscaled%matched .eqv. tbase%matched),'scaling leaves assignments')
    call require(maxval(abs(tscaled%principal_cosines(2:,:)- &
        tbase%principal_cosines(2:,:))) < 2.0e-12_dp,'principal angles unit-invariant')
    call require(maxval(abs(tscaled%similarity(2:,:)-tbase%similarity(2:,:))) < 2.0e-12_dp, &
        'worst cosine unit-invariant')
    do f = 1, 3
        do i = 1, 3
            call require(maxval(abs(tscaled%basis_a(f,i,:)*conjg(d(i))- &
                tbase%basis_a(f,i,:))) < 2.0e-11_dp,'aligned A basis scales correctly')
            call require(maxval(abs(tscaled%basis_b(f,i,:)*conjg(d(i))- &
                tbase%basis_b(f,i,:))) < 2.0e-11_dp,'aligned B basis scales correctly')
        end do
    end do
    tlegacy = track_canonical_subspaces(scaled,min_similarity=0.8_dp)
    call require(.not. all(tlegacy%matched(2,:)), &
        'legacy Euclidean matching breaks under extreme channel scaling')

    invalid = base
    deallocate(invalid%covariance_a)
    tbase = track_canonical_subspaces(invalid,covariance_metric=.true.)
    call require(tbase%status == canonical_subspace_invalid_input,'missing covariances rejected')
    invalid = base
    invalid%covariance_a(1,1,2) = cmplx(0.5_dp,0.2_dp,dp)
    tbase = track_canonical_subspaces(invalid,covariance_metric=.true.)
    call require(tbase%status == canonical_subspace_invalid_input,'non-Hermitian metric rejected')
    invalid = base
    invalid%covariance_a(1,1,1) = (-1.0_dp,0.0_dp)
    tbase = track_canonical_subspaces(invalid,covariance_metric=.true.)
    call require(tbase%status == canonical_subspace_invalid_input,'nonpositive metric rejected')
    invalid = base
    invalid%covariance_a(1,1,1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    tbase = track_canonical_subspaces(invalid,covariance_metric=.true.)
    call require(tbase%status == canonical_subspace_invalid_input,'NaN metric rejected')
    invalid = base
    invalid%valid(1) = .false.
    invalid%covariance_a(1,:,:) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    tbase = track_canonical_subspaces(invalid,covariance_metric=.true.)
    call require(tbase%status == canonical_subspace_success,'invalid masked bin skips metric')
    call require(all(tbase%segment_start(2,:)),'invalid masked metric causes restart')

    allocate(tf%frequencies(3),tf%times(2),tf%valid(3,2),tf%squared_modes(3,2,2), &
        tf%coefficients_a(3,2,3,2),tf%coefficients_b(3,2,3,2), &
        tf%covariance_a(3,2,3,3),tf%covariance_b(3,2,3,3))
    tf%frequencies = base%frequencies
    tf%times = [0.0_dp,1.0_dp]
    tf%valid = .true.
    do k = 1, 2
        tf%squared_modes(:,k,:) = base%squared_modes
        tf%coefficients_a(:,k,:,:) = base%coefficients_a
        tf%coefficients_b(:,k,:,:) = base%coefficients_b
        tf%covariance_a(:,k,:,:) = base%covariance_a
        tf%covariance_b(:,k,:,:) = base%covariance_b
    end do
    tf%valid(2,2) = .false.
    tf%covariance_a(2,2,:,:) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    ttf = track_canonical_subspaces(tf,covariance_metric=.true.)
    call require(ttf%status == canonical_subspace_success,'time-frequency metric accepted')
    call require(ttf%covariance_metric_used,'TF metric mode recorded')
    call require(all(ttf%matched(2,1,:)),'first window matched')
    call require(all(ttf%track_id(2,2,:) == 0),'second window invalid mask')
    call require(all(ttf%segment_start(3,2,:)),'independent second-window restart')

    ! Verify exported conditional spectral metrics use exactly the same
    ! diagonal shrinkage, trace ridge and Schur complement as the eigenproblem.
    allocate(spectrum%frequencies(3),spectrum%power(3,5,5))
    spectrum%frequencies = base%frequencies
    spectrum%status = 0
    input = (0.0_dp,0.0_dp)
    do i = 1, 5
        input(i,i) = cmplx(1.0_dp+0.1_dp*i,0.0_dp,dp)
        do j = i+1, 5
            input(i,j) = cmplx(0.04_dp*i*j,0.02_dp*(j-i),dp)
        end do
    end do
    cs = matmul(input,transpose(conjg(input)))
    do f = 1, 3
        spectrum%power(f,:,:) = cs*(1.0_dp+0.06_dp*f)
    end do
    result = canonical_coherence_matrix(spectrum,[1,2],[3,4],conditioned_on=[5], &
        ridge=0.04_dp,shrinkage=0.15_dp,return_vectors=.true.,return_covariances=.true.)
    call require(result%status == 0 .and. all(result%valid),'conditional estimator valid')
    call require(all(shape(result%covariance_a)==[3,2,2]),'A covariance allocation')
    call require(all(shape(result%covariance_b)==[3,2,2]),'B covariance allocation')
    do f = 1, 3
        mean_power = sum([(real(spectrum%power(f,i,i),dp),i=1,5)])/5.0_dp
        reg = 0.85_dp*spectrum%power(f,:,:)
        do i = 1, 5
            reg(i,i) = cmplx(real(spectrum%power(f,i,i),dp)+0.04_dp*mean_power,0.0_dp,dp)
        end do
        do j = 1, 2
            do i = 1, 2
                a(i,j) = reg(i,j)-reg(i,5)*reg(5,j)/reg(5,5)
                b(i,j) = reg(i+2,j+2)-reg(i+2,5)*reg(5,j+2)/reg(5,5)
            end do
        end do
        call require(maxval(abs(result%covariance_a(f,:,:)-a)) < 1.0e-12_dp, &
            'returned conditional A covariance matches Schur complement')
        call require(maxval(abs(result%covariance_b(f,:,:)-b)) < 1.0e-12_dp, &
            'returned conditional B covariance matches Schur complement')
    end do
    tbase = track_canonical_subspaces(result,covariance_metric=.true.)
    call require(tbase%status == canonical_subspace_success,'estimated metrics accepted')
    result = canonical_coherence_matrix(spectrum,[1,2],[3,4])
    call require(.not. allocated(result%covariance_a),'default metrics unallocated')
    result = canonical_coherence_matrix(spectrum,[1,2],[3,4],return_covariances=.true.)
    call require(allocated(result%covariance_a) .and. .not. allocated(result%coefficients_a), &
        'metrics without vectors supported')
    bands(:,1) = [0.1_dp,0.8_dp]
    bands(:,2) = [1.1_dp,1.8_dp]
    band = band_canonical_coherence_matrix(spectrum,bands,[1,2],[3,4],return_vectors=.true., &
        return_covariances=.true.)
    call require(all(shape(band%covariance_a)==[2,2,2]),'band conditional spectra')
    allocate(spec_tf%frequencies(3),spec_tf%times(2),spec_tf%power(3,2,5,5))
    spec_tf%frequencies = spectrum%frequencies
    spec_tf%times = [0.0_dp,1.0_dp]
    spec_tf%status = 0
    do k = 1, 2
        spec_tf%power(:,k,:,:) = spectrum%power
    end do
    tf_result = canonical_coherence_tf_matrix(spec_tf,[1,2],[3,4],return_vectors=.true., &
        return_covariances=.true.)
    call require(all(shape(tf_result%covariance_a)==[3,2,2,2]),'TF conditional spectra')
    ttf = track_canonical_subspaces(tf_result,covariance_metric=.true.)
    call require(ttf%status == canonical_subspace_success,'TF estimated metrics accepted')
    tf_result = band_canonical_coherence_tf_matrix(spec_tf,bands,[1,2],[3,4], &
        return_vectors=.true.,return_covariances=.true.)
    call require(all(shape(tf_result%covariance_a)==[2,2,2,2]),'TF band conditional spectra')
    print *, 'PASS covariance-metric canonical tracking'

contains
    subroutine require(condition,label)
        logical, intent(in) :: condition !! required test assertion
        character(*), intent(in) :: label !! descriptive assertion label
        if (condition) return
        print *, 'FAIL: ',label
        error stop 1
    end subroutine require
end program test_signal_canonical_covmetric_tracking
