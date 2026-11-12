! SPDX-License-Identifier: MIT
program test_signal_canonical_coherence
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_group_coherence, only : group_coherence_matrix, group_coherence_result
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, &
        canonical_coherence_tf_result, canonical_coherence_matrix, &
        canonical_coherence_tf_matrix, band_canonical_coherence_matrix, &
        band_canonical_coherence_tf_matrix, multitaper_canonical_coherence, &
        multitaper_tf_canonical_coherence, multitaper_band_canonical_coherence, &
        multitaper_tf_band_canonical_coherence
    use scifort_signal_partial_coherence, only : partial_coherence_invalid_input, &
        partial_coherence_singular
    implicit none

    type(multitaper_matrix_result) :: s, v
    type(multitaper_spectral_matrix_spectrogram_result) :: st
    type(canonical_coherence_result) :: r, ri, band
    type(canonical_coherence_tf_result) :: rt, bt
    type(group_coherence_result) :: group
    complex(dp) :: z(96,5), mix(5,5), cov(5,5)
    real(dp) :: x(96,5), ang, bands(2,2), product_modes
    integer :: i, j, f, t

    mix = cmplx(0.0_dp,0.0_dp,dp)
    do i = 1, 5
        mix(i,i) = cmplx(1.1_dp+0.18_dp*i,0.0_dp,dp)
        do j = i+1, 5
            mix(i,j) = cmplx(0.12_dp*i*j,-0.03_dp*(i+j),dp)
        end do
    end do
    cov = matmul(mix,transpose(conjg(mix)))
    allocate(s%frequencies(4),s%power(4,5,5),s%cv_shrinkage(4))
    s%frequencies = [-1.0_dp,-0.4_dp,0.5_dp,1.2_dp]
    s%cv_shrinkage = [0.0_dp,0.1_dp,0.3_dp,0.8_dp]
    s%effective_tapers = 5.0_dp
    s%status = 0
    do f = 1, 4
        s%power(f,:,:) = cov*(1.0_dp+0.15_dp*f)
    end do
    r = canonical_coherence_matrix(s,[1,3],[2,4],conditioned_on=[5],ridge=0.02_dp)
    call require(r%status == 0 .and. all(r%valid),'conditional modes valid')
    call require(all(shape(r%squared_modes)==[4,2]),'mode output shape')
    do f = 1, 4
        call require(all(r%squared_modes(f,:) >= 0.0_dp) .and. &
            all(r%squared_modes(f,:) <= 1.0_dp),'mode bounds')
        call require(r%squared_modes(f,1) >= r%squared_modes(f,2),'descending eigenvalues')
        product_modes = 1.0_dp-product(1.0_dp-r%squared_modes(f,:))
        call near(product_modes,r%coherence(f),1.0e-12_dp,'Wilks product identity')
    end do
    ri = canonical_coherence_matrix(s,[4,2],[3,1],conditioned_on=[5],ridge=0.02_dp)
    call require(ri%status == 0,'permuted groups valid')
    call near(maxval(abs(ri%squared_modes-r%squared_modes)),0.0_dp,1.0e-12_dp, &
        'permutation and A/B symmetry')
    group = group_coherence_matrix(s,[1,3],[2,4],conditioned_on=[5],ridge=0.02_dp)
    call near(maxval(abs(group%coherence-r%coherence)),0.0_dp,1.0e-14_dp, &
        'same group aggregate')
    r = canonical_coherence_matrix(s,[1,3],[2,4],conditioned_on=[5],cv_shrinkage=.true.)
    call require(r%status == 0,'CV shrinkage')
    call near(maxval(abs(r%shrinkage_used-s%cv_shrinkage)),0.0_dp,0.0_dp,'CV coefficients')
    r = canonical_coherence_matrix(s,[1,3],[2,4],conditioned_on=[5],auto_shrinkage=.true.)
    call require(r%status == 0 .and. all(r%valid),'automatic shrinkage')
    ! Precisely known conditional roots, with nontrivial complex phase.
    ! An independent fifth channel does not change the roots.
    v = s
    v%power = cmplx(0.0_dp,0.0_dp,dp)
    do f = 1, 4
        do i = 1, 5
            v%power(f,i,i) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        v%power(f,1,3) = cmplx(0.0_dp,0.9_dp,dp)
        v%power(f,3,1) = conjg(v%power(f,1,3))
        v%power(f,2,4) = cmplx(0.4_dp,0.0_dp,dp)
        v%power(f,4,2) = conjg(v%power(f,2,4))
    end do
    r = canonical_coherence_matrix(v,[1,2],[3,4],conditioned_on=[5])
    call require(r%status == 0,'known complex-root matrix valid')
    call near(maxval(abs(r%squared_modes(:,1)-0.81_dp)),0.0_dp,1.0e-12_dp, &
        'first analytic root')
    call near(maxval(abs(r%squared_modes(:,2)-0.16_dp)),0.0_dp,1.0e-12_dp, &
        'second analytic root')
    v%power(:,2,4) = cmplx(0.0_dp,0.0_dp,dp)
    v%power(:,4,2) = cmplx(0.0_dp,0.0_dp,dp)
    r = canonical_coherence_matrix(v,[1,2],[3,4],conditioned_on=[5])
    call near(maxval(abs(r%squared_modes(:,2))),0.0_dp,1.0e-14_dp, &
        'zero canonical mode')
    v%power(:,1,3) = cmplx(0.0_dp,0.999999_dp,dp)
    v%power(:,3,1) = conjg(v%power(:,1,3))
    r = canonical_coherence_matrix(v,[1,2],[3,4],conditioned_on=[5])
    call require(r%status == 0,'nearly perfect mode matrix valid')
    call near(maxval(abs(r%squared_modes(:,1)-0.999999_dp**2)), &
        0.0_dp,2.0e-12_dp,'nearly perfect mode')
    r = canonical_coherence_matrix(s,[1,3],[2,4],conditioned_on=[5],shrinkage=1.0_dp)
    call near(maxval(abs(r%squared_modes)),0.0_dp,1.0e-14_dp,'full shrinkage independence')
    r = canonical_coherence_matrix(s,[1],[2],conditioned_on=[3,4,5],ridge=0.03_dp)
    group = group_coherence_matrix(s,[1],[2],conditioned_on=[3,4,5],ridge=0.03_dp)
    call near(maxval(abs(r%squared_modes(:,1)-group%coherence)),0.0_dp,1.0e-12_dp, &
        'singleton/partial coherence')
    r = canonical_coherence_matrix(s,[1],[2,3,4,5],ridge=0.03_dp)
    group = group_coherence_matrix(s,[1],[2,3,4,5],ridge=0.03_dp)
    call near(maxval(abs(r%squared_modes(:,1)-group%coherence)),0.0_dp,1.0e-12_dp, &
        'singleton/multiple coherence')

    allocate(st%frequencies(4),st%times(3),st%power(4,3,5,5),st%cv_shrinkage(4,3))
    st%frequencies = s%frequencies
    st%times = [0.1_dp,0.2_dp,0.3_dp]
    st%status = 0
    do t = 1, 3
        st%power(:,t,:,:) = s%power*real(t,dp)
        st%cv_shrinkage(:,t) = s%cv_shrinkage
    end do
    rt = canonical_coherence_tf_matrix(st,[1,3],[2,4],conditioned_on=[5],ridge=0.02_dp)
    call require(rt%status == 0 .and. all(rt%valid),'time-frequency valid')
    call require(all(shape(rt%squared_modes)==[4,3,2]),'time-frequency mode shape')
    ri = canonical_coherence_matrix(s,[1,3],[2,4],conditioned_on=[5],ridge=0.02_dp)
    do t = 1, 3
        call near(maxval(abs(rt%squared_modes(:,t,:)-ri%squared_modes)), &
            0.0_dp,1.0e-12_dp,'scale invariant time-frequency')
    end do
    bands(:,1) = [-0.8_dp,0.1_dp]
    bands(:,2) = [0.2_dp,1.0_dp]
    band = band_canonical_coherence_matrix(s,bands,[1,3],[2,4],conditioned_on=[5], &
        ridge=0.03_dp,shrinkage=0.1_dp)
    bt = band_canonical_coherence_tf_matrix(st,bands,[1,3],[2,4],conditioned_on=[5], &
        ridge=0.03_dp,shrinkage=0.1_dp)
    call require(band%status == 0 .and. bt%status == 0,'integrated band variants')
    do t = 1, 3
        call near(maxval(abs(bt%squared_modes(:,t,:)-band%squared_modes)), &
            0.0_dp,1.0e-12_dp,'band time-frequency consistency')
    end do

    r = canonical_coherence_matrix(s,[1,1],[2])
    call require(r%status == partial_coherence_invalid_input,'duplicate indices')
    call require(size(r%squared_modes)==0,'invalid result empty')
    r = canonical_coherence_matrix(s,[1],[2],conditioned_on=[1])
    call require(r%status == partial_coherence_invalid_input,'overlap rejected')
    r = canonical_coherence_matrix(s,[1],[2],shrinkage=1.1_dp)
    call require(r%status == partial_coherence_invalid_input,'shrinkage validation')
    r = canonical_coherence_matrix(s,[1],[2],auto_shrinkage=.true.,cv_shrinkage=.true.)
    call require(r%status == partial_coherence_invalid_input,'conflicting shrinkage')
    v = s
    v%power(2,:,:) = cmplx(0.0_dp,0.0_dp,dp)
    r = canonical_coherence_matrix(v,[1,3],[2,4],conditioned_on=[5])
    call require(r%status == partial_coherence_singular,'singular CSD status')
    call require(.not. r%valid(2) .and. r%valid(1),'cell singular mask')

    do i = 1, 96
        ang = 0.13_dp*real(i,dp)
        do j = 1, 5
            x(i,j) = sin(ang*j)+0.3_dp*cos(0.27_dp*i+0.4_dp*j)
            z(i,j) = cmplx(x(i,j),0.2_dp*sin(0.19_dp*i+0.7_dp*j),dp)
        end do
    end do
    r = multitaper_canonical_coherence(x,3.0_dp,[1,3],[2,4],conditioned_on=[5], &
        kmax=6,ridge=0.04_dp,shrinkage=0.2_dp)
    call require(r%status == 0 .and. all(r%valid),'real multitaper API')
    r = multitaper_canonical_coherence(z,3.0_dp,[1,3],[2,4],conditioned_on=[5], &
        kmax=6,ridge=0.04_dp,weighting='adaptive',cv_shrinkage=.true.)
    call require(r%status == 0 .and. all(r%valid),'complex adaptive CV API')
    rt = multitaper_tf_canonical_coherence(x,3.0_dp,64,[1,3],[2,4],conditioned_on=[5], &
        kmax=6,ridge=0.04_dp,shrinkage=0.2_dp)
    call require(rt%status == 0 .and. all(rt%valid),'real time-frequency API')
    rt = multitaper_tf_canonical_coherence(z,3.0_dp,64,[1,3],[2,4],conditioned_on=[5], &
        kmax=6,ridge=0.04_dp,weighting='eigen',auto_shrinkage=.true.)
    call require(rt%status == 0 .and. all(rt%valid),'complex TF eigen auto API')
    bands(:,1) = [0.04_dp,0.18_dp]
    bands(:,2) = [0.22_dp,0.40_dp]
    band = multitaper_band_canonical_coherence(x,3.0_dp,bands,[1,3],[2,4], &
        conditioned_on=[5],kmax=6,ridge=0.04_dp,shrinkage=0.2_dp)
    call require(band%status == 0 .and. all(band%valid),'real band API')
    band = multitaper_band_canonical_coherence(z,3.0_dp,bands,[1,3],[2,4], &
        conditioned_on=[5],kmax=6,ridge=0.04_dp,shrinkage=0.2_dp)
    call require(band%status == 0 .and. all(band%valid),'complex band API')
    bt = multitaper_tf_band_canonical_coherence(x,3.0_dp,64,bands,[1,3],[2,4], &
        conditioned_on=[5],kmax=6,ridge=0.04_dp,shrinkage=0.2_dp)
    call require(bt%status == 0 .and. all(bt%valid),'real time-frequency band API')
    bt = multitaper_tf_band_canonical_coherence(z,3.0_dp,64,bands,[1,3],[2,4], &
        conditioned_on=[5],kmax=6,ridge=0.04_dp,shrinkage=0.2_dp)
    call require(bt%status == 0 .and. all(bt%valid),'complex time-frequency band API')

    print *, 'PASS canonical spectral coherence native regression'
contains
    subroutine require(condition,label)
        logical, intent(in) :: condition !! condition that must hold
        character(*), intent(in) :: label !! test-case description
        if (condition) return
        print *, 'FAIL:', label
        error stop 1
    end subroutine require
    subroutine near(actual,expected,tolerance,label)
        real(dp), intent(in) :: actual !! computed scalar
        real(dp), intent(in) :: expected !! target scalar
        real(dp), intent(in) :: tolerance !! maximum absolute error
        character(*), intent(in) :: label !! test-case description
        if (abs(actual-expected) <= tolerance) return
        print *, 'FAIL:',label,actual,expected
        error stop 1
    end subroutine near
end program test_signal_canonical_coherence
