! SPDX-License-Identifier: MIT
program test_signal_canonical_vectors
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, &
        canonical_coherence_tf_result, canonical_coherence_matrix, canonical_coherence_tf_matrix, &
        band_canonical_coherence_matrix, band_canonical_coherence_tf_matrix, &
        multitaper_canonical_coherence, multitaper_tf_canonical_coherence, &
        multitaper_band_canonical_coherence, multitaper_tf_band_canonical_coherence
    use scifort_signal_partial_coherence, only : partial_coherence_invalid_input, &
        partial_coherence_singular
    implicit none
    type(multitaper_matrix_result) :: s, s_zero
    type(multitaper_spectral_matrix_spectrogram_result) :: st
    type(canonical_coherence_result) :: r, r2, band
    type(canonical_coherence_tf_result) :: tf, tf_band
    complex(dp) :: mixing(6,6), covariance(6,6), reg(6,6)
    complex(dp) :: laa(3,3), lbb(2,2), lab(3,2)
    complex(dp) :: actual_a(3,2), actual_b(2,2), cross(2,2), norms_a(2,2), norms_b(2,2)
    complex(dp) :: unit(2,2)
    real(dp) :: x(128,5), bands(2,2)
    complex(dp) :: z(128,5)
    integer :: f, i, j, k, t

    mixing = (0.0_dp,0.0_dp)
    do i = 1, 6
        mixing(i,i) = cmplx(1.3_dp+0.13_dp*i,0.0_dp,dp)
        do j = i+1, 6
            mixing(i,j) = cmplx(0.07_dp*i*j,0.09_dp*(j-i),dp)
        end do
    end do
    covariance = matmul(mixing,transpose(conjg(mixing)))
    allocate(s%frequencies(4),s%power(4,6,6),s%cv_shrinkage(4))
    s%frequencies = [0.0_dp,0.2_dp,0.4_dp,0.6_dp]
    s%cv_shrinkage = [0.0_dp,0.07_dp,0.3_dp,0.85_dp]
    s%effective_tapers = 6.0_dp
    s%status = 0
    do f = 1, 4
        s%power(f,:,:) = covariance*(0.9_dp+0.11_dp*f)
    end do
    r = canonical_coherence_matrix(s,[1,3,5],[2,4],conditioned_on=[6], &
        ridge=0.03_dp,shrinkage=0.08_dp,return_vectors=.true.)
    call require(r%status == 0 .and. all(r%valid),'conditional vectors valid')
    call require(all(shape(r%coefficients_a)==[4,3,2]),'A coefficients shape')
    call require(all(shape(r%coefficients_b)==[4,2,2]),'B coefficients shape')
    call require(all(shape(r%directions_stable)==[4,2]),'direction flags shape')
    call require(all(r%directions_stable),'isolated nonzero roots')
    ! Conditions are checked using scalar arrays below to remain portable.
    do f = 1, 4
        reg = (1.0_dp-0.08_dp)*s%power(f,:,:)
        do i = 1, 6
            reg(i,i) = cmplx(real(s%power(f,i,i),dp)+0.03_dp* &
                sum([(real(s%power(f,k,k),dp),k=1,6)])/6.0_dp,0.0_dp,dp)
        end do
        do j = 1, 3
            do i = 1, 3
                laa(i,j) = reg(idx_a(i),idx_a(j))- &
                    reg(idx_a(i),6)*reg(6,idx_a(j))/reg(6,6)
            end do
        end do
        do j = 1, 2
            do i = 1, 2
                lbb(i,j) = reg(idx_b(i),idx_b(j))- &
                    reg(idx_b(i),6)*reg(6,idx_b(j))/reg(6,6)
            end do
        end do
        do j = 1, 2
            do i = 1, 3
                lab(i,j) = reg(idx_a(i),idx_b(j))- &
                    reg(idx_a(i),6)*reg(6,idx_b(j))/reg(6,6)
            end do
        end do
        actual_a = r%coefficients_a(f,:,:)
        actual_b = r%coefficients_b(f,:,:)
        norms_a = matmul(transpose(conjg(actual_a)),matmul(laa,actual_a))
        norms_b = matmul(transpose(conjg(actual_b)),matmul(lbb,actual_b))
        cross = matmul(transpose(conjg(actual_a)),matmul(lab,actual_b))
        unit = (0.0_dp,0.0_dp)
        do k = 1, 2
            unit(k,k) = (1.0_dp,0.0_dp)
        end do
        call close_matrix(norms_a,unit,3.0e-11_dp,'group A unit covariance')
        call close_matrix(norms_b,unit,3.0e-11_dp,'group B unit covariance')
        unit = (0.0_dp,0.0_dp)
        do k = 1, 2
            unit(k,k) = cmplx(sqrt(r%squared_modes(f,k)),0.0_dp,dp)
        end do
        call close_matrix(cross,unit,3.0e-11_dp,'diagonal positive cross covariance')
        do k = 1, 2
            i = maxloc(abs(actual_a(:,k)),dim=1)
            call require(abs(aimag(actual_a(i,k))) < 1.0e-12_dp, 'real phase anchor')
            call require(real(actual_a(i,k),dp) > 0.0_dp,'positive phase anchor')
        end do
    end do
    r2 = canonical_coherence_matrix(s,[1,3,5],[2,4],conditioned_on=[6], &
        ridge=0.03_dp,shrinkage=0.08_dp)
    call require(.not. allocated(r2%coefficients_a),'default remains allocation-free')
    call require(.not. allocated(r2%directions_stable),'no default flags allocation')
    call require(maxval(abs(r%squared_modes-r2%squared_modes)) < 1.0e-12_dp, &
        'vector request leaves strengths unchanged')
    r2 = canonical_coherence_matrix(s,[2,4],[1,3,5],conditioned_on=[6], &
        ridge=0.03_dp,shrinkage=0.08_dp,return_vectors=.true.)
    call require(r2%status == 0 .and. all(r2%valid),'reversed rectangular groups')
    call require(all(shape(r2%coefficients_a)==[4,2,2]),'reversed A shape')
    call require(all(shape(r2%coefficients_b)==[4,3,2]),'reversed B shape')
    call require(maxval(abs(r%squared_modes-r2%squared_modes)) < 1.0e-12_dp, &
        'swapped strengths')
    r = canonical_coherence_matrix(s,[1,3],[2,4],conditioned_on=[5,6], &
        cv_shrinkage=.true.,return_vectors=.true.)
    call require(r%status == 0 .and. all(r%valid),'CV vectors')
    r = canonical_coherence_matrix(s,[1,3],[2,4],conditioned_on=[5,6], &
        auto_shrinkage=.true.,return_vectors=.true.)
    call require(r%status == 0 .and. all(r%valid),'automatic shrinkage vectors')
    r = canonical_coherence_matrix(s,[1,1],[2],return_vectors=.true.)
    call require(r%status == partial_coherence_invalid_input,'vector invalid input')
    call require(.not. allocated(r%coefficients_a),'invalid input unallocated')
    s_zero = s
    s_zero%power = (0.0_dp,0.0_dp)
    r = canonical_coherence_matrix(s_zero,[1],[2],return_vectors=.true.)
    call require(r%status == partial_coherence_singular .and. .not. any(r%valid), &
        'singular vector cells')
    call require(all(ieee_is_nan(real(r%coefficients_a,dp))), 'singular A NaNs')
    call require(all(ieee_is_nan(real(r%coefficients_b,dp))), 'singular B NaNs')

    ! Identity covariance: each canonical mode is zero; directions are not unique.
    s_zero%power = (0.0_dp,0.0_dp)
    do i = 1, 6
        s_zero%power(:,i,i) = (1.0_dp,0.0_dp)
    end do
    r = canonical_coherence_matrix(s_zero,[1,2],[3,4],return_vectors=.true.)
    call require(r%status == 0 .and. .not. any(r%directions_stable), &
        'zero repeated modes have unstable orientations')
    call require(maxval(abs(r%squared_modes)) <= 1.0e-14_dp,'zero modes')
    call require(maxval(abs(r%coefficients_a)) <= 1.0_dp+1.0e-12_dp,'zero mode coefficients')

    ! Repeated nonzero roots yield a valid normalized basis, but individual
    ! directions are not identifiable. Mark those directions as unstable.
    s_zero%power(:,1,3) = cmplx(0.6_dp,0.0_dp,dp)
    s_zero%power(:,3,1) = cmplx(0.6_dp,0.0_dp,dp)
    s_zero%power(:,2,4) = cmplx(0.6_dp,0.0_dp,dp)
    s_zero%power(:,4,2) = cmplx(0.6_dp,0.0_dp,dp)
    r = canonical_coherence_matrix(s_zero,[1,2],[3,4],return_vectors=.true.)
    call require(r%status == 0 .and. .not. any(r%directions_stable), &
        'repeated nonzero roots flag unstable directions')
    call require(maxval(abs(r%squared_modes-0.36_dp)) < 1.0e-12_dp, &
        'repeated nonzero root values')
    s_zero%power(:,2,4) = cmplx(0.6000001_dp,0.0_dp,dp)
    s_zero%power(:,4,2) = cmplx(0.6000001_dp,0.0_dp,dp)
    r = canonical_coherence_matrix(s_zero,[1,2],[3,4],return_vectors=.true.)
    call require(r%status == 0 .and. .not. any(r%directions_stable), &
        'near-tied roots flagged')
    ! One positive and one null root exercise the arbitrary null-space basis.
    s_zero%power(:,2,4) = (0.0_dp,0.0_dp)
    s_zero%power(:,4,2) = (0.0_dp,0.0_dp)
    r = canonical_coherence_matrix(s_zero,[1,2],[3,4],return_vectors=.true.)
    call require(r%status == 0 .and. all(r%directions_stable(:,1)) .and. &
        .not. any(r%directions_stable(:,2)), 'mixed isolated and null roots')
    call require(maxval(abs(r%squared_modes(:,2))) < 1.0e-12_dp, &
        'zero root alongside positive root')

    allocate(st%frequencies(4),st%times(2),st%power(4,2,6,6))
    st%status = 0
    st%frequencies = s%frequencies
    st%times = [0.0_dp,1.0_dp]
    do t = 1, 2
        st%power(:,t,:,:) = s%power*real(t,dp)
    end do
    tf = canonical_coherence_tf_matrix(st,[1,3,5],[2,4],conditioned_on=[6], &
        ridge=0.03_dp,shrinkage=0.08_dp,return_vectors=.true.)
    call require(tf%status == 0 .and. all(tf%valid),'TF vectors')
    call require(all(shape(tf%coefficients_a)==[4,2,3,2]),'TF A shape')
    call require(all(shape(tf%coefficients_b)==[4,2,2,2]),'TF B shape')
    bands(:,1) = [0.05_dp,0.35_dp]
    bands(:,2) = [0.25_dp,0.45_dp]
    band = band_canonical_coherence_matrix(s,bands,[1,3],[2,4], &
        return_vectors=.true.,ridge=0.01_dp)
    tf_band = band_canonical_coherence_tf_matrix(st,bands,[1,3],[2,4], &
        return_vectors=.true.,ridge=0.01_dp)
    call require(band%status == 0 .and. tf_band%status == 0,'band vectors')
    call require(all(shape(tf_band%coefficients_a)==[2,2,2,2]),'TF band coefficient shape')
    call require(maxval(abs(band%squared_modes-tf_band%squared_modes(:,1,:))) &
        < 1.0e-11_dp,'band/TF consistency')
    do i = 1, 128
        do j = 1, 5
            x(i,j) = sin(0.03_dp*i*j)+0.2_dp*cos(0.12_dp*i+0.2_dp*j)
            z(i,j) = cmplx(x(i,j),0.2_dp*sin(0.2_dp*i*j),dp)
        end do
    end do
    r = multitaper_canonical_coherence(x,3.0_dp,[1,2],[3,4],kmax=5, &
        ridge=0.05_dp,return_vectors=.true.)
    call require(r%status == 0 .and. all(r%valid),'direct real vectors')
    r = multitaper_canonical_coherence(z,3.0_dp,[1,2],[3,4],kmax=5, &
        ridge=0.05_dp,return_vectors=.true.)
    call require(r%status == 0 .and. all(r%valid),'direct complex vectors')
    tf = multitaper_tf_canonical_coherence(x,3.0_dp,64,[1,2],[3,4],kmax=5, &
        ridge=0.05_dp,return_vectors=.true.)
    call require(tf%status == 0 .and. all(tf%valid),'direct TF real vectors')
    tf = multitaper_tf_canonical_coherence(z,3.0_dp,64,[1,2],[3,4],kmax=5, &
        ridge=0.05_dp,return_vectors=.true.)
    call require(tf%status == 0 .and. all(tf%valid),'direct TF complex vectors')
    band = multitaper_band_canonical_coherence(x,3.0_dp,bands,[1,2],[3,4], &
        kmax=5,ridge=0.05_dp,return_vectors=.true.)
    call require(band%status == 0 .and. all(band%valid),'direct band real vectors')
    band = multitaper_band_canonical_coherence(z,3.0_dp,bands,[1,2],[3,4], &
        kmax=5,ridge=0.05_dp,return_vectors=.true.)
    call require(band%status == 0 .and. all(band%valid),'direct band complex vectors')
    tf_band = multitaper_tf_band_canonical_coherence(x,3.0_dp,64,bands,[1,2],[3,4], &
        kmax=5,ridge=0.05_dp,return_vectors=.true.)
    call require(tf_band%status == 0 .and. all(tf_band%valid),'direct TF band real vectors')
    tf_band = multitaper_tf_band_canonical_coherence(z,3.0_dp,64,bands,[1,2],[3,4], &
        kmax=5,ridge=0.05_dp,return_vectors=.true.)
    call require(tf_band%status == 0 .and. all(tf_band%valid),'direct TF band complex vectors')
    print *, 'PASS canonical coefficient vectors'

contains
    pure integer function idx_a(i) result(idx)
        integer, intent(in) :: i !! one-based channel index into A
        integer, parameter :: a(3) = [1,3,5]
        idx = a(i)
    end function idx_a
    pure integer function idx_b(i) result(idx)
        integer, intent(in) :: i !! one-based channel index into B
        integer, parameter :: b(2) = [2,4]
        idx = b(i)
    end function idx_b
    subroutine require(cond,label)
        logical, intent(in) :: cond !! passing condition
        character(*), intent(in) :: label !! diagnostic label
        if (cond) return
        print *, 'FAIL:',label
        error stop 1
    end subroutine require
    subroutine close_matrix(actual,expected,tolerance,label)
        complex(dp), intent(in) :: actual(:,:) !! computed matrix
        complex(dp), intent(in) :: expected(:,:) !! target matrix
        real(dp), intent(in) :: tolerance !! absolute discrepancy allowed
        character(*), intent(in) :: label !! descriptive label
        if (maxval(abs(actual-expected)) <= tolerance) return
        print *, 'FAIL:',label,maxval(abs(actual-expected))
        error stop 1
    end subroutine close_matrix
end program test_signal_canonical_vectors
