! SPDX-License-Identifier: MIT
program test_signal_group_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_group_coherence, only : group_coherence_result, &
        group_coherence_tf_result, group_coherence_matrix, group_coherence_tf_matrix, &
        band_group_coherence_matrix, band_group_coherence_tf_matrix, &
        multitaper_group_coherence, multitaper_tf_group_coherence, &
        multitaper_band_group_coherence, multitaper_tf_band_group_coherence
    use scifort_signal_multiple_coherence, only : multiple_coherence_matrix, multiple_coherence_result
    use scifort_signal_partial_coherence, only : partial_coherence_matrix, partial_coherence_result, &
        partial_coherence_invalid_input, partial_coherence_singular
    implicit none

    type(multitaper_matrix_result) :: s, v, integrated
    type(multitaper_spectral_matrix_spectrogram_result) :: st
    type(group_coherence_result) :: result, counterpart, a, b
    type(group_coherence_tf_result) :: time_result, bt
    type(partial_coherence_result) :: partial
    type(multiple_coherence_result) :: multiple
    real(dp) :: bands(2,2), x(96,4), angle
    complex(dp) :: matrix(4,4), basis(4,4), z(96,4)
    integer :: i, j, k, t

    basis = (0.0_dp,0.0_dp)
    do i = 1, 4
        basis(i,i) = cmplx(1.3_dp+0.1_dp*i,0.0_dp,dp)
    end do
    basis(1,2) = cmplx(0.4_dp,0.1_dp,dp)
    basis(1,3) = cmplx(0.2_dp,-0.5_dp,dp)
    basis(2,3) = cmplx(0.7_dp,0.2_dp,dp)
    basis(2,4) = cmplx(-0.3_dp,0.2_dp,dp)
    basis(3,4) = cmplx(0.45_dp,-0.15_dp,dp)
    matrix = matmul(basis,transpose(conjg(basis)))
    allocate(s%frequencies(4),s%power(4,4,4),s%cv_shrinkage(4))
    s%frequencies = [-1.6_dp,-0.25_dp,0.7_dp,1.7_dp]
    s%status = 0
    s%effective_tapers = 5.0_dp
    s%cv_shrinkage = [0.0_dp,0.15_dp,0.45_dp,0.8_dp]
    do k = 1, 4
        s%power(k,:,:) = matrix*(1.0_dp+0.06_dp*k)
    end do
    result = group_coherence_matrix(s,[1],[2],conditioned_on=[3,4],ridge=0.05_dp)
    partial = partial_coherence_matrix(s,ridge=0.05_dp)
    call require(result%status == 0,'group coherence success')
    call near(maxval(abs(result%coherence-partial%coherence(:,1,2))),0.0_dp,1.0e-13_dp, &
        'partial coherence singleton/group identity')
    result = group_coherence_matrix(s,[1],[2,3,4],ridge=0.02_dp)
    multiple = multiple_coherence_matrix(s,ridge=0.02_dp)
    call near(maxval(abs(result%coherence-multiple%coherence(:,1))),0.0_dp,1.0e-13_dp, &
        'multiple coherence singleton/group identity')
    result = group_coherence_matrix(s,[1,3],[2],conditioned_on=[4],ridge=0.03_dp,shrinkage=0.2_dp)
    counterpart = group_coherence_matrix(s,[2],[3,1],conditioned_on=[4],ridge=0.03_dp,shrinkage=0.2_dp)
    call require(result%status == 0 .and. counterpart%status == 0,'symmetric block coherence')
    call near(maxval(abs(result%coherence-counterpart%coherence)),0.0_dp,1.0e-13_dp, &
        'invariance to group labels and channel ordering')
    call require(all(result%log_information >= 0.0_dp),'information is nonnegative')
    call near(maxval(abs(result%coherence-(1.0_dp-exp(-result%log_information)))), &
        0.0_dp,5.0e-14_dp,'log information transform')
    result = group_coherence_matrix(s,[1,3],[2],conditioned_on=[4],auto_shrinkage=.true.)
    call require(result%status == 0,'automatic diagonal shrinkage')
    call require(all(result%shrinkage_used >= 0.0_dp) .and. &
        all(result%shrinkage_used <= 1.0_dp),'automatic shrinkage bounds')
    result = group_coherence_matrix(s,[1,3],[2],conditioned_on=[4],cv_shrinkage=.true.)
    call require(result%status == 0,'CV shrinkage path')
    call near(maxval(abs(result%shrinkage_used-s%cv_shrinkage)),0.0_dp,0.0_dp, &
        'uses stored CV coefficients')

    allocate(st%frequencies(4),st%times(3),st%power(4,3,4,4),st%cv_shrinkage(4,3))
    st%frequencies = s%frequencies
    st%times = [0.1_dp,0.2_dp,0.3_dp]
    st%status = 0
    st%effective_tapers = 5.0_dp
    do t = 1, 3
        st%power(:,t,:,:) = s%power*real(t,dp)
        st%cv_shrinkage(:,t) = s%cv_shrinkage/real(t,dp)
    end do
    time_result = group_coherence_tf_matrix(st,[1,3],[2],conditioned_on=[4], &
        ridge=0.03_dp,cv_shrinkage=.true.)
    call require(time_result%status == 0,'time-frequency group coherence')
    do t = 1, 3
        v = s
        v%power = st%power(:,t,:,:)
        v%cv_shrinkage = st%cv_shrinkage(:,t)
        result = group_coherence_matrix(v,[1,3],[2],conditioned_on=[4], &
            ridge=0.03_dp,cv_shrinkage=.true.)
        call near(maxval(abs(time_result%coherence(:,t)-result%coherence)), &
            0.0_dp,5.0e-14_dp,'TF vs static')
    end do

    bands(:,1) = [-1.2_dp,0.4_dp]
    bands(:,2) = [0.2_dp,1.5_dp]
    a = band_group_coherence_matrix(s,bands,[1,3],[2],conditioned_on=[4], &
        shrinkage=0.15_dp,ridge=0.01_dp)
    integrated = integrate_reference(s,bands)
    b = group_coherence_matrix(integrated,[1,3],[2],conditioned_on=[4], &
        shrinkage=0.15_dp,ridge=0.01_dp)
    call require(a%status == 0,'band-integrated group coherence')
    call near(maxval(abs(a%coherence-b%coherence)),0.0_dp,5.0e-14_dp, &
        'band group vs integrated CSD coherence')
    bt = band_group_coherence_tf_matrix(st,bands,[1,3],[2],conditioned_on=[4], &
        shrinkage=0.15_dp,ridge=0.01_dp)
    call require(bt%status == 0 .and. all(bt%valid),'TF band group coherence')
    do t = 1, 3
        call near(maxval(abs(bt%coherence(:,t)-a%coherence)),0.0_dp,1.0e-13_dp, &
            'scale invariance of band group coherence')
    end do

    result = group_coherence_matrix(s,[1,1],[2])
    call require(result%status == partial_coherence_invalid_input,'reject duplicate group indices')
    result = group_coherence_matrix(s,[1],[2],conditioned_on=[1])
    call require(result%status == partial_coherence_invalid_input,'reject overlapping indices')
    result = group_coherence_matrix(s,[1],[5])
    call require(result%status == partial_coherence_invalid_input,'reject out of range')
    result = group_coherence_matrix(s,[integer::],[2])
    call require(result%status == partial_coherence_invalid_input,'reject empty group')
    result = group_coherence_matrix(s,[1],[2],shrinkage=-0.01_dp)
    call require(result%status == partial_coherence_invalid_input,'reject negative shrinkage')
    result = group_coherence_matrix(s,[1],[2],auto_shrinkage=.true.,cv_shrinkage=.true.)
    call require(result%status == partial_coherence_invalid_input,'reject conflicting shrinkage')
    result = band_group_coherence_matrix(s,bands,[1],[2],conditioned_on=[3], &
        shrinkage=0.2_dp)
    call require(result%status == 0,'band routine handles conditioned groups')
    v = s
    v%power(1,1,2) = cmplx(0.6_dp,0.7_dp,dp)
    result = group_coherence_matrix(v,[1],[2])
    call require(result%status == partial_coherence_invalid_input,'reject non-Hermitian spectrum')
    v = s
    v%power(2,:,:) = (0.0_dp,0.0_dp)
    result = group_coherence_matrix(v,[1],[2])
    call require(result%status == partial_coherence_singular,'report singular spectral cell')
    call require(.not. result%valid(2) .and. result%valid(1),'per-cell validity masks')

    ! Direct real/complex entry points retain the multitaper options.
    do k = 1, 96
        angle = 0.09_dp*real(k,dp)
        do j = 1, 4
            x(k,j) = sin(angle*j)+0.3_dp*cos(0.21_dp*k-0.7_dp*j)
            z(k,j) = cmplx(x(k,j),0.2_dp*cos(0.17_dp*k+0.3_dp*j),dp)
        end do
    end do
    result = multitaper_group_coherence(x,3.0_dp,[1],[2,3],conditioned_on=[4], &
        kmax=6,nfft=128,ridge=0.02_dp,weighting='eigen',auto_shrinkage=.true.)
    call require(result%status == 0 .and. all(result%valid),'real DPSS auto group coherence')
    result = multitaper_group_coherence(z,3.0_dp,[1],[2,3],conditioned_on=[4], &
        kmax=6,nfft=128,ridge=0.02_dp,weighting='adaptive',cv_shrinkage=.true.)
    call require(result%status == 0 .and. all(result%valid),'complex adaptive CV group coherence')
    time_result = multitaper_tf_group_coherence(z,3.0_dp,64,[1],[2,3], &
        conditioned_on=[4],noverlap=32,kmax=6,ridge=0.03_dp,shrinkage=0.2_dp)
    call require(time_result%status == 0 .and. size(time_result%times)>1, &
        'complex DPSS TF group coherence')
    bands(:,1) = [0.05_dp,0.2_dp]
    bands(:,2) = [0.2_dp,0.4_dp]
    result = multitaper_band_group_coherence(x,3.0_dp,bands,[1],[2,3], &
        conditioned_on=[4],kmax=6,ridge=0.02_dp,shrinkage=0.1_dp)
    call require(result%status == 0 .and. all(result%valid),'real DPSS band group coherence')
    bt = multitaper_tf_band_group_coherence(z,3.0_dp,64,bands,[1],[2,3], &
        conditioned_on=[4],kmax=6,ridge=0.03_dp,shrinkage=0.1_dp)
    call require(bt%status == 0 .and. all(bt%valid),'complex DPSS band TF group coherence')
    print *, 'PASS group coherence native regressions'
contains
    function integrate_reference(source,limits) result(out)
        use scifort_signal_band_coherence, only : integrate_spectral_bands
        type(multitaper_matrix_result), intent(in) :: source !! input cross spectra
        real(dp), intent(in) :: limits(:,:) !! integration bands
        type(multitaper_matrix_result) :: out
        out = integrate_spectral_bands(source,limits)
    end function integrate_reference

    subroutine require(ok,message)
        logical, intent(in) :: ok !! assertion success
        character(*), intent(in) :: message !! descriptive assertion label
        if (.not. ok) then
            print *, 'FAIL: ',message
            error stop 1
        end if
    end subroutine require

    subroutine near(actual,expected,tolerance,message)
        real(dp), intent(in) :: actual !! computed difference
        real(dp), intent(in) :: expected !! reference zero
        real(dp), intent(in) :: tolerance !! maximum absolute error
        character(*), intent(in) :: message !! assertion description
        if (.not. (abs(actual-expected)<=tolerance)) then
            print *, 'FAIL: ',message, actual, expected, tolerance
            error stop 1
        end if
    end subroutine near
end program test_signal_group_coherence
