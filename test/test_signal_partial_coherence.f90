! SPDX-License-Identifier: MIT
program test_signal_partial_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, multitaper_spectral_matrix
    use scifort_signal_partial_coherence, only : partial_coherence_result, partial_coherence_matrix, &
        multitaper_partial_coherence, partial_coherence_success, &
        partial_coherence_singular, partial_coherence_invalid_input
    implicit none
    type(multitaper_matrix_result) :: mat
    type(partial_coherence_result) :: pc, loaded, scaled, from_data
    real(dp) :: data(65,3), t
    complex(dp) :: cdata(65,3)
    integer :: i, j

    allocate(mat%frequencies(3), mat%power(3,3,3))
    mat%status = 0
    mat%frequencies = [0.0_dp,0.25_dp,0.5_dp]
    mat%power = (0.0_dp,0.0_dp)
    do i = 1,3
        mat%power(1,i,i) = (1.0_dp,0.0_dp)
        mat%power(2,i,i) = (1.0_dp,0.0_dp)
        mat%power(3,i,i) = (1.0_dp,0.0_dp)
    end do
    ! Chain graph: X1 - X2 - X3. Conditional 1-3 coherence is zero.
    mat%power(1,1,2) = (0.5_dp,0.0_dp)
    mat%power(1,2,1) = conjg(mat%power(1,1,2))
    mat%power(1,2,3) = (0.5_dp,0.0_dp)
    mat%power(1,3,2) = conjg(mat%power(1,2,3))
    mat%power(1,1,3) = (0.25_dp,0.0_dp)
    mat%power(1,3,1) = conjg(mat%power(1,1,3))
    ! Complex cross-spectrum, positive-definite (via small off-diagonals).
    mat%power(2,1,2) = cmplx(0.1_dp,0.23_dp,dp)
    mat%power(2,2,1) = conjg(mat%power(2,1,2))
    mat%power(2,2,3) = cmplx(0.20_dp,-0.11_dp,dp)
    mat%power(2,3,2) = conjg(mat%power(2,2,3))
    mat%power(2,1,3) = cmplx(0.14_dp,0.04_dp,dp)
    mat%power(2,3,1) = conjg(mat%power(2,1,3))
    ! Singular all-ones rank-one matrix.
    mat%power(3,:,:) = (1.0_dp,0.0_dp)

    pc = partial_coherence_matrix(mat)
    if (pc%status /= partial_coherence_singular) error stop 'expected singular bin'
    if (.not. all(pc%valid(:2))) error stop 'lost valid frequency'
    if (pc%valid(3)) error stop 'singular matrix reported valid'
    if (.not. ieee_is_nan(pc%coherence(3,1,2))) error stop 'missing NaN'
    if (abs(pc%coherence(1,1,3)) > 1.0e-25_dp) error stop 'chain conditional independence'
    if (abs(pc%coherence(1,1,2) - 0.2_dp) > 1.0e-13_dp) error stop 'chain partial coherence'
    if (maxval(abs(pc%coherence(2,:,:)-transpose(pc%coherence(2,:,:)))) > 1.0e-13_dp) &
        error stop 'squared-coherence asymmetry'

    loaded = partial_coherence_matrix(mat, ridge=0.02_dp)
    if (loaded%status /= partial_coherence_success) error stop 'loaded singular matrix'
    if (.not. all(loaded%valid)) error stop 'invalid loaded bin'
    if (any(loaded%coherence < 0.0_dp) .or. any(loaded%coherence > 1.0_dp)) &
        error stop 'coherence outside [0,1]'
    do i = 1,3
        if (maxval(abs(loaded%coherence(:,i,i)-1.0_dp)) > 1.0e-14_dp) error stop 'diagonal'
    end do
    mat%power = 1000.0_dp*mat%power
    scaled = partial_coherence_matrix(mat, ridge=0.02_dp)
    if (scaled%status /= partial_coherence_success) error stop 'scale test failed'
    if (maxval(abs(scaled%coherence-loaded%coherence)) > 1.0e-12_dp) error stop 'non-invariant loading'
    mat%power = mat%power/1000.0_dp
    pc = partial_coherence_matrix(mat, ridge=-1.0_dp)
    if (pc%status /= partial_coherence_invalid_input) error stop 'accepted negative ridge'
    mat%power(1,1,2) = cmplx(0.0_dp,0.5_dp,dp)
    pc = partial_coherence_matrix(mat)
    if (pc%status /= partial_coherence_invalid_input) error stop 'accepted non-Hermitian matrix'

    do i = 1,65
        t = real(i,dp)
        data(i,1) = sin(0.213_dp*t)+0.16_dp*cos(0.547_dp*t)
        data(i,2) = 0.8_dp*data(i,1)+0.6_dp*sin(0.47_dp*t)
        data(i,3) = 0.4_dp*data(i,1)+0.5_dp*data(i,2)+0.36_dp*cos(0.39_dp*t)
    end do
    from_data = multitaper_partial_coherence(data,4.0_dp,kmax=7,nfft=65,ridge=0.05_dp)
    if (from_data%status /= partial_coherence_success) error stop 'real signals failed'
    mat = multitaper_spectral_matrix(data,4.0_dp,kmax=7,nfft=65)
    loaded = partial_coherence_matrix(mat, ridge=0.05_dp)
    if (maxval(abs(from_data%coherence-loaded%coherence)) > 1.0e-13_dp) &
        error stop 'convenience interface discrepancy'
    from_data = multitaper_partial_coherence(transpose(data),4.0_dp, &
        kmax=7,nfft=65,axis=2,ridge=0.05_dp)
    if (maxval(abs(from_data%coherence-loaded%coherence)) > 1.0e-13_dp) &
        error stop 'axis option discrepancy'
    do j=1,3
        do i=1,65
            cdata(i,j) = cmplx(data(i,j),sin(real(i*j,dp)*0.191_dp),dp)
        end do
    end do
    from_data = multitaper_partial_coherence(cdata,4.0_dp,kmax=7,nfft=65,ridge=0.05_dp)
    if (from_data%status /= partial_coherence_success) error stop 'complex signals failed'
    mat = multitaper_spectral_matrix(cdata,4.0_dp,kmax=7,nfft=65)
    loaded = partial_coherence_matrix(mat, ridge=0.05_dp)
    if (maxval(abs(from_data%coherence-loaded%coherence)) > 1.0e-13_dp) &
        error stop 'complex convenience interface discrepancy'
    print '(a)', 'PASS partial coherence native tests'
end program test_signal_partial_coherence
