! SPDX-License-Identifier: MIT
program test_signal_partial_shrinkage_cv
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_multitaper, only : multitaper_matrix_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix_spectrogram
    use scifort_signal_partial_coherence, only : partial_coherence_matrix, &
        partial_coherence_tf_matrix, multitaper_partial_coherence, multitaper_tf_partial_coherence, &
        partial_coherence_result, partial_coherence_tf_result, &
        partial_coherence_success, partial_coherence_invalid_input
    implicit none

    type(multitaper_matrix_result) :: spectral, regular
    type(multitaper_spectral_matrix_spectrogram_result) :: tf, regular_tf
    type(partial_coherence_result) :: result1, result2
    type(partial_coherence_tf_result) :: resulttf1, resulttf2
    real(dp) :: x(96,4), time
    complex(dp) :: z(96,4)
    character(8) :: mode
    integer :: i, ch, m

    do ch = 1, 4
        do i = 1, 96
            time = real(i,dp)
            x(i,ch) = sin(0.08_dp*time+0.5_dp*ch) + 0.3_dp*cos(time*(0.11_dp+0.015_dp*ch)) + &
                0.18_dp*sin(time*(0.33_dp+0.004_dp*ch*ch))
            z(i,ch) = cmplx(x(i,ch),0.26_dp*cos(time*(0.09_dp+0.014_dp*ch)) + &
                0.17_dp*sin(time*0.23_dp+real(ch,dp)),dp)
        end do
    end do

    do m = 1, 3
        select case (m)
        case (1)
            mode = 'equal'
        case (2)
            mode = 'eigen'
        case (3)
            mode = 'adaptive'
        end select
        regular = multitaper_spectral_matrix(x,3.5_dp,kmax=6,nfft=128,weighting=trim(mode))
        spectral = multitaper_spectral_matrix(x,3.5_dp,kmax=6,nfft=128,weighting=trim(mode), &
            cv_shrinkage=.true.)
        if (spectral%status /= spectral_status_success) error stop 'real spectral status'
        if (.not. allocated(spectral%cv_shrinkage)) error stop 'missing cv coefficients'
        if (size(spectral%cv_shrinkage) /= 65) error stop 'cv shape'
        if (minval(spectral%cv_shrinkage) < 0.0_dp .or. &
            maxval(spectral%cv_shrinkage) > 1.0_dp) error stop 'coefficient range'
        if (maxval(abs(spectral%power-regular%power)) > 3.0e-11_dp) error stop 'real spectrum changed'

        result1 = multitaper_partial_coherence(x,3.5_dp,kmax=6,nfft=128,weighting=trim(mode), &
            cv_shrinkage=.true.,ridge=0.025_dp,return_coherency=.true.)
        result2 = partial_coherence_matrix(spectral,cv_shrinkage=.true.,ridge=0.025_dp, &
            return_coherency=.true.)
        call check_static(result1,result2)

        regular = multitaper_spectral_matrix(transpose(z),3.5_dp,kmax=6,nfft=129,axis=2, &
            weighting=trim(mode))
        spectral = multitaper_spectral_matrix(transpose(z),3.5_dp,kmax=6,nfft=129,axis=2, &
            weighting=trim(mode),cv_shrinkage=.true.)
        if (maxval(abs(spectral%power-regular%power)) > 3.0e-11_dp) error stop 'complex spectrum changed'
        result1 = multitaper_partial_coherence(transpose(z),3.5_dp,kmax=6,nfft=129,axis=2, &
            weighting=trim(mode),cv_shrinkage=.true.,ridge=0.025_dp,return_coherency=.true.)
        result2 = partial_coherence_matrix(spectral,cv_shrinkage=.true.,ridge=0.025_dp, &
            return_coherency=.true.)
        call check_static(result1,result2)

        regular_tf = multitaper_spectral_matrix_spectrogram(x,2.5_dp,32,noverlap=16, &
            kmax=6,nfft=48,weighting=trim(mode))
        tf = multitaper_spectral_matrix_spectrogram(x,2.5_dp,32,noverlap=16, &
            kmax=6,nfft=48,weighting=trim(mode),cv_shrinkage=.true.)
        if (tf%status /= spectral_status_success) error stop 'real TF spectral status'
        if (size(tf%cv_shrinkage,1) /= 25 .or. size(tf%cv_shrinkage,2) /= 5) &
            error stop 'TF coefficient dimensions'
        if (minval(tf%cv_shrinkage) < 0.0_dp .or. &
            maxval(tf%cv_shrinkage) > 1.0_dp) error stop 'TF coefficient range'
        if (maxval(abs(tf%power-regular_tf%power)) > 1.0e-12_dp) error stop 'real TF spectrum changed'
        resulttf1 = multitaper_tf_partial_coherence(x,2.5_dp,32,noverlap=16,kmax=6, &
            nfft=48,weighting=trim(mode),cv_shrinkage=.true.,ridge=0.025_dp,return_coherency=.true.)
        resulttf2 = partial_coherence_tf_matrix(tf,cv_shrinkage=.true.,ridge=0.025_dp, &
            return_coherency=.true.)
        call check_tf(resulttf1,resulttf2)

        regular_tf = multitaper_spectral_matrix_spectrogram(transpose(z),2.5_dp,32,noverlap=16, &
            kmax=6,nfft=48,weighting=trim(mode),axis=2)
        tf = multitaper_spectral_matrix_spectrogram(transpose(z),2.5_dp,32,noverlap=16, &
            kmax=6,nfft=48,weighting=trim(mode),axis=2,cv_shrinkage=.true.)
        if (maxval(abs(tf%power-regular_tf%power)) > 1.0e-12_dp) error stop 'complex TF spectrum changed'
        resulttf1 = multitaper_tf_partial_coherence(transpose(z),2.5_dp,32,noverlap=16, &
            kmax=6,nfft=48,weighting=trim(mode),axis=2,cv_shrinkage=.true.,ridge=0.025_dp, &
            return_coherency=.true.)
        resulttf2 = partial_coherence_tf_matrix(tf,cv_shrinkage=.true.,ridge=0.025_dp, &
            return_coherency=.true.)
        call check_tf(resulttf1,resulttf2)
    end do

    spectral = multitaper_spectral_matrix(x,3.5_dp,kmax=6)
    result1 = partial_coherence_matrix(spectral,cv_shrinkage=.true.)
    if (result1%status /= partial_coherence_invalid_input) error stop 'missing CV matrix accepted'
    tf = multitaper_spectral_matrix_spectrogram(x,2.5_dp,32,kmax=6)
    resulttf1 = partial_coherence_tf_matrix(tf,cv_shrinkage=.true.)
    if (resulttf1%status /= partial_coherence_invalid_input) error stop 'missing CV TF accepted'
    spectral = multitaper_spectral_matrix(x,3.5_dp,kmax=6,cv_shrinkage=.true.)
    result1 = partial_coherence_matrix(spectral,cv_shrinkage=.true.,shrinkage=0.2_dp)
    if (result1%status /= partial_coherence_invalid_input) error stop 'CV and manual accepted'
    result1 = partial_coherence_matrix(spectral,cv_shrinkage=.true.,auto_shrinkage=.true.,effective_tapers=6.0_dp)
    if (result1%status /= partial_coherence_invalid_input) error stop 'CV and Gaussian accepted'
    result1 = multitaper_partial_coherence(x,3.5_dp,kmax=6,cv_shrinkage=.true.,shrinkage=0.2_dp)
    if (result1%status /= partial_coherence_invalid_input) error stop 'signal CV and manual accepted'
    result1 = multitaper_partial_coherence(x,3.5_dp,kmax=6,cv_shrinkage=.true.,auto_shrinkage=.true.)
    if (result1%status /= partial_coherence_invalid_input) error stop 'signal CV and Gaussian accepted'
    spectral = multitaper_spectral_matrix(x,3.5_dp,kmax=2,cv_shrinkage=.true.)
    if (spectral%status == spectral_status_success) error stop 'accepted 2-taper CV'
    spectral = multitaper_spectral_matrix(x(:,1:1),3.5_dp,kmax=6,cv_shrinkage=.true.)
    if (spectral%status == spectral_status_success) error stop 'accepted single-channel CV'

    spectral = multitaper_spectral_matrix(x,3.5_dp,kmax=6,cv_shrinkage=.true.)
    spectral%cv_shrinkage(1) = ieee_value(0.0_dp,ieee_quiet_nan)
    result1 = partial_coherence_matrix(spectral,cv_shrinkage=.true.)
    if (result1%status /= partial_coherence_invalid_input) error stop 'accepted NaN CV'
    spectral%cv_shrinkage(1) = -0.01_dp
    result1 = partial_coherence_matrix(spectral,cv_shrinkage=.true.)
    if (result1%status /= partial_coherence_invalid_input) error stop 'accepted negative CV'
    ! Corrupted metadata and incompatible time-frequency modes must fail.
    spectral = multitaper_spectral_matrix(x,3.5_dp,kmax=6,cv_shrinkage=.true.)
    spectral%cv_shrinkage(1) = 1.01_dp
    result1 = partial_coherence_matrix(spectral,cv_shrinkage=.true.)
    if (result1%status /= partial_coherence_invalid_input) error stop 'accepted CV > 1'
    deallocate(spectral%cv_shrinkage)
    allocate(spectral%cv_shrinkage(2))
    spectral%cv_shrinkage = 0.2_dp
    result1 = partial_coherence_matrix(spectral,cv_shrinkage=.true.)
    if (result1%status /= partial_coherence_invalid_input) error stop 'accepted wrong CV shape'

    tf = multitaper_spectral_matrix_spectrogram(x,2.5_dp,32,kmax=6,cv_shrinkage=.true.)
    resulttf1 = partial_coherence_tf_matrix(tf,cv_shrinkage=.true.,shrinkage=0.2_dp)
    if (resulttf1%status /= partial_coherence_invalid_input) error stop 'TF CV and manual accepted'
    resulttf1 = multitaper_tf_partial_coherence(x,2.5_dp,32,kmax=6,cv_shrinkage=.true., &
        auto_shrinkage=.true.)
    if (resulttf1%status /= partial_coherence_invalid_input) error stop 'TF CV and Gaussian accepted'
    tf%cv_shrinkage(1,1) = 1.5_dp
    resulttf1 = partial_coherence_tf_matrix(tf,cv_shrinkage=.true.)
    if (resulttf1%status /= partial_coherence_invalid_input) error stop 'accepted TF alpha > 1'
    deallocate(tf%cv_shrinkage)
    allocate(tf%cv_shrinkage(1,1))
    tf%cv_shrinkage = 0.5_dp
    resulttf1 = partial_coherence_tf_matrix(tf,cv_shrinkage=.true.)
    if (resulttf1%status /= partial_coherence_invalid_input) error stop 'accepted wrong TF CV shape'

    ! All-zero spectra have undefined inverse and must remain explicitly singular.
    spectral = multitaper_spectral_matrix(0.0_dp*x,3.5_dp,kmax=6,cv_shrinkage=.true.)
    if (spectral%status /= spectral_status_success) error stop 'zero CV spectral status'
    result1 = partial_coherence_matrix(spectral,cv_shrinkage=.true.)
    if (any(result1%valid)) error stop 'zero spectral inverse accepted'
    if (any(result1%shrinkage_used == result1%shrinkage_used)) &
        error stop 'zero spectrum shrinkage should be NaN'

    print '(a)', 'PASS LOO diagonal shrinkage for equal/eigen/adaptive, real/complex, static/TF'

contains
    subroutine check_static(a,b)
        type(partial_coherence_result), intent(in) :: a !! high-level output
        type(partial_coherence_result), intent(in) :: b !! direct output
        if (a%status /= partial_coherence_success .or. b%status /= partial_coherence_success) &
            error stop 'static CV status'
        if (.not. all(a%valid) .or. .not. all(b%valid)) error stop 'static CV singular'
        if (maxval(abs(a%shrinkage_used-b%shrinkage_used)) > 1.0e-13_dp) &
            error stop 'static CV coefficient mismatch'
        if (maxval(abs(a%coherency-b%coherency)) > 1.0e-12_dp) &
            error stop 'static CV coherency mismatch'
        if (maxval(abs(a%precision-b%precision)) > 1.0e-10_dp) &
            error stop 'static CV precision mismatch'
    end subroutine check_static

    subroutine check_tf(a,b)
        type(partial_coherence_tf_result), intent(in) :: a !! high-level output
        type(partial_coherence_tf_result), intent(in) :: b !! direct output
        if (a%status /= partial_coherence_success .or. b%status /= partial_coherence_success) &
            error stop 'TF CV status'
        if (.not. all(a%valid) .or. .not. all(b%valid)) error stop 'TF CV singular'
        if (maxval(abs(a%shrinkage_used-b%shrinkage_used)) > 1.0e-13_dp) &
            error stop 'TF CV coefficient mismatch'
        if (maxval(abs(a%coherency-b%coherency)) > 1.0e-12_dp) &
            error stop 'TF CV coherency mismatch'
        if (maxval(abs(a%precision-b%precision)) > 1.0e-10_dp) &
            error stop 'TF CV precision mismatch'
    end subroutine check_tf
end program test_signal_partial_shrinkage_cv
