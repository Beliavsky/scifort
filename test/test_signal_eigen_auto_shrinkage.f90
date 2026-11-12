! SPDX-License-Identifier: MIT
program test_signal_eigen_auto_shrinkage
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram, multitaper_effective_tapers
    use scifort_signal_partial_coherence, only : partial_coherence_matrix, &
        partial_coherence_tf_matrix, multitaper_partial_coherence, &
        multitaper_tf_partial_coherence, partial_coherence_result, &
        partial_coherence_tf_result, partial_coherence_success, partial_coherence_invalid_input
    implicit none

    type(dpss_result) :: tapers
    type(multitaper_matrix_result) :: m
    type(multitaper_spectral_matrix_spectrogram_result) :: sm
    type(partial_coherence_result) :: measured, direct
    type(partial_coherence_tf_result) :: measured_tf, direct_tf
    real(dp) :: x(96,4), w(3), neff, t
    complex(dp) :: z(96,4)
    integer :: i, ch

    if (multitaper_effective_tapers([1.0_dp,1.0_dp,1.0_dp,1.0_dp]) /= 4.0_dp) &
        error stop 'equal effective count'
    if (abs(multitaper_effective_tapers([1.0_dp,2.0_dp,0.0_dp])-1.8_dp) > 1.0e-14_dp) &
        error stop 'unequal effective count'
    if (abs(multitaper_effective_tapers([huge(1.0_dp)/4,huge(1.0_dp)/2,0.0_dp])-1.8_dp) &
        > 1.0e-14_dp) error stop 'effective count overflow'
    if (multitaper_effective_tapers([0.0_dp,0.0_dp]) /= 0.0_dp) error stop 'zero weights'
    if (multitaper_effective_tapers([1.0_dp,-1.0_dp]) /= 0.0_dp) error stop 'negative weights'
    if (multitaper_effective_tapers([real(dp)::]) /= 0.0_dp) error stop 'empty weights'
    w = [1.0_dp,1.0_dp,ieee_value(0.0_dp,ieee_quiet_nan)]
    if (multitaper_effective_tapers(w) /= 0.0_dp) error stop 'NaN weights'

    do ch = 1,4
        do i = 1,96
            t = real(i,dp)
            x(i,ch) = sin(t*(0.035_dp+0.014_dp*ch)+0.15_dp*ch) + &
                0.35_dp*cos(t*(0.23_dp+0.011_dp*ch)+0.04_dp*ch) + &
                0.13_dp*sin(t*0.58_dp + real(ch*ch,dp))
            z(i,ch) = cmplx(x(i,ch), &
                0.25_dp*sin(t*(0.12_dp+0.03_dp*ch))+0.11_dp*cos(t*0.29_dp*ch),dp)
        end do
    end do

    tapers = dpss(96,2.5_dp,kmax=6,sym=.false.,norm='l2')
    neff = multitaper_effective_tapers(tapers%ratios)
    if (neff <= 1.0_dp .or. neff >= 6.0_dp) error stop 'expected nonuniform eigen weights'
    m = multitaper_spectral_matrix(x,2.5_dp,kmax=6,weighting='eigen',nfft=128)
    if (abs(m%effective_tapers-neff) > 1.0e-13_dp) error stop 'real metadata'
    measured = multitaper_partial_coherence(x,2.5_dp,kmax=6,nfft=128, &
        weighting='eigen',auto_shrinkage=.true.,ridge=0.025_dp,return_coherency=.true.)
    direct = partial_coherence_matrix(m,auto_shrinkage=.true.,effective_tapers=neff, &
        ridge=0.025_dp,return_coherency=.true.)
    call compare_static(measured,direct)

    m = multitaper_spectral_matrix(transpose(z),2.5_dp,kmax=6,axis=2, &
        weighting='eigen',nfft=129)
    measured = multitaper_partial_coherence(transpose(z),2.5_dp,kmax=6,axis=2, &
        nfft=129,weighting='eigen',auto_shrinkage=.true.,ridge=0.025_dp, &
        return_coherency=.true.)
    direct = partial_coherence_matrix(m,auto_shrinkage=.true.,effective_tapers=neff, &
        ridge=0.025_dp,return_coherency=.true.)
    call compare_static(measured,direct)

    tapers = dpss(32,2.5_dp,kmax=6,sym=.false.,norm='l2')
    neff = multitaper_effective_tapers(tapers%ratios)
    sm = multitaper_spectral_matrix_spectrogram(x,2.5_dp,32,noverlap=16, &
        kmax=6,weighting='eigen',nfft=48)
    if (abs(sm%effective_tapers-neff) > 1.0e-13_dp) error stop 'TF metadata'
    measured_tf = multitaper_tf_partial_coherence(x,2.5_dp,32,noverlap=16, &
        kmax=6,nfft=48,weighting='eigen',auto_shrinkage=.true., &
        ridge=0.025_dp,return_coherency=.true.)
    direct_tf = partial_coherence_tf_matrix(sm,auto_shrinkage=.true.,effective_tapers=neff, &
        ridge=0.025_dp,return_coherency=.true.)
    call compare_tf(measured_tf,direct_tf)

    sm = multitaper_spectral_matrix_spectrogram(transpose(z),2.5_dp,32, &
        noverlap=16,kmax=6,axis=2,weighting='eigen',nfft=48)
    measured_tf = multitaper_tf_partial_coherence(transpose(z),2.5_dp,32, &
        noverlap=16,kmax=6,axis=2,nfft=48,weighting='eigen', &
        auto_shrinkage=.true.,ridge=0.025_dp,return_coherency=.true.)
    direct_tf = partial_coherence_tf_matrix(sm,auto_shrinkage=.true.,effective_tapers=neff, &
        ridge=0.025_dp,return_coherency=.true.)
    call compare_tf(measured_tf,direct_tf)

    m = multitaper_spectral_matrix(x,2.5_dp,kmax=1,weighting='eigen')
    if (m%effective_tapers /= 1.0_dp) error stop 'single taper count'
    measured = multitaper_partial_coherence(x,2.5_dp,kmax=1,weighting='eigen', &
        auto_shrinkage=.true.)
    if (measured%status /= partial_coherence_success) error stop 'single taper status'
    if (maxval(abs(measured%shrinkage_used-1.0_dp)) /= 0.0_dp) &
        error stop 'single taper shrinkage'
    m = multitaper_spectral_matrix(x,2.5_dp,kmax=6)
    if (m%effective_tapers /= 6.0_dp) error stop 'equal metadata'
    m = multitaper_spectral_matrix(x,2.5_dp,kmax=6,weighting='adaptive')
    if (m%effective_tapers /= 0.0_dp) error stop 'adaptive metadata'
    measured = multitaper_partial_coherence(x,2.5_dp,kmax=6,weighting='adaptive', &
        auto_shrinkage=.true.)
    if (measured%status /= partial_coherence_invalid_input) error stop 'adaptive accepted'
    measured_tf = multitaper_tf_partial_coherence(z,2.5_dp,32,kmax=6, &
        weighting='adaptive',auto_shrinkage=.true.)
    if (measured_tf%status /= partial_coherence_invalid_input) error stop 'TF adaptive accepted'
    m = multitaper_spectral_matrix(x,-1.0_dp,kmax=6,weighting='eigen')
    if (m%effective_tapers /= 0.0_dp) error stop 'invalid metadata'
    print '(a)', 'PASS eigen-weighted auto shrinkage and effective taper count'

contains
    subroutine compare_static(a,b)
        type(partial_coherence_result), intent(in) :: a !! high-level signal estimator
        type(partial_coherence_result), intent(in) :: b !! direct-matrix estimator

        if (a%status /= partial_coherence_success .or. b%status /= partial_coherence_success) &
            error stop 'static eigen status'
        if (.not. all(a%valid) .or. .not. all(b%valid)) error stop 'static eigen singular'
        if (maxval(abs(a%shrinkage_used-b%shrinkage_used)) > 1.0e-14_dp) &
            error stop 'static eigen alpha mismatch'
        if (maxval(abs(a%coherency-b%coherency)) > 1.0e-13_dp) &
            error stop 'static eigen coherency mismatch'
        if (maxval(abs(a%precision-b%precision)) > 1.0e-12_dp) &
            error stop 'static eigen precision mismatch'
    end subroutine compare_static

    subroutine compare_tf(a,b)
        type(partial_coherence_tf_result), intent(in) :: a !! high-level TF estimator
        type(partial_coherence_tf_result), intent(in) :: b !! direct-matrix TF estimator

        if (a%status /= partial_coherence_success .or. b%status /= partial_coherence_success) &
            error stop 'TF eigen status'
        if (.not. all(a%valid) .or. .not. all(b%valid)) error stop 'TF eigen singular'
        if (maxval(abs(a%shrinkage_used-b%shrinkage_used)) > 1.0e-14_dp) &
            error stop 'TF eigen alpha mismatch'
        if (maxval(abs(a%coherency-b%coherency)) > 1.0e-13_dp) &
            error stop 'TF eigen coherency mismatch'
        if (maxval(abs(a%precision-b%precision)) > 1.0e-12_dp) &
            error stop 'TF eigen precision mismatch'
    end subroutine compare_tf
end program test_signal_eigen_auto_shrinkage
