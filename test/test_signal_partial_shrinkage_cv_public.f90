! SPDX-License-Identifier: MIT
program test_signal_partial_shrinkage_cv_public
    use scifort_signal, only : dp, multitaper_spectral_matrix, multitaper_spectral_matrix_spectrogram, &
        multitaper_partial_coherence, multitaper_tf_partial_coherence, &
        partial_coherence_matrix, partial_coherence_tf_matrix, &
        multitaper_matrix_result, multitaper_spectral_matrix_spectrogram_result, &
        partial_coherence_result, partial_coherence_tf_result, partial_coherence_success
    implicit none
    real(dp) :: x(64,3), t
    integer :: i,j
    type(multitaper_matrix_result) :: spectrum
    type(multitaper_spectral_matrix_spectrogram_result) :: moving
    type(partial_coherence_result) :: a,b
    type(partial_coherence_tf_result) :: c,d

    do j = 1,3
        do i = 1,64
            t = real(i,dp)
            x(i,j) = sin(t*(0.07_dp+0.03_dp*j))+0.5_dp*cos(t*0.17_dp+0.2_dp*j)
        end do
    end do
    spectrum = multitaper_spectral_matrix(x,2.5_dp,kmax=5,weighting='adaptive', &
        cv_shrinkage=.true.)
    a = multitaper_partial_coherence(x,2.5_dp,kmax=5,weighting='adaptive', &
        cv_shrinkage=.true.,return_coherency=.true.,ridge=0.02_dp)
    b = partial_coherence_matrix(spectrum,cv_shrinkage=.true.,return_coherency=.true.,ridge=0.02_dp)
    if (a%status /= partial_coherence_success .or. b%status /= partial_coherence_success) &
        error stop 'public static CV status'
    if (maxval(abs(a%shrinkage_used-b%shrinkage_used)) > 1.0e-12_dp) &
        error stop 'public static shrinkage mismatch'
    if (maxval(abs(a%coherency-b%coherency)) > 1.0e-12_dp) &
        error stop 'public static coherency mismatch'

    moving = multitaper_spectral_matrix_spectrogram(x,2.5_dp,32,kmax=5, &
        weighting='adaptive',cv_shrinkage=.true.)
    c = multitaper_tf_partial_coherence(x,2.5_dp,32,kmax=5, &
        weighting='adaptive',cv_shrinkage=.true.,ridge=0.02_dp)
    d = partial_coherence_tf_matrix(moving,cv_shrinkage=.true.,ridge=0.02_dp)
    if (c%status /= partial_coherence_success .or. d%status /= partial_coherence_success) &
        error stop 'public TF CV status'
    if (maxval(abs(c%shrinkage_used-d%shrinkage_used)) > 1.0e-12_dp) &
        error stop 'public TF shrinkage mismatch'
    print '(a)', 'PASS public LOO cross-validated partial coherence'
end program test_signal_partial_shrinkage_cv_public
