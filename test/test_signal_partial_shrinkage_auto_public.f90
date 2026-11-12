! SPDX-License-Identifier: MIT
! Verify the top-level scifort_signal re-exports automatic shrinkage keywords.
program test_signal_partial_shrinkage_auto_public
    use scifort_signal, only : dp, multitaper_partial_coherence, &
        multitaper_tf_partial_coherence, partial_coherence_matrix, &
        multitaper_spectral_matrix, multitaper_matrix_result, multitaper_effective_tapers, &
        partial_coherence_result, partial_coherence_tf_result, &
        partial_coherence_success
    implicit none
    real(dp) :: x(64,3), theta
    type(multitaper_matrix_result) :: spectrum
    type(partial_coherence_result) :: a, b
    type(partial_coherence_tf_result) :: t
    integer :: i, ch

    do ch = 1,3
        do i = 1,64
            theta = real(i,dp)
            x(i,ch) = sin(0.1_dp*theta+real(ch,dp))+0.2_dp*cos(0.31_dp*theta*real(ch,dp))
        end do
    end do
    spectrum = multitaper_spectral_matrix(x,2.5_dp,kmax=4)
    a = partial_coherence_matrix(spectrum,auto_shrinkage=.true.,effective_tapers=4.0_dp, &
        ridge=0.03_dp,return_coherency=.true.)
    b = multitaper_partial_coherence(x,2.5_dp,kmax=4,auto_shrinkage=.true., &
        ridge=0.03_dp,return_coherency=.true.)
    if (a%status /= partial_coherence_success .or. b%status /= partial_coherence_success) &
        error stop 'public automatic shrinkage status'
    if (maxval(abs(a%coherency-b%coherency))>1.0e-13_dp) error stop 'public shrinkage mismatch'
    if (maxval(abs(a%shrinkage_used-b%shrinkage_used))>1.0e-13_dp) &
        error stop 'public shrinkage coefficients mismatch'
    t = multitaper_tf_partial_coherence(x,2.5_dp,32,noverlap=16,kmax=4, &
        auto_shrinkage=.true.,ridge=0.03_dp)
    if (t%status /= partial_coherence_success) error stop 'public time-frequency automatic status'
    if (size(t%shrinkage_used,2) /= 3) error stop 'public TF time grid wrong'
    spectrum = multitaper_spectral_matrix(x,2.5_dp,kmax=6,weighting='eigen')
    if (spectrum%effective_tapers <= 1.0_dp .or. spectrum%effective_tapers >= 6.0_dp) &
        error stop 'public eigen effective count'
    if (abs(multitaper_effective_tapers([1.0_dp,2.0_dp,0.0_dp])-1.8_dp) > 1.0e-14_dp) &
        error stop 'public weighted sample count'
    a = partial_coherence_matrix(spectrum,auto_shrinkage=.true., &
        effective_tapers=spectrum%effective_tapers,ridge=0.03_dp)
    b = multitaper_partial_coherence(x,2.5_dp,kmax=6,weighting='eigen', &
        auto_shrinkage=.true.,ridge=0.03_dp)
    if (a%status /= partial_coherence_success .or. b%status /= partial_coherence_success) &
        error stop 'public eigen automatic status'
    if (maxval(abs(a%coherence-b%coherence)) > 1.0e-13_dp) &
        error stop 'public eigen coherence mismatch'
    t = multitaper_tf_partial_coherence(x,2.5_dp,32,noverlap=16,kmax=6, &
        weighting='eigen',auto_shrinkage=.true.,ridge=0.03_dp)
    if (t%status /= partial_coherence_success .or. .not. all(t%valid)) &
        error stop 'public TF eigen status'

    print *, 'PASS public automatic spectral shrinkage API'
end program test_signal_partial_shrinkage_auto_public
