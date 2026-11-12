! SPDX-License-Identifier: MIT
! Check the new multiple-coherence surface via the top-level scifort_signal module.
program test_signal_multiple_coherence_public
    use scifort_signal, only : dp, multiple_coherence_matrix, multiple_coherence_tf_matrix, &
        multitaper_multiple_coherence, multitaper_tf_multiple_coherence, &
        multiple_coherence_result, multiple_coherence_tf_result, &
        multitaper_matrix_result, multitaper_spectral_matrix_spectrogram_result, &
        partial_coherence_success
    implicit none
    type(multitaper_matrix_result) :: s
    type(multitaper_spectral_matrix_spectrogram_result) :: st
    type(multiple_coherence_result) :: r
    type(multiple_coherence_tf_result) :: rt
    real(dp) :: x(80,2)
    complex(dp) :: z(80,2)
    integer :: j

    allocate(s%frequencies(1),s%power(1,2,2))
    s%status = 0
    s%frequencies = [0.1_dp]
    s%power(1,:,:) = reshape([cmplx(2.0_dp,0.0_dp,dp), &
        cmplx(0.5_dp,-0.3_dp,dp),cmplx(0.5_dp,0.3_dp,dp), &
        cmplx(1.0_dp,0.0_dp,dp)],[2,2])
    r = multiple_coherence_matrix(s)
    if (r%status /= partial_coherence_success) error stop 'public static matrix'
    if (abs(r%coherence(1,1)-0.17_dp) > 2.0e-13_dp) error stop 'public value'

    allocate(st%frequencies(1),st%times(1),st%power(1,1,2,2))
    st%status = 0
    st%frequencies = s%frequencies
    st%times = [0.0_dp]
    st%power(1,1,:,:) = s%power(1,:,:)
    rt = multiple_coherence_tf_matrix(st)
    if (rt%status /= partial_coherence_success) error stop 'public TF matrix'
    if (abs(rt%coherence(1,1,2)-r%coherence(1,2)) > 2.0e-13_dp) error stop 'TF value'

    do j=1,80
        x(j,1) = cos(0.15_dp*real(j,dp))
        x(j,2) = sin(0.13_dp*real(j,dp))+0.3_dp*x(j,1)
        z(j,1) = cmplx(x(j,1),0.2_dp*x(j,2),dp)
        z(j,2) = cmplx(x(j,2),0.3_dp*x(j,1),dp)
    end do
    r = multitaper_multiple_coherence(x,3.0_dp,kmax=4,ridge=0.04_dp)
    if (r%status /= partial_coherence_success) error stop 'public real direct'
    rt = multitaper_tf_multiple_coherence(z,3.0_dp,40, &
        kmax=4,ridge=0.04_dp,weighting='eigen')
    if (rt%status /= partial_coherence_success) error stop 'public complex TF direct'
    print *, 'PASS public multiple coherence'
end program test_signal_multiple_coherence_public
