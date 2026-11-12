! SPDX-License-Identifier: MIT
program test_signal_partial_coherence_public
    use scifort_signal, only : dp, multitaper_partial_coherence, &
        partial_coherence_result, partial_coherence_success, partial_coherence_matrix, &
        multitaper_spectral_matrix, multitaper_matrix_result
    implicit none
    real(dp) :: x(48,3)
    complex(dp) :: z(48,3)
    type(partial_coherence_result) :: real_result, complex_result, direct
    type(multitaper_matrix_result) :: spectral
    integer :: i

    do i=1,48
        x(i,1) = sin(0.17_dp*real(i,dp))
        x(i,2) = cos(0.32_dp*real(i,dp)) + x(i,1)*0.2_dp
        x(i,3) = sin(0.46_dp*real(i,dp)) + x(i,2)*0.3_dp
    end do
    z = cmplx(x,0.1_dp,dp)
    z(:,2) = z(:,2) + cmplx(0.0_dp,x(:,3),dp)
    real_result = multitaper_partial_coherence(x,3.5_dp,kmax=6,ridge=0.05_dp)
    if (real_result%status /= partial_coherence_success) error stop 'public real'
    if (size(real_result%coherence,2) /= 3) error stop 'public output shape'
    complex_result = multitaper_partial_coherence(z,3.5_dp,kmax=6,ridge=0.05_dp)
    if (complex_result%status /= partial_coherence_success) error stop 'public complex'
    spectral = multitaper_spectral_matrix(z,3.5_dp,kmax=6)
    direct = partial_coherence_matrix(spectral,ridge=0.05_dp)
    if (maxval(abs(complex_result%coherence-direct%coherence)) > 5.0e-13_dp) &
        error stop 'public matrix adapter'
    print '(a)', 'PASS public partial-coherence interface'
end program test_signal_partial_coherence_public
