! SPDX-License-Identifier: MIT
program test_signal_tf_partial_coherence_public
    use scifort_signal, only : dp, partial_coherence_tf_result, &
        partial_coherence_tf_matrix, multitaper_tf_partial_coherence, &
        partial_coherence_success, multitaper_spectral_matrix_spectrogram, &
        multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(partial_coherence_tf_result) :: a, b
    type(multitaper_spectral_matrix_spectrogram_result) :: matrix
    real(dp) :: x(64,3)
    complex(dp) :: y(64,3)
    integer :: i, j

    do j=1,3
        do i=1,64
            x(i,j)=sin(real(i*j,dp)*0.21_dp)+cos(real(i+j,dp)*0.43_dp)
            y(i,j)=cmplx(x(i,j),sin(real(i+j,dp)*0.12_dp),dp)
        end do
    end do
    a = multitaper_tf_partial_coherence(x,2.5_dp,32,kmax=5,ridge=0.05_dp)
    if (a%status /= partial_coherence_success) error stop 'real public interface'
    matrix = multitaper_spectral_matrix_spectrogram(x,2.5_dp,32,kmax=5)
    b = partial_coherence_tf_matrix(matrix,ridge=0.05_dp)
    if (maxval(abs(a%coherence-b%coherence)) > 1.0e-14_dp) error stop 'public matrix interface'
    a = multitaper_tf_partial_coherence(y,2.5_dp,32,kmax=5,ridge=0.05_dp)
    if (a%status /= partial_coherence_success) error stop 'complex public interface'
    print '(a)', 'PASS public multitaper time-frequency partial coherence'
end program test_signal_tf_partial_coherence_public
