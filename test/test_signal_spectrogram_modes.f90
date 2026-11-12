! SPDX-License-Identifier: MIT
! Native API, invalid options and phase-unwrapping tests.
program test_signal_spectrogram_modes
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectrogram_real => spectrogram, &
        spectrogram_result, spectral_status_success, spectral_status_invalid_input, &
        spectral_unwrap_phase
    use scifort_signal_spectral_complex, only : spectrogram_complex
    implicit none

    integer :: i
    real(dp) :: x(33), phase(3,2), pi
    complex(dp) :: z(33)
    type(spectrogram_result) :: a, b, c

    pi = acos(-1.0_dp)
    do i = 1, size(x)
        x(i) = sin(0.32_dp*real(i,dp)) + cos(0.91_dp*real(i,dp))
        z(i) = cmplx(x(i),0.7_dp*cos(0.37_dp*real(i,dp)),dp)
    end do
    a = spectrogram_real(x,nperseg=9,noverlap=2)
    b = spectrogram_real(x,nperseg=9,noverlap=2,mode='psd')
    if (a%status /= spectral_status_success) error stop 'PSD regression: status'
    if (maxval(abs(a%power-b%power)) > 1.0e-13_dp) error stop 'PSD regression: output'
    if (allocated(a%values)) error stop 'PSD unexpectedly has complex output'

    c = spectrogram_real(x,nperseg=9,noverlap=2,mode='complex')
    if (c%status /= spectral_status_success) error stop 'complex mode: status'
    if (.not. allocated(c%values)) error stop 'complex mode: missing coefficients'
    if (size(c%values,1) /= 5) error stop 'real one-sided frequency count'
    if (size(c%power) /= 0) error stop 'complex mode returned real power array'
    b = spectrogram_real(x,nperseg=9,noverlap=2,mode='magnitude')
    if (maxval(abs(b%power - abs(c%values))) > 1.0e-13_dp) error stop 'magnitude mismatch'

    c = spectrogram_complex(z,nperseg=9,noverlap=2,mode='complex')
    if (c%status /= spectral_status_success) error stop 'complex input: status'
    if (size(c%values,1) /= 9) error stop 'complex input: frequency bins'
    b = spectrogram_complex(z,nperseg=9,noverlap=2,mode='angle')
    if (maxval(abs(b%power - atan2(aimag(c%values),real(c%values,dp)))) > 1.0e-13_dp) &
        error stop 'angle output mismatch'
    b = spectrogram_complex(z,nperseg=9,noverlap=2,mode='nonsense')
    if (b%status /= spectral_status_invalid_input) error stop 'invalid mode accepted'
    if (size(b%power) /= 0) error stop 'invalid mode: nonempty output'
    b = spectrogram_real(x,mode='nonsense')
    if (b%status /= spectral_status_invalid_input) error stop 'real invalid mode accepted'

    phase(:,1) = [0.9_dp*pi,-0.9_dp*pi,-0.8_dp*pi]
    phase(:,2) = [-0.9_dp*pi,0.9_dp*pi,0.8_dp*pi]
    call spectral_unwrap_phase(phase)
    if (maxval(abs(phase(:,1) - [0.9_dp*pi,1.1_dp*pi,1.2_dp*pi])) > 1.0e-12_dp) &
        error stop 'positive phase unwrap'
    if (maxval(abs(phase(:,2) - [-0.9_dp*pi,-1.1_dp*pi,-1.2_dp*pi])) > 1.0e-12_dp) &
        error stop 'negative phase unwrap'
    print *, 'PASS: native legacy spectrogram modes'
end program test_signal_spectrogram_modes
