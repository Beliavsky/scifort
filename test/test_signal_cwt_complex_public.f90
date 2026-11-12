! SPDX-License-Identifier: MIT
! Verify complex CWT overloads are exported by the public signal module.
program test_signal_cwt_complex_public
    use scifort_signal, only : dp, cwt_morlet2, cwt_complex_wavelet, &
        morlet2_wavelet, cwt_peak_success, &
        find_peaks_cwt_complex, find_peaks_cwt_magnitude
    implicit none
    real(dp) :: x(14)
    complex(dp) :: z(14)
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), wave(:)
    integer, allocatable :: peaks(:), from_coefficients(:)
    integer :: i, code

    do i = 1, size(x)
        x(i) = sin(0.13_dp*i)
        z(i) = cmplx(x(i),cos(0.7_dp*i),dp)
    end do
    a = cwt_morlet2(x,[1.0_dp,2.0_dp,3.5_dp],status=code)
    if (code /= cwt_peak_success) error stop 'public real CWT'
    b = cwt_complex_wavelet(x,[1.0_dp,2.0_dp,3.5_dp],morlet2_wavelet,status=code)
    if (code /= cwt_peak_success) error stop 'public callback CWT'
    if (maxval(abs(a-b)) > 1.0e-13_dp) error stop 'public callback values'
    c = cwt_morlet2(z,[0.75_dp,2.5_dp],omega0=4.0_dp,status=code)
    if (code /= cwt_peak_success .or. any(shape(c) /= [2,14])) error stop 'public complex CWT'
    peaks = find_peaks_cwt_complex(x,[1.0_dp,2.0_dp,3.5_dp],status=code)
    if (code /= cwt_peak_success) error stop 'public real CWT peaks'
    from_coefficients = find_peaks_cwt_magnitude(a,[1.0_dp,2.0_dp,3.5_dp],status=code)
    if (code /= cwt_peak_success) error stop 'public magnitude CWT peaks'
    if (size(peaks) /= size(from_coefficients)) error stop 'public real peak count'
    if (any(peaks /= from_coefficients)) error stop 'public real peak positions'
    peaks = find_peaks_cwt_complex(z,[0.75_dp,2.5_dp],omega0=4.0_dp,status=code)
    if (code /= cwt_peak_success) error stop 'public complex CWT peaks'
    from_coefficients = find_peaks_cwt_magnitude(c,[0.75_dp,2.5_dp],status=code)
    if (code /= cwt_peak_success) error stop 'public complex magnitude peaks'
    if (size(peaks) /= size(from_coefficients)) error stop 'public complex peak count'
    if (any(peaks /= from_coefficients)) error stop 'public complex peak positions'
    wave = morlet2_wavelet(9.0_dp,2.0_dp)
    if (size(wave) /= 9) error stop 'public Morlet helper'
    print '(a)', 'PASS complex CWT public interfaces'
end program test_signal_cwt_complex_public
