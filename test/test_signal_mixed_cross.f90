! SPDX-License-Identifier: MIT
! Native regression for mixed real/complex cross-spectral analysis.
program test_signal_mixed_cross
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_result, cross_spectral_result, &
        spectral_status_success, spectral_status_invalid_input
    use scifort_signal_spectral_complex, only : csd_complex, coherence_complex, &
        csd_real_complex, csd_complex_real, coherence_real_complex, coherence_complex_real
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result, &
        short_fft_success, short_fft_invalid_input
    implicit none
    real(dp) :: x(87), xs(41)
    complex(dp) :: z(87), zs(41)
    complex(dp), allocatable :: cross(:,:), target(:,:)
    type(cross_spectral_result) :: c, cc
    type(spectral_result) :: q, qq
    type(short_time_fft) :: plan
    type(short_time_fft_result) :: sx, sz
    real(dp), parameter :: tol = 2.0e-11_dp
    integer :: i, status

    do i = 1, size(x)
        x(i) = 0.7_dp*cos(0.19_dp*i) + sin(0.43_dp*i) + 0.013_dp*i
        z(i) = cmplx(0.4_dp*sin(0.23_dp*i) + x(i), &
            cos(0.37_dp*i) - 0.15_dp*sin(0.09_dp*i), dp)
    end do
    xs = x(2:82:2)
    zs = z(2:82:2)

    c = csd_real_complex(x,z,nperseg=15,noverlap=7,nfft=19, &
        window='gaussian',window_parameter=3.5_dp,average='median', &
        scaling='spectrum',detrend='linear',return_onesided=.true.)
    cc = csd_complex(cmplx(x,0.0_dp,dp),z,nperseg=15,noverlap=7,nfft=19, &
        window='gaussian',window_parameter=3.5_dp,average='median', &
        scaling='spectrum',detrend='linear')
    if (c%status /= spectral_status_success .or. cc%status /= spectral_status_success) &
        error stop 'real-complex CSD status'
    if (size(c%power) /= 19 .or. size(cc%power) /= 19) error stop 'mixed CSD must be two-sided'
    if (maxval(abs(c%power-cc%power)) > tol) error stop 'real-complex CSD mismatch'
    if (maxval(abs(c%frequencies-cc%frequencies)) > tol) error stop 'CSD frequencies'

    c = csd_complex_real(z,x,nperseg=14,noverlap=5,nfft=16,detrend='none')
    cc = csd_complex(z,cmplx(x,0.0_dp,dp),nperseg=14,noverlap=5,nfft=16,detrend='none')
    if (c%status /= spectral_status_success) error stop 'complex-real CSD status'
    if (maxval(abs(c%power-cc%power)) > tol) error stop 'complex-real CSD mismatch'

    ! Call directly with a noncontiguous section (not a prepacked copy).
    c = csd_real_complex(x(2:82:2),z(2:82:2),nperseg=11,noverlap=4,nfft=17)
    cc = csd_real_complex(xs,zs,nperseg=11,noverlap=4,nfft=17)
    if (c%status /= spectral_status_success) error stop 'strided CSD status'
    if (maxval(abs(c%power-cc%power)) > tol) error stop 'strided CSD mismatch'

    q = coherence_real_complex(xs,zs,nperseg=12,noverlap=3,nfft=15)
    qq = coherence_complex(cmplx(xs,0.0_dp,dp),zs,nperseg=12,noverlap=3,nfft=15)
    if (q%status /= spectral_status_success .or. size(q%power) /= 15) &
        error stop 'real-complex coherence status'
    if (maxval(abs(q%power-qq%power)) > tol) error stop 'real-complex coherence mismatch'
    q = coherence_complex_real(zs,xs,nperseg=11,noverlap=4,nfft=14, &
        window='hann',detrend='linear')
    qq = coherence_complex(zs,cmplx(xs,0.0_dp,dp),nperseg=11,noverlap=4, &
        nfft=14,window='hann',detrend='linear')
    if (q%status /= spectral_status_success) error stop 'complex-real coherence status'
    if (maxval(abs(q%power-qq%power)) > tol) error stop 'complex-real coherence mismatch'

    c = csd_real_complex(xs,zs,nperseg=9,noverlap=5,nfft=13, &
        detrend_callback=center_complex)
    cc = csd_complex(cmplx(xs,0.0_dp,dp),zs,nperseg=9,noverlap=5,nfft=13, &
        detrend_callback=center_complex)
    if (c%status /= spectral_status_success) error stop 'callback status'
    if (maxval(abs(c%power-cc%power)) > tol) error stop 'callback mismatch'
    c = csd_real_complex(xs,zs,nperseg=9,noverlap=5,nfft=13, &
        detrend='constant',detrend_callback=center_complex)
    if (c%status /= spectral_status_invalid_input) error stop 'callback conflict accepted'
    if (size(c%power) /= 0) error stop 'invalid CSD returned bins'
    c = csd_complex_real(zs,xs,nperseg=9,noverlap=9)
    if (c%status /= spectral_status_invalid_input) error stop 'noverlap >= nperseg accepted'

    call plan%initialize_window('hann',9,3,24.0_dp,mfft=13, &
        fft_mode='centered',scale_to='psd',phase_shift=2)
    if (plan%status /= short_fft_success) error stop 'plan initialization'
    call plan%cross_spectrogram(xs,zs,cross,status,padding='odd')
    if (status /= short_fft_success) error stop 'real-complex cross spectrogram status'
    call plan%stft(xs,sx,padding='odd')
    call plan%stft(zs,sz,padding='odd')
    if (sx%status /= short_fft_success .or. sz%status /= short_fft_success) error stop 'STFT status'
    target = sx%values * conjg(sz%values)
    if (any(shape(cross) /= shape(target))) error stop 'cross shape'
    if (maxval(abs(cross-target)) > tol) error stop 'real-complex cross mismatch'
    call plan%cross_spectrogram(x(2:82:2),z(2:82:2),cross,status,padding='odd')
    if (status /= short_fft_success) error stop 'strided mixed cross spectrum status'
    if (maxval(abs(cross-target)) > tol) error stop 'strided cross mismatch'
    call plan%cross_spectrogram(zs,xs,cross,status,padding='zeros')
    if (status /= short_fft_success) error stop 'complex-real cross spectrogram status'
    call plan%stft(zs,sz,padding='zeros')
    call plan%stft(xs,sx,padding='zeros')
    target = sz%values * conjg(sx%values)
    if (maxval(abs(cross-target)) > tol) error stop 'complex-real cross mismatch'

    call plan%cross_spectrogram(x(1:20),z(1:21),cross,status)
    if (status /= short_fft_invalid_input) error stop 'unequal signal lengths accepted'
    if (size(cross) /= 0) error stop 'failure did not return empty output'
    call plan%initialize_window('hann',9,3,24.0_dp,mfft=13,fft_mode='onesided')
    call plan%cross_spectrogram(xs,zs,cross,status)
    if (status /= short_fft_invalid_input) error stop 'mixed one-sided cross spectrum accepted'
    if (size(cross) /= 0) error stop 'invalid one-sided spectrum not empty'

    print '(a)', 'test_signal_mixed_cross passed'
contains
    subroutine center_complex(a)
        complex(dp), intent(inout) :: a(:) !! single unwindowed spectral segment
        a = a - sum(a)/real(size(a),dp)
    end subroutine center_complex
end program test_signal_mixed_cross
