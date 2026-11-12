! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_time_frequency
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_psd, multitaper_spectrogram, &
        multitaper_spectrogram_result, multitaper_spectrogram_channels_result, &
        multitaper_harmonic_ftest, multitaper_harmonic_result
    use scifort_signal_spectral, only : spectral_result, spectral_status_success, spectral_status_invalid_input
    implicit none
    integer, parameter :: n = 160, segment_length = 40
    real(dp), parameter :: nw = 3.5_dp, fs = 100.0_dp
    real(dp) :: x(n), y(n), data(n,2), phase, pi
    complex(dp) :: z(n), zdata(n,2), tone(n)
    complex(dp), parameter :: amplitude = (1.1_dp,-0.6_dp)
    type(multitaper_spectrogram_result) :: st, st2
    type(multitaper_spectrogram_channels_result) :: channels, channels2
    type(multitaper_harmonic_result) :: ft, fc, zero, bad
    type(spectral_result) :: psd
    integer :: i

    pi = acos(-1.0_dp)
    do i = 1, n
        phase = 2.0_dp*pi*real(i-1,dp)/real(n,dp)
        x(i) = sin(15.0_dp*phase) + 0.2_dp*cos(37.0_dp*phase) + 0.05_dp*cos(9.0_dp*phase)
        y(i) = cos(23.0_dp*phase) + 0.35_dp*x(i)
        z(i) = cmplx(x(i),y(i),dp)
        tone(i) = amplitude*exp(cmplx(0.0_dp,27.0_dp*phase,dp))
    end do
    data(:,1) = x
    data(:,2) = y
    zdata(:,1) = z
    zdata(:,2) = tone

    st = multitaper_spectrogram(x,nw,segment_length,noverlap=20,kmax=5,fs=fs,nfft=64, &
        detrend=.false.,weighting='adaptive')
    if (st%status /= spectral_status_success) error stop 'real spectrogram status'
    if (size(st%times) /= 7 .or. size(st%power,1) /= 33) error stop 'spectrogram shape'
    if (abs(st%times(1)-0.2_dp) > 1.0e-12_dp) error stop 'spectrogram first time'
    if (abs(st%times(7)-1.4_dp) > 1.0e-12_dp) error stop 'spectrogram last time'
    if (any(st%power < 0.0_dp)) error stop 'negative spectral density'
    psd = multitaper_psd(x(1:segment_length),nw,kmax=5,fs=fs,nfft=64, &
        detrend=.false.,weighting='adaptive')
    if (psd%status /= spectral_status_success) error stop 'reference PSD status'
    if (maxval(abs(st%power(:,1)-psd%power)) > 2.0e-12_dp) &
        error stop 'spectrogram adaptive PSD parity'
    psd = multitaper_psd(x(41:80),nw,kmax=5,fs=fs,nfft=64,detrend=.false.,weighting='eigen')
    st2 = multitaper_spectrogram(x,nw,segment_length,noverlap=0,kmax=5,fs=fs,nfft=64, &
        detrend=.false.,weighting='eigen')
    if (size(st2%times) /= 4) error stop 'no-overlap count'
    if (maxval(abs(st2%power(:,2)-psd%power)) > 2.0e-12_dp) &
        error stop 'spectrogram eigen PSD parity'
    channels = multitaper_spectrogram(data,nw,segment_length,noverlap=20,kmax=5,fs=fs, &
        nfft=64,detrend=.false.,weighting='adaptive')
    channels2 = multitaper_spectrogram(transpose(data),nw,segment_length,noverlap=20,kmax=5,fs=fs, &
        nfft=64,detrend=.false.,weighting='adaptive',axis=2)
    if (channels%status /= spectral_status_success) error stop 'channel PSD status'
    if (maxval(abs(channels%power-channels2%power)) > 3.0e-12_dp) error stop 'sample axis mismatch'
    if (maxval(abs(channels%power(:,:,1)-st%power)) > 3.0e-12_dp) error stop '1-D channel parity'
    st = multitaper_spectrogram(z,nw,segment_length,noverlap=20,kmax=5,fs=fs,nfft=65, &
        detrend=.true.,weighting='adaptive')
    channels = multitaper_spectrogram(zdata,nw,segment_length,noverlap=20,kmax=5,fs=fs,nfft=65, &
        detrend=.true.,weighting='adaptive')
    if (st%status /= spectral_status_success) error stop 'complex spectrogram status'
    if (size(st%power,1) /= 65) error stop 'complex two-sided FFT bins'
    if (maxval(abs(channels%power(:,:,1)-st%power)) > 3.0e-12_dp) error stop 'complex channel parity'
    if (st%frequencies(65) >= 0.0_dp) error stop 'negative FFT frequencies'

    ft = multitaper_harmonic_ftest(x,nw,kmax=5,fs=fs,nfft=n)
    if (ft%status /= spectral_status_success) error stop 'real F status'
    if (size(ft%pvalue) /= n/2+1) error stop 'real F frequency bins'
    if (any(ft%pvalue < 0.0_dp) .or. any(ft%pvalue > 1.0_dp)) error stop 'invalid probabilities'
    if (ft%pvalue(16) >= 1.0e-5_dp) error stop 'real line not detected'
    fc = multitaper_harmonic_ftest(tone,nw,kmax=5,fs=fs,nfft=n)
    if (fc%status /= spectral_status_success) error stop 'complex F status'
    if (size(fc%pvalue) /= n) error stop 'complex F frequency bins'
    if (fc%pvalue(28) >= 1.0e-8_dp) error stop 'complex line not detected'
    if (abs(fc%amplitude(28)-amplitude) > 1.0e-11_dp) error stop 'harmonic amplitude'
    if (fc%frequencies(n) >= 0.0_dp) error stop 'complex negative frequency'
    zero = multitaper_harmonic_ftest(0.0_dp*x,nw,kmax=5)
    if (zero%status /= spectral_status_success) error stop 'zero signal F status'
    if (any(zero%f_statistic /= 0.0_dp) .or. any(zero%pvalue /= 1.0_dp)) error stop 'zero signal F values'

    st = multitaper_spectrogram(x,nw,nperseg=0)
    if (st%status /= spectral_status_invalid_input) error stop 'short segment accepted'
    st = multitaper_spectrogram(x,nw,nperseg=40,noverlap=40)
    if (st%status /= spectral_status_invalid_input) error stop 'full overlap accepted'
    st = multitaper_spectrogram(x,nw,nperseg=40,nfft=39)
    if (st%status /= spectral_status_invalid_input) error stop 'short nfft accepted'
    st = multitaper_spectrogram(x,nw,nperseg=40,weighting='unknown')
    if (st%status /= spectral_status_invalid_input) error stop 'unknown weighting accepted'
    channels = multitaper_spectrogram(data,nw,40,axis=3)
    if (channels%status /= spectral_status_invalid_input) error stop 'invalid axis accepted'
    bad = multitaper_harmonic_ftest(x,nw,kmax=1)
    if (bad%status /= spectral_status_invalid_input) error stop 'one taper F accepted'
    bad = multitaper_harmonic_ftest(x,nw,nfft=n-1)
    if (bad%status /= spectral_status_invalid_input) error stop 'short F FFT accepted'
    print '(a)', 'Multitaper spectrogram and harmonic F-test passed'
end program test_signal_multitaper_time_frequency
