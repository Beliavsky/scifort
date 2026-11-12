! SPDX-License-Identifier: MIT
! Self-contained complex/two-sided spectral analysis regressions.
program test_signal_spectral_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal
    implicit none
    complex(dp) :: x(12), y(12)
    real(dp) :: rate
    type(spectral_result) :: p, q
    type(cross_spectral_result) :: c, d
    type(stft_result) :: z
    type(spectrogram_result) :: s
    integer :: i

    rate = 12.0_dp
    x = [(cmplx(cos(2.0_dp*acos(-1.0_dp)*real(i-1,dp)/12.0_dp), &
                     sin(2.0_dp*acos(-1.0_dp)*real(i-1,dp)/12.0_dp),dp),i=1,12)]
    y = cmplx(0.5_dp,1.25_dp,dp)*x
    p = periodogram(x,fs=rate,window='boxcar',detrend='none',scaling='spectrum')
    if (p%status /= spectral_status_success) error stop 'complex periodogram status'
    if (size(p%power) /= 12) error stop 'complex spectrum must be two-sided'
    if (maxval(abs(p%power - [(merge(1.0_dp,0.0_dp,i==2),i=1,12)])) > 1.0e-12_dp) &
        error stop 'complex sinusoid periodogram'
    if (abs(p%frequencies(7)+6.0_dp)>1.0e-14_dp) error stop 'negative Nyquist'
    if (abs(p%frequencies(8)+5.0_dp)>1.0e-14_dp) error stop 'fft frequency order'
    q = welch(x,fs=rate,window='boxcar',nperseg=12,noverlap=0,detrend='none',scaling='spectrum')
    if (maxval(abs(p%power - q%power))>1.0e-12_dp) error stop 'complex Welch'
    c = csd(x,y,fs=rate,window='boxcar',nperseg=12,noverlap=0,detrend='none',scaling='spectrum')
    if (abs(c%power(2)-cmplx(0.5_dp,1.25_dp,dp))>1.0e-12_dp) error stop 'cross-spectrum sign'
    d = csd(y,x,fs=rate,window='boxcar',nperseg=12,noverlap=0,detrend='none',scaling='spectrum')
    if (maxval(abs(c%power-conjg(d%power)))>1.0e-12_dp) error stop 'CSD conjugacy'
    q = coherence(x,y,fs=rate,nperseg=12,noverlap=0,detrend='none')
    if (abs(q%power(2)-1.0_dp)>1.0e-12_dp) error stop 'complex coherence'
    z = stft(x,fs=rate,window='boxcar',nperseg=12,noverlap=0,nfft=12, &
        boundary='none',padded=.false.,detrend='none')
    if (z%status/=spectral_status_success) error stop 'complex STFT status'
    if (size(z%frequencies)/=12 .or. size(z%times)/=1) error stop 'complex STFT shape'
    if (abs(z%values(2,1)-1.0_dp)>1.0e-12_dp) error stop 'complex STFT coefficient'
    s = spectrogram(x,fs=rate,window='boxcar',nperseg=12,noverlap=0,nfft=12,detrend='none')
    if (s%status/=spectral_status_success .or. size(s%power,1)/=12) error stop 'spectrogram shape'
    if (abs(s%power(2,1)-1.0_dp)>1.0e-12_dp) error stop 'spectrogram power normalization'
    p = periodogram(x,return_onesided=.true.,detrend='none')
    if (size(p%power)/=12) error stop 'complex must ignore one-sided request'
    x(3) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    p = welch(x)
    if (p%status /= spectral_status_invalid_input) error stop 'nonfinite complex input'
    print '(a)', 'test_signal_spectral_complex passed'
end program test_signal_spectral_complex
