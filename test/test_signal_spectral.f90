! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_spectral
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_signal
    implicit none

    type(spectral_result) :: p, q
    type(cross_spectral_result) :: c, d
    type(stft_result) :: t
    type(spectrogram_result) :: s
    real(dp), allocatable :: w(:)
    real(dp) :: a(8), b(8), sinusoid(8), many(35), pi
    integer :: i, status

    pi = acos(-1.0_dp)
    a = [(real(i-1,dp),i=1,8)]
    b = a**2 / 8.0_dp
    sinusoid = [(sin(2.0_dp*pi*real(i-1,dp)/8.0_dp),i=1,8)]
    many = [(sin(2.0_dp*pi*real(i-1,dp)/7.0_dp) + &
        0.1_dp * real(mod(i,5),dp),i=1,35)]

    w = get_window('hann',8,status=status)
    if (status /= spectral_status_success) error stop 'hann status'
    call close_real(w,[0.0_dp,0.1464466094067262_dp,0.5_dp, &
        0.8535533905932737_dp,1.0_dp,0.8535533905932737_dp, &
        0.5_dp,0.1464466094067262_dp])
    w = get_window('boxcar',4)
    call close_real(w,[1.0_dp,1.0_dp,1.0_dp,1.0_dp])
    w = get_window('hann',1)
    call close_real(w,[1.0_dp])
    w = get_window('bad-name',8,status=status)
    if (status /= spectral_status_invalid_input .or. size(w) /= 0) error stop 'bad window'

    p = periodogram(sinusoid, fs=8.0_dp, window='boxcar', detrend='none')
    if (p%status /= spectral_status_success) error stop 'periodogram status'
    call close_real(p%frequencies,[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp])
    call close_real(p%power,[0.0_dp,0.5_dp,0.0_dp,0.0_dp,0.0_dp])
    p = periodogram(a, fs=8.0_dp, nfft=8)
    call close_real(p%power(1:3),[0.0_dp,3.414213562373095_dp,1.0_dp])
    p = periodogram(many,fs=1.0_dp,window='boxcar',detrend='none',scaling='spectrum')
    if (abs(sum(p%power) - sum(many**2)/real(size(many),dp)) > 3.0e-12_dp) &
        error stop 'non-power-two Parseval'
    p = periodogram(a,fs=1.0_dp,nfft=4)
    q = periodogram(a(1:4),fs=1.0_dp)
    call close_real(p%power,q%power)

    p = welch(a,fs=8.0_dp,nperseg=4,noverlap=2,nfft=8)
    call close_real(p%power(1:3),[1.0_dp/12.0_dp,0.2047588984322123_dp,5.0_dp/24.0_dp])
    c = csd(a,a,fs=8.0_dp,nperseg=4,noverlap=2,nfft=8)
    call close_complex(c%power,cmplx(p%power,0.0_dp,dp))
    c = csd(a,b,fs=8.0_dp,nperseg=4,noverlap=2,nfft=8)
    d = csd(b,a,fs=8.0_dp,nperseg=4,noverlap=2,nfft=8)
    call close_complex(c%power,conjg(d%power))
    q = coherence(a,a,fs=8.0_dp,nperseg=4,noverlap=2,nfft=8)
    call close_real(q%power,spread(1.0_dp,1,size(q%power)))

    p = welch(many,fs=1.0_dp,nperseg=7,noverlap=3,nfft=13,average='median')
    if (p%status /= spectral_status_success) error stop 'median status'
    if (any(.not.(p%power >= 0.0_dp))) error stop 'median nonnegative'

    t = stft(a,fs=8.0_dp,nperseg=4,noverlap=2,nfft=8)
    if (t%status /= spectral_status_success) error stop 'stft status'
    call close_real(t%times,[0.0_dp,0.25_dp,0.5_dp,0.75_dp,1.0_dp])
    call close_real(t%frequencies,[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp])
    if (abs(t%values(1,2)-2.0_dp) > 1.0e-12_dp) error stop 'stft amplitude'
    t = stft(a,nperseg=3,noverlap=1,nfft=5)
    if (abs(t%times(1)) > 1.0e-12_dp) error stop 'odd stft boundary time'
    s = spectrogram(a,fs=8.0_dp,nperseg=4,noverlap=2,nfft=8)
    if (s%status /= spectral_status_success) error stop 'spectrogram status'
    call close_real(s%times,[0.25_dp,0.5_dp,0.75_dp])
    call close_real(s%power(1,:),spread(0.09375_dp,1,3))

    p = periodogram(a,fs=-1.0_dp)
    if (p%status /= spectral_status_invalid_input) error stop 'negative fs status'
    p = welch(a,nperseg=4,noverlap=4)
    if (p%status /= spectral_status_invalid_input) error stop 'overlap status'
    p = welch(a,nperseg=8,nfft=4)
    if (p%status /= spectral_status_invalid_input) error stop 'short FFT status'
    p = periodogram(a,scaling='bad')
    if (p%status /= spectral_status_invalid_input) error stop 'bad scaling'
    t = stft(a,boundary='wrap')
    if (t%status /= spectral_status_invalid_input) error stop 'unsupported extension'
    s = spectrogram(a,window='bad')
    if (s%status /= spectral_status_invalid_input) error stop 'unsupported window'
    b(3) = ieee_value(0.0_dp,ieee_quiet_nan)
    c = csd(a,b)
    if (c%status /= spectral_status_invalid_input) error stop 'NaN signal'
    p = coherence(spread(0.0_dp,1,8),spread(0.0_dp,1,8),nperseg=4)
    if (p%status /= spectral_status_success .or. .not. all(ieee_is_nan(p%power))) &
        error stop 'undefined coherence'
    print '(a)', 'test_signal_spectral passed'
contains
    subroutine close_real(actual,expected)
        real(dp), intent(in) :: actual(:) !! actual values
        real(dp), intent(in) :: expected(:) !! reference values
        if (size(actual) /= size(expected)) error stop 'spectral real length mismatch'
        if (any(abs(actual-expected)>2.0e-10_dp)) then
            print *, 'actual ', actual
            print *, 'expected ', expected
            error stop 'spectral real values mismatch'
        end if
    end subroutine close_real
    subroutine close_complex(actual,expected)
        complex(dp), intent(in) :: actual(:) !! actual spectrum
        complex(dp), intent(in) :: expected(:) !! reference spectrum
        if (size(actual) /= size(expected)) error stop 'spectral complex length mismatch'
        if (any(abs(actual-expected)>2.0e-10_dp)) error stop 'spectral complex values mismatch'
    end subroutine close_complex
end program test_signal_spectral
