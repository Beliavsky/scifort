program test_signal_short_time_fft
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result, &
        short_time_spectrogram_result, short_fft_success, short_fft_invalid_input, short_fft_noninvertible
    implicit none
    type(short_time_fft) :: plan, from_dual, tight
    type(short_time_fft_result) :: s
    type(short_time_spectrogram_result) :: spec
    real(dp) :: x(35), w(9), tight_reference(9)
    real(dp), allocatable :: y(:)
    complex(dp) :: z(35)
    complex(dp), allocatable :: zrec(:), cross(:,:)
    integer :: i, status
    real(dp), parameter :: tol = 3.0e-12_dp

    do i = 1,size(x)
        x(i) = cos(0.18_dp*real(i,dp)) + 0.3_dp*sin(0.37_dp*real(i,dp))
        z(i) = cmplx(x(i),0.2_dp*cos(0.13_dp*real(i,dp)),dp)
    end do
    call plan%initialize_window('hann', 8, 3, 32.0_dp, mfft=16)
    if (plan%status /= short_fft_success) error stop 'init'
    call plan%stft(x,s)
    if (s%status /= short_fft_success) error stop 'stft'
    if (size(s%values,1) /= 9 .or. size(s%values,2) /= plan%p_num(size(x))) error stop 'shape'
    if (s%p0 /= plan%p_min) error stop 'slice origin'
    call plan%istft(s%values,y,status,k1=size(x))
    if (status /= short_fft_success) error stop 'istft'
    if (maxval(abs(y-x)) > tol) error stop 'real roundtrip'
    call plan%spectrogram(x,spec)
    if (spec%status /= short_fft_success) error stop 'spectrogram'
    if (maxval(abs(spec%power-abs(s%values)**2)) > tol) error stop 'spectral power'
    call plan%cross_spectrogram(x,2.0_dp*x,cross,status)
    if (status /= short_fft_success) error stop 'cross status'
    if (maxval(abs(cross-2.0_dp*abs(s%values)**2)) > tol) error stop 'cross spectrum'

    call plan%initialize_window('hann', 9, 3, 10.0_dp, mfft=11, &
        fft_mode='centered', scale_to='psd', phase_shift=2)
    if (plan%status /= short_fft_success) error stop 'complex init'
    call plan%stft(z,s)
    if (s%status /= short_fft_success) error stop 'complex stft'
    call plan%istft(s%values,zrec,status,k1=size(z))
    if (status /= short_fft_success) error stop 'complex istft'
    if (maxval(abs(zrec-z)) > tol) error stop 'complex roundtrip'

    call plan%initialize_window('hann', 8, 2, 12.0_dp,fft_mode='onesided2X',scale_to='psd')
    call plan%stft(x,s)
    if (s%status /= short_fft_success) error stop 'onesided2X'
    call plan%istft(s%values,y,status,k1=size(x))
    if (status /= short_fft_success) error stop 'onesided2X inverse'
    if (maxval(abs(y-x)) > tol) error stop 'onesided2X roundtrip'

    ! Building from the canonical dual recovers the original analysis window.
    call plan%initialize_window('hann', 9, 3, 10.0_dp, mfft=12, fft_mode='twosided')
    if (plan%status /= short_fft_success) error stop 'dual reference setup'
    call from_dual%initialize_from_dual(plan%dual_win,3,10.0_dp,mfft=12,fft_mode='twosided')
    if (from_dual%status /= short_fft_success) error stop 'from dual status'
    if (maxval(abs(plan%win-from_dual%win)) > tol) error stop 'from dual analysis window'
    if (maxval(abs(plan%dual_win-from_dual%dual_win)) > tol) error stop 'from dual synthesis window'
    call from_dual%stft(z,s)
    if (s%status /= short_fft_success) error stop 'from dual FFT'
    call from_dual%istft(s%values,zrec,status,k1=size(z))
    if (status /= short_fft_success .or. maxval(abs(zrec-z)) > tol) error stop 'from dual roundtrip'

    ! A tight analysis window has a dual differing by just mfft for unitary scaling.
    tight_reference = [ (exp(-((real(i-1,dp)-4.0_dp)/2.0_dp)**2),i=1,9) ]
    call tight%initialize_equal_dual(tight_reference,3,10.0_dp,mfft=12, &
        fft_mode='twosided',scale_to='unitary')
    if (tight%status /= short_fft_success) error stop 'tight window status'
    if (maxval(abs(tight%dual_win - 12.0_dp*tight%win)) > tol) error stop 'tight dual identity'
    call tight%stft(z,s)
    if (s%status /= short_fft_success) error stop 'tight FFT'
    if (abs(sum(abs(s%values)**2)-sum(abs(z)**2)) > tol*real(size(z),dp)) error stop 'unitary Parseval'
    call tight%istft(s%values,zrec,status,k1=size(z))
    if (status /= short_fft_success .or. maxval(abs(zrec-z)) > tol) error stop 'tight roundtrip'

    ! Complex-valued analysis windows have conjugated forward analysis and
    ! an un-conjugated synthesis dual, as in ShortTimeFFT.
    call tight%initialize_equal_dual(cmplx(tight_reference,0.12_dp*tight_reference,dp), &
        3,10.0_dp,mfft=12,fft_mode='centered')
    if (tight%status /= short_fft_success) error stop 'complex tight status'
    call tight%stft(z,s)
    call tight%istft(s%values,zrec,status,k1=size(z))
    if (status /= short_fft_success .or. maxval(abs(zrec-z)) > tol) error stop 'complex tight roundtrip'

    w = 0.0_dp
    w(1) = 1.0_dp
    call plan%initialize(w, 3, 1.0_dp)
    if (plan%status /= short_fft_noninvertible) error stop 'singular window'
    call plan%initialize_window('hann',8,3,-1.0_dp)
    if (plan%status /= short_fft_invalid_input) error stop 'negative sampling rate'
    call plan%initialize_window('hann',8,3,1.0_dp,mfft=7)
    if (plan%status /= short_fft_invalid_input) error stop 'short fft'
    print *, 'test_signal_short_time_fft: all checks passed'
end program test_signal_short_time_fft
