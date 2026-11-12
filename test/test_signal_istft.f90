! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_istft
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : stft, stft_result, get_window
    use scifort_signal_spectral_complex, only : stft_complex
    use scifort_signal_istft, only : istft, check_cola, check_nola, &
        istft_status_success, istft_status_invalid_input, istft_status_noninvertible
    implicit none

    integer, parameter :: ns = 48
    real(dp) :: x(ns)
    complex(dp) :: xc(ns)
    real(dp), allocatable :: yr(:), times(:), hann(:), symmetry(:)
    complex(dp), allocatable :: yc(:)
    type(stft_result) :: result
    integer :: i, status

    x = [(sin(real(i,dp)*0.17_dp) + 0.2_dp*cos(real(i,dp)*0.9_dp), i=1,ns)]
    xc = cmplx(x, [(0.12_dp*cos(real(i,dp)*0.43_dp), i=1,ns)], dp)

    result = stft(x, fs=8.0_dp, window='hann', nperseg=8, noverlap=4, nfft=16)
    if (result%status /= 0) error stop 'real stft'
    call istft(result%values, yr, status, fs=8.0_dp, window='hann', &
        nperseg=8, noverlap=4, nfft=16, times=times)
    if (status /= istft_status_success) error stop 'real inverse status'
    if (size(yr) /= ns .or. size(times) /= ns) error stop 'real inverse dimensions'
    if (maxval(abs(yr-x)) > 2.0e-12_dp) error stop 'real inverse mismatch'
    if (abs(times(ns)-real(ns-1,dp)/8.0_dp) > 1.0e-14_dp) error stop 'inverse time grid'
    ! Defaults must infer an even unpadded FFT length from the input spectrum.
    result = stft(x, window='hamming', nperseg=8, noverlap=6)
    call istft(result%values, yr, status, window='hamming', noverlap=6)
    if (status /= istft_status_success .or. size(yr) /= ns) error stop 'default inverse shape'
    if (maxval(abs(yr-x)) > 3.0e-12_dp) error stop 'default inverse values'

    ! Odd-length, full-spectrum input must recover complex-valued data.
    result = stft_complex(xc, window='boxcar', nperseg=7, noverlap=4, nfft=7)
    if (result%status /= 0) error stop 'complex stft status'
    call istft(result%values, yc, status, window='boxcar', nperseg=7, &
        noverlap=4, nfft=7)
    if (status /= istft_status_success) error stop 'complex inverse status'
    if (size(yc) < ns) error stop 'complex inverse shape'
    if (maxval(abs(yc(:ns)-xc)) > 4.0e-12_dp) error stop 'complex inverse mismatch'
    if (size(yc) > ns) then
        if (maxval(abs(yc(ns+1:))) > 4.0e-12_dp) error stop 'complex inverse padding'
    end if
    result = stft(x, window='hann', nperseg=7, noverlap=4, nfft=7)
    call istft(result%values, yr, status, window='hann', nperseg=7, noverlap=4)
    if (status /= istft_status_success) error stop 'odd real inverse'
    if (maxval(abs(yr(:ns)-x)) > 2.0e-12_dp) error stop 'odd real values'

    result = stft(x, window='boxcar', nperseg=8, noverlap=0, boundary='none', padded=.false.)
    call istft(result%values, yr, status, window='boxcar', nperseg=8, &
        noverlap=0, boundary=.false.)
    if (status /= istft_status_success .or. size(yr) /= ns) error stop 'no boundary'
    if (maxval(abs(yr-x)) > 2.0e-12_dp) error stop 'no boundary samples'
    result = stft(x, window='hann', nperseg=8, noverlap=4, boundary='none', padded=.false.)
    call istft(result%values, yr, status, window='hann', nperseg=8, noverlap=4, boundary=.false.)
    if (status /= istft_status_noninvertible) error stop 'NOLA status'
    if (maxval(abs(yr(2:)-x(2:))) > 2.0e-12_dp) error stop 'NOLA partial output'

    if (.not. check_cola('hann', 8, 4, status=status)) error stop 'Hann COLA'
    if (status /= istft_status_success) error stop 'Hann check status'
    if (.not. check_nola('hann', 8, 4)) error stop 'Hann NOLA'
    if (check_cola('hann', 8, 0)) error stop 'Hann no-overlap COLA'
    if (check_nola('hann', 8, 0)) error stop 'Hann no-overlap NOLA'
    if (.not. check_cola('boxcar', 8, 0)) error stop 'boxcar COLA'
    if (.not. check_nola('boxcar', 8, 0)) error stop 'boxcar NOLA'
    if (check_cola('hann', 8, 2)) error stop 'non-COLA Hann'
    hann = get_window('hann', 8)
    if (.not. check_cola(hann, 8, 4)) error stop 'explicit COLA'
    if (.not. check_nola(hann, 8, 4)) error stop 'explicit NOLA'
    symmetry = [(0.5_dp-0.5_dp*cos(2.0_dp*acos(-1.0_dp)*real(i-1,dp)/7.0_dp),i=1,8)]
    if (check_cola(symmetry, 8, 4)) error stop 'symmetric Hann COLA'
    if (.not. check_nola(symmetry, 8, 4)) error stop 'symmetric Hann NOLA'
    if (check_nola([1.0_dp,0.0_dp,1.0_dp,0.0_dp],4,2)) error stop 'checkerboard NOLA'
    if (check_cola('invalid',8,4,status=status)) error stop 'unknown window'
    if (status /= istft_status_invalid_input) error stop 'unknown window status'
    if (check_nola('hann',8,-1,status=status)) error stop 'invalid overlap'
    if (status /= istft_status_invalid_input) error stop 'invalid overlap status'

    result = stft(x, window='hann', nperseg=8, noverlap=4)
    call istft(result%values, yr, status, nperseg=8, noverlap=4, fs=-2.0_dp)
    if (status /= istft_status_invalid_input .or. size(yr) /= 0) error stop 'invalid sampling rate'
    call istft(result%values, yr, status, nperseg=8, nfft=5)
    if (status /= istft_status_invalid_input) error stop 'invalid fft size'
    call istft(result%values, yr, status, nperseg=8, input_onesided=.false.)
    if (status /= istft_status_invalid_input) error stop 'invalid real spectrum type'
    call istft(result%values, yr, status, nperseg=8, window='hann', window_values=hann)
    if (status /= istft_status_invalid_input) error stop 'duplicate window inputs'
    call istft(result%values(:,1:0), yr, status)
    if (status /= istft_status_invalid_input) error stop 'empty stft'
    call istft(result%values, yr, status, nperseg=8, noverlap=4, window_values=hann)
    if (status /= istft_status_success) error stop 'explicit inverse window'
    if (maxval(abs(yr-x)) > 2.0e-12_dp) error stop 'explicit inverse values'

    print '(a)', 'test_signal_istft passed'
end program test_signal_istft
