! SPDX-License-Identifier: MIT
! Verify multidimensional ShortTimeFFT and automatic REML via umbrella module.
program test_signal_short_time_fft_nd_public
    use scifort_signal, only : dp, short_time_fft, short_time_fft_result_3d, &
        short_time_fft_result_4d, short_time_spectrogram_result_3d, &
        short_time_spectrogram_result_4d, short_fft_success, &
        whittaker_result, whittaker_henderson, whittaker_status_success
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result_3d) :: s2
    type(short_time_fft_result_4d) :: s3
    type(short_time_spectrogram_result_3d) :: power2
    type(short_time_spectrogram_result_4d) :: power3
    type(whittaker_result) :: fit
    real(dp) :: x(2,24), y(2,3,24)
    real(dp), allocatable :: rec(:,:), rec3(:,:,:)
    integer :: i, j, k, status

    do i=1,24
        do j=1,2
            x(j,i) = sin(0.15_dp*i)+0.02_dp*j*i
        end do
        do k=1,3
            do j=1,2
                y(j,k,i) = sin(0.15_dp*i)+0.02_dp*j*i+0.1_dp*k
            end do
        end do
    end do
    call plan%initialize_window('hann',8,3,16.0_dp,mfft=16)
    if (plan%status /= short_fft_success) error stop 'public 2D plan'
    call plan%stft(x,s2)
    if (s2%status /= short_fft_success) error stop 'public 2D STFT'
    call plan%spectrogram(x,power2)
    if (power2%status /= short_fft_success) error stop 'public 2D spectrogram'
    call plan%istft(s2%values,rec,status,k1=24)
    if (status /= short_fft_success) error stop 'public 2D inverse'
    if (maxval(abs(rec-x)) > 1.e-12_dp) error stop 'public 2D roundtrip'
    call plan%stft(y,s3)
    if (s3%status /= short_fft_success) error stop 'public 3D STFT'
    call plan%spectrogram(y,power3)
    if (power3%status /= short_fft_success) error stop 'public 3D spectrogram'
    call plan%istft(s3%values,rec3,status,k1=24)
    if (status /= short_fft_success) error stop 'public 3D inverse'
    if (maxval(abs(rec3-y)) > 1.e-12_dp) error stop 'public 3D roundtrip'
    fit = whittaker_henderson(x(1,:))
    if (fit%status /= whittaker_status_success) error stop 'public WH REML'
    if (fit%lamb <= 0.0_dp) error stop 'public WH penalty'
    print *, 'test_signal_short_time_fft_nd_public passed'
end program test_signal_short_time_fft_nd_public
