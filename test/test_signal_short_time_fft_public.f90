! SPDX-License-Identifier: MIT
program test_signal_short_time_fft_public
    use scifort_signal, only : dp, short_time_fft, short_time_fft_result, &
        short_time_spectrogram_result, short_fft_success
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result) :: spectrum
    type(short_time_spectrogram_result) :: power
    real(dp), allocatable :: data(:)
    integer :: status

    call plan%initialize_window('hann',8,2,16.0_dp,mfft=16)
    if (plan%status /= short_fft_success) error stop 'public init'
    call plan%stft([1.0_dp,2.0_dp,3.0_dp,2.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp],spectrum)
    if (spectrum%status /= short_fft_success) error stop 'public stft'
    call plan%istft(spectrum%values,data,status,k1=8)
    if (status /= short_fft_success) error stop 'public istft'
    if (maxval(abs(data-[1.0_dp,2.0_dp,3.0_dp,2.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]))>1.e-12_dp) &
        error stop 'public roundtrip'
    call plan%spectrogram(data,power)
    if (power%status /= short_fft_success) error stop 'public spectrogram'
    print *, 'public ShortTimeFFT smoke test passed'
end program test_signal_short_time_fft_public
