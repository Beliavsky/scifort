! SPDX-License-Identifier: MIT
! Verify the maximum supported STFT ranks through the public umbrella module.
program test_signal_short_time_fft_high_rank_public
    use scifort_signal, only : dp, short_time_fft, short_time_fft_result_7d, &
        short_time_fft_result_15d, short_time_spectrogram_result_7d, &
        short_time_spectrogram_result_15d, short_fft_success
    implicit none
    type(short_time_fft) :: sft
    type(short_time_fft_result_7d) :: frame6
    type(short_time_fft_result_15d) :: frame14
    type(short_time_spectrogram_result_7d) :: power6
    type(short_time_spectrogram_result_15d) :: power14
    real(dp) :: x6(2,1,1,1,1,24)
    real(dp), allocatable :: y6(:,:,:,:,:,:)
    real(dp) :: x14(2,1,1,1,1,1,1,1,1,1,1,1,1,24)
    real(dp), allocatable :: y14(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    integer :: status,i

    call sft%initialize_window('hamming',7,2,16.0_dp,mfft=12)
    if (sft%status /= short_fft_success) error stop 'high-rank public plan'
    x6 = reshape([(sin(0.07_dp*real(i,dp)),i=1,size(x6))],shape(x6))
    call sft%stft(x6,frame6,axis=6)
    if (frame6%status /= short_fft_success) error stop 'rank6 public STFT'
    call sft%istft(frame6%values,y6,status,k1=24,f_axis=6)
    if (status /= short_fft_success) error stop 'rank6 public iSTFT'
    if (maxval(abs(x6-y6)) > 2.e-12_dp) error stop 'rank6 public round trip'
    call sft%spectrogram(x6,power6,axis=6)
    if (power6%status /= short_fft_success) error stop 'rank6 public power'
    x14 = reshape([(sin(0.07_dp*real(i,dp)),i=1,size(x14))],shape(x14))
    call sft%stft(x14,frame14,axis=14)
    if (frame14%status /= short_fft_success) error stop 'rank14 public STFT'
    call sft%istft(frame14%values,y14,status,k1=24,f_axis=14)
    if (status /= short_fft_success) error stop 'rank14 public iSTFT'
    if (maxval(abs(x14-y14)) > 2.e-12_dp) error stop 'rank14 public round trip'
    call sft%spectrogram(x14,power14,axis=14)
    if (power14%status /= short_fft_success) error stop 'rank14 public power'
    print *, 'test_signal_short_time_fft_high_rank_public: all checks passed'
end program test_signal_short_time_fft_high_rank_public
