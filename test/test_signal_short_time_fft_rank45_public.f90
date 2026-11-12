! SPDX-License-Identifier: MIT
! Verify rank-four and rank-five STFT exports from the public API.
program test_signal_short_time_fft_rank45_public
    use scifort_signal, only : dp, short_time_fft, short_time_fft_result_5d, &
        short_time_fft_result_6d, short_time_spectrogram_result_5d, &
        short_time_spectrogram_result_6d, short_fft_success
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result_5d) :: spec4
    type(short_time_fft_result_6d) :: spec5
    type(short_time_spectrogram_result_5d) :: p4
    type(short_time_spectrogram_result_6d) :: p5
    real(dp) :: x4(2,2,2,24), x5(2,2,2,2,24)
    real(dp), allocatable :: y4(:,:,:,:), y5(:,:,:,:,:)
    integer :: i,j,k,l,m,status

    do l=1,24
        do k=1,2
            do j=1,2
                do i=1,2
                    x4(i,j,k,l) = sin(0.11_dp*l+0.09_dp*i+0.04_dp*j+0.07_dp*k)
                    do m=1,2
                        x5(i,j,k,m,l) = x4(i,j,k,l)+0.02_dp*m
                    end do
                end do
            end do
        end do
    end do
    call plan%initialize_window('hann',8,3,16.0_dp,mfft=16)
    if (plan%status /= short_fft_success) error stop 'rank45 public plan'
    call plan%stft(x4,spec4)
    if (spec4%status /= short_fft_success) error stop 'rank4 public STFT'
    call plan%istft(spec4%values,y4,status,k1=24)
    if (status /= short_fft_success) error stop 'rank4 public inverse'
    if (maxval(abs(x4-y4)) > 1.e-12_dp) error stop 'rank4 public samples'
    call plan%spectrogram(x4,p4)
    if (p4%status /= short_fft_success) error stop 'rank4 public power'
    call plan%stft(x5,spec5)
    if (spec5%status /= short_fft_success) error stop 'rank5 public STFT'
    call plan%istft(spec5%values,y5,status,k1=24)
    if (status /= short_fft_success) error stop 'rank5 public inverse'
    if (maxval(abs(x5-y5)) > 1.e-12_dp) error stop 'rank5 public samples'
    call plan%spectrogram(x5,p5)
    if (p5%status /= short_fft_success) error stop 'rank5 public power'
    print *, 'test_signal_short_time_fft_rank45_public passed'
end program test_signal_short_time_fft_rank45_public
