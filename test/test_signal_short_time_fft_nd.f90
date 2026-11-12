! SPDX-License-Identifier: MIT
! Rank-two and rank-three ShortTimeFFT regression, including arbitrary inverse axes.
program test_signal_short_time_fft_nd
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result_3d, &
        short_time_fft_result_4d, short_time_spectrogram_result_3d, &
        short_time_spectrogram_result_4d, short_fft_success, short_fft_invalid_input
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result_3d) :: s2
    type(short_time_fft_result_4d) :: s3
    type(short_time_spectrogram_result_3d) :: power2
    type(short_time_spectrogram_result_4d) :: power3
    real(dp) :: x2(35,3), x3(3,35,2), yrev(3,35)
    real(dp), allocatable :: y2(:,:), y3(:,:,:)
    complex(dp) :: z2(35,3), z3(3,35,2)
    complex(dp), allocatable :: zy2(:,:), zy3(:,:,:), cross2(:,:,:), cross3(:,:,:,:)
    complex(dp), allocatable :: reordered3(:,:,:), reordered4(:,:,:,:)
    integer :: i, j, k, st
    real(dp), parameter :: tol = 5.0e-12_dp

    do j=1,3
        do i=1,35
            x2(i,j) = sin(0.13_dp*i+0.07_dp*j)+0.2_dp*cos(0.19_dp*i-0.1_dp*j)
            z2(i,j) = cmplx(x2(i,j),0.17_dp*sin(0.11_dp*i+0.03_dp*j),dp)
        end do
    end do
    yrev = transpose(x2)
    do k=1,2
        do j=1,35
            do i=1,3
                x3(i,j,k) = sin(0.13_dp*i+0.07_dp*j+0.09_dp*k) + &
                    0.2_dp*cos(0.19_dp*i-0.1_dp*j+0.05_dp*k)
                z3(i,j,k) = cmplx(x3(i,j,k),0.17_dp*sin(0.11_dp*i+0.03_dp*j-0.02_dp*k),dp)
            end do
        end do
    end do

    call plan%initialize_window('hann',9,3,20.0_dp,mfft=13,fft_mode='onesided2X',scale_to='psd')
    if (plan%status /= short_fft_success) error stop 'nd plan'
    call plan%stft(x2,s2,axis=1,padding='edge')
    if (s2%status /= short_fft_success) error stop '2D real STFT axis 1'
    if (any(shape(s2%values) /= [7,3,plan%p_num(35)])) error stop '2D spectrum dimensions'
    call plan%istft(s2%values,y2,st,k1=35,f_axis=1)
    if (st /= short_fft_success) error stop '2D inverse axis 1'
    if (any(shape(y2) /= shape(x2)) .or. maxval(abs(y2-x2)) > tol) error stop '2D roundtrip axis 1'
    call plan%cross_spectrogram(x2,2.0_dp*x2,cross2,st,axis=1,padding='edge')
    if (st /= short_fft_success) error stop '2D cross'
    if (maxval(abs(cross2-2.0_dp*abs(s2%values)**2)) > tol) error stop '2D cross values'
    call plan%spectrogram(x2,power2,axis=1,padding='edge')
    if (power2%status /= short_fft_success) error stop '2D spectrogram'
    if (maxval(abs(power2%power-abs(s2%values)**2)) > tol) error stop '2D power'
    call plan%stft(yrev,s2,axis=-1,padding='zeros')
    if (s2%status /= short_fft_success) error stop '2D default axis'
    if (any(shape(s2%values) /= [3,7,plan%p_num(35)])) error stop '2D default shape'
    call plan%istft(s2%values,y2,st,k1=35)
    if (st /= short_fft_success) error stop '2D default inverse'
    if (maxval(abs(y2-yrev)) > tol) error stop '2D default roundtrip'
    call plan%stft(x2,s2,axis=0)
    if (s2%status /= short_fft_invalid_input) error stop 'invalid forward axis'
    call plan%stft(yrev,s2,axis=2)
    call plan%istft(s2%values,y2,st,k1=35,f_axis=3)
    if (st /= short_fft_invalid_input) error stop 'invalid inverse axis'
    call plan%istft(s2%values,y2,st,k1=35,f_axis=2,t_axis=2)
    if (st /= short_fft_invalid_input) error stop 'duplicate inverse axes'

    call plan%stft(x3,s3,axis=2,padding='even')
    if (s3%status /= short_fft_success) error stop '3D real STFT'
    if (any(shape(s3%values) /= [3,7,2,plan%p_num(35)])) error stop '3D spectrum dimensions'
    call plan%istft(s3%values,y3,st,k1=35,f_axis=2)
    if (st /= short_fft_success) error stop '3D inverse'
    if (any(shape(y3) /= shape(x3)) .or. maxval(abs(y3-x3)) > tol) error stop '3D roundtrip'
    call plan%stft_detrend(x3,s3,'constant',axis=2,padding='odd')
    if (s3%status /= short_fft_success) error stop '3D detrend'
    call plan%stft(x3,s3,axis=2,padding='even')
    call plan%spectrogram(x3,power3,axis=2,padding='even')
    if (power3%status /= short_fft_success) error stop '3D spectrogram'
    if (maxval(abs(power3%power-abs(s3%values)**2)) > tol) error stop '3D power'
    call plan%cross_spectrogram(x3,2.0_dp*x3,cross3,st,axis=2,padding='even')
    if (st /= short_fft_success) error stop '3D cross status'
    if (maxval(abs(cross3-2.0_dp*abs(s3%values)**2)) > tol) error stop '3D cross values'

    ! A nonstandard time axis (1) and frequency axis (4) should preserve the
    ! ordering of the other dimensions after inverse reconstruction.
    allocate(reordered4(size(s3%values,4),3,2,size(s3%values,2)))
    do k=1,size(s3%values,4)
        do j=1,size(s3%values,2)
            reordered4(k,:,:,j) = s3%values(:,j,:,k)
        end do
    end do
    call plan%istft(reordered4,y3,st,k1=35,f_axis=4,t_axis=1)
    if (st /= short_fft_success) error stop 'nonstandard 4D inverse'
    if (any(shape(y3) /= [3,2,35])) error stop 'nonstandard 4D inverse shape'
    do k=1,2
        do j=1,35
            if (maxval(abs(y3(:,k,j)-x3(:,j,k))) > tol) error stop 'nonstandard 4D inverse values'
        end do
    end do

    call plan%initialize_window('hann',8,2,20.0_dp,mfft=16,fft_mode='centered')
    if (plan%status /= short_fft_success) error stop 'nd complex plan'
    call plan%stft(z2,s2,axis=1,padding='odd')
    if (s2%status /= short_fft_success) error stop 'complex 2D forward'
    call plan%istft(s2%values,zy2,st,k1=35,f_axis=1)
    if (st /= short_fft_success .or. maxval(abs(zy2-z2)) > tol) error stop 'complex 2D inverse'
    call plan%stft_detrend(z2,s2,'linear',axis=1,padding='even')
    if (s2%status /= short_fft_success) error stop 'complex 2D detrending'
    call plan%stft(z2,s2,axis=1,padding='odd')
    allocate(reordered3(size(s2%values,3),3,size(s2%values,1)))
    do k=1,size(s2%values,3)
        do j=1,size(s2%values,1)
            reordered3(k,:,j) = s2%values(j,:,k)
        end do
    end do
    call plan%istft(reordered3,zy2,st,k1=35,f_axis=3,t_axis=1)
    if (st /= short_fft_success) error stop 'nonstandard 3D inverse'
    if (any(shape(zy2) /= [3,35])) error stop 'nonstandard 3D inverse shape'
    if (maxval(abs(zy2-transpose(z2))) > tol) error stop 'nonstandard 3D inverse values'

    call plan%stft(z3,s3,axis=2,padding='edge')
    if (s3%status /= short_fft_success) error stop 'complex 3D forward'
    call plan%istft(s3%values,zy3,st,k1=35,f_axis=2)
    if (st /= short_fft_success) error stop 'complex 3D inverse'
    if (maxval(abs(zy3-z3)) > tol) error stop 'complex 3D roundtrip'
    call plan%stft_detrend(z3,s3,'linear',axis=2)
    if (s3%status /= short_fft_success) error stop 'complex 3D detrending'
    print *, 'test_signal_short_time_fft_nd: all checks passed'
end program test_signal_short_time_fft_nd
