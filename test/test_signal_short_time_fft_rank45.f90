! SPDX-License-Identifier: MIT
! High-rank ShortTimeFFT: full transforms, arbitrary axes, and invalid inputs.
program test_signal_short_time_fft_rank45
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, &
        short_time_fft_result_5d, short_time_fft_result_6d, &
        short_time_spectrogram_result_5d, short_time_spectrogram_result_6d, &
        short_fft_success, short_fft_invalid_input
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result_5d) :: s4
    type(short_time_fft_result_6d) :: s5
    type(short_time_spectrogram_result_5d) :: power4
    type(short_time_spectrogram_result_6d) :: power5
    real(dp) :: x4(2,35,3,2), x5(2,3,35,2,2)
    complex(dp) :: z4(2,35,3,2), z5(2,3,35,2,2)
    real(dp), allocatable :: y4(:,:,:,:), y5(:,:,:,:,:)
    complex(dp), allocatable :: w4(:,:,:,:), w5(:,:,:,:,:)
    complex(dp), allocatable :: cross4(:,:,:,:,:), cross5(:,:,:,:,:,:)
    complex(dp), allocatable :: rotated4(:,:,:,:,:), rotated5(:,:,:,:,:,:)
    integer :: a,b,c,d,e,f,g,status
    real(dp), parameter :: tol = 5.0e-12_dp

    do d=1,2
        do c=1,3
            do b=1,35
                do a=1,2
                    x4(a,b,c,d) = sin(0.11_dp*a+0.09_dp*b+0.07_dp*c+0.05_dp*d) + &
                        0.2_dp*cos(0.12_dp*a-0.19_dp*b+0.13_dp*d)
                    z4(a,b,c,d) = cmplx(x4(a,b,c,d), &
                        0.14_dp*sin(0.08_dp*a-0.05_dp*b+0.03_dp*c),dp)
                end do
            end do
        end do
    end do
    do e=1,2
        do d=1,2
            do c=1,35
                do b=1,3
                    do a=1,2
                        x5(a,b,c,d,e) = sin(0.11_dp*a+0.09_dp*b+0.07_dp*c+0.05_dp*d+0.03_dp*e) + &
                            0.2_dp*cos(0.12_dp*a-0.19_dp*c+0.13_dp*d)
                        z5(a,b,c,d,e) = cmplx(x5(a,b,c,d,e), &
                            0.14_dp*sin(0.08_dp*a-0.05_dp*c+0.03_dp*e),dp)
                    end do
                end do
            end do
        end do
    end do

    call plan%initialize_window('hann',9,3,20.0_dp,mfft=13, &
        fft_mode='onesided2X',scale_to='psd')
    if (plan%status /= short_fft_success) error stop 'high-rank plan'
    call plan%stft(x4,s4,axis=2,padding='even')
    if (s4%status /= short_fft_success) error stop 'rank4 real stft'
    if (any(shape(s4%values) /= [2,7,3,2,plan%p_num(35)])) error stop 'rank4 shape'
    call plan%istft(s4%values,y4,status,k1=35,f_axis=2)
    if (status /= short_fft_success) error stop 'rank4 real inverse'
    if (any(shape(y4) /= shape(x4))) error stop 'rank4 recovered shape'
    if (maxval(abs(y4-x4)) > tol) error stop 'rank4 real roundtrip'
    call plan%spectrogram(x4,power4,axis=2,padding='even')
    if (power4%status /= short_fft_success) error stop 'rank4 spectrum status'
    if (maxval(abs(power4%power-abs(s4%values)**2)) > tol) error stop 'rank4 spectrum values'
    call plan%cross_spectrogram(x4,2*x4,cross4,status,axis=2,padding='even')
    if (status /= short_fft_success) error stop 'rank4 cross status'
    if (maxval(abs(cross4-2*abs(s4%values)**2)) > tol) error stop 'rank4 cross values'
    call plan%stft_detrend(x4,s4,'linear',axis=2,padding='odd')
    if (s4%status /= short_fft_success) error stop 'rank4 detrend'

    call plan%stft(x5,s5,axis=3,padding='edge')
    if (s5%status /= short_fft_success) error stop 'rank5 real stft'
    if (any(shape(s5%values) /= [2,3,7,2,2,plan%p_num(35)])) error stop 'rank5 shape'
    call plan%istft(s5%values,y5,status,k1=35,f_axis=3)
    if (status /= short_fft_success) error stop 'rank5 real inverse'
    if (any(shape(y5) /= shape(x5))) error stop 'rank5 recovered shape'
    if (maxval(abs(y5-x5)) > tol) error stop 'rank5 real roundtrip'
    call plan%spectrogram(x5,power5,axis=3,padding='edge')
    if (power5%status /= short_fft_success) error stop 'rank5 spectrum status'
    if (maxval(abs(power5%power-abs(s5%values)**2)) > tol) error stop 'rank5 spectrum values'
    call plan%cross_spectrogram(x5,2*x5,cross5,status,axis=3,padding='edge')
    if (status /= short_fft_success) error stop 'rank5 cross status'
    if (maxval(abs(cross5-2*abs(s5%values)**2)) > tol) error stop 'rank5 cross values'
    call plan%stft_detrend(x5,s5,'constant',axis=3,padding='odd')
    if (s5%status /= short_fft_success) error stop 'rank5 detrend'
    call plan%stft(x4,s4,axis=0)
    if (s4%status /= short_fft_invalid_input) error stop 'rank4 axis zero'
    call plan%stft(x5,s5,axis=-6)
    if (s5%status /= short_fft_invalid_input) error stop 'rank5 axis invalid'
    call plan%stft(x5,s5,axis=3,padding='even')
    call plan%istft(s5%values,y5,status,k1=35,f_axis=1,t_axis=1)
    if (status /= short_fft_invalid_input) error stop 'rank5 duplicate inverse axes'

    ! Permute axes: time axis moves to 1; frequency becomes the final axis.
    ! The inverse must retain the ordering of remaining batch dimensions.
    call plan%stft(x4,s4,axis=2,padding='even')
    allocate(rotated4(size(s4%values,5),2,3,2,7))
    do e=1,size(s4%values,5)
        do f=1,7
            do d=1,2
                do c=1,3
                    do a=1,2
                        rotated4(e,a,c,d,f) = s4%values(a,f,c,d,e)
                    end do
                end do
            end do
        end do
    end do
    call plan%istft(rotated4,y4,status,k1=35,f_axis=5,t_axis=1)
    if (status /= short_fft_success) error stop 'rank4 permuted inverse'
    if (any(shape(y4) /= [2,3,2,35])) error stop 'rank4 permuted shape'
    do d=1,2
        do c=1,3
            do b=1,35
                if (maxval(abs(y4(:,c,d,b)-x4(:,b,c,d))) > tol) error stop 'rank4 permuted values'
            end do
        end do
    end do

    call plan%stft(x5,s5,axis=3,padding='even')
    allocate(rotated5(size(s5%values,6),2,3,2,2,7))
    do g=1,size(s5%values,6)
        do f=1,7
            do e=1,2
                do d=1,2
                    do b=1,3
                        do a=1,2
                            rotated5(g,a,b,d,e,f) = s5%values(a,b,f,d,e,g)
                        end do
                    end do
                end do
            end do
        end do
    end do
    call plan%istft(rotated5,y5,status,k1=35,f_axis=6,t_axis=1)
    if (status /= short_fft_success) error stop 'rank5 permuted inverse'
    if (any(shape(y5) /= [2,3,2,2,35])) error stop 'rank5 permuted shape'
    do e=1,2
        do d=1,2
            do c=1,35
                if (maxval(abs(y5(:,:,d,e,c)-x5(:,:,c,d,e))) > tol) error stop 'rank5 permuted values'
            end do
        end do
    end do

    call plan%initialize_window('hann',8,2,20.0_dp,mfft=16,fft_mode='centered')
    if (plan%status /= short_fft_success) error stop 'complex plan'
    call plan%stft(z4,s4,axis=2,padding='odd')
    call plan%istft(s4%values,w4,status,k1=35,f_axis=2)
    if (status /= short_fft_success) error stop 'complex rank4 inverse'
    if (maxval(abs(w4-z4)) > tol) error stop 'complex rank4 roundtrip'
    call plan%stft_detrend(z4,s4,'constant',axis=2,padding='zeros')
    if (s4%status /= short_fft_success) error stop 'complex rank4 detrend'
    call plan%stft(z5,s5,axis=-3,padding='odd')
    call plan%istft(s5%values,w5,status,k1=35,f_axis=3)
    if (status /= short_fft_success) error stop 'complex rank5 inverse'
    if (maxval(abs(w5-z5)) > tol) error stop 'complex rank5 roundtrip'
    call plan%stft_detrend(z5,s5,'linear',axis=3,padding='even')
    if (s5%status /= short_fft_success) error stop 'complex rank5 detrend'
    call plan%spectrogram(z5,power5,axis=3,padding='odd')
    if (power5%status /= short_fft_success) error stop 'complex rank5 spectrogram'
    call plan%cross_spectrogram(z4,z4,cross4,status,axis=2,padding='odd')
    if (status /= short_fft_success) error stop 'complex rank4 cross'
    print *, 'test_signal_short_time_fft_rank45: all checks passed'
end program test_signal_short_time_fft_rank45
