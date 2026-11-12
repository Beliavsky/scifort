! SPDX-License-Identifier: MIT
! Round-trip, spectrogram, and invalid-axis regression for input ranks 6--14.
program test_signal_short_time_fft_high_rank
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft
    implicit none
    call test_rank6()
    call test_rank7()
    call test_rank8()
    call test_rank9()
    call test_rank10()
    call test_rank11()
    call test_rank12()
    call test_rank13()
    call test_rank14()
    print *, 'test_signal_short_time_fft_high_rank: all checks passed'
contains
    subroutine test_rank6()
        type(short_time_fft) :: sft
        type(short_time_fft_result_7d) :: frame, reference
        type(short_time_spectrogram_result_7d) :: power
        integer :: dims(6), expected(7), axis, j, k, status
        integer :: src_shape(7), perm_shape(7)
        complex(dp), allocatable :: reordered(:,:,:,:,:,:,:)
        real(dp), allocatable :: x(:,:,:,:,:,:), y(:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:), w(:,:,:,:,:,:)
        complex(dp), allocatable :: cross(:,:,:,:,:,:,:)
        do axis=1,6
            dims = 1
            dims(axis) = 27
            j=1
            if (j == axis) j=2
            dims(j)=2
            j=6
            if (j == axis) j=5
            dims(j)=2
            allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)),z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
            x = reshape([(sin(0.11_dp*real(k,dp))+ &
                0.12_dp*cos(0.07_dp*real(k,dp)),k=1,product(dims))],dims)
            z = cmplx(x,0.2_dp*reshape([(sin(0.09_dp*real(k,dp)), &
                k=1,product(dims))],dims),dp)
            call sft%initialize_window('hann',9,3,21.0_dp,mfft=16, &
                fft_mode='onesided2X',scale_to='psd')
            if (sft%status /= short_fft_success) error stop 'high-rank real plan'
            call sft%stft(x,frame,axis=axis)
            if (frame%status /= short_fft_success) error stop 'high-rank real STFT'
            expected(1:6) = dims
            expected(axis) = sft%frequency_count()
            expected(7) = sft%p_num(27)
            if (any(shape(frame%values) /= expected)) error stop 'high-rank real shape'
            call sft%istft(frame%values,y,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank real iSTFT'
            if (any(shape(y) /= dims)) error stop 'high-rank real round-trip shape'
            if (maxval(abs(x-y)) > 2.0e-12_dp) error stop 'high-rank real round-trip'
            if (axis == 6) then
                src_shape = shape(frame%values)
                perm_shape(1) = src_shape(7)
                perm_shape(2:7) = src_shape(1:6)
                reordered = reshape(frame%values,perm_shape, &
                    order=[(j,j=2,7),1])
                call sft%istft(reordered,y,status,k1=27,f_axis=-1,t_axis=1)
                if (status /= short_fft_success) &
                    error stop 'high-rank reordered inverse status'
                if (maxval(abs(x-y)) > 2.0e-12_dp) &
                    error stop 'high-rank reordered inverse values'
            end if
            call sft%spectrogram(x,power,axis=axis)
            if (power%status /= short_fft_success) error stop 'high-rank real power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank power values'
            call sft%cross_spectrogram(x,2*x,cross,status,axis=axis)
            if (status /= short_fft_success) error stop 'high-rank cross status'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank cross values'
            call sft%stft_detrend(x,reference,'linear',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank real detrend'
            if (axis == 6) then
                call sft%stft(x,reference,axis=-1)
                if (reference%status /= short_fft_success) error stop 'high-rank negative axis'
                if (maxval(abs(reference%values-frame%values)) > 2.0e-12_dp) &
                    error stop 'high-rank negative axis values'
            end if
            call sft%initialize_window('hamming',9,3,21.0_dp,mfft=13,fft_mode='centered')
            if (sft%status /= short_fft_success) error stop 'high-rank complex plan'
            call sft%stft(z,frame,axis=axis,padding='even')
            if (frame%status /= short_fft_success) error stop 'high-rank complex STFT'
            call sft%istft(frame%values,w,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank complex iSTFT'
            if (any(shape(w) /= dims)) error stop 'high-rank complex round-trip shape'
            if (maxval(abs(z-w)) > 2.0e-12_dp) error stop 'high-rank complex round-trip'
            if (axis == 6) then
                src_shape = shape(frame%values)
                perm_shape(1) = src_shape(7)
                perm_shape(2:7) = src_shape(1:6)
                reordered = reshape(frame%values,perm_shape, &
                    order=[(j,j=2,7),1])
                call sft%istft(reordered,w,status,k1=27,f_axis=-1,t_axis=1)
                if (status /= short_fft_success) &
                    error stop 'high-rank reordered complex inverse status'
                if (maxval(abs(z-w)) > 2.0e-12_dp) &
                    error stop 'high-rank reordered complex inverse values'
            end if
            call sft%spectrogram(z,power,axis=axis,padding='even')
            if (power%status /= short_fft_success) error stop 'high-rank complex power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex power values'
            call sft%cross_spectrogram(z,2*z,cross,status,axis=axis,padding='even')
            if (status /= short_fft_success) error stop 'high-rank complex cross'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex cross values'
            call sft%stft_detrend(z,reference,'constant',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank complex detrend'
            deallocate(x,z)
        end do
    end subroutine test_rank6

    subroutine test_rank7()
        type(short_time_fft) :: sft
        type(short_time_fft_result_8d) :: frame, reference
        type(short_time_spectrogram_result_8d) :: power
        integer :: dims(7), expected(8), axis, j, k, status
        real(dp), allocatable :: x(:,:,:,:,:,:,:), y(:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:), w(:,:,:,:,:,:,:)
        complex(dp), allocatable :: cross(:,:,:,:,:,:,:,:)
        do axis=1,7
            dims = 1
            dims(axis) = 27
            j=1
            if (j == axis) j=2
            dims(j)=2
            j=7
            if (j == axis) j=6
            dims(j)=2
            allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)),z(dims(1),dims(2),dims(3), &
            dims(4),dims(5),dims(6),dims(7)))
            x = reshape([(sin(0.11_dp*real(k,dp))+ &
                0.12_dp*cos(0.07_dp*real(k,dp)),k=1,product(dims))],dims)
            z = cmplx(x,0.2_dp*reshape([(sin(0.09_dp*real(k,dp)), &
                k=1,product(dims))],dims),dp)
            call sft%initialize_window('hann',9,3,21.0_dp,mfft=16, &
                fft_mode='onesided2X',scale_to='psd')
            if (sft%status /= short_fft_success) error stop 'high-rank real plan'
            call sft%stft(x,frame,axis=axis)
            if (frame%status /= short_fft_success) error stop 'high-rank real STFT'
            expected(1:7) = dims
            expected(axis) = sft%frequency_count()
            expected(8) = sft%p_num(27)
            if (any(shape(frame%values) /= expected)) error stop 'high-rank real shape'
            call sft%istft(frame%values,y,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank real iSTFT'
            if (any(shape(y) /= dims)) error stop 'high-rank real round-trip shape'
            if (maxval(abs(x-y)) > 2.0e-12_dp) error stop 'high-rank real round-trip'
            call sft%spectrogram(x,power,axis=axis)
            if (power%status /= short_fft_success) error stop 'high-rank real power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank power values'
            call sft%cross_spectrogram(x,2*x,cross,status,axis=axis)
            if (status /= short_fft_success) error stop 'high-rank cross status'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank cross values'
            call sft%stft_detrend(x,reference,'linear',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank real detrend'
            if (axis == 7) then
                call sft%stft(x,reference,axis=-1)
                if (reference%status /= short_fft_success) error stop 'high-rank negative axis'
                if (maxval(abs(reference%values-frame%values)) > 2.0e-12_dp) &
                    error stop 'high-rank negative axis values'
            end if
            call sft%initialize_window('hamming',9,3,21.0_dp,mfft=13,fft_mode='centered')
            if (sft%status /= short_fft_success) error stop 'high-rank complex plan'
            call sft%stft(z,frame,axis=axis,padding='even')
            if (frame%status /= short_fft_success) error stop 'high-rank complex STFT'
            call sft%istft(frame%values,w,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank complex iSTFT'
            if (any(shape(w) /= dims)) error stop 'high-rank complex round-trip shape'
            if (maxval(abs(z-w)) > 2.0e-12_dp) error stop 'high-rank complex round-trip'
            call sft%spectrogram(z,power,axis=axis,padding='even')
            if (power%status /= short_fft_success) error stop 'high-rank complex power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex power values'
            call sft%cross_spectrogram(z,2*z,cross,status,axis=axis,padding='even')
            if (status /= short_fft_success) error stop 'high-rank complex cross'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex cross values'
            call sft%stft_detrend(z,reference,'constant',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank complex detrend'
            deallocate(x,z)
        end do
    end subroutine test_rank7

    subroutine test_rank8()
        type(short_time_fft) :: sft
        type(short_time_fft_result_9d) :: frame, reference
        type(short_time_spectrogram_result_9d) :: power
        integer :: dims(8), expected(9), axis, j, k, status
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:), y(:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:), w(:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: cross(:,:,:,:,:,:,:,:,:)
        do axis=1,8
            dims = 1
            dims(axis) = 27
            j=1
            if (j == axis) j=2
            dims(j)=2
            j=8
            if (j == axis) j=7
            dims(j)=2
            allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8)),z(dims(1),dims(2), &
            dims(3),dims(4),dims(5),dims(6),dims(7),dims(8)))
            x = reshape([(sin(0.11_dp*real(k,dp))+ &
                0.12_dp*cos(0.07_dp*real(k,dp)),k=1,product(dims))],dims)
            z = cmplx(x,0.2_dp*reshape([(sin(0.09_dp*real(k,dp)), &
                k=1,product(dims))],dims),dp)
            call sft%initialize_window('hann',9,3,21.0_dp,mfft=16, &
                fft_mode='onesided2X',scale_to='psd')
            if (sft%status /= short_fft_success) error stop 'high-rank real plan'
            call sft%stft(x,frame,axis=axis)
            if (frame%status /= short_fft_success) error stop 'high-rank real STFT'
            expected(1:8) = dims
            expected(axis) = sft%frequency_count()
            expected(9) = sft%p_num(27)
            if (any(shape(frame%values) /= expected)) error stop 'high-rank real shape'
            call sft%istft(frame%values,y,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank real iSTFT'
            if (any(shape(y) /= dims)) error stop 'high-rank real round-trip shape'
            if (maxval(abs(x-y)) > 2.0e-12_dp) error stop 'high-rank real round-trip'
            call sft%spectrogram(x,power,axis=axis)
            if (power%status /= short_fft_success) error stop 'high-rank real power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank power values'
            call sft%cross_spectrogram(x,2*x,cross,status,axis=axis)
            if (status /= short_fft_success) error stop 'high-rank cross status'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank cross values'
            call sft%stft_detrend(x,reference,'linear',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank real detrend'
            if (axis == 8) then
                call sft%stft(x,reference,axis=-1)
                if (reference%status /= short_fft_success) error stop 'high-rank negative axis'
                if (maxval(abs(reference%values-frame%values)) > 2.0e-12_dp) &
                    error stop 'high-rank negative axis values'
            end if
            call sft%initialize_window('hamming',9,3,21.0_dp,mfft=13,fft_mode='centered')
            if (sft%status /= short_fft_success) error stop 'high-rank complex plan'
            call sft%stft(z,frame,axis=axis,padding='even')
            if (frame%status /= short_fft_success) error stop 'high-rank complex STFT'
            call sft%istft(frame%values,w,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank complex iSTFT'
            if (any(shape(w) /= dims)) error stop 'high-rank complex round-trip shape'
            if (maxval(abs(z-w)) > 2.0e-12_dp) error stop 'high-rank complex round-trip'
            call sft%spectrogram(z,power,axis=axis,padding='even')
            if (power%status /= short_fft_success) error stop 'high-rank complex power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex power values'
            call sft%cross_spectrogram(z,2*z,cross,status,axis=axis,padding='even')
            if (status /= short_fft_success) error stop 'high-rank complex cross'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex cross values'
            call sft%stft_detrend(z,reference,'constant',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank complex detrend'
            deallocate(x,z)
        end do
    end subroutine test_rank8

    subroutine test_rank9()
        type(short_time_fft) :: sft
        type(short_time_fft_result_10d) :: frame, reference
        type(short_time_spectrogram_result_10d) :: power
        integer :: dims(9), expected(10), axis, j, k, status
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:), y(:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:), w(:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: cross(:,:,:,:,:,:,:,:,:,:)
        do axis=1,9
            dims = 1
            dims(axis) = 27
            j=1
            if (j == axis) j=2
            dims(j)=2
            j=9
            if (j == axis) j=8
            dims(j)=2
            allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9)),z(dims(1), &
            dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9)))
            x = reshape([(sin(0.11_dp*real(k,dp))+ &
                0.12_dp*cos(0.07_dp*real(k,dp)),k=1,product(dims))],dims)
            z = cmplx(x,0.2_dp*reshape([(sin(0.09_dp*real(k,dp)), &
                k=1,product(dims))],dims),dp)
            call sft%initialize_window('hann',9,3,21.0_dp,mfft=16, &
                fft_mode='onesided2X',scale_to='psd')
            if (sft%status /= short_fft_success) error stop 'high-rank real plan'
            call sft%stft(x,frame,axis=axis)
            if (frame%status /= short_fft_success) error stop 'high-rank real STFT'
            expected(1:9) = dims
            expected(axis) = sft%frequency_count()
            expected(10) = sft%p_num(27)
            if (any(shape(frame%values) /= expected)) error stop 'high-rank real shape'
            call sft%istft(frame%values,y,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank real iSTFT'
            if (any(shape(y) /= dims)) error stop 'high-rank real round-trip shape'
            if (maxval(abs(x-y)) > 2.0e-12_dp) error stop 'high-rank real round-trip'
            call sft%spectrogram(x,power,axis=axis)
            if (power%status /= short_fft_success) error stop 'high-rank real power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank power values'
            call sft%cross_spectrogram(x,2*x,cross,status,axis=axis)
            if (status /= short_fft_success) error stop 'high-rank cross status'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank cross values'
            call sft%stft_detrend(x,reference,'linear',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank real detrend'
            if (axis == 9) then
                call sft%stft(x,reference,axis=-1)
                if (reference%status /= short_fft_success) error stop 'high-rank negative axis'
                if (maxval(abs(reference%values-frame%values)) > 2.0e-12_dp) &
                    error stop 'high-rank negative axis values'
            end if
            call sft%initialize_window('hamming',9,3,21.0_dp,mfft=13,fft_mode='centered')
            if (sft%status /= short_fft_success) error stop 'high-rank complex plan'
            call sft%stft(z,frame,axis=axis,padding='even')
            if (frame%status /= short_fft_success) error stop 'high-rank complex STFT'
            call sft%istft(frame%values,w,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank complex iSTFT'
            if (any(shape(w) /= dims)) error stop 'high-rank complex round-trip shape'
            if (maxval(abs(z-w)) > 2.0e-12_dp) error stop 'high-rank complex round-trip'
            call sft%spectrogram(z,power,axis=axis,padding='even')
            if (power%status /= short_fft_success) error stop 'high-rank complex power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex power values'
            call sft%cross_spectrogram(z,2*z,cross,status,axis=axis,padding='even')
            if (status /= short_fft_success) error stop 'high-rank complex cross'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex cross values'
            call sft%stft_detrend(z,reference,'constant',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank complex detrend'
            deallocate(x,z)
        end do
    end subroutine test_rank9

    subroutine test_rank10()
        type(short_time_fft) :: sft
        type(short_time_fft_result_11d) :: frame, reference
        type(short_time_spectrogram_result_11d) :: power
        integer :: dims(10), expected(11), axis, j, k, status
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:), y(:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:), w(:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: cross(:,:,:,:,:,:,:,:,:,:,:)
        do axis=1,10
            dims = 1
            dims(axis) = 27
            j=1
            if (j == axis) j=2
            dims(j)=2
            j=10
            if (j == axis) j=9
            dims(j)=2
            allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10)), &
            z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10)))
            x = reshape([(sin(0.11_dp*real(k,dp))+ &
                0.12_dp*cos(0.07_dp*real(k,dp)),k=1,product(dims))],dims)
            z = cmplx(x,0.2_dp*reshape([(sin(0.09_dp*real(k,dp)), &
                k=1,product(dims))],dims),dp)
            call sft%initialize_window('hann',9,3,21.0_dp,mfft=16, &
                fft_mode='onesided2X',scale_to='psd')
            if (sft%status /= short_fft_success) error stop 'high-rank real plan'
            call sft%stft(x,frame,axis=axis)
            if (frame%status /= short_fft_success) error stop 'high-rank real STFT'
            expected(1:10) = dims
            expected(axis) = sft%frequency_count()
            expected(11) = sft%p_num(27)
            if (any(shape(frame%values) /= expected)) error stop 'high-rank real shape'
            call sft%istft(frame%values,y,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank real iSTFT'
            if (any(shape(y) /= dims)) error stop 'high-rank real round-trip shape'
            if (maxval(abs(x-y)) > 2.0e-12_dp) error stop 'high-rank real round-trip'
            call sft%spectrogram(x,power,axis=axis)
            if (power%status /= short_fft_success) error stop 'high-rank real power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank power values'
            call sft%cross_spectrogram(x,2*x,cross,status,axis=axis)
            if (status /= short_fft_success) error stop 'high-rank cross status'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank cross values'
            call sft%stft_detrend(x,reference,'linear',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank real detrend'
            if (axis == 10) then
                call sft%stft(x,reference,axis=-1)
                if (reference%status /= short_fft_success) error stop 'high-rank negative axis'
                if (maxval(abs(reference%values-frame%values)) > 2.0e-12_dp) &
                    error stop 'high-rank negative axis values'
            end if
            call sft%initialize_window('hamming',9,3,21.0_dp,mfft=13,fft_mode='centered')
            if (sft%status /= short_fft_success) error stop 'high-rank complex plan'
            call sft%stft(z,frame,axis=axis,padding='even')
            if (frame%status /= short_fft_success) error stop 'high-rank complex STFT'
            call sft%istft(frame%values,w,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank complex iSTFT'
            if (any(shape(w) /= dims)) error stop 'high-rank complex round-trip shape'
            if (maxval(abs(z-w)) > 2.0e-12_dp) error stop 'high-rank complex round-trip'
            call sft%spectrogram(z,power,axis=axis,padding='even')
            if (power%status /= short_fft_success) error stop 'high-rank complex power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex power values'
            call sft%cross_spectrogram(z,2*z,cross,status,axis=axis,padding='even')
            if (status /= short_fft_success) error stop 'high-rank complex cross'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex cross values'
            call sft%stft_detrend(z,reference,'constant',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank complex detrend'
            deallocate(x,z)
        end do
    end subroutine test_rank10

    subroutine test_rank11()
        type(short_time_fft) :: sft
        type(short_time_fft_result_12d) :: frame, reference
        type(short_time_spectrogram_result_12d) :: power
        integer :: dims(11), expected(12), axis, j, k, status
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:,:), y(:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:), w(:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: cross(:,:,:,:,:,:,:,:,:,:,:,:)
        do axis=1,11
            dims = 1
            dims(axis) = 27
            j=1
            if (j == axis) j=2
            dims(j)=2
            j=11
            if (j == axis) j=10
            dims(j)=2
            allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10), &
            dims(11)),z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11)))
            x = reshape([(sin(0.11_dp*real(k,dp))+ &
                0.12_dp*cos(0.07_dp*real(k,dp)),k=1,product(dims))],dims)
            z = cmplx(x,0.2_dp*reshape([(sin(0.09_dp*real(k,dp)), &
                k=1,product(dims))],dims),dp)
            call sft%initialize_window('hann',9,3,21.0_dp,mfft=16, &
                fft_mode='onesided2X',scale_to='psd')
            if (sft%status /= short_fft_success) error stop 'high-rank real plan'
            call sft%stft(x,frame,axis=axis)
            if (frame%status /= short_fft_success) error stop 'high-rank real STFT'
            expected(1:11) = dims
            expected(axis) = sft%frequency_count()
            expected(12) = sft%p_num(27)
            if (any(shape(frame%values) /= expected)) error stop 'high-rank real shape'
            call sft%istft(frame%values,y,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank real iSTFT'
            if (any(shape(y) /= dims)) error stop 'high-rank real round-trip shape'
            if (maxval(abs(x-y)) > 2.0e-12_dp) error stop 'high-rank real round-trip'
            call sft%spectrogram(x,power,axis=axis)
            if (power%status /= short_fft_success) error stop 'high-rank real power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank power values'
            call sft%cross_spectrogram(x,2*x,cross,status,axis=axis)
            if (status /= short_fft_success) error stop 'high-rank cross status'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank cross values'
            call sft%stft_detrend(x,reference,'linear',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank real detrend'
            if (axis == 11) then
                call sft%stft(x,reference,axis=-1)
                if (reference%status /= short_fft_success) error stop 'high-rank negative axis'
                if (maxval(abs(reference%values-frame%values)) > 2.0e-12_dp) &
                    error stop 'high-rank negative axis values'
            end if
            call sft%initialize_window('hamming',9,3,21.0_dp,mfft=13,fft_mode='centered')
            if (sft%status /= short_fft_success) error stop 'high-rank complex plan'
            call sft%stft(z,frame,axis=axis,padding='even')
            if (frame%status /= short_fft_success) error stop 'high-rank complex STFT'
            call sft%istft(frame%values,w,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank complex iSTFT'
            if (any(shape(w) /= dims)) error stop 'high-rank complex round-trip shape'
            if (maxval(abs(z-w)) > 2.0e-12_dp) error stop 'high-rank complex round-trip'
            call sft%spectrogram(z,power,axis=axis,padding='even')
            if (power%status /= short_fft_success) error stop 'high-rank complex power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex power values'
            call sft%cross_spectrogram(z,2*z,cross,status,axis=axis,padding='even')
            if (status /= short_fft_success) error stop 'high-rank complex cross'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex cross values'
            call sft%stft_detrend(z,reference,'constant',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank complex detrend'
            deallocate(x,z)
        end do
    end subroutine test_rank11

    subroutine test_rank12()
        type(short_time_fft) :: sft
        type(short_time_fft_result_13d) :: frame, reference
        type(short_time_spectrogram_result_13d) :: power
        integer :: dims(12), expected(13), axis, j, k, status
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:,:,:), y(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:), w(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: cross(:,:,:,:,:,:,:,:,:,:,:,:,:)
        do axis=1,12
            dims = 1
            dims(axis) = 27
            j=1
            if (j == axis) j=2
            dims(j)=2
            j=12
            if (j == axis) j=11
            dims(j)=2
            allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10), &
            dims(11),dims(12)),z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9), &
            dims(10),dims(11),dims(12)))
            x = reshape([(sin(0.11_dp*real(k,dp))+ &
                0.12_dp*cos(0.07_dp*real(k,dp)),k=1,product(dims))],dims)
            z = cmplx(x,0.2_dp*reshape([(sin(0.09_dp*real(k,dp)), &
                k=1,product(dims))],dims),dp)
            call sft%initialize_window('hann',9,3,21.0_dp,mfft=16, &
                fft_mode='onesided2X',scale_to='psd')
            if (sft%status /= short_fft_success) error stop 'high-rank real plan'
            call sft%stft(x,frame,axis=axis)
            if (frame%status /= short_fft_success) error stop 'high-rank real STFT'
            expected(1:12) = dims
            expected(axis) = sft%frequency_count()
            expected(13) = sft%p_num(27)
            if (any(shape(frame%values) /= expected)) error stop 'high-rank real shape'
            call sft%istft(frame%values,y,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank real iSTFT'
            if (any(shape(y) /= dims)) error stop 'high-rank real round-trip shape'
            if (maxval(abs(x-y)) > 2.0e-12_dp) error stop 'high-rank real round-trip'
            call sft%spectrogram(x,power,axis=axis)
            if (power%status /= short_fft_success) error stop 'high-rank real power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank power values'
            call sft%cross_spectrogram(x,2*x,cross,status,axis=axis)
            if (status /= short_fft_success) error stop 'high-rank cross status'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank cross values'
            call sft%stft_detrend(x,reference,'linear',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank real detrend'
            if (axis == 12) then
                call sft%stft(x,reference,axis=-1)
                if (reference%status /= short_fft_success) error stop 'high-rank negative axis'
                if (maxval(abs(reference%values-frame%values)) > 2.0e-12_dp) &
                    error stop 'high-rank negative axis values'
            end if
            call sft%initialize_window('hamming',9,3,21.0_dp,mfft=13,fft_mode='centered')
            if (sft%status /= short_fft_success) error stop 'high-rank complex plan'
            call sft%stft(z,frame,axis=axis,padding='even')
            if (frame%status /= short_fft_success) error stop 'high-rank complex STFT'
            call sft%istft(frame%values,w,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank complex iSTFT'
            if (any(shape(w) /= dims)) error stop 'high-rank complex round-trip shape'
            if (maxval(abs(z-w)) > 2.0e-12_dp) error stop 'high-rank complex round-trip'
            call sft%spectrogram(z,power,axis=axis,padding='even')
            if (power%status /= short_fft_success) error stop 'high-rank complex power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex power values'
            call sft%cross_spectrogram(z,2*z,cross,status,axis=axis,padding='even')
            if (status /= short_fft_success) error stop 'high-rank complex cross'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex cross values'
            call sft%stft_detrend(z,reference,'constant',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank complex detrend'
            deallocate(x,z)
        end do
    end subroutine test_rank12

    subroutine test_rank13()
        type(short_time_fft) :: sft
        type(short_time_fft_result_14d) :: frame, reference
        type(short_time_spectrogram_result_14d) :: power
        integer :: dims(13), expected(14), axis, j, k, status
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:,:,:,:), y(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:,:), w(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: cross(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        do axis=1,13
            dims = 1
            dims(axis) = 27
            j=1
            if (j == axis) j=2
            dims(j)=2
            j=13
            if (j == axis) j=12
            dims(j)=2
            allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10), &
            dims(11),dims(12),dims(13)),z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8), &
            dims(9),dims(10),dims(11),dims(12),dims(13)))
            x = reshape([(sin(0.11_dp*real(k,dp))+ &
                0.12_dp*cos(0.07_dp*real(k,dp)),k=1,product(dims))],dims)
            z = cmplx(x,0.2_dp*reshape([(sin(0.09_dp*real(k,dp)), &
                k=1,product(dims))],dims),dp)
            call sft%initialize_window('hann',9,3,21.0_dp,mfft=16, &
                fft_mode='onesided2X',scale_to='psd')
            if (sft%status /= short_fft_success) error stop 'high-rank real plan'
            call sft%stft(x,frame,axis=axis)
            if (frame%status /= short_fft_success) error stop 'high-rank real STFT'
            expected(1:13) = dims
            expected(axis) = sft%frequency_count()
            expected(14) = sft%p_num(27)
            if (any(shape(frame%values) /= expected)) error stop 'high-rank real shape'
            call sft%istft(frame%values,y,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank real iSTFT'
            if (any(shape(y) /= dims)) error stop 'high-rank real round-trip shape'
            if (maxval(abs(x-y)) > 2.0e-12_dp) error stop 'high-rank real round-trip'
            call sft%spectrogram(x,power,axis=axis)
            if (power%status /= short_fft_success) error stop 'high-rank real power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank power values'
            call sft%cross_spectrogram(x,2*x,cross,status,axis=axis)
            if (status /= short_fft_success) error stop 'high-rank cross status'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank cross values'
            call sft%stft_detrend(x,reference,'linear',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank real detrend'
            if (axis == 13) then
                call sft%stft(x,reference,axis=-1)
                if (reference%status /= short_fft_success) error stop 'high-rank negative axis'
                if (maxval(abs(reference%values-frame%values)) > 2.0e-12_dp) &
                    error stop 'high-rank negative axis values'
            end if
            call sft%initialize_window('hamming',9,3,21.0_dp,mfft=13,fft_mode='centered')
            if (sft%status /= short_fft_success) error stop 'high-rank complex plan'
            call sft%stft(z,frame,axis=axis,padding='even')
            if (frame%status /= short_fft_success) error stop 'high-rank complex STFT'
            call sft%istft(frame%values,w,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank complex iSTFT'
            if (any(shape(w) /= dims)) error stop 'high-rank complex round-trip shape'
            if (maxval(abs(z-w)) > 2.0e-12_dp) error stop 'high-rank complex round-trip'
            call sft%spectrogram(z,power,axis=axis,padding='even')
            if (power%status /= short_fft_success) error stop 'high-rank complex power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex power values'
            call sft%cross_spectrogram(z,2*z,cross,status,axis=axis,padding='even')
            if (status /= short_fft_success) error stop 'high-rank complex cross'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex cross values'
            call sft%stft_detrend(z,reference,'constant',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank complex detrend'
            deallocate(x,z)
        end do
    end subroutine test_rank13

    subroutine test_rank14()
        type(short_time_fft) :: sft
        type(short_time_fft_result_15d) :: frame, reference
        type(short_time_spectrogram_result_15d) :: power
        integer :: dims(14), expected(15), axis, j, k, status
        integer :: src_shape(15), perm_shape(15)
        complex(dp), allocatable :: reordered(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:), y(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:,:,:), w(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: cross(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        do axis=1,14
            dims = 1
            dims(axis) = 27
            j=1
            if (j == axis) j=2
            dims(j)=2
            j=14
            if (j == axis) j=13
            dims(j)=2
            allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10), &
            dims(11),dims(12),dims(13),dims(14)),z(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7), &
            dims(8),dims(9),dims(10),dims(11),dims(12),dims(13),dims(14)))
            x = reshape([(sin(0.11_dp*real(k,dp))+ &
                0.12_dp*cos(0.07_dp*real(k,dp)),k=1,product(dims))],dims)
            z = cmplx(x,0.2_dp*reshape([(sin(0.09_dp*real(k,dp)), &
                k=1,product(dims))],dims),dp)
            call sft%initialize_window('hann',9,3,21.0_dp,mfft=16, &
                fft_mode='onesided2X',scale_to='psd')
            if (sft%status /= short_fft_success) error stop 'high-rank real plan'
            call sft%stft(x,frame,axis=axis)
            if (frame%status /= short_fft_success) error stop 'high-rank real STFT'
            expected(1:14) = dims
            expected(axis) = sft%frequency_count()
            expected(15) = sft%p_num(27)
            if (any(shape(frame%values) /= expected)) error stop 'high-rank real shape'
            call sft%istft(frame%values,y,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank real iSTFT'
            if (any(shape(y) /= dims)) error stop 'high-rank real round-trip shape'
            if (maxval(abs(x-y)) > 2.0e-12_dp) error stop 'high-rank real round-trip'
            if (axis == 14) then
                src_shape = shape(frame%values)
                perm_shape(1) = src_shape(15)
                perm_shape(2:15) = src_shape(1:14)
                reordered = reshape(frame%values,perm_shape, &
                    order=[(j,j=2,15),1])
                call sft%istft(reordered,y,status,k1=27,f_axis=-1,t_axis=1)
                if (status /= short_fft_success) &
                    error stop 'high-rank reordered inverse status'
                if (maxval(abs(x-y)) > 2.0e-12_dp) &
                    error stop 'high-rank reordered inverse values'
            end if
            call sft%spectrogram(x,power,axis=axis)
            if (power%status /= short_fft_success) error stop 'high-rank real power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank power values'
            call sft%cross_spectrogram(x,2*x,cross,status,axis=axis)
            if (status /= short_fft_success) error stop 'high-rank cross status'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank cross values'
            call sft%stft_detrend(x,reference,'linear',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank real detrend'
            if (axis == 14) then
                call sft%stft(x,reference,axis=-1)
                if (reference%status /= short_fft_success) error stop 'high-rank negative axis'
                if (maxval(abs(reference%values-frame%values)) > 2.0e-12_dp) &
                    error stop 'high-rank negative axis values'
            end if
            call sft%initialize_window('hamming',9,3,21.0_dp,mfft=13,fft_mode='centered')
            if (sft%status /= short_fft_success) error stop 'high-rank complex plan'
            call sft%stft(z,frame,axis=axis,padding='even')
            if (frame%status /= short_fft_success) error stop 'high-rank complex STFT'
            call sft%istft(frame%values,w,status,k1=27,f_axis=axis)
            if (status /= short_fft_success) error stop 'high-rank complex iSTFT'
            if (any(shape(w) /= dims)) error stop 'high-rank complex round-trip shape'
            if (maxval(abs(z-w)) > 2.0e-12_dp) error stop 'high-rank complex round-trip'
            if (axis == 14) then
                src_shape = shape(frame%values)
                perm_shape(1) = src_shape(15)
                perm_shape(2:15) = src_shape(1:14)
                reordered = reshape(frame%values,perm_shape, &
                    order=[(j,j=2,15),1])
                call sft%istft(reordered,w,status,k1=27,f_axis=-1,t_axis=1)
                if (status /= short_fft_success) &
                    error stop 'high-rank reordered complex inverse status'
                if (maxval(abs(z-w)) > 2.0e-12_dp) &
                    error stop 'high-rank reordered complex inverse values'
            end if
            call sft%spectrogram(z,power,axis=axis,padding='even')
            if (power%status /= short_fft_success) error stop 'high-rank complex power'
            if (maxval(abs(power%power-abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex power values'
            call sft%cross_spectrogram(z,2*z,cross,status,axis=axis,padding='even')
            if (status /= short_fft_success) error stop 'high-rank complex cross'
            if (maxval(abs(cross-2*abs(frame%values)**2)) > 2.0e-12_dp) &
                error stop 'high-rank complex cross values'
            call sft%stft_detrend(z,reference,'constant',axis=axis)
            if (reference%status /= short_fft_success) error stop 'high-rank complex detrend'
            deallocate(x,z)
        end do
    end subroutine test_rank14

end program test_signal_short_time_fft_high_rank
