! SPDX-License-Identifier: MIT
! Callable ShortTimeFFT detrending for real/complex one-dimensional signals.
program test_signal_stft_callable_detrend
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result, &
        short_time_fft_result_3d, short_time_fft_result_4d, &
        short_fft_success, short_fft_invalid_input
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result) :: out, ref
    type(short_time_fft_result_3d) :: nd2
    type(short_time_fft_result_4d) :: nd3
    real(dp) :: x(58), unchanged(58), xr2(2,58), xr3(2,58,3)
    complex(dp) :: z(58), zc2(58,2), zc3(2,3,58)
    integer :: i, j, k, calls

    do i=1,size(x)
        x(i) = sin(0.15_dp*i)+0.2_dp*cos(0.4_dp*i)+0.003_dp*i*i
        z(i) = cmplx(x(i),0.1_dp*cos(0.3_dp*i),dp)
    end do
    unchanged = x
    call plan%initialize_window('hann',9,3,20.0_dp,mfft=15,fft_mode='twosided')
    if (plan%status /= short_fft_success) error stop 'initialization'
    call plan%stft_detrend(x,out,detrend=identity_real)
    if (out%status /= short_fft_success) error stop 'real callback'
    call plan%stft(x,ref)
    if (ref%status /= short_fft_success) error stop 'STFT'
    if (maxval(abs(out%values-ref%values)) > 1.0e-12_dp) error stop 'identity callback mismatch'
    if (any(x /= unchanged)) error stop 'callback changed input'

    call plan%stft_detrend(z,out,detrend=identity_complex, padding='odd')
    if (out%status /= short_fft_success) error stop 'complex callback'
    call plan%stft(z,ref,padding='odd')
    if (maxval(abs(out%values-ref%values)) > 1.0e-12_dp) error stop 'complex identity mismatch'

    call plan%stft_detrend(x(2:58:2),out,detrend=nonlinear_real,padding='even')
    if (out%status /= short_fft_success) error stop 'strided input'
    if (size(out%values,2) /= plan%p_num(29)) error stop 'strided shape'

    call plan%stft_detrend(x,out,detrend=bad_real)
    if (out%status /= short_fft_invalid_input) error stop 'nonfinite real callback'
    if (allocated(out%values) .or. allocated(out%times) .or. &
        allocated(out%frequencies)) error stop 'partial bad real callback'
    call plan%stft_detrend(z,out,detrend=bad_complex)
    if (out%status /= short_fft_invalid_input) error stop 'nonfinite complex callback'
    if (allocated(out%values)) error stop 'partial bad complex callback'

    ! Callable detrending along arbitrary pre-existing rank-2/3 axes.
    do j=1,2
        xr2(j,:) = real(j,dp)*x
        zc2(:,j) = z * real(j,dp)
        do k=1,3
            xr3(j,:,k) = x + 0.2_dp*real(j+k,dp)
            zc3(j,k,:) = z + cmplx(0.01_dp*j,0.02_dp*k,dp)
        end do
    end do
    call plan%stft_detrend(xr2,nd2,detrend=nonlinear_real,axis=-1)
    if (nd2%status /= short_fft_success) error stop 'rank-2 real callback'
    do j=1,2
        call plan%stft_detrend(xr2(j,:),ref,detrend=nonlinear_real)
        if (maxval(abs(nd2%values(j,:,:)-ref%values)) > 1.0e-12_dp) &
            error stop 'rank-2 real channel'
    end do
    call plan%stft_detrend(zc2,nd2,detrend=nonlinear_complex,axis=1)
    if (nd2%status /= short_fft_success) error stop 'rank-2 complex callback'
    do j=1,2
        call plan%stft_detrend(zc2(:,j),ref,detrend=nonlinear_complex)
        if (maxval(abs(nd2%values(:,j,:)-ref%values)) > 1.0e-12_dp) &
            error stop 'rank-2 complex channel'
    end do
    call plan%stft_detrend(xr3,nd3,detrend=nonlinear_real,axis=-2,padding='even')
    if (nd3%status /= short_fft_success) error stop 'rank-3 real callback'
    do k=1,3
        do j=1,2
            call plan%stft_detrend(xr3(j,:,k),ref,detrend=nonlinear_real,padding='even')
            if (maxval(abs(nd3%values(j,:,k,:)-ref%values)) > 1.0e-12_dp) &
                error stop 'rank-3 real channel'
        end do
    end do
    call plan%stft_detrend(zc3,nd3,detrend=nonlinear_complex)
    if (nd3%status /= short_fft_success) error stop 'rank-3 complex callback'
    do k=1,3
        do j=1,2
            call plan%stft_detrend(zc3(j,k,:),ref,detrend=nonlinear_complex)
            if (maxval(abs(nd3%values(j,k,:,:)-ref%values)) > 1.0e-12_dp) &
                error stop 'rank-3 complex channel'
        end do
    end do
    calls=0
    call plan%stft_detrend(xr2,nd2,detrend=bad_second_real)
    if (nd2%status /= short_fft_invalid_input) error stop 'rank-2 invalid callback'
    if (allocated(nd2%values)) error stop 'partial rank-2 result'

    call plan%stft_detrend(x,out,detrend='linear')
    if (out%status /= short_fft_success) error stop 'existing string API'
    call plan%stft_detrend(x,out,detrend=nonlinear_real)
    if (out%status /= short_fft_success) error stop 'callback recovery after invalid'
    do k=1,size(out%values,2)
        if (any(abs(out%values(:,k)) > huge(1.0_dp))) error stop 'nonfinite coefficient'
    end do
    print *, 'test_signal_stft_callable_detrend: all checks passed'
contains
    subroutine identity_real(samples)
        real(dp), intent(inout) :: samples(:) !! caller-provided real time slice
    end subroutine identity_real
    subroutine identity_complex(samples)
        complex(dp), intent(inout) :: samples(:) !! caller-provided complex time slice
    end subroutine identity_complex
    subroutine nonlinear_real(samples)
        real(dp), intent(inout) :: samples(:) !! real time slice
        samples = samples - sum(samples)/real(size(samples),dp) + 0.07_dp*samples**2
    end subroutine nonlinear_real
    subroutine nonlinear_complex(samples)
        complex(dp), intent(inout) :: samples(:) !! complex time slice
        samples = samples - sum(samples)/real(size(samples),dp) + 0.07_dp*samples**2 &
            + cmplx(0.0_dp,0.04_dp,dp)*samples
    end subroutine nonlinear_complex
    subroutine bad_second_real(samples)
        real(dp), intent(inout) :: samples(:) !! second callback makes invalid result
        calls=calls+1
        if (calls==2) samples(1)=ieee_value(0.0_dp,ieee_quiet_nan)
    end subroutine bad_second_real
    subroutine bad_real(samples)
        real(dp), intent(inout) :: samples(:) !! real time slice
        samples(1) = ieee_value(0.0_dp,ieee_quiet_nan)
    end subroutine bad_real
    subroutine bad_complex(samples)
        complex(dp), intent(inout) :: samples(:) !! complex time slice
        samples(1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    end subroutine bad_complex
end program test_signal_stft_callable_detrend
