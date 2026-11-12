! SPDX-License-Identifier: MIT
! Callback frequency weights for Fourier resampling, including ND/strided APIs.
program test_signal_resample_callable_window
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : resample, resample_into, &
        fourier_resample_status_success, fourier_resample_status_invalid_input
    implicit none
    real(dp), parameter :: tol = 6.0e-12_dp
    real(dp) :: x(7), x8(8), x2(7,3), out2(11,3), ref2(11,3)
    real(dp) :: x3(2,7,2), out3(2,11,2), buf(2,22,2), snapshot(2,22,2)
    real(dp) :: x15(1,1,1,1,1,1,1,7,1,1,1,1,1,1,1)
    real(dp) :: out15(1,1,1,1,1,1,1,11,1,1,1,1,1,1,1)
    real(dp), allocatable :: y(:), weights(:), t(:), rr2(:,:), rr3(:,:,:)
    complex(dp) :: z(8), z2(8,2), y2(5,2), f2(8,2), refc(5,2)
    complex(dp) :: z3(2,8,2), outc3(2,5,2)
    complex(dp) :: z15(1,1,8,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: outc15(1,1,5,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp), allocatable :: yc(:), cweights(:), cr3(:,:,:)
    integer :: i, j, k, st, calls, expected_n
    logical :: reject

    x = [1.0_dp,2.0_dp,-1.0_dp,4.0_dp,0.0_dp,1.0_dp,3.0_dp]
    x8 = [x,0.25_dp]
    do i = 1,8
        z(i) = cmplx(cos(real(i,dp)),sin(0.5_dp*real(i,dp)),dp)
    end do
    calls = 0
    expected_n = 7
    reject = .false.
    weights = real_response(fftfreq(7))
    calls = 0
    y = resample(x,11,window_function=real_response,status=st)
    if (st /= fourier_resample_status_success .or. calls /= 1) error stop 'scalar callback called once'
    call check(y,resample(x,11,frequency_window=weights))
    ! Even FFT grid must contain negative Nyquist frequency, not positive.
    expected_n = 8
    calls = 0
    y = resample(x8,13,window_function=real_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'even scalar callback'
    call check(y,resample(x8,13,frequency_window=real_response(fftfreq(8))))
    expected_n = 7
    do j=1,3
        x2(:,j) = x*real(j,dp)
    end do
    calls = 0
    out2 = -1.0_dp
    call resample_into(x2,out2,11,axis=1,window_function=real_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'multichannel callback once'
    do j=1,3
        ref2(:,j) = resample(x2(:,j),11,frequency_window=weights)
    end do
    if (maxval(abs(out2-ref2)) > tol) error stop 'multichannel results'
    calls = 0
    rr2 = resample(x2,11,axis=1,window_function=real_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'rank2 callback once'
    if (maxval(abs(rr2-ref2)) > tol) error stop 'rank2 callback values'
    do i=1,2
        do k=1,2
            x3(i,:,k) = x*real(i+k,dp)
        end do
    end do
    calls = 0
    rr3 = resample(x3,11,axis=2,window_function=real_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'rank3 callback once'
    calls = 0
    call resample_into(x3,out3,11,axis=2,window_function=real_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'rank3 into callback once'
    if (maxval(abs(out3-rr3)) > tol) error stop 'rank3 into callback values'

    ! Caller-output arrays of arbitrary rank and noncontiguous slices.
    x15(1,1,1,1,1,1,1,:,1,1,1,1,1,1,1) = x
    calls = 0
    call resample_into(x15,out15,11,axis=8,window_function=real_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'rank15 callback once'
    call check(out15(1,1,1,1,1,1,1,:,1,1,1,1,1,1,1),ref2(:,1))
    buf = -37.0_dp
    calls = 0
    call resample_into(x3(:,::1,:),buf(:,2:22:2,:),11,axis=2, &
        window_function=real_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'strided callback once'
    if (maxval(abs(buf(:,2:22:2,:)-rr3)) > tol) error stop 'strided resample values'
    if (any(buf(:,1:21:2,:) /= -37.0_dp)) error stop 'strided gaps changed'

    expected_n = 8
    cweights = complex_response(fftfreq(8))
    calls = 0
    yc = resample(z,5,complex_window_function=complex_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'complex callback once'
    if (maxval(abs(yc-resample(z,5,complex_window=cweights))) > tol) error stop 'complex callback values'
    do j=1,2
        z2(:,j) = z*real(j,dp)
        f2(:,j) = z*real(j,dp)
    end do
    calls = 0
    call resample_into(z2,y2,5,axis=1,complex_window_function=complex_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'complex multichannel callback once'
    do j=1,2
        refc(:,j) = resample(z2(:,j),5,complex_window=cweights)
    end do
    if (maxval(abs(y2-refc)) > tol) error stop 'complex multichannel results'
    do i=1,2
        do k=1,2
            z3(i,:,k) = z*real(i+k,dp)
        end do
    end do
    calls = 0
    cr3 = resample(z3,5,axis=2,complex_window_function=complex_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'complex rank3 callback once'
    calls = 0
    call resample_into(z3,outc3,5,axis=2,complex_window_function=complex_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'complex rank3 into callback once'
    if (maxval(abs(cr3-outc3)) > tol) error stop 'complex rank3 results'
    z15(1,1,:,1,1,1,1,1,1,1,1,1,1,1,1) = z
    calls = 0
    call resample_into(z15,outc15,5,axis=3,complex_window_function=complex_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'complex rank15 callback once'
    if (maxval(abs(outc15(1,1,:,1,1,1,1,1,1,1,1,1,1,1,1)- &
        resample(z,5,complex_window=cweights))) > tol) error stop 'complex rank15 values'
    calls = 0
    call resample_into(f2,y2,5,axis=1,domain='freq', &
        complex_window_function=complex_response,status=st)
    if (st /= 0 .or. calls /= 1) error stop 'frequency domain callback'
    do j=1,2
        refc(:,j) = resample(f2(:,j),5,domain='freq',complex_window=cweights)
    end do
    if (maxval(abs(y2-refc)) > tol) error stop 'frequency domain values'

    ! Mutually exclusive windows and nonfinite callback responses are rejected.
    expected_n = 7
    calls = 0
    y = resample(x,11,window_function=real_response,frequency_window=weights,status=st)
    if (st /= fourier_resample_status_invalid_input .or. size(y) /= 0 .or. calls /= 0) &
        error stop 'reject competing scalar windows'
    calls = 0
    snapshot = buf
    call resample_into(x3,buf(:,2:22:2,:),11,axis=2, &
        window='hann',window_function=real_response,status=st)
    if (st /= fourier_resample_status_invalid_input .or. calls /= 0) error stop 'conflicting ND windows'
    if (any(snapshot /= buf)) error stop 'competing windows changed output'
    reject = .true.
    calls = 0
    call resample_into(x3,buf(:,2:22:2,:),11,axis=2,window_function=real_response,status=st)
    if (st /= fourier_resample_status_invalid_input .or. calls /= 1) error stop 'nonfinite real window'
    if (any(snapshot /= buf)) error stop 'nonfinite window changed output'
    calls = 0
    call resample_into(x3,buf(:,2:22:2,:),11,axis=5,window_function=real_response,status=st)
    if (st /= fourier_resample_status_invalid_input .or. calls /= 0) error stop 'invalid axis callback'
    if (any(snapshot /= buf)) error stop 'invalid axis changed output'
    expected_n = 8
    calls = 0
    call resample_into(z2,y2,5,complex_window_function=complex_response,status=st)
    if (st /= fourier_resample_status_invalid_input .or. calls /= 1) error stop 'nonfinite complex window'
    reject = .false.
    expected_n = 7
    calls = 0
    y = resample(x,11,window_function=real_response,t=[1.0_dp,1.25_dp], &
        output_times=t,status=st)
    if (st /= fourier_resample_status_success .or. calls /= 1) error stop 'sample times callback'
    if (size(t) /= 11) error stop 'callback sample times length'
    if (abs(t(11)-(1.0_dp+10.0_dp*7.0_dp*0.25_dp/11.0_dp)) > tol) &
        error stop 'callback sample times wrong'
    print '(a)', 'test_signal_resample_callable_window passed'

contains
    function fftfreq(n) result(f)
        integer, intent(in) :: n !! transform length
        real(dp) :: f(n)
        integer :: l
        do l=0,n-1
            if (l <= (n-1)/2) then
                f(l+1) = real(l,dp)/real(n,dp)
            else
                f(l+1) = real(l-n,dp)/real(n,dp)
            end if
        end do
    end function fftfreq

    function real_response(f) result(w)
        real(dp), intent(in) :: f(:) !! FFT-order frequencies
        real(dp) :: w(size(f))
        integer :: n
        n = size(f)
        calls = calls+1
        if (n /= expected_n) error stop 'wrong callback frequency count'
        if (any(abs(f-fftfreq(n)) > 2.0e-15_dp)) error stop 'wrong callback frequency ordering'
        w = 0.75_dp + 0.25_dp*sin(11.0_dp*f) + 0.3_dp*f*f
        if (reject) w(1) = ieee_value(0.0_dp,ieee_quiet_nan)
    end function real_response

    function complex_response(f) result(w)
        real(dp), intent(in) :: f(:) !! FFT-order frequencies
        complex(dp) :: w(size(f))
        integer :: n
        n = size(f)
        calls = calls+1
        if (n /= expected_n) error stop 'wrong complex callback length'
        if (any(abs(f-fftfreq(n)) > 2.0e-15_dp)) error stop 'wrong complex callback ordering'
        w = cmplx(1.0_dp+0.2_dp*cos(8.0_dp*f),0.4_dp*sin(7.0_dp*f)-0.1_dp*f,dp)
        if (reject) w(1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    end function complex_response

    subroutine check(a,b)
        real(dp), intent(in) :: a(:) !! computed samples
        real(dp), intent(in) :: b(:) !! reference samples
        if (size(a) /= size(b)) error stop 'shape mismatch'
        if (maxval(abs(a-b)) > tol) error stop 'resample values mismatch'
    end subroutine check
end program test_signal_resample_callable_window
