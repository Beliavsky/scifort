! SPDX-License-Identifier: MIT
! Verify rank-independent overlap-add including real FFT blocking.
program test_signal_oaconvolve_assumed_rank
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : fftconvolve_nd, oaconvolve_nd, oaconvolve
    use scifort_signal_convolve, only : convolve, convolution_status_success, &
        convolution_status_invalid_input
    implicit none
    real(dp) :: a(513), b(7), y(519), ref(519), same(513), valid(507)
    real(dp) :: r2(257,3), k2(5,1), out2(261,3), fft2(261,3)
    real(dp) :: r2last(2,257), k2last(1,7), out2last(2,263), fft2last(2,263)
    real(dp) :: r3(259,2,9), k3(5,1,3), out3(263,2,11), fft3(263,2,11)
    real(dp) :: r5(257,2,2,1,1), k5(5,1,1,1,1), out5(261,2,2,1,1)
    real(dp) :: fft5(261,2,2,1,1)
    real(dp) :: r15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,513)
    real(dp) :: k15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,7)
    real(dp) :: out15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,519)
    real(dp) :: fft15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,519)
    real(dp) :: noncontig(261,4), wrong(518), empty(0), empty_out(0)
    complex(dp) :: ca(513), cb(7), cy(519), cf(519)
    integer :: i, j, code

    a = [(sin(0.08_dp*i)+0.05_dp*i/513.0_dp,i=1,size(a))]
    b = [(cos(0.37_dp*i),i=1,size(b))]
    call oaconvolve_nd(a,b,y,status=code)
    if (code /= convolution_status_success) error stop '1d OA failed'
    ref = convolve(a,b)
    if (maxval(abs(y-ref)) > 3.e-12_dp) error stop '1d OA direct mismatch'
    if (maxval(abs(y-oaconvolve(a,b))) > 3.e-12_dp) error stop '1d legacy mismatch'
    call oaconvolve_nd(a,b,same,mode='same',status=code)
    if (code /= convolution_status_success) error stop 'same failed'
    if (maxval(abs(same-ref(4:516))) > 3.e-12_dp) error stop 'same mismatch'
    call oaconvolve_nd(b,a,valid,mode='valid',status=code)
    if (code /= convolution_status_success) error stop 'reversed valid failed'
    if (maxval(abs(valid-ref(7:513))) > 3.e-12_dp) error stop 'reversed valid mismatch'
    call oaconvolve_nd(b,a,y,mode='full',status=code)
    if (code /= convolution_status_success) error stop 'reversed OA failed'
    if (maxval(abs(y-ref)) > 3.e-12_dp) error stop 'reversed OA mismatch'

    do j=1,3
        r2(:,j) = a(:257)*real(j,dp)
    end do
    k2(:,1) = b(:5)
    call oaconvolve_nd(r2,k2,out2,axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'rank2 axis1 failed'
    call fftconvolve_nd(r2,k2,fft2,axes=[1],status=code)
    if (maxval(abs(out2-fft2)) > 3.e-12_dp) error stop 'rank2 axis1 mismatch'
    do j=1,2
        r2last(j,:) = sin(0.01_dp*j*[(real(i,dp),i=1,257)])
    end do
    k2last(1,:) = b
    call oaconvolve_nd(r2last,k2last,out2last,axes=[2],status=code)
    if (code /= convolution_status_success) error stop 'rank2 last axis failed'
    call fftconvolve_nd(r2last,k2last,fft2last,axes=[2],status=code)
    if (maxval(abs(out2last-fft2last)) > 3.e-12_dp) error stop 'rank2 last axis mismatch'

    r3 = reshape([(sin(0.05_dp*i),i=1,size(r3))],shape(r3))
    k3 = reshape([(cos(0.04_dp*i),i=1,size(k3))],shape(k3))
    call oaconvolve_nd(r3,k3,out3,axes=[1,3],status=code)
    if (code /= convolution_status_success) error stop 'rank3 multiple axes failed'
    call fftconvolve_nd(r3,k3,fft3,axes=[1,3],status=code)
    if (maxval(abs(out3-fft3)) > 2.e-11_dp) error stop 'rank3 multiple axes mismatch'

    r5 = reshape([(sin(0.03_dp*i),i=1,size(r5))],shape(r5))
    k5 = reshape([(cos(0.21_dp*i),i=1,size(k5))],shape(k5))
    call oaconvolve_nd(r5,k5,out5,axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'rank5 failed'
    call fftconvolve_nd(r5,k5,fft5,axes=[1],status=code)
    if (maxval(abs(out5-fft5)) > 4.e-12_dp) error stop 'rank5 mismatch'

    r15 = reshape(a,shape(r15))
    k15 = reshape(b,shape(k15))
    call oaconvolve_nd(r15,k15,out15,axes=[15],status=code)
    if (code /= convolution_status_success) error stop 'rank15 axis15 failed'
    call fftconvolve_nd(r15,k15,fft15,axes=[15],status=code)
    if (maxval(abs(out15-fft15)) > 3.e-12_dp) error stop 'rank15 axis15 mismatch'

    ca = cmplx(a,0.13_dp*a,dp)
    cb = cmplx(b,-0.23_dp*b,dp)
    call oaconvolve_nd(ca,cb,cy,status=code)
    if (code /= convolution_status_success) error stop 'complex/complex failed'
    call fftconvolve_nd(ca,cb,cf,status=code)
    if (maxval(abs(cy-cf)) > 5.e-12_dp) error stop 'complex/complex mismatch'
    call oaconvolve_nd(a,cb,cy,status=code)
    if (code /= convolution_status_success) error stop 'real/complex failed'
    call fftconvolve_nd(a,cb,cf,status=code)
    if (maxval(abs(cy-cf)) > 5.e-12_dp) error stop 'real/complex mismatch'
    call oaconvolve_nd(ca,b,cy,status=code)
    if (code /= convolution_status_success) error stop 'complex/real failed'
    call fftconvolve_nd(ca,b,cf,status=code)
    if (maxval(abs(cy-cf)) > 5.e-12_dp) error stop 'complex/real mismatch'

    noncontig = -99.0_dp
    call oaconvolve_nd(r2(:,1:3:2),k2,noncontig(:,1:3:2),axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'noncontiguous output failed'
    if (any(noncontig(:,2) /= -99.0_dp)) error stop 'noncontiguous overwrite'
    if (maxval(abs(noncontig(:,1)-fft2(:,1))) > 3.e-12_dp) &
        error stop 'noncontiguous output mismatch'

    call oaconvolve_nd(a,b,wrong,status=code)
    if (code /= convolution_status_invalid_input) error stop 'wrong shape not rejected'
    if (.not. all(ieee_is_nan(wrong))) error stop 'wrong shape not NaN'
    call oaconvolve_nd(a,b,y,axes=[1,1],status=code)
    if (code /= convolution_status_invalid_input) error stop 'duplicate axis not rejected'
    call oaconvolve_nd(a,b,y,axes=[2],status=code)
    if (code /= convolution_status_invalid_input) error stop 'out-of-range axis'
    call oaconvolve_nd(a,b,y,mode='bad',status=code)
    if (code /= convolution_status_invalid_input) error stop 'invalid mode not rejected'
    call oaconvolve_nd(a,b,y,axes=[integer::],status=code)
    if (code /= convolution_status_invalid_input) error stop 'empty axes not rejected'
    call oaconvolve_nd(empty,b,empty_out,status=code)
    if (code /= convolution_status_invalid_input) error stop 'empty input not rejected'
    block
        real(dp) :: nan_a(513)
        nan_a = a
        nan_a(20) = ieee_value(0.0_dp,ieee_quiet_nan)
        call oaconvolve_nd(nan_a,b,y,status=code)
        if (code /= convolution_status_invalid_input) error stop 'NaN not rejected'
        if (.not. all(ieee_is_nan(y))) error stop 'NaN not propagated to output'
    end block
    print '(a)', 'test_signal_oaconvolve_assumed_rank passed'
end program test_signal_oaconvolve_assumed_rank
