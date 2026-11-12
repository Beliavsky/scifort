! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multirate_nd
    use scifort_signal, only : dp, resample, resample_poly, upfirdn, decimate, &
        resample_status_success, resample_status_invalid_input, &
        fourier_resample_status_success, fourier_resample_status_invalid_input
    implicit none
    real(dp) :: x(12,4), xr(5,12,3), h(3)
    complex(dp) :: z(12,4), hz(3)
    real(dp), allocatable :: a(:,:), b(:,:,:), v(:), times(:)
    complex(dp), allocatable :: c(:,:), d(:,:), t(:)
    integer :: i, j, k, st
    real(dp), parameter :: tol = 4.0e-12_dp

    h = [0.2_dp, 0.6_dp, 0.2_dp]
    hz = cmplx(h,[0.1_dp,-0.2_dp,0.3_dp],dp)
    do j = 1, 4
        do i = 1, 12
            x(i,j) = sin(0.17_dp*real(i*j,dp)) + 0.03_dp*real(i+j,dp)
            z(i,j) = cmplx(x(i,j), cos(0.2_dp*real(i+2*j,dp)),dp)
        end do
    end do
    do k = 1, 3
        do j = 1, 12
            do i = 1, 5
                xr(i,j,k) = sin(real(i+j+k,dp)*0.12_dp)
            end do
        end do
    end do

    ! A Fourier resampling axis changes exactly one array extent.
    a = resample(x,9,axis=1,status=st)
    if (st/=fourier_resample_status_success .or. any(shape(a)/=[9,4])) error stop 'resample rank2 shape'
    do j=1,4
        v=resample(x(:,j),9)
        if (maxval(abs(a(:,j)-v))>tol) error stop 'resample rank2 axis1'
    end do
    c = resample(z,7,window='hann',axis=-1,status=st)
    if (st/=fourier_resample_status_success .or. any(shape(c)/=[12,7])) error stop 'resample rank2 complex shape'
    do i=1,12
        t=resample(z(i,:),7,window='hann')
        if (maxval(abs(c(i,:)-t))>tol) error stop 'resample rank2 axis2'
    end do
    b=resample(xr,8,axis=2,t=[0.0_dp,0.1_dp],output_times=times,status=st)
    if (st/=fourier_resample_status_success .or. any(shape(b)/=[5,8,3])) error stop 'resample rank3'
    if (size(times)/=8) error stop 'Fourier resample times'
    if (abs(times(8)-1.05_dp)>tol) error stop 'Fourier resample time endpoint'
    do k=1,3
        do i=1,5
            v=resample(xr(i,:,k),8)
            if (maxval(abs(b(i,:,k)-v))>tol) error stop 'resample rank3 axis2'
        end do
    end do

    ! Polyphase resampling handles real/complex arrays and first/last axes.
    a=resample_poly(x,3,2,taps=h,axis=1,status=st)
    if (st/=resample_status_success .or. any(shape(a)/=[18,4])) error stop 'poly rank2 shape'
    do j=1,4
        v=resample_poly(x(:,j),3,2,taps=h)
        if (maxval(abs(a(:,j)-v))>tol) error stop 'poly rank2 axis1'
    end do
    c=resample_poly(z,2,3,taps=h,padtype='mean',axis=2,status=st)
    if (st/=resample_status_success .or. any(shape(c)/=[12,3])) error stop 'poly complex shape'
    do i=1,12
        t=resample_poly(z(i,:),2,3,taps=h,padtype='mean')
        if (maxval(abs(c(i,:)-t))>tol) error stop 'poly rank2 complex'
    end do
    b=resample_poly(xr,2,3,taps=h,axis=-1,status=st)
    if (st/=resample_status_success .or. any(shape(b)/=[5,12,2])) error stop 'poly rank3 shape'

    ! Explicit four real/complex FIR/signal combinations.
    a=upfirdn(h,x,up=2,down=3,axis=1,status=st)
    if (st/=resample_status_success) error stop 'upfirdn rank2 status'
    do j=1,4
        v=upfirdn(h,x(:,j),up=2,down=3)
        if (maxval(abs(a(:,j)-v))>tol) error stop 'upfirdn rank2'
    end do
    c=upfirdn(h,z,up=2,down=3,axis=2,status=st)
    if (st/=resample_status_success) error stop 'upfirdn rc'
    d=upfirdn(hz,x,up=2,down=3,axis=2,status=st)
    if (st/=resample_status_success) error stop 'upfirdn cr'
    d=upfirdn(hz,z,up=2,down=3,axis=2,status=st)
    if (st/=resample_status_success) error stop 'upfirdn cc'
    b=upfirdn(h,xr,up=2,down=3,axis=2,status=st)
    if (st/=resample_status_success) error stop 'upfirdn rank3'

    a=decimate(x,2,n=12,ftype='fir',axis=1,status=st)
    if (st/=resample_status_success .or. any(shape(a)/=[6,4])) error stop 'decimate rank2 shape'
    do j=1,4
        v=decimate(x(:,j),2,n=12,ftype='fir')
        if (maxval(abs(a(:,j)-v))>tol) error stop 'decimate rank2'
    end do
    c=decimate(z,2,n=12,ftype='fir',zero_phase=.false.,axis=-1,status=st)
    if (st/=resample_status_success .or. any(shape(c)/=[12,2])) error stop 'decimate complex shape'
    b=decimate(xr,2,n=12,ftype='fir',axis=2,status=st)
    if (st/=resample_status_success .or. any(shape(b)/=[5,6,3])) error stop 'decimate rank3 shape'

    ! Invalid axes and parameters have defined zero-extent results.
    a=resample(x,4,axis=0,status=st)
    if (st/=fourier_resample_status_invalid_input .or. size(a)/=0) error stop 'bad Fourier axis'
    a=resample_poly(x,3,0,axis=1,status=st)
    if (st/=resample_status_invalid_input .or. size(a)/=0) error stop 'bad poly factor'
    a=upfirdn(h,x,axis=-3,status=st)
    if (st/=resample_status_invalid_input .or. size(a)/=0) error stop 'bad upfirdn axis'
    a=decimate(x,0,axis=1,status=st)
    if (st/=resample_status_invalid_input .or. size(a)/=0) error stop 'bad decimation'

    print '(a)', 'test_signal_multirate_nd passed'
end program test_signal_multirate_nd
