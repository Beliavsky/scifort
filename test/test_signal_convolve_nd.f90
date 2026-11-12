! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_convolve_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, convolve_nd, correlate_nd, convolution_shape_nd, &
        convolve_nd_success, convolve_nd_invalid_input, convolve2d, correlate2d
    implicit none
    integer :: i, j, t, status
    integer, allocatable :: dims(:)
    real(dp) :: x(5,4), w(2,3), z(9,8), zr(9,8), large(10,8), target(4,3)
    real(dp) :: a1(3), b1(2), y1(4)
    real(dp), allocatable :: big(:,:), kernel(:,:), got(:,:), reference(:,:)
    real(dp) :: xa(3,1,2), wa(1,2,1), y3(3,2,2)
    real(dp) :: x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: w15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: y15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    complex(dp) :: cx(5,4), cw(2,3), cy(6,6), cy2(6,6), cy3(6,6)
    integer :: cases

    cases = 0
    do j=1,4
        do i=1,5
            x(i,j) = sin(real(i*7-j*5,dp))
            cx(i,j) = cmplx(x(i,j),0.1_dp*real(i-2*j,dp),dp)
        end do
    end do
    do j=1,3
        do i=1,2
            w(i,j) = cos(real(i*3+j*7,dp))
            cw(i,j) = cmplx(w(i,j),0.2_dp*real(i+j,dp),dp)
        end do
    end do
    do t=1,3
        select case(t)
        case(1)
            dims = convolution_shape_nd(shape(x),shape(w),'full',status=status)
        case(2)
            dims = convolution_shape_nd(shape(x),shape(w),'same',status=status)
        case(3)
            dims = convolution_shape_nd(shape(x),shape(w),'valid',status=status)
        end select
        if (status/=convolve_nd_success) error stop 'layout'
        if (size(dims)/=2) error stop 'layout rank'
        select case(t)
        case(1)
            call convolve_nd(x,w,z(:6,:6),mode='full',method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'convolve direct full'
            if (maxval(abs(z(:6,:6)-convolve2d(x,w,'full')))>2.e-12_dp) error stop 'convolve full'
            call correlate_nd(x,w,zr(:6,:6),mode='full',method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'correlate direct full'
            if (maxval(abs(zr(:6,:6)-correlate2d(x,w,'full')))>2.e-12_dp) error stop 'correlate full'
            call correlate_nd(cx,cw,cy,method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'complex correlate'
            call correlate_nd(cx,cw,cy2,method='fft',status=status)
            if (status/=convolve_nd_success) error stop 'complex correlate fft'
            if (maxval(abs(cy-cy2))>2.e-11_dp) error stop 'complex correlation difference'
            call convolve_nd(cx,cw,cy3,method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'complex convolve'
            call convolve_nd(cx,cw,cy2,method='fft',status=status)
            if (status/=convolve_nd_success) error stop 'complex convolve fft'
            if (maxval(abs(cy3-cy2))>2.e-11_dp) error stop 'complex convolution difference'
            call correlate_nd(x,cw,cy2,method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'mixed correlation'
            call convolve_nd(cx,w,cy3,method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'mixed convolution'
        case(2)
            call convolve_nd(x,w,z(:5,:4),mode='same',method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'convolve same'
            if (maxval(abs(z(:5,:4)-convolve2d(x,w,'same')))>2.e-12_dp) error stop 'convolve same values'
            call correlate_nd(x,w,zr(:5,:4),mode='same',method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'correlate same'
            ! For even second operands correlate2d uses a different center.
            call correlate_nd(x,w,z(:5,:4),mode='same',method='fft',status=status)
            if (status/=convolve_nd_success) error stop 'correlate fft same'
            if (maxval(abs(zr(:5,:4)-z(:5,:4)))>2.e-11_dp) error stop 'correlate same direct/fft'
        case(3)
            call convolve_nd(x,w,z(:4,:2),mode='valid',method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'convolve valid'
            if (maxval(abs(z(:4,:2)-convolve2d(x,w,'valid')))>2.e-12_dp) error stop 'convolve valid values'
            call correlate_nd(x,w,zr(:4,:2),mode='valid',method='direct',status=status)
            if (status/=convolve_nd_success) error stop 'correlate valid'
            if (maxval(abs(zr(:4,:2)-correlate2d(x,w,'valid')))>2.e-12_dp) error stop 'correlate valid values'
        end select
        cases=cases+1
    end do

    ! 1D: both operands reversed in valid mode, and auto dispatch.
    a1 = [1.0_dp,2.0_dp,-1.0_dp]
    b1 = [0.5_dp,3.0_dp]
    call convolve_nd(a1,b1,y1,'full','auto',status=status)
    if (status/=convolve_nd_success) error stop 'rank1 auto'
    if (maxval(abs(y1-[0.5_dp,4.0_dp,5.5_dp,-3.0_dp]))>1.e-12_dp) error stop 'rank1 convolution'
    call correlate_nd(a1,b1,y1,'full','direct',status=status)
    if (status/=convolve_nd_success) error stop 'rank1 correlate'
    if (maxval(abs(y1-[3.0_dp,6.5_dp,-2.0_dp,-0.5_dp]))>1.e-12_dp) error stop 'rank1 correlation'
    cases=cases+2

    ! Wrong shape, nan input, duplicate axes, and unsupported methods must
    ! not modify the caller's result (including strided sections).
    a1(2) = ieee_value(0.0_dp,ieee_quiet_nan)
    y1 = 12.0_dp
    call convolve_nd(a1,b1,y1,'full','direct',status=status)
    if (status/=convolve_nd_invalid_input .or. any(y1/=12.0_dp)) &
        error stop 'nan must reject without writing'
    a1(2) = 2.0_dp
    cases=cases+1

    ! Selected axes, singleton broadcasting and asymmetric nonconvolved shape.
    xa = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[3,1,2])
    wa = reshape([2.0_dp,-1.0_dp],[1,2,1])
    call convolve_nd(xa,wa,y3,axes=[2],method='direct',status=status)
    if (status/=convolve_nd_success) error stop 'axis2 direct'
    call convolve_nd(xa,wa,y3,axes=[2],method='fft',status=status)
    if (status/=convolve_nd_success) error stop 'axis2 fft'
    call correlate_nd(xa,wa,y3,axes=[2],method='direct',status=status)
    if (status/=convolve_nd_success) error stop 'axis2 correlation'
    cases=cases+3

    ! Rank-15 matrix with two nontrivial axes, and mixed length orders.
    x15 = 0.0_dp
    w15 = 0.0_dp
    x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=1.0_dp
    x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)=3.0_dp
    w15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=2.0_dp
    w15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)=-1.0_dp
    call correlate_nd(x15,w15,y15,method='direct',status=status)
    if (status/=convolve_nd_success) error stop 'rank15 direct'
    if (abs(y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)+1.0_dp)>1.e-12_dp) &
        error stop 'rank15 check'
    call correlate_nd(x15,w15,y15,method='fft',status=status)
    if (status/=convolve_nd_success) error stop 'rank15 fft'
    cases=cases+2

    ! Noncontiguous output, inputs and alias-free copyback from assumed-rank.
    large=0.0_dp
    large(1:9:2,1:7:2)=x
    call convolve_nd(large(1:9:2,1:7:2),w,z(:6,:6),method='direct',status=status)
    if (status/=convolve_nd_success) error stop 'noncontiguous inputs'
    if (maxval(abs(z(:6,:6)-convolve2d(x,w)))>1.e-12_dp) error stop 'strided values'
    zr=-9.0_dp
    call convolve_nd(x,w,zr(1:7:2,1:5:2),mode='valid',method='direct',status=status)
    if (status/=convolve_nd_invalid_input) error stop 'strided output invalid shape'
    if (any(zr/=-9.0_dp)) error stop 'strided error changed output'
    zr=-9.0_dp
    call convolve_nd(x,w,zr(1:7:2,1:3:2),mode='valid',method='direct',status=status)
    if (status/=convolve_nd_success) error stop 'strided output valid'
    if (maxval(abs(zr(1:7:2,1:3:2)-convolve2d(x,w,'valid')))>1.e-12_dp) &
        error stop 'strided output mismatch'
    cases=cases+3

    ! Automatic dispatch should select FFT when the direct work is large.
    allocate(big(140,100),kernel(7,5),got(146,104),reference(146,104))
    do j=1,100
        do i=1,140
            big(i,j)=cos(real(3*i-5*j,dp)/21.0_dp)
        end do
    end do
    do j=1,5
        do i=1,7
            kernel(i,j)=sin(real(i*j,dp)/10.0_dp)
        end do
    end do
    call convolve_nd(big,kernel,got,method='auto',status=status)
    if (status/=convolve_nd_success) error stop 'auto large convolution'
    call convolve_nd(big,kernel,reference,method='fft',status=status)
    if (status/=convolve_nd_success) error stop 'fft large convolution'
    if (maxval(abs(got-reference))>2.e-11_dp) error stop 'auto fft convolution'
    call correlate_nd(big,kernel,got,method='auto',status=status)
    if (status/=convolve_nd_success) error stop 'auto large correlation'
    call correlate_nd(big,kernel,reference,method='fft',status=status)
    if (status/=convolve_nd_success) error stop 'fft large correlation'
    if (maxval(abs(got-reference))>2.e-11_dp) error stop 'auto fft correlation'
    cases=cases+2

    target=17.0_dp
    call convolve_nd(x,w,target,mode='full',status=status)
    if (status/=convolve_nd_invalid_input .or. any(target/=17.0_dp)) &
        error stop 'invalid output overwritten'
    call convolve_nd(x,w,z(:6,:6),mode='bogus',status=status)
    if (status/=convolve_nd_invalid_input) error stop 'invalid mode'
    call correlate_nd(x,w,z(:6,:6),axes=[1,1],status=status)
    if (status/=convolve_nd_invalid_input) error stop 'duplicate axis'
    call correlate_nd(x,w,z(:6,:6),method='unknown',status=status)
    if (status/=convolve_nd_invalid_input) error stop 'invalid method'
    dims=convolution_shape_nd(shape(x),shape(w),mode='valid',axes=[1,1],status=status)
    if (status/=convolve_nd_invalid_input .or. size(dims)/=0) error stop 'invalid helper'
    cases=cases+5
    print '(a,i0,a)', 'PASS rank-independent convolution/correlation ', cases, ' native cases'
end program test_signal_convolve_nd
