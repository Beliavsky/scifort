! SPDX-License-Identifier: MIT
! Checked regression for Gustafsson initial states through the public API.
program test_signal_filtfilt_gust
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : filtfilt, filtfilt_gust, &
        filtfilt_status_success, filtfilt_status_invalid_input
    implicit none
    real(dp) :: b(3), a(3), xs(60), a2(1), b2(1)
    real(dp) :: m2(20,7), t3(8,12,10)
    real(dp), allocatable :: real_y(:), real_y2(:), ym(:,:), yt(:,:,:)
    complex(dp) :: cm2(13,17), ct3(7,11,8), cb(2), ca(2)
    complex(dp), allocatable :: cy(:), cym(:,:), cyt(:,:,:)
    integer :: i,j,k, status

    b = [0.1_dp,0.2_dp,0.1_dp]
    a = [1.0_dp,-0.4_dp,0.2_dp]
    do i=1,size(xs)
        xs(i) = sin(real(i,dp)*0.1_dp) + 0.25_dp*cos(real(i,dp)*0.3_dp)
    end do
    real_y = filtfilt(b,a,xs, method='gust', status=status)
    if (status /= filtfilt_status_success) error stop 'rank-1 Gustafsson status'
    if (maxval(abs(real_y(1:5)-[ &
        0.0765282272491676_dp, &
        0.0997889956879121_dp, &
        0.1137093903354660_dp, &
        0.1206475564666222_dp, &
        0.1243310375151487_dp ])) > 2.0e-12_dp) error stop 'SciPy reference'
    real_y2 = filtfilt_gust(b,a,xs,status=status)
    if (status /= filtfilt_status_success) error stop 'direct Gustafsson status'
    if (maxval(abs(real_y-real_y2)) > 5.0e-14_dp) error stop 'direct Gustafsson differs'
    real_y2 = filtfilt(b,a,xs,method='gust',irlen=5,padlen=-10,status=status)
    if (status /= filtfilt_status_success) error stop 'irlen with ignored padlen'
    if (.not. all(abs(real_y2) < 1.0e4_dp)) error stop 'irlen result not finite'
    real_y = filtfilt(b,a,xs(1:60:2),method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'strided signal status'
    real_y2 = filtfilt_gust(b,a,xs(1:60:2))
    if (maxval(abs(real_y-real_y2)) > 5.0e-13_dp) error stop 'strided data differs'
    real_y2 = filtfilt(b,a,xs,padlen=9,status=status)
    if (status /= filtfilt_status_success) error stop 'default pad status'
    if (maxval(abs(real_y2-filtfilt(b,a,xs,method='pad',padlen=9))) > 1.e-14_dp) &
        error stop 'pad method changed'

    real_y2 = filtfilt(b,a,xs,method='gust',irlen=0,status=status)
    if (status /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(real_y2))) &
        error stop 'invalid irlen accepted'
    real_y2 = filtfilt(b,a,xs,method='unknown',status=status)
    if (status /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(real_y2))) &
        error stop 'invalid method accepted'
    real_y2 = filtfilt(b,[0.0_dp,0.0_dp,0.0_dp],xs,method='gust',status=status)
    if (status /= filtfilt_status_invalid_input) error stop 'bad denominator accepted'
    b2 = [0.25_dp]
    a2 = [0.5_dp]
    real_y2 = filtfilt(b2,a2,xs(1:1),method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'zeroth-order status'
    if (abs(real_y2(1) - 0.25_dp * xs(1)) > 1.0e-14_dp) error stop 'zeroth-order result'

    do j=1,size(m2,2)
        do i=1,size(m2,1)
            m2(i,j) = sin(0.1_dp*i*j)
        end do
    end do
    ym = filtfilt(b,a,m2,axis=1,method='gust',irlen=3,status=status)
    if (status /= filtfilt_status_success) error stop 'real matrix axis 1 status'
    do j=1,size(m2,2)
        real_y = filtfilt_gust(b,a,m2(:,j),irlen=3)
        if (maxval(abs(ym(:,j)-real_y)) > 3.e-11_dp) error stop 'real matrix axis 1 values'
    end do
    ym = filtfilt(b,a,m2,axis=2,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'real matrix axis 2 status'
    do i=1,size(m2,1)
        real_y = filtfilt_gust(b,a,m2(i,:))
        if (maxval(abs(ym(i,:)-real_y)) > 3.e-11_dp) error stop 'real matrix axis 2 values'
    end do
    block
        real(dp) :: empty(16,0)
        real(dp), allocatable :: ye(:,:)
        ye = filtfilt(b,a,empty,axis=1,method='gust',status=status)
        if (status /= filtfilt_status_success .or. size(ye) /= 0) &
            error stop 'empty outer channels status'
        ye = filtfilt(b,[0._dp,0._dp,0._dp],empty,axis=1,method='gust',status=status)
        if (status /= filtfilt_status_invalid_input) &
            error stop 'empty channels invalid denominator'
    end block
    ym = filtfilt(b,a,m2,axis=3,method='gust',status=status)
    if (status /= filtfilt_status_invalid_input) error stop 'invalid matrix axis'

    do k=1,size(t3,3)
        do j=1,size(t3,2)
            do i=1,size(t3,1)
                t3(i,j,k) = sin(0.04_dp*i*j+k*0.2_dp)
            end do
        end do
    end do
    do k=1,3
        yt = filtfilt(b,a,t3,axis=k,method='gust',irlen=3,status=status)
        if (status /= filtfilt_status_success) error stop 'real rank3 status'
        select case (k)
        case (1)
            real_y = filtfilt_gust(b,a,t3(:,3,4),irlen=3)
            if (maxval(abs(yt(:,3,4)-real_y)) > 2.e-11_dp) error stop 'rank3 axis1'
        case (2)
            real_y = filtfilt_gust(b,a,t3(2,:,5),irlen=3)
            if (maxval(abs(yt(2,:,5)-real_y)) > 2.e-11_dp) error stop 'rank3 axis2'
        case (3)
            real_y = filtfilt_gust(b,a,t3(6,7,:),irlen=3)
            if (maxval(abs(yt(6,7,:)-real_y)) > 2.e-11_dp) error stop 'rank3 axis3'
        end select
    end do

    do j=1,size(cm2,2)
        do i=1,size(cm2,1)
            cm2(i,j) = cmplx(cos(0.05_dp*i*j),sin(0.13_dp*i+0.2_dp*j),dp)
        end do
    end do
    cb = cmplx([0.2_dp, 0.1_dp],[0.1_dp,-0.05_dp],dp)
    ca = cmplx([1.0_dp,-0.4_dp],[0.0_dp,0.2_dp],dp)
    cym = filtfilt(cb,ca,cm2,axis=1,method='gust',irlen=5,status=status)
    if (status /= filtfilt_status_success) error stop 'complex matrix axis 1 status'
    cy = filtfilt_gust(cb,ca,cm2(:,4),irlen=5,status=status)
    if (status /= filtfilt_status_success) error stop 'complex direct status'
    if (maxval(abs(cym(:,4)-cy)) > 2.e-11_dp) error stop 'complex matrix values'
    do k=1,size(ct3,3)
        do j=1,size(ct3,2)
            do i=1,size(ct3,1)
                ct3(i,j,k) = cmplx(sin(0.1_dp*i*j+k),cos(0.02_dp*i+j+k),dp)
            end do
        end do
    end do
    cyt = filtfilt(cb,ca,ct3,axis=3,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'complex rank3 status'
    cy = filtfilt_gust(cb,ca,ct3(2,3,:))
    if (maxval(abs(cyt(2,3,:)-cy)) > 2.e-11_dp) error stop 'complex rank3 values'
    cy = filtfilt_gust(b,a,cm2(:,1))
    if (maxval(abs(cy - filtfilt_gust(cmplx(b,0._dp,dp), &
        cmplx(a,0._dp,dp),cm2(:,1)))) > 5.e-13_dp) error stop 'mixed coefficients'
    print *, 'PASS: filtfilt Gustafsson public API tests'
end program test_signal_filtfilt_gust
