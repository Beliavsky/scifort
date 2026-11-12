! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_filtfilt_axes
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal, only : dp, filtfilt, sosfiltfilt, &
        filtfilt_status_success, filtfilt_status_invalid_input, filtfilt_status_no_steady_state
    implicit none
    real(dp) :: b(3), a(3), sos(1,6), x2(55,4), x3(3,50,5)
    complex(dp) :: cb(3), ca(3), csos(1,6), cx2(55,4), cx3(3,50,5)
    real(dp), allocatable :: y2(:,:), y3(:,:,:), r(:)
    complex(dp), allocatable :: cy2(:,:), cy3(:,:,:), cr(:)
    character(len=8) :: mode(4)
    integer :: i,j,k,m,st

    b=[0.4_dp,0.25_dp,0.1_dp]
    a=[1.0_dp,-0.2_dp,0.04_dp]
    sos(1,:)=[b(1),b(2),b(3),a(1),a(2),a(3)]
    cb=cmplx(b,0.0_dp,dp)
    ca=cmplx(a,0.0_dp,dp)
    cb(2)=cb(2)+cmplx(0.0_dp,0.03_dp,dp)
    ca(3)=ca(3)+cmplx(0.0_dp,0.02_dp,dp)
    csos(1,:)=[cb(1),cb(2),cb(3),ca(1),ca(2),ca(3)]
    do j=1,4
        do i=1,55
            x2(i,j)=sin(0.1_dp*i)+cos(0.07_dp*j*i)
        end do
    end do
    do k=1,5
        do j=1,50
            do i=1,3
                x3(i,j,k)=sin(0.09_dp*(i+j+k))+0.02_dp*j
            end do
        end do
    end do
    cx2=cmplx(x2,0.2_dp*x2,dp)
    cx3=cmplx(x3,-0.15_dp*x3,dp)
    mode=['odd     ','even    ','constant','none    ']
    do m=1,4
        y2=filtfilt(b,a,x2,axis=1,padtype=trim(mode(m)),status=st)
        if (st/=filtfilt_status_success) error stop 'rank2 filtfilt'
        do j=1,4
            r=filtfilt(b,a,x2(:,j),padtype=trim(mode(m)))
            if (any(abs(y2(:,j)-r)>1.e-13_dp)) error stop 'rank2 lane'
        end do
        if (any(abs(y2-sosfiltfilt(sos,x2,axis=1,padtype=trim(mode(m))))>1.e-13_dp)) &
            error stop 'rank2 SOS reverse'
        cy2=filtfilt(cb,ca,cx2,axis=1,padtype=trim(mode(m)),status=st)
        if (st/=filtfilt_status_success) error stop 'rank2 complex'
        do j=1,4
            cr=filtfilt(cb,ca,cx2(:,j),padtype=trim(mode(m)))
            if (any(abs(cy2(:,j)-cr)>1.e-13_dp)) error stop 'rank2 complex lane'
        end do
        if (any(abs(cy2-sosfiltfilt(csos,cx2,axis=1,padtype=trim(mode(m))))>1.e-13_dp)) &
            error stop 'complex SOS rank2'
    end do
    y2=filtfilt(b,a,transpose(x2))
    do j=1,4
        r=filtfilt(b,a,x2(:,j))
        if (any(abs(y2(j,:)-r)>1.e-13_dp)) error stop 'rank2 default'
    end do
    do k=1,3
        y3=filtfilt(b,a,x3,axis=k,padlen=0,status=st)
        if (st/=filtfilt_status_success) error stop 'rank3 filtfilt'
        if (any(abs(y3-sosfiltfilt(sos,x3,axis=k,padlen=0))>1.e-13_dp)) &
            error stop 'rank3 sosfiltfilt'
        cy3=filtfilt(cb,ca,cx3,axis=k,padlen=0,status=st)
        if (st/=filtfilt_status_success) error stop 'rank3 complex'
        if (any(abs(cy3-sosfiltfilt(csos,cx3,axis=k,padlen=0))>1.e-13_dp)) &
            error stop 'rank3 complex sosfiltfilt'
    end do
    y3=sosfiltfilt(sos,x3)
    if (any(abs(y3-sosfiltfilt(sos,x3,axis=3))>1.e-13_dp)) error stop 'default axis3'
    y2=filtfilt(b,a,x2,axis=0,status=st)
    if (st/=filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(y2))) error stop 'bad axis'
    y3=sosfiltfilt(sos,x3,axis=3,padlen=-1,status=st)
    if (st/=filtfilt_status_invalid_input) error stop 'bad padlen'
    y3=filtfilt(b,a,x3,axis=1,status=st)
    if (st/=filtfilt_status_invalid_input) error stop 'short lane'
    y2=filtfilt([1.0_dp],[1.0_dp,-1.0_dp],x2,axis=1,status=st)
    if (st/=filtfilt_status_no_steady_state) error stop 'unit pole'
    y2=sosfiltfilt(sos,x2,axis=1,padtype='no-such-mode',status=st)
    if (st/=filtfilt_status_invalid_input) error stop 'bad padtype'
    print '(a)', 'test_signal_filtfilt_axes passed'
end program
