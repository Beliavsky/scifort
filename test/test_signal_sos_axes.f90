! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_sos_axes
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal, only : dp, sosfilt, sos_status_success, sos_status_invalid_input
    implicit none
    real(dp) :: sos(2,6), x2(35,5), x3(4,35,3), zi2(2,2,5)
    complex(dp) :: csos(2,6), cx2(35,5), cx3(4,35,3)
    real(dp), allocatable :: y2(:,:), y3(:,:,:), split(:,:), test(:), zf2(:,:,:), zf3(:,:,:,:)
    real(dp), allocatable :: znext(:,:,:)
    complex(dp), allocatable :: cy2(:,:), cy3(:,:,:), czf3(:,:,:,:)
    integer :: i, j, k, status

    sos(1,:) = [0.2_dp,0.3_dp,0.1_dp,1.0_dp,-0.3_dp,0.12_dp]
    sos(2,:) = [0.8_dp,0.1_dp,0.2_dp,1.0_dp,-0.2_dp,0.04_dp]
    do j=1,5
        do i=1,35
            x2(i,j) = sin(real(i+j,dp)*0.17_dp)+0.01_dp*i
        end do
    end do
    do k=1,3
        do j=1,35
            do i=1,4
                x3(i,j,k) = cos(0.11_dp*(i+j+k)) + 0.05_dp*j
            end do
        end do
    end do
    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = csos(1,2) + cmplx(0.0_dp,0.04_dp,dp)
    csos(2,5) = csos(2,5) + cmplx(0.0_dp,0.01_dp,dp)
    cx2 = cmplx(x2,0.1_dp*x2,dp)
    cx3 = cmplx(x3,0.15_dp*x3,dp)

    y2 = sosfilt(sos,x2,axis=1,zf=zf2,status=status)
    if (status /= sos_status_success) error stop 'rank2 sosfilt axis1'
    if (any(shape(zf2) /= [2,2,5])) error stop 'rank2 sosfilt state shape'
    do j=1,5
        test = sosfilt(sos,x2(:,j))
        if (any(abs(y2(:,j)-test)>1.e-13_dp)) error stop 'rank2 sosfilt lane'
    end do
    zi2 = 0.025_dp
    y2 = sosfilt(sos,x2,zi=zi2,axis=1,zf=zf2,status=status)
    if (status /= sos_status_success) error stop 'rank2 sosfilt zi'
    do j=1,5
        test = sosfilt(sos,x2(:,j),zi=zi2(:,:,j))
        if (any(abs(y2(:,j)-test)>1.e-13_dp)) error stop 'rank2 sosfilt zi lane'
    end do
    split = sosfilt(sos,x2(1:13,:),zi=zi2,axis=1,zf=znext,status=status)
    if (status /= sos_status_success) error stop 'first block'
    y2 = sosfilt(sos,x2(14:,:),zi=znext,axis=1,zf=zf2,status=status)
    if (status /= sos_status_success) error stop 'second block'
    ! Every full lane must match two independent filtered blocks.
    do j=1,5
        test = sosfilt(sos,x2(:,j),zi=zi2(:,:,j))
        if (any(abs([split(:,j),y2(:,j)]-test)>1.e-13_dp)) error stop 'streaming axis1'
    end do

    y2 = sosfilt(sos,transpose(x2))
    if (any(shape(y2) /= [5,35])) error stop 'default rank2 axis'
    do j=1,5
        test = sosfilt(sos,x2(:,j))
        if (any(abs(y2(j,:)-test)>1.e-13_dp)) error stop 'default rank2 lane'
    end do
    cy2 = sosfilt(csos,cx2,axis=1)
    do j=1,5
        if (any(abs(cy2(:,j)-sosfilt(csos,cx2(:,j)))>1.e-13_dp)) &
            error stop 'rank2 complex'
    end do

    do k=1,3
        y3 = sosfilt(sos,x3,axis=k,zf=zf3,status=status)
        if (status /= sos_status_success) error stop 'rank3 sosfilt'
        if (size(zf3,k+1) /= 2) error stop 'rank3 state axis'
        select case (k)
        case (1)
            do j=1,3
                do i=1,35
                    test=sosfilt(sos,x3(:,i,j))
                    if (any(abs(y3(:,i,j)-test)>1.e-13_dp)) error stop 'rank3 axis1'
                end do
            end do
        case (2)
            do j=1,3
                do i=1,4
                    test=sosfilt(sos,x3(i,:,j))
                    if (any(abs(y3(i,:,j)-test)>1.e-13_dp)) error stop 'rank3 axis2'
                end do
            end do
        case (3)
            do j=1,35
                do i=1,4
                    test=sosfilt(sos,x3(i,j,:))
                    if (any(abs(y3(i,j,:)-test)>1.e-13_dp)) error stop 'rank3 axis3'
                end do
            end do
        end select
        cy3 = sosfilt(csos,cx3,axis=k,zf=czf3,status=status)
        if (status /= sos_status_success) error stop 'rank3 complex'
        if (any(shape(czf3) /= [2, merge(2,4,k==1), &
            merge(2,35,k==2),merge(2,3,k==3)])) error stop 'complex final-state shape'
    end do
    y3 = sosfilt(sos,x3)
    if (any(abs(y3-sosfilt(sos,x3,axis=3))>1.e-13_dp)) error stop 'rank3 default axis'

    y2=sosfilt(sos,x2,axis=0,status=status)
    if (status/=sos_status_invalid_input .or. .not. all(ieee_is_nan(y2))) error stop 'bad axis'
    y2=sosfilt(sos,x2,axis=1,zi=zi2(:1,:,:),zf=zf2,status=status)
    if (status/=sos_status_invalid_input) error stop 'bad zi shape'
    if (size(zf2) /= 0) error stop 'bad zi zf state'
    y3=sosfilt(sos,x3,axis=4,status=status)
    if (status/=sos_status_invalid_input .or. .not. all(ieee_is_nan(y3))) error stop 'rank3 bad axis'
    print '(a)', 'test_signal_sos_axes passed'
end program
