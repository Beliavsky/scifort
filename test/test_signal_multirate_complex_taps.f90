! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multirate_complex_taps
    use scifort_signal, only : dp, resample_poly, upfirdn, &
        resample_status_success, resample_status_invalid_input
    implicit none
    real(dp) :: x(9), a(9,3), b(2,9,3), hreal(3), samples(3), isolated(9)
    complex(dp) :: taps(3), z(9), complex_cval
    real(dp), allocatable :: r(:), fallback(:)
    complex(dp), allocatable :: y(:), yc(:), yrank2(:,:), yrank3(:,:,:), anti(:)
    integer :: i, j, k, st
    real(dp), parameter :: tol = 3.0e-12_dp

    taps = cmplx([0.15_dp,0.7_dp,0.15_dp], [0.3_dp,-0.2_dp,0.1_dp],dp)
    hreal = [0.15_dp,0.7_dp,0.15_dp]
    complex_cval = cmplx(0.7_dp,-0.25_dp,dp)
    do i = 1,9
        x(i) = sin(0.24_dp*real(i,dp))+0.2_dp*real(i,dp)
        z(i) = cmplx(x(i),cos(0.39_dp*real(i,dp)),dp)
        do j=1,3
            a(i,j) = x(i)+0.11_dp*real(j,dp)
        end do
    end do
    do k=1,3
        do j=1,9
            do i=1,2
                b(i,j,k) = a(j,k) + real(i,dp)*0.07_dp
            end do
        end do
    end do

    y=resample_poly(x,3,2,taps=taps,status=st)
    if (st/=resample_status_success .or. size(y)/=14) error stop 'complex taps real input'
    yc=resample_poly(cmplx(x,0.0_dp,dp),3,2,taps=taps,status=st)
    if (st/=resample_status_success) error stop 'complex taps complex input'
    if (maxval(abs(y-yc))>tol) error stop 'real/complex input consistency'
    r=resample_poly(x,3,2,taps=hreal,status=st)
    if (st/=resample_status_success) error stop 'original real taps'
    if (maxval(abs(real(y,dp)-r))>1.0_dp) error stop 'complex filter sanity'

    do k=1,3
        select case(k)
        case (1)
            y=resample_poly(z,3,2,taps=taps,padtype='median',status=st)
        case (2)
            y=resample_poly(z,3,2,taps=taps,padtype='minimum',status=st)
        case (3)
            y=resample_poly(z,3,2,taps=taps,padtype='maximum',status=st)
        end select
        if (st/=resample_status_success .or. size(y)/=14) error stop 'statistical padtype'
    end do
    y=resample_poly(x,3,2,taps=taps,padtype='constant',cval=complex_cval,status=st)
    if (st/=resample_status_success) error stop 'complex constant padding'
    yrank2=resample_poly(a,3,2,taps=taps,axis=1,status=st)
    if (st/=resample_status_success .or. any(shape(yrank2)/=[14,3])) error stop 'rank2 complex taps'
    do j=1,3
        y=resample_poly(a(:,j),3,2,taps=taps)
        if (maxval(abs(y-yrank2(:,j)))>tol) error stop 'rank2 complex tap axis1'
    end do
    yrank2=resample_poly(a,3,2,taps=taps,axis=-1,padtype='median',status=st)
    if (st/=resample_status_success .or. any(shape(yrank2)/=[9,5])) error stop 'rank2 axis negative'
    do i=1,9
        y=resample_poly(a(i,:),3,2,taps=taps,padtype='median')
        if (maxval(abs(y-yrank2(i,:)))>tol) error stop 'rank2 complex tap axis2'
    end do
    yrank3=resample_poly(b,2,3,taps=taps,axis=2,padtype='antireflect',status=st)
    if (st/=resample_status_success .or. any(shape(yrank3)/=[2,6,3])) error stop 'rank3 complex taps'
    do k=1,3
        do i=1,2
            y=resample_poly(b(i,:,k),2,3,taps=taps,padtype='antireflect')
            if (maxval(abs(y-yrank3(i,:,k)))>tol) error stop 'rank3 complex tap axis2'
        end do
    end do

    samples=[1.0_dp,3.0_dp,6.0_dp]
    block
        real(dp) :: hot(15)
        hot=0.0_dp
        hot(11)=1.0_dp
        r=upfirdn(hot,samples,mode='antireflect',status=st)
        if (st/=resample_status_success) error stop 'antireflect status'
        if (abs(r(10)+1.0_dp)>tol .or. abs(r(14)-9.0_dp)>tol) error stop 'antireflect values'
        r=upfirdn(hot,samples,mode='antisymmetric',status=st)
        if (st/=resample_status_success) error stop 'antisymmetric status'
        if (abs(r(10)+1.0_dp)>tol .or. abs(r(14)+6.0_dp)>tol) error stop 'antisymmetric values'
    end block
    yrank2=resample_poly(a,2,3,taps=taps,axis=0,status=st)
    if (st/=resample_status_invalid_input .or. size(yrank2)/=0) error stop 'invalid axis'
    y=resample_poly(z,2,3,taps=taps,padtype='badmode',status=st)
    if (st/=resample_status_invalid_input .or. size(y)/=0) error stop 'invalid mode'
    y=resample_poly(x,2,3,taps=taps(1:0),status=st)
    if (st/=resample_status_invalid_input .or. size(y)/=0) error stop 'empty complex taps'
    y=resample_poly(x,2,3,taps=taps,padtype='mean',cval=complex_cval,status=st)
    if (st/=resample_status_invalid_input .or. size(y)/=0) error stop 'irrelevant cval'
    print '(a)', 'test_signal_multirate_complex_taps passed'
end program test_signal_multirate_complex_taps
