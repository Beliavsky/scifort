! SPDX-License-Identifier: MIT
! Public assumed-rank spline filtering and coordinate evaluation regression.
program test_signal_spline_nd
    use scifort_signal, only : dp, spline_filter_nd, spline_filter1d, &
        cspline_nd, qspline_nd, map_coordinates_nd, &
        spline_status_success, spline_status_invalid_input
    implicit none
    real(dp) :: x(7,8,4), coeff(7,8,4), oneaxis(7,8,4), twice(7,8,4)
    real(dp) :: scratch(7,8,4), gridvals(224), coord(3,224), query(3,5), vals(5), rawvals(5)
    real(dp) :: flat15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: output15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: coord15(15,3), vals15(3), vals2(3), coords2(2,3), sx(2,3), sc(2,3)
    real(dp) :: sub(4,4,4), copy_sub(4,4,4)
    complex(dp) :: z(7,8,4), cz(7,8,4), zvals(5)
    integer :: i,j,k,p,order,status

    do k=1,4
        do j=1,8
            do i=1,7
                x(i,j,k)=sin(0.17_dp*i+0.41_dp*j) + cos(0.51_dp*k-0.12_dp*j)
            end do
        end do
    end do
    z=cmplx(x,0.3_dp*x,dp)
    p=0
    do k=1,4
        do j=1,8
            do i=1,7
                p=p+1
                coord(:,p)=real([i,j,k],dp)
            end do
        end do
    end do
    query(:,1)=[1.0_dp,1.0_dp,1.0_dp]
    query(:,2)=[4.3_dp,2.5_dp,3.1_dp]
    query(:,3)=[7.0_dp,8.0_dp,4.0_dp]
    query(:,4)=[-3.25_dp,10.4_dp,1.9_dp]
    query(:,5)=[7.5_dp,-4.2_dp,5.8_dp]

    do order=0,5
        call spline_filter_nd(x,coeff,order,status)
        if (status/=spline_status_success) error stop 'real nd status'
        call map_coordinates_nd(coeff,coord,gridvals,order,status)
        if (status/=spline_status_success) error stop 'grid map status'
        if (maxval(abs(gridvals-reshape(x,[size(x)])))>3.e-12_dp) &
            error stop 'exact knots reconstruction'

        call spline_filter1d(x,oneaxis,axis=1,order=order,status=status)
        if (status/=spline_status_success) error stop 'axis 1 status'
        call spline_filter1d(oneaxis,twice,axis=2,order=order,status=status)
        if (status/=spline_status_success) error stop 'axis 2 status'
        call spline_filter1d(twice,scratch,axis=3,order=order,status=status)
        if (status/=spline_status_success) error stop 'axis 3 status'
        if (maxval(abs(coeff-scratch))>1.e-12_dp) error stop 'all-axis vs one-axis'

        call spline_filter_nd(z,cz,order,status)
        if (status/=spline_status_success) error stop 'complex nd status'
        if (maxval(abs(cz-cmplx(coeff,0.3_dp*coeff,dp)))>1.e-12_dp) &
            error stop 'complex coefficient linearity'
        call map_coordinates_nd(cz,query,zvals,order,status)
        if (status/=spline_status_success) error stop 'complex evaluation'
        call map_coordinates_nd(coeff,query,vals,order,status)
        if (status/=spline_status_success) error stop 'real evaluation'
        if (maxval(abs(zvals-cmplx(vals,0.3_dp*vals,dp)))>1.e-12_dp) &
            error stop 'complex evaluation linearity'
        call map_coordinates_nd(x,query,rawvals,order,status,prefilter=.true.)
        if (status/=spline_status_success) error stop 'prefilter evaluation'
        if (maxval(abs(vals-rawvals))>1.e-12_dp) error stop 'prefilter mismatch'
    end do

    call cspline_nd(x,scratch,status)
    call spline_filter_nd(x,coeff,order=3,status=status)
    if (maxval(abs(scratch-coeff))>1.e-12_dp) error stop 'cubic alias'
    call qspline_nd(x,scratch,status)
    call spline_filter_nd(x,coeff,order=2,status=status)
    if (maxval(abs(scratch-coeff))>1.e-12_dp) error stop 'quadratic alias'

    copy_sub=x(1:7:2,1:8:2,:)
    call spline_filter_nd(copy_sub,sub,order=3,status=status)
    if (status/=spline_status_success) error stop 'strided input'
    call spline_filter_nd(x(1:7:2,1:8:2,:),copy_sub,order=3,status=status)
    if (maxval(abs(sub-copy_sub))>1.e-12_dp) error stop 'strided input vs copy'
    call spline_filter_nd(x(1:7:2,1:8:2,:),scratch(1:7:2,1:8:2,:),order=3,status=status)
    if (maxval(abs(sub-scratch(1:7:2,1:8:2,:)))>1.e-12_dp) &
        error stop 'strided output'

    call map_coordinates_nd(x,query,vals,3,status,mode='constant',prefilter=.true.,cval=-12.5_dp)
    if (status/=spline_status_success) error stop 'constant mode status'
    if (abs(vals(4)+12.5_dp)>1.e-14_dp .or. abs(vals(5)+12.5_dp)>1.e-14_dp) &
        error stop 'constant out of bounds'
    call map_coordinates_nd(x,query,vals,order=9,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid order'
    call map_coordinates_nd(x,query(1:2,:),vals,order=3,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid coords rank'
    call spline_filter_nd(x,coeff,order=6,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid filter order'
    call spline_filter1d(x,coeff,axis=4,order=3,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid axis'
    call spline_filter_nd(x(1:3,:,:),coeff,order=3,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid output shape'

    sx=reshape([(real(i,dp),i=1,6)],[2,3])
    flat15=reshape(sx,shape(flat15))
    call cspline_nd(flat15,output15,status)
    if (status/=spline_status_success) error stop 'rank-15 coefficient status'
    call cspline_nd(sx,sc,status)
    if (abs(sum(output15)-sum(sc))>1.e-12_dp) error stop 'rank-15 coefficient sum'
    coords2(:,1)=[1.0_dp,1.0_dp]
    coords2(:,2)=[1.75_dp,2.6_dp]
    coords2(:,3)=[-4.75_dp,5.4_dp]
    coord15=1.0_dp
    coord15(1:2,:)=coords2
    call map_coordinates_nd(output15,coord15,vals15,3,status)
    if (status/=spline_status_success) error stop 'rank-15 map status'
    call map_coordinates_nd(sc,coords2,vals2,3,status)
    if (maxval(abs(vals15-vals2))>1.e-12_dp) error stop 'rank-15 vs rank-2'

    ! Every intermediate assumed rank 6-14 must use the same numerical kernel.
    block
        real(dp) :: data(2,3,1,1,1,1), coefficients(2,3,1,1,1,1)
        real(dp) :: q(6,3), v(3)
        data=reshape(sx,shape(data))
        q=1.0_dp
        q(1:2,:)=coords2
        call spline_filter_nd(data,coefficients,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 6 coefficient status'
        call map_coordinates_nd(coefficients,q,v,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 6 evaluation status'
        if (maxval(abs(v-vals2))>1.e-12_dp) error stop 'rank 6 mismatch'
    end block
    block
        real(dp) :: data(2,3,1,1,1,1,1), coefficients(2,3,1,1,1,1,1)
        real(dp) :: q(7,3), v(3)
        data=reshape(sx,shape(data))
        q=1.0_dp
        q(1:2,:)=coords2
        call spline_filter_nd(data,coefficients,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 7 coefficient status'
        call map_coordinates_nd(coefficients,q,v,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 7 evaluation status'
        if (maxval(abs(v-vals2))>1.e-12_dp) error stop 'rank 7 mismatch'
    end block
    block
        real(dp) :: data(2,3,1,1,1,1,1,1), coefficients(2,3,1,1,1,1,1,1)
        real(dp) :: q(8,3), v(3)
        data=reshape(sx,shape(data))
        q=1.0_dp
        q(1:2,:)=coords2
        call spline_filter_nd(data,coefficients,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 8 coefficient status'
        call map_coordinates_nd(coefficients,q,v,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 8 evaluation status'
        if (maxval(abs(v-vals2))>1.e-12_dp) error stop 'rank 8 mismatch'
    end block
    block
        real(dp) :: data(2,3,1,1,1,1,1,1,1), coefficients(2,3,1,1,1,1,1,1,1)
        real(dp) :: q(9,3), v(3)
        data=reshape(sx,shape(data))
        q=1.0_dp
        q(1:2,:)=coords2
        call spline_filter_nd(data,coefficients,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 9 coefficient status'
        call map_coordinates_nd(coefficients,q,v,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 9 evaluation status'
        if (maxval(abs(v-vals2))>1.e-12_dp) error stop 'rank 9 mismatch'
    end block
    block
        real(dp) :: data(2,3,1,1,1,1,1,1,1,1), coefficients(2,3,1,1,1,1,1,1,1,1)
        real(dp) :: q(10,3), v(3)
        data=reshape(sx,shape(data))
        q=1.0_dp
        q(1:2,:)=coords2
        call spline_filter_nd(data,coefficients,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 10 coefficient status'
        call map_coordinates_nd(coefficients,q,v,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 10 evaluation status'
        if (maxval(abs(v-vals2))>1.e-12_dp) error stop 'rank 10 mismatch'
    end block
    block
        real(dp) :: data(2,3,1,1,1,1,1,1,1,1,1), coefficients(2,3,1,1,1,1,1,1,1,1,1)
        real(dp) :: q(11,3), v(3)
        data=reshape(sx,shape(data))
        q=1.0_dp
        q(1:2,:)=coords2
        call spline_filter_nd(data,coefficients,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 11 coefficient status'
        call map_coordinates_nd(coefficients,q,v,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 11 evaluation status'
        if (maxval(abs(v-vals2))>1.e-12_dp) error stop 'rank 11 mismatch'
    end block
    block
        real(dp) :: data(2,3,1,1,1,1,1,1,1,1,1,1), coefficients(2,3,1,1,1,1,1,1,1,1,1,1)
        real(dp) :: q(12,3), v(3)
        data=reshape(sx,shape(data))
        q=1.0_dp
        q(1:2,:)=coords2
        call spline_filter_nd(data,coefficients,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 12 coefficient status'
        call map_coordinates_nd(coefficients,q,v,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 12 evaluation status'
        if (maxval(abs(v-vals2))>1.e-12_dp) error stop 'rank 12 mismatch'
    end block
    block
        real(dp) :: data(2,3,1,1,1,1,1,1,1,1,1,1,1), coefficients(2,3,1,1,1,1,1,1,1,1,1,1,1)
        real(dp) :: q(13,3), v(3)
        data=reshape(sx,shape(data))
        q=1.0_dp
        q(1:2,:)=coords2
        call spline_filter_nd(data,coefficients,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 13 coefficient status'
        call map_coordinates_nd(coefficients,q,v,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 13 evaluation status'
        if (maxval(abs(v-vals2))>1.e-12_dp) error stop 'rank 13 mismatch'
    end block
    block
        real(dp) :: data(2,3,1,1,1,1,1,1,1,1,1,1,1,1), coefficients(2,3,1,1,1,1,1,1,1,1,1,1,1,1)
        real(dp) :: q(14,3), v(3)
        data=reshape(sx,shape(data))
        q=1.0_dp
        q(1:2,:)=coords2
        call spline_filter_nd(data,coefficients,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 14 coefficient status'
        call map_coordinates_nd(coefficients,q,v,order=3,status=status)
        if (status/=spline_status_success) error stop 'rank 14 evaluation status'
        if (maxval(abs(v-vals2))>1.e-12_dp) error stop 'rank 14 mismatch'
    end block

    print '(a)', 'PASS: rank-15 spline filtering and coordinate evaluation'
end program test_signal_spline_nd
