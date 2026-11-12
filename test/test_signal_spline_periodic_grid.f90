! SPDX-License-Identifier: MIT
! Periodic IIR prefilter, explicit boundary modes and shaped coordinate grids.
program test_signal_spline_periodic_grid
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_filter_nd, spline_filter1d, &
        map_coordinates_nd, map_coordinates_grid_nd
    use scifort_signal_spline, only : spline_status_success, spline_status_invalid_input
    implicit none
    real(dp) :: x(11), cf(11), af(11), rec(11), coords1(1,11)
    real(dp) :: input(7,9), coord(2,5,4), out(5,4), flat_vals(20), &
        out_strided(10,8), coord_strided(2,10,8), rank15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: bad_coord(3,5,4), z_coord(15,3), rank15_out(3)
    complex(dp) :: z(7,9), z_out(5,4), z_flat(20)
    integer :: i, j, degree, status, n

    do i=1,11
        x(i)=sin(0.27_dp*i)+0.08_dp*i
        coords1(1,i)=real(i,dp)
    end do
    do n=2,11
        do degree=2,5
            call spline_filter_nd(x(:n),cf(:n),order=degree,status=status,mode='grid-wrap')
            if (status/=spline_status_success) error stop 'periodic filter status'
            call spline_filter1d(x(:n),af(:n),1,degree,status,mode='grid-wrap')
            if (status/=spline_status_success) error stop 'periodic axis filter status'
            if (maxval(abs(af(:n)-cf(:n)))>1.e-13_dp) &
                error stop 'periodic one axis vs all axes'
            call map_coordinates_nd(cf(:n),coords1(:,:n),rec(:n),degree,status,mode='grid-wrap')
            if (status/=spline_status_success) error stop 'periodic reconstruction status'
            if (maxval(abs(rec(:n)-x(:n)))>1.e-11_dp) &
                error stop 'periodic reconstruction'
            call spline_filter1d(x(:n),af(:n),1,degree,status,mode='reflect')
            if (status/=spline_status_success) error stop 'reflect filter status'
            call map_coordinates_nd(af(:n),coords1(:,:n),rec(:n),degree,status,mode='reflect')
            if (maxval(abs(rec(:n)-x(:n)))>1.e-11_dp) &
                error stop 'exact reflect reconstruction'
            call spline_filter1d(x(:n),cf(:n),1,degree,status,mode='grid-mirror')
            if (maxval(abs(cf(:n)-af(:n)))>1.e-13_dp) &
                error stop 'grid-mirror filter alias'
            call spline_filter1d(x(:n),cf(:n),1,degree,status,mode='nearest')
            if (maxval(abs(cf(:n)-af(:n)))>1.e-13_dp) &
                error stop 'nearest spline filter convention'
        end do
    end do
    call spline_filter_nd(x,cf,order=3,status=status,mode='invalid')
    if (status/=spline_status_invalid_input) error stop 'invalid spline filter mode'

    do j=1,9
        do i=1,7
            input(i,j)=sin(.2_dp*i+.45_dp*j)+.02_dp*i*j
        end do
    end do
    z=cmplx(input,.3_dp*input+0.1_dp,dp)
    do j=1,4
        do i=1,5
            coord(:,i,j)=[0.8_dp+.9_dp*i+.22_dp*j,1.0_dp+.42_dp*i+.95_dp*j]
        end do
    end do
    call map_coordinates_grid_nd(input,coord,out,order=5,mode='grid-wrap',status=status)
    if (status/=spline_status_success) error stop 'real shaped grid status'
    call map_coordinates_nd(input,reshape(coord,[2,20]),flat_vals,order=5, &
        mode='grid-wrap',prefilter=.true.,status=status)
    if (maxval(abs(out-reshape(flat_vals,[5,4])))>2.e-12_dp) &
        error stop 'real grid vs flat coordinates'
    call map_coordinates_grid_nd(input,coord,out_strided(1:10:2,1:8:2), &
        order=5,mode='grid-wrap',status=status)
    if (status/=spline_status_success) error stop 'strided grid status'
    if (maxval(abs(out-out_strided(1:10:2,1:8:2)))>2.e-12_dp) &
        error stop 'strided grid output'
    coord_strided=0.0_dp
    coord_strided(:,1:10:2,1:8:2)=coord
    call map_coordinates_grid_nd(input,coord_strided(:,1:10:2,1:8:2), &
        out,order=5,mode='grid-wrap',status=status)
    if (status/=spline_status_success) error stop 'strided coordinates status'
    if (maxval(abs(out-reshape(flat_vals,[5,4])))>2.e-12_dp) &
        error stop 'strided coordinates values'
    call map_coordinates_grid_nd(z,coord,z_out,order=3,mode='grid-constant', &
        cval=(-2.0_dp,1.0_dp),status=status)
    if (status/=spline_status_success) error stop 'complex shaped grid status'
    call map_coordinates_nd(z,reshape(coord,[2,20]),z_flat,order=3, &
        mode='grid-constant',prefilter=.true.,cval=(-2.0_dp,1.0_dp),status=status)
    if (maxval(abs(z_out-reshape(z_flat,[5,4])))>2.e-12_dp) &
        error stop 'complex grid vs flat coordinates'

    rank15=1.0_dp
    rank15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=2.0_dp
    z_coord=1.0_dp
    z_coord(1,:)=[1.0_dp,1.5_dp,2.0_dp]
    z_coord(2,:)=[1.0_dp,2.5_dp,3.0_dp]
    call map_coordinates_grid_nd(rank15,z_coord,rank15_out,order=3, &
        mode='mirror',status=status)
    if (status/=spline_status_success) error stop 'rank15 grid status'
    if (any(abs(rank15_out)>3.0_dp)) error stop 'rank15 grid magnitude'

    bad_coord=1.0_dp
    call map_coordinates_grid_nd(input,bad_coord,out,order=3,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid grid coord components'
    call map_coordinates_grid_nd(input,coord,out_strided,order=3,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid grid shape'
    call map_coordinates_grid_nd(input,coord,out,order=9,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid grid order'
    coord(1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call map_coordinates_grid_nd(input,coord,out,status=status)
    if (status/=spline_status_invalid_input) error stop 'nonfinite grid coordinates'
    print '(a)', 'PASS: periodic prefilter and shaped coordinate grids'
end program test_signal_spline_periodic_grid
