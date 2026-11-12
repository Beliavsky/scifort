! SPDX-License-Identifier: MIT
! Boundary modes and affine/shift transforms using rank-independent splines.
program test_signal_spline_modes_transform
    use scifort_signal, only : dp, map_coordinates_nd, affine_transform_nd, &
        shift_nd, spline_filter_nd, spline_status_success, spline_status_invalid_input, &
        spline_nd_status_too_expensive
    implicit none
    character(len=16), parameter :: modes(8) = [character(len=16) :: &
        'mirror','constant','reflect','grid-mirror','nearest', &
        'wrap','grid-wrap','grid-constant']
    real(dp) :: x(5,6),y(4,7),mat(2,2),offs(2),coord(2,28),v(28)
    real(dp) :: strided(8,14),sh(5,6),id(2,2),ay(5,6)
    real(dp) :: h15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: out15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: c15(15,15),offset15(15),displ(2)
    real(dp) :: five_d(2,2,2,2,2), coords5(5,1), five_out(1)
    complex(dp) :: z(5,6),w(4,7),vc(28),shc(5,6),shca(5,6),c15input(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: c15output(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: i,j,k,m,order,status

    do j=1,6
        do i=1,5
            x(i,j)=sin(0.3_dp*i+0.27_dp*j)+cos(0.17_dp*j-0.22_dp*i)
        end do
    end do
    z=cmplx(x,0.21_dp+0.41_dp*x,dp)
    mat(1,:)=[0.83_dp,0.18_dp]
    mat(2,:)=[-0.11_dp,0.92_dp]
    offs=[-0.7_dp,0.31_dp]
    id=0.0_dp
    id(1,1)=1.0_dp
    id(2,2)=1.0_dp
    displ=[0.28_dp,-1.22_dp]

    do j=1,7
        do i=1,4
            k=i+4*(j-1)
            coord(:,k)=1.0_dp+offs+matmul(mat,real([i-1,j-1],dp))
        end do
    end do
    do m=1,size(modes)
        do order=0,5
            call affine_transform_nd(x,mat,offs,y,order,status,mode=trim(modes(m)))
            if (status/=spline_status_success) error stop 'real affine transform status'
            call map_coordinates_nd(x,coord,v,order,status,mode=trim(modes(m)))
            if (status/=spline_status_success) error stop 'real coordinate status'
            if (maxval(abs(reshape(y,[28])-v))>3.e-12_dp) error stop 'affine vs map'

            call affine_transform_nd(z,mat,offs,w,order,status,mode=trim(modes(m)))
            if (status/=spline_status_success) error stop 'complex affine status'
            call map_coordinates_nd(z,coord,vc,order,status,mode=trim(modes(m)))
            if (maxval(abs(reshape(w,[28])-vc))>3.e-12_dp) &
                error stop 'complex affine vs map'
            call affine_transform_nd(x,mat,offs,strided(1:8:2,1:14:2), &
                order,status,mode=trim(modes(m)))
            if (status/=spline_status_success) error stop 'strided affine status'
            if (maxval(abs(y-strided(1:8:2,1:14:2)))>3.e-12_dp) &
                error stop 'strided affine output'

            call shift_nd(x,displ,sh,order,status,mode=trim(modes(m)))
            if (status/=spline_status_success) error stop 'shift status'
            call affine_transform_nd(x,id,-displ,ay,order,status,mode=trim(modes(m)))
            if (status/=spline_status_success) error stop 'shift affine status'
            if (maxval(abs(sh-ay))>3.e-12_dp) error stop 'shift vs affine'
            call shift_nd(z,displ,shc,order,status,mode=trim(modes(m)))
            if (status/=spline_status_success) error stop 'complex shift status'
            call affine_transform_nd(z,id,-displ,shca,order,status,mode=trim(modes(m)))
            if (status/=spline_status_success) error stop 'complex shift affine status'
            if (maxval(abs(shc-shca))>3.e-12_dp) &
                error stop 'complex shift vs affine'
        end do
    end do

    ! Validate prefilter=True on all boundary modes at integer knots.
    do m=1,8
        call affine_transform_nd(x,id,[0.0_dp,0.0_dp],ay,order=3,status=status, &
            mode=trim(modes(m)),prefilter=.true.)
        if (status/=spline_status_success) error stop 'prefilter affine status'
        if (maxval(abs(x-ay))>3.e-12_dp) error stop 'prefilter exact knots'
    end do

    call affine_transform_nd(x,id,[0.0_dp,0.0_dp],ay,order=3,status=status, &
        mode='nearest',prefilter=.true.)
    if (status/=spline_status_success) error stop 'nearest prepadding status'
    call affine_transform_nd(x,id,[0.0_dp,0.0_dp],ay,order=3,status=status, &
        mode='grid-constant',prefilter=.true.)
    if (status/=spline_status_success) error stop 'grid constant prepadding status'
    call affine_transform_nd(x,id,[0.0_dp],ay,order=1,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid affine offset'
    call affine_transform_nd(x,mat(:,1:1),offs,ay,order=1,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid affine matrix'
    call affine_transform_nd(x,id,offs,ay,order=7,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid affine order'
    call affine_transform_nd(x,id,offs,ay,mode='invalid',status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid mode'
    call shift_nd(x,[1.0_dp],ay,status=status)
    if (status/=spline_status_invalid_input) error stop 'invalid shift shape'

    h15=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],shape(h15))
    c15input=cmplx(h15,0.28_dp*h15,dp)
    c15=0.0_dp
    offset15=0.0_dp
    do i=1,15
        c15(i,i)=1.0_dp
    end do
    call affine_transform_nd(h15,c15,offset15,out15,order=3, &
        prefilter=.true.,status=status,mode='grid-wrap')
    if (status/=spline_status_success) error stop 'rank-15 affine status'
    if (maxval(abs(h15-out15))>1.e-11_dp) error stop 'rank-15 reconstruction'
    do m=1,2
        if (m==1) then
            call affine_transform_nd(h15,c15,offset15,out15,order=5, &
                prefilter=.true.,status=status,mode='nearest')
        else
            call affine_transform_nd(h15,c15,offset15,out15,order=5, &
                prefilter=.true.,status=status,mode='grid-constant',cval=-1.25_dp)
        end if
        if (status/=spline_status_success) error stop 'rank-15 prepad status'
        if (maxval(abs(h15-out15))>5.e-11_dp) error stop 'rank-15 prepad knots'
    end do
    call shift_nd(c15input,offset15,c15output,order=1,status=status)
    if (status/=spline_status_success) error stop 'rank-15 complex shift status'
    if (maxval(abs(c15input-c15output))>1.e-13_dp) &
        error stop 'rank-15 complex shift values'

    five_d=1.0_dp
    coords5=1.0_dp
    call map_coordinates_nd(five_d,coords5,five_out,order=5, &
        prefilter=.true.,status=status,mode='grid-constant')
    if (status/=spline_nd_status_too_expensive) &
        error stop 'explicit prepad allocation guard'

    print '(a)', 'PASS: spline boundary modes and affine/shift transforms'
end program test_signal_spline_modes_transform
