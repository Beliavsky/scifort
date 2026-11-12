! SPDX-License-Identifier: MIT
! Rank-independent SciFort zoom and rotate regression.
program test_signal_zoom_rotate_nd
    use scifort_signal, only : dp, zoom_nd, rotate_nd, zoom_shape_nd, &
        rotate_shape_nd, spline_status_success, spline_status_invalid_input, &
        affine_transform_nd
    implicit none
    real(dp) :: x(5,6), x3(5,3,4), full(14,12)
    real(dp) :: y(7,4), yscalar(10,12)
    real(dp) :: h15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        out15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: z15(15)
    complex(dp) :: cx(5,6), cy(7,4), cz15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: co15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: turn15(3,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp), allocatable :: rr(:,:,:), rot(:,:)
    integer, allocatable :: dims(:)
    integer :: i,j,k,status,ang

    do j=1,6
        do i=1,5
            x(i,j)=sin(0.3_dp*i+0.4_dp*j)+0.2_dp*i
        end do
    end do
    cx=cmplx(x,0.27_dp*x+0.13_dp,dp)
    dims=zoom_shape_nd(shape(x),[1.4_dp,0.7_dp],status)
    if(status/=spline_status_success .or. any(dims/=[7,4])) error stop 'zoom shape'
    call zoom_nd(x,[1.4_dp,0.7_dp],y,order=3,mode='mirror', &
        prefilter=.true.,status=status)
    if(status/=spline_status_success) error stop 'real zoom'
    call zoom_nd(cx,[1.4_dp,0.7_dp],cy,order=3,mode='mirror', &
        prefilter=.true.,status=status)
    if(status/=spline_status_success) error stop 'complex zoom'
    if(maxval(abs(real(cy,dp)-y))>1.e-12_dp) error stop 'complex zoom real component'
    call zoom_nd(x,2.0_dp,yscalar,order=1,status=status)
    if(status/=spline_status_success) error stop 'scalar zoom'
    full=0.0_dp
    call zoom_nd(x,[1.4_dp,0.7_dp],full(1:14:2,1:8:2), &
        order=3,prefilter=.true.,status=status)
    if(status/=spline_status_success) error stop 'noncontiguous zoom'
    if(maxval(abs(full(1:14:2,1:8:2)-y))>1.e-12_dp) &
        error stop 'noncontiguous zoom result'

    dims=zoom_shape_nd([5],0.5_dp,status)
    if(status/=spline_status_success .or. any(dims/=[2])) error stop 'ties to even (2.5)'
    dims=zoom_shape_nd([7],0.5_dp,status)
    if(status/=spline_status_success .or. any(dims/=[4])) error stop 'ties to even (3.5)'
    dims=zoom_shape_nd([5,7],[-1.0_dp,1.0_dp],status)
    if(status/=spline_status_invalid_input) error stop 'negative zoom factor'
    call zoom_nd(x,[1.4_dp],y,status=status)
    if(status/=spline_status_invalid_input) error stop 'zoom rank mismatch'
    call zoom_nd(x,[1.4_dp,0.7_dp],yscalar,status=status)
    if(status/=spline_status_invalid_input) error stop 'zoom output mismatch'

    ! The input-first-axis/last-axis plane is selected with 1-based axes.
    do k=1,4
        do j=1,3
            do i=1,5
                x3(i,j,k)=sin(0.23_dp*i+0.4_dp*j-0.3_dp*k)
            end do
        end do
    end do
    do ang=-90,90,45
        dims=rotate_shape_nd(shape(x3),real(ang,dp),axes=[3,1],status=status)
        if(status/=spline_status_success) error stop 'rotate shape'
        allocate(rr(dims(1),dims(2),dims(3)))
        call rotate_nd(x3,real(ang,dp),rr,axes=[3,1],order=1,status=status)
        if(status/=spline_status_success) error stop 'rotate 3d'
        deallocate(rr)
    end do
    dims=rotate_shape_nd(shape(x),90.0_dp,status=status)
    if(any(dims/=[6,5])) error stop 'right angle shape'
    dims=rotate_shape_nd(shape(x),32.0_dp,reshape=.false.,status=status)
    if(any(dims/=[5,6])) error stop 'fixed rotate shape'
    dims=rotate_shape_nd(shape(x),32.0_dp,axes=[1,1],status=status)
    if(status/=spline_status_invalid_input) error stop 'repeated axes'
    allocate(rot(6,5))
    call rotate_nd(x,90.0_dp,rot,order=0,status=status,prefilter=.false.)
    if(status/=spline_status_success) error stop 'rotate right angle'
    do j=1,5
        do i=1,6
            if(abs(rot(i,j)-x(j,7-i))>1.e-12_dp) &
                error stop 'ninety degree orientation'
        end do
    end do
    call rotate_nd(x,90.0_dp,y,status=status)
    if(status/=spline_status_invalid_input) error stop 'rotate output mismatch'

    h15=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],shape(h15))
    cz15=cmplx(h15,h15*.2_dp,dp)
    z15=1.0_dp
    dims=zoom_shape_nd(shape(h15),z15,status)
    if(status/=spline_status_success .or. any(dims/=shape(h15))) &
        error stop 'rank15 zoom shape'
    call zoom_nd(h15,z15,out15,order=1,status=status)
    if(status/=spline_status_success .or. maxval(abs(out15-h15))>1.e-12_dp) &
        error stop 'rank15 real zoom'
    call zoom_nd(cz15,1.0_dp,co15,order=1,status=status)
    if(status/=spline_status_success .or. maxval(abs(co15-cz15))>1.e-12_dp) &
        error stop 'rank15 complex scalar zoom'
    dims=rotate_shape_nd(shape(h15),0.0_dp,axes=[1,2],status=status)
    if(any(dims/=shape(h15))) error stop 'rank15 rotation shape'
    call rotate_nd(cz15,0.0_dp,co15,axes=[1,2],order=1,status=status)
    if(status/=spline_status_success .or. maxval(abs(co15-cz15))>1.e-12_dp) &
        error stop 'rank15 complex rotation'
    dims=rotate_shape_nd(shape(cz15),90.0_dp,axes=[2,1],status=status)
    if(status/=spline_status_success .or. any(dims/=shape(turn15))) &
        error stop 'rank15 right-angle shape'
    call rotate_nd(cz15,90.0_dp,turn15,axes=[2,1],order=1,status=status)
    if(status/=spline_status_success) error stop 'rank15 right-angle rotation'
    if(maxval(abs(reshape(turn15,[6])-cmplx([5.0_dp,3.0_dp,1.0_dp, &
        6.0_dp,4.0_dp,2.0_dp],.2_dp*[5.0_dp,3.0_dp,1.0_dp, &
        6.0_dp,4.0_dp,2.0_dp],dp)))>1.e-12_dp) &
        error stop 'rank15 right-angle contents'

    print '(a)', 'PASS: zoom_nd and rotate_nd ranks 2, 3, 15'
end program test_signal_zoom_rotate_nd
