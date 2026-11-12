! SPDX-License-Identifier: MIT
! Deterministic binary reference output for scipy.ndimage.zoom/rotate.
program zoom_rotate_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : zoom_nd, rotate_nd, zoom_shape_nd, rotate_shape_nd
    use scifort_signal_spline, only : spline_status_success
    implicit none
    character(len=16), parameter :: modes(7) = [character(len=16) :: &
        'mirror','constant','nearest','reflect','grid-wrap','grid-constant','wrap']
    real(dp), parameter :: factors(2)=[1.4_dp,0.7_dp]
    integer, parameter :: orders(4)=[0,1,3,5], angles(4)=[0,23,90,-37]
    real(dp) :: x(5,6),xx(5,3,4),cvalue
    complex(dp) :: cx(5,6),cx3(5,3,4),ccvalue
    real(dp), allocatable :: y(:,:), rot(:,:), rot3(:,:,:)
    complex(dp), allocatable :: yc(:,:), rotc(:,:), rotc3(:,:,:)
    integer, allocatable :: dims(:)
    integer :: i,j,k,l,m,o,g,f,ang,rs,status,u
    logical :: pre,grid,grow

    cvalue=0.35_dp
    ccvalue=cmplx(0.35_dp,0.21_dp,dp)
    do j=1,6
        do i=1,5
            x(i,j)=sin(.31_dp*i+.27_dp*j)+.13_dp*i-.04_dp*j
        end do
    end do
    cx=cmplx(x,.27_dp*x-.19_dp,dp)
    do k=1,4
        do j=1,3
            do i=1,5
                xx(i,j,k)=sin(.31_dp*i+.17_dp*j-.29_dp*k)+.13_dp*i
            end do
        end do
    end do
    cx3=cmplx(xx,.27_dp*xx-.19_dp,dp)
    open(newunit=u,file='zoom_rotate_results.bin',access='stream',&
         form='unformatted',status='replace',action='write')
    dims=zoom_shape_nd(shape(x),factors,status)
    if(status/=spline_status_success) error stop 'zoom shape'
    allocate(y(dims(1),dims(2)),yc(dims(1),dims(2)))
    do m=1,size(modes)
        do o=1,size(orders)
            do g=0,1
                grid=(g==1)
                do f=0,1
                    pre=(f==1)
                    call zoom_nd(x,factors,y,order=orders(o),status=status,&
                        grid_mode=grid,mode=trim(modes(m)),prefilter=pre,cval=cvalue)
                    if(status/=spline_status_success) error stop 'real zoom'
                    write(u) y
                    call zoom_nd(cx,factors,yc,order=orders(o),status=status,&
                        grid_mode=grid,mode=trim(modes(m)),prefilter=pre,cval=ccvalue)
                    if(status/=spline_status_success) error stop 'complex zoom'
                    write(u) yc
                end do
            end do
        end do
    end do
    deallocate(y,yc)

    do m=1,size(modes)
        do o=1,size(orders)
            do l=1,size(angles)
                do rs=0,1
                    grow=(rs==1)
                    dims=rotate_shape_nd(shape(x),real(angles(l),dp),reshape=grow,status=status)
                    if(status/=spline_status_success) error stop 'rotate shape'
                    allocate(rot(dims(1),dims(2)),rotc(dims(1),dims(2)))
                    do f=0,1
                        pre=(f==1)
                        call rotate_nd(x,real(angles(l),dp),rot,reshape=grow,&
                            order=orders(o),status=status,mode=trim(modes(m)), &
                            prefilter=pre,cval=cvalue)
                        if(status/=spline_status_success) error stop 'real rotate'
                        write(u) rot
                        call rotate_nd(cx,real(angles(l),dp),rotc,reshape=grow,&
                            order=orders(o),status=status,mode=trim(modes(m)), &
                            prefilter=pre,cval=ccvalue)
                        if(status/=spline_status_success) error stop 'complex rotate'
                        write(u) rotc
                    end do
                    deallocate(rot,rotc)
                end do
            end do
        end do
    end do

    do m=1,2
        do o=1,3
            do l=1,3
                dims=rotate_shape_nd(shape(xx),real(angles(l),dp),axes=[1,3],&
                    reshape=(m==1),status=status)
                if(status/=spline_status_success) error stop '3d shape'
                allocate(rot3(dims(1),dims(2),dims(3)), &
                    rotc3(dims(1),dims(2),dims(3)))
                call rotate_nd(xx,real(angles(l),dp),rot3,axes=[1,3], &
                    reshape=(m==1),order=orders(o),status=status, &
                    mode='mirror',prefilter=.true.)
                if(status/=spline_status_success) error stop '3d rotate'
                write(u) rot3
                call rotate_nd(cx3,real(angles(l),dp),rotc3,axes=[1,3], &
                    reshape=(m==1),order=orders(o),status=status, &
                    mode='mirror',prefilter=.true.)
                if(status/=spline_status_success) error stop '3d complex rotate'
                write(u) rotc3
                deallocate(rot3,rotc3)
            end do
        end do
    end do
    close(u)
    print '(a)', 'PASS: zoom/rotate reference driver'
end program zoom_rotate_reference_driver
