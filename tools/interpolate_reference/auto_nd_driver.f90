! SPDX-License-Identifier: MIT
! Standard-input reference driver for native 4D/5D Delaunay and LinearND.
program auto_nd_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_scattered, only : linear_nd_interpolator, complex_linear_nd_interpolator_type, &
        delaunay_nd
    implicit none
    integer :: n,d,nchan,nq,rescale,i,j,code
    integer, allocatable :: cells(:,:)
    real(dp), allocatable :: sites(:,:),queries(:,:),vr(:,:),vi(:,:)
    complex(dp), allocatable :: values(:,:),result(:,:)
    type(complex_linear_nd_interpolator_type) :: interp
    read(*,*) n,d,nchan,nq,rescale
    allocate(sites(n,d),values(n,nchan),queries(nq,d),vr(n,nchan),vi(n,nchan))
    do i=1,n
        read(*,*) sites(i,:)
    end do
    do i=1,n
        read(*,*) vr(i,:),vi(i,:)
    end do
    values=cmplx(vr,vi,dp)
    do i=1,nq
        read(*,*) queries(i,:)
    end do
    call delaunay_nd(sites,cells,code)
    write(*,'(i0,1x,i0)') code, size(cells,2)
    if (code/=0) stop
    do i=1,size(cells,2)
        write(*,'(*(i0,1x))') cells(:,i)
    end do
    interp=linear_nd_interpolator(sites,values,fill_value=cmplx(-99.0_dp,7.0_dp,dp), &
                                   rescale=(rescale==1),status=code)
    write(*,'(i0)') code
    if (code/=0) stop
    result=interp%evaluate(queries,status=code)
    write(*,'(i0)') code
    do i=1,nq
        write(*,'(*(es24.16e3,1x))') real(result(i,:)),aimag(result(i,:))
    end do
end program auto_nd_driver
