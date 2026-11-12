! SPDX-License-Identifier: MIT
! Stream protocol for supplied Delaunay simplices in any geometric dimension.
program simplex_nd_driver
    use scifort_kinds, only: dp
    use scifort_interpolate_scattered, only: linear_nd_interpolator, &
        linear_nd_interpolator_type, complex_linear_nd_interpolator_type
    implicit none
    integer :: n,d,nch,ns,nq,i,status
    integer, allocatable :: simplices(:,:)
    real(dp), allocatable :: x(:,:),yreal(:,:),yimag(:,:),xi(:,:),vreal(:,:)
    complex(dp), allocatable :: yc(:,:),vc(:,:)
    type(linear_nd_interpolator_type) :: p
    type(complex_linear_nd_interpolator_type) :: pc
    read(*,*) n,d,nch,ns,nq
    allocate(x(n,d),yreal(n,nch),yimag(n,nch),xi(nq,d), &
             simplices(d+1,ns),yc(n,nch))
    do i=1,n
        read(*,*) x(i,:)
    end do
    do i=1,n
        read(*,*) yreal(i,:),yimag(i,:)
    end do
    do i=1,ns
        read(*,*) simplices(:,i)
    end do
    do i=1,nq
        read(*,*) xi(i,:)
    end do
    p=linear_nd_interpolator(x,yreal,simplices=simplices,rescale=.true., &
                             fill_value=-91.0_dp,status=status)
    write(*,'(i0)') status
    if (status/=0) stop
    vreal=p%evaluate(xi,status=status)
    write(*,'(i0)') status
    if (status/=0) stop
    yc=cmplx(yreal,yimag,dp)
    pc=linear_nd_interpolator(x,yc,simplices=simplices, &
                              fill_value=cmplx(-91.0_dp,17.0_dp,dp),status=status)
    write(*,'(i0)') status
    if (status/=0) stop
    vc=pc%evaluate(xi,status=status)
    write(*,'(i0)') status
    if (status/=0) stop
    do i=1,nq
        write(*,'(*(es26.17e3,1x))') vreal(i,:), real(vc(i,:),dp),aimag(vc(i,:))
    end do
end program simplex_nd_driver
