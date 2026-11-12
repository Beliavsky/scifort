! SPDX-License-Identifier: MIT
! CPU-time benchmark for rank-generic supplied-simplex lookup.
program bench_interpolate_simplex_bvh
    use scifort_kinds, only : dp
    use scifort_interpolate_scattered, only : linear_nd_interpolator, linear_nd_interpolator_type
    implicit none
    integer :: n,d,ns,nq,i,code,ntimes
    integer, allocatable :: simplices(:,:)
    real(dp), allocatable :: points(:,:),data(:),query(:,:),values(:,:)
    real(dp) :: t0,t1,checksum
    type(linear_nd_interpolator_type) :: obj
    read(*,*) n,d,ns,nq,ntimes
    allocate(points(n,d),data(n),simplices(d+1,ns),query(nq,d))
    do i=1,n
        read(*,*) points(i,:)
    end do
    read(*,*) data
    do i=1,ns
        read(*,*) simplices(:,i)
    end do
    do i=1,nq
        read(*,*) query(i,:)
    end do
    obj=linear_nd_interpolator(points,data,simplices=simplices,fill_value=-99.0_dp,status=code)
    if (code/=0) error stop 'benchmark interpolation construction failure'
    checksum=0.0_dp
    call cpu_time(t0)
    do i=1,ntimes
        values=obj%evaluate(query,status=code)
        if (code/=0) error stop 'benchmark interpolation evaluation failure'
        checksum=checksum+sum(values)
    end do
    call cpu_time(t1)
    print '(f16.9,1x,es24.16e3)',t1-t0,checksum
end program
