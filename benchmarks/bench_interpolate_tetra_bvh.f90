! SPDX-License-Identifier: MIT
! Repeatable end-to-end benchmark of 3D LinearNDInterpolator construction
! and batch evaluation.  Compile against the current or predecessor source.
program bench_interpolate_tetra_bvh
    use scifort_kinds, only: dp
    use scifort_interpolate_scattered, only: linear_nd_interpolator, linear_nd_interpolator_type
    implicit none
    integer :: nsite, nquery, i, nseed, status, ntetra
    integer, allocatable :: seed(:), simplex(:,:)
    real(dp), allocatable :: x(:,:), q(:,:), y(:,:), z(:,:)
    real(dp) :: t0,t1,t2,checksum,err
    type(linear_nd_interpolator_type) :: interpolator
    character(len=40) :: arg

    nsite=160
    nquery=12000
    if (command_argument_count()>=1) then
        call get_command_argument(1,arg)
        read(arg,*) nsite
    end if
    if (command_argument_count()>=2) then
        call get_command_argument(2,arg)
        read(arg,*) nquery
    end if
    if (nsite<4 .or. nquery<1) error stop 'invalid benchmark dimensions'
    allocate(x(nsite,3),q(nquery,3),y(nsite,3))
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed=20261031
    call random_seed(put=seed)
    call random_number(x)
    x(:,1)=1.5_dp*x(:,1)
    x(:,2)=0.9_dp*x(:,2)
    x(:,3)=0.6_dp*x(:,3)
    y(:,1)=sin(x(:,1)+2.0_dp*x(:,2))+x(:,3)**2
    y(:,2)=cos(3.0_dp*x(:,3)-x(:,1))+x(:,2)**2
    y(:,3)=x(:,1)*x(:,2)-exp(-x(:,3))
    call random_number(q)
    q(:,1)=1.5_dp*(-0.03_dp+1.06_dp*q(:,1))
    q(:,2)=0.9_dp*(-0.03_dp+1.06_dp*q(:,2))
    q(:,3)=0.6_dp*(-0.03_dp+1.06_dp*q(:,3))
    do i=1,min(8,nquery)
        q(i,:)=x(i,:)
    end do
    call cpu_time(t0)
    interpolator=linear_nd_interpolator(x,y,fill_value=-12.0_dp,status=status)
    call cpu_time(t1)
    if (status/=0 .or. .not.interpolator%valid()) error stop 'bench construction'
    simplex=interpolator%simplices()
    ntetra=size(simplex,2)
    z=interpolator%evaluate(q,status=status)
    call cpu_time(t2)
    if (status/=0) error stop 'bench evaluation'
    err=0.0_dp
    do i=1,min(8,nquery)
        err=max(err,maxval(abs(z(i,:)-y(i,:))))
    end do
    if (err>2.0e-10_dp) error stop 'bench sample-site interpolation'
    checksum=sum(z)
    write(*,'(a,i0)') 'tetrahedra=',ntetra
    write(*,'(a,es24.15)') 'construct_s=',t1-t0
    write(*,'(a,es24.15)') 'evaluate_s=',t2-t1
    write(*,'(a,es24.15)') 'checksum=',checksum
    write(*,'(a,i0)') 'filled=',count(z(:,1)==-12.0_dp)
    write(*,'(a,es24.15)') 'vertex_error=',err
end program bench_interpolate_tetra_bvh
