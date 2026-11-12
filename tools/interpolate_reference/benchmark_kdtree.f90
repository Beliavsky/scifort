! SPDX-License-Identifier: MIT
! Deterministic single-run benchmark for 2D and higher-dimensional nearest searches.
program benchmark_kdtree
    use scifort_kinds, only: dp
    use scifort_interpolate_kdtree, only: kd_tree_type, kd_tree
    implicit none
    type(kd_tree_type) :: index
    integer :: n, nq, dim, i, j, best, result, visits, total_visits, sz
    integer, allocatable :: seed(:), expected(:)
    real(dp), allocatable :: points(:,:), queries(:,:)
    real(dp) :: t0,t1,t2,t3,dmin,d2,build_time,kd_time,linear_time
    character(len=32) :: argument
    n=12000
    nq=1500
    dim=2
    if (command_argument_count()>=1) then
        call get_command_argument(1,argument)
        read(argument,*) n
    end if
    if (command_argument_count()>=2) then
        call get_command_argument(2,argument)
        read(argument,*) nq
    end if
    if (command_argument_count()>=3) then
        call get_command_argument(3,argument)
        read(argument,*) dim
    end if
    if (n<1 .or. nq<1 .or. dim<1) error stop 'positive sizes required'
    allocate(points(n,dim),queries(nq,dim),expected(nq))
    call random_seed(size=sz)
    allocate(seed(sz))
    seed=20261031
    call random_seed(put=seed)
    call random_number(points)
    call random_number(queries)
    call cpu_time(t0)
    index=kd_tree(points)
    call cpu_time(t1)
    total_visits=0
    do i=1,nq
        result=index%nearest(queries(i,:),visited=visits)
        total_visits=total_visits+visits
        if (result<1 .or. result>n) error stop 'invalid nearest site'
    end do
    call cpu_time(t2)
    do i=1,nq
        best=1
        dmin=huge(1.0_dp)
        do j=1,n
            d2=sum((points(j,:)-queries(i,:))**2)
            if (d2<dmin) then
                best=j
                dmin=d2
            end if
        end do
        expected(i)=best
    end do
    call cpu_time(t3)
    do i=1,nq
        result=index%nearest(queries(i,:))
        if (expected(i)/=result) error stop 'benchmark correctness: nearest index differs'
    end do
    build_time=t1-t0
    kd_time=t2-t1
    linear_time=t3-t2
    write(*,'(a,i0)') 'sites=',n
    write(*,'(a,i0)') 'queries=',nq
    write(*,'(a,i0)') 'dimensions=',dim
    write(*,'(a,es14.5)') 'build seconds=',build_time
    write(*,'(a,es14.5)') 'kdtree query seconds=',kd_time
    write(*,'(a,es14.5)') 'brute force query seconds=',linear_time
    write(*,'(a,f10.2)') 'mean exact sites visited=',real(total_visits,dp)/real(nq,dp)
    write(*,'(a,f12.6)') 'mean site fraction visited=',real(total_visits,dp)/(real(nq,dp)*real(n,dp))
end program benchmark_kdtree
