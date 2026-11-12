! SPDX-License-Identifier: MIT
! Standalone machine-readable comparison protocol for SciPy KDTree.
program kdtree_periodic_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_spatial, only : kd_tree, kd_tree_type, kdtree_success
    implicit none
    integer :: n,m,d,i,k,code
    integer, allocatable :: ids(:),ball(:),pairs(:,:),rp(:),columns(:),row(:),col(:)
    integer(int64) :: count
    real(dp) :: r,p
    real(dp), allocatable :: xa(:,:),xb(:,:),box(:),q(:),dist(:),values(:)
    character(len=32) :: ptext
    type(kd_tree_type) :: a,b
    read(*,*) n,m,d,r,ptext
    if (trim(ptext)=='inf') then
        p=ieee_value(0.0_dp,ieee_positive_inf)
    else
        read(ptext,*) p
    end if
    allocate(xa(n,d),xb(m,d),box(d),q(d))
    read(*,*) box
    read(*,*) q
    do i=1,n
        read(*,*) xa(i,:)
    end do
    do i=1,m
        read(*,*) xb(i,:)
    end do
    a=kd_tree(xa,boxsize=box,status=code)
    if (code/=kdtree_success) error stop 'bad first tree'
    b=kd_tree(xb,boxsize=box,status=code)
    if (code/=kdtree_success) error stop 'bad second tree'
    k=min(6,n)
    call a%query(q,ids,k=k,p=p,distance=dist,status=code)
    if (code/=kdtree_success) error stop 'bad query'
    print '(a,1x,i0)', 'Q',k
    do i=1,k
        print '(i0,1x,es27.18e3)', ids(i),dist(i)
    end do
    call a%query_ball_point(q,r,ball,p=p,status=code)
    if (code/=kdtree_success) error stop 'bad ball_point'
    print '(a,1x,i0)', 'B',size(ball)
    do i=1,size(ball)
        print '(i0)', ball(i)
    end do
    call a%query_pairs(r,pairs,p=p,status=code)
    if (code/=kdtree_success) error stop 'bad query_pairs'
    print '(a,1x,i0)', 'P',size(pairs,1)
    do i=1,size(pairs,1)
        print '(i0,1x,i0)', pairs(i,:)
    end do
    call a%query_ball_tree(b,r,rp,columns,p=p,status=code)
    if (code/=kdtree_success) error stop 'bad query_ball_tree'
    print '(a,1x,i0)', 'R',size(rp)
    do i=1,size(rp)
        print '(i0)', rp(i)
    end do
    print '(a,1x,i0)', 'J',size(columns)
    do i=1,size(columns)
        print '(i0)', columns(i)
    end do
    count=a%count_neighbors(b,r,p=p,status=code)
    if (code/=kdtree_success) error stop 'bad count_neighbors'
    print '(a,1x,i0)', 'C',count
    call a%sparse_distance_matrix(b,r,row,col,values,p=p,status=code)
    if (code/=kdtree_success) error stop 'bad sparse_distance_matrix'
    print '(a,1x,i0)', 'S',size(row)
    do i=1,size(row)
        print '(i0,1x,i0,1x,es27.18e3)', row(i),col(i),values(i)
    end do
end program kdtree_periodic_driver
