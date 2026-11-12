! SPDX-License-Identifier: MIT
! Line-oriented reference driver for independent scipy.spatial.cKDTree checks.
program kdtree_pairs_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_spatial, only : kd_tree, kd_tree_type, kdtree_success
    implicit none
    type(kd_tree_type) :: a, b
    integer :: n, m, dim, i, code
    real(dp) :: radius
    real(dp), allocatable :: xa(:,:), xb(:,:), distance(:)
    integer, allocatable :: pairs(:,:), row_ptr(:), neighbors(:), rows(:), columns(:)
    integer(int64) :: neighbor_count

    read (*,*) n, m, dim, radius
    allocate(xa(n,dim),xb(m,dim))
    do i=1,n
        read (*,*) xa(i,:)
    end do
    do i=1,m
        read (*,*) xb(i,:)
    end do
    a=kd_tree(xa,status=code)
    if (code/=kdtree_success) error stop 'first tree invalid'
    b=kd_tree(xb,status=code)
    if (code/=kdtree_success) error stop 'second tree invalid'
    call a%query_pairs(radius,pairs,status=code)
    if (code/=kdtree_success) error stop 'query_pairs invalid'
    print '(a,1x,i0)', 'P',size(pairs,1)
    do i=1,size(pairs,1)
        print '(i0,1x,i0)', pairs(i,1),pairs(i,2)
    end do
    call a%query_ball_tree(b,radius,row_ptr,neighbors,status=code)
    if (code/=kdtree_success) error stop 'query_ball_tree invalid'
    print '(a,1x,i0)', 'R',size(row_ptr)
    do i=1,size(row_ptr)
        print '(i0)',row_ptr(i)
    end do
    print '(a,1x,i0)', 'J',size(neighbors)
    do i=1,size(neighbors)
        print '(i0)',neighbors(i)
    end do
    neighbor_count=a%count_neighbors(b,radius,status=code)
    if (code/=kdtree_success) error stop 'count_neighbors invalid'
    print '(a,1x,i0)', 'C',neighbor_count
    call a%sparse_distance_matrix(b,radius,rows,columns,distance,status=code)
    if (code/=kdtree_success) error stop 'sparse_distance_matrix invalid'
    print '(a,1x,i0)', 'S',size(distance)
    do i=1,size(distance)
        print '(i0,1x,i0,1x,es25.17e3)', rows(i),columns(i),distance(i)
    end do
end program kdtree_pairs_driver
