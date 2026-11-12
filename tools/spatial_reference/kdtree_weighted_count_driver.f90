! SPDX-License-Identifier: MIT
! Text I/O test adapter for scipy.spatial cKDTree.count_neighbors comparison.
program kdtree_weighted_count_driver
    use, intrinsic :: iso_fortran_env, only: int64
    use scifort_kinds, only: dp
    use scifort_spatial, only: kd_tree, kd_tree_type
    implicit none
    type(kd_tree_type) :: a,b
    real(dp), allocatable :: x(:,:),y(:,:),box(:),radii(:),wa(:),wb(:),weighted(:)
    integer(int64), allocatable :: counts(:)
    real(dp) :: p,scalar_result
    integer(int64) :: scalar_count
    integer :: n,m,d,nr,i,code,flag,visited
    logical :: cumulative
    read(*,*) n,m,d,nr,p,flag
    cumulative=flag/=0
    allocate(x(n,d),y(m,d),box(d),radii(nr),wa(n),wb(m))
    read(*,*) box
    if (nr>0) read(*,*) radii
    do i=1,n
        read(*,*) x(i,:)
    end do
    do i=1,m
        read(*,*) y(i,:)
    end do
    read(*,*) wa
    read(*,*) wb
    a=kd_tree(x,boxsize=box,status=code)
    if (code/=0) error stop 'failed to build A'
    b=kd_tree(y,boxsize=box,status=code)
    if (code/=0) error stop 'failed to build B'
    counts=a%count_neighbors(b,radii,cumulative=cumulative,p=p,status=code,visited=visited)
    if (code/=0) error stop 'integer count failed'
    print '(a,1x,i0)', 'CI',size(counts)
    do i=1,size(counts)
        print '(i0)', counts(i)
    end do
    weighted=a%weighted_count_neighbors(b,radii,weights_self=wa,weights_other=wb, &
        cumulative=cumulative,p=p,status=code,visited=visited)
    if (code/=0) error stop 'weighted count failed'
    print '(a,1x,i0)', 'WW',size(weighted)
    do i=1,size(weighted)
        print '(es26.17e3)', weighted(i)
    end do
    weighted=a%weighted_count_neighbors(b,radii,weights_self=wa, &
        cumulative=cumulative,p=p,status=code)
    if (code/=0) error stop 'one-sided count A failed'
    print '(a,1x,i0)', 'WA',size(weighted)
    do i=1,size(weighted)
        print '(es26.17e3)', weighted(i)
    end do
    weighted=a%weighted_count_neighbors(b,radii,weights_other=wb, &
        cumulative=cumulative,p=p,status=code)
    if (code/=0) error stop 'one-sided count B failed'
    print '(a,1x,i0)', 'WB',size(weighted)
    do i=1,size(weighted)
        print '(es26.17e3)', weighted(i)
    end do
    weighted=a%weighted_count_neighbors(b,radii,cumulative=cumulative,p=p,status=code)
    if (code/=0) error stop 'unit weight count failed'
    print '(a,1x,i0)', 'WU',size(weighted)
    do i=1,size(weighted)
        print '(es26.17e3)', weighted(i)
    end do
    if (nr>0) then
        scalar_count=a%count_neighbors(b,radii(1),p=p,status=code)
        if (code/=0) error stop 'scalar count failed'
        print '(a,1x,i0)', 'SC',scalar_count
        scalar_result=a%weighted_count_neighbors(b,radii(1), &
            weights_self=wa,weights_other=wb,p=p,status=code)
        if (code/=0) error stop 'scalar weighted count failed'
        print '(a,1x,es26.17e3)', 'SW',scalar_result
    end if
end program kdtree_weighted_count_driver
