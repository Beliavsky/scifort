! SPDX-License-Identifier: MIT
! Standalone cluster reference driver for SciPy numerical oracle comparisons.
program cluster_hierarchy_driver
    use scifort_kinds, only : dp
    use scifort_cluster_hierarchy, only : linkage, fcluster_distance, &
        fcluster_maxclust, cophenetic_distances, leaves_list, cluster_success
    implicit none
    character(len=32) :: method, metric
    integer :: n, p, status, i, t
    real(dp) :: threshold
    real(dp), allocatable :: x(:,:), z(:,:), d(:)
    integer, allocatable :: flat(:), maxflat(:), leaf_order(:)

    read(*,*) n,p,method,metric
    read(*,*) threshold,t
    allocate(x(n,p))
    do i=1,n
        read(*,*) x(i,:)
    end do
    z=linkage(x,trim(method),trim(metric),status)
    if (status/=cluster_success) then
        print *, 'INVALID'
        stop
    end if
    d=cophenetic_distances(z,status)
    if (status/=cluster_success) error stop 'cophenetic failure'
    flat=fcluster_distance(z,threshold,status)
    if (status/=cluster_success) error stop 'flat failure'
    maxflat=fcluster_maxclust(z,t,status)
    if (status/=cluster_success) error stop 'maxflat failure'
    leaf_order=leaves_list(z,status)
    if (status/=cluster_success) error stop 'leaf failure'
    do i=1,n-1
        print '(2(i0,1x),es25.16,1x,i0)', int(z(i,1)),int(z(i,2)), &
            z(i,3),int(z(i,4))
    end do
    print '(*(es25.16,1x))', d
    print '(*(i0,1x))', flat
    print '(*(i0,1x))', maxflat
    print '(*(i0,1x))', leaf_order
end program cluster_hierarchy_driver
