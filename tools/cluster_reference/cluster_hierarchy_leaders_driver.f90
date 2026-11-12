! SPDX-License-Identifier: MIT
! One-based linkage and leaf labels in, one-based cluster leader IDs out.
program cluster_hierarchy_leaders_driver
    use scifort_kinds, only : dp
    use scifort_cluster_hierarchy, only : leaders
    implicit none
    integer :: n, i, k, status
    real(dp), allocatable :: z(:,:)
    integer, allocatable :: labels(:), leader_pairs(:,:)

    read(*,*) n
    allocate(z(n-1,4),labels(n))
    do k=1,n-1
        read(*,*) z(k,:)
    end do
    read(*,*) labels
    leader_pairs = leaders(z,labels,status)
    print '(i0,1x,i0)', status,size(leader_pairs,1)
    if (status /= 0) stop
    do i=1,size(leader_pairs,1)
        print '(i0,1x,i0)', leader_pairs(i,:)
    end do
end program cluster_hierarchy_leaders_driver
