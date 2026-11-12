! SPDX-License-Identifier: MIT
! Read a SciPy-style hierarchy converted to one-based IDs, then run the
! optimal-leaf-ordering solver on the original condensed distances.
program cluster_hierarchy_ordering_driver
    use scifort_kinds, only : dp
    use scifort_cluster_hierarchy, only : optimal_leaf_ordering, leaves_list
    implicit none
    integer :: n, i, k, status
    real(dp), allocatable :: y(:), z(:,:), ordered(:,:)
    integer, allocatable :: leaves(:)

    read(*,*) n
    allocate(y(n*(n-1)/2),z(n-1,4))
    do i=1,size(y)
        read(*,*) y(i)
    end do
    do k=1,n-1
        read(*,*) z(k,:)
    end do
    ordered = optimal_leaf_ordering(z,y,status)
    print '(i0,1x,i0)', status,size(ordered,1)
    if (status /= 0) stop
    do k=1,n-1
        print '(4(es25.16,1x))', ordered(k,:)
    end do
    leaves = leaves_list(ordered,status)
    print '(*(i0,1x))', leaves
end program cluster_hierarchy_ordering_driver
