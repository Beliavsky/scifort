! SPDX-License-Identifier: MIT
! Streaming numerical reference driver for compare_cluster_hierarchy_layout_scipy.py
program cluster_hierarchy_layout_reference_driver
    use scifort_kinds, only : dp
    use scifort_cluster_hierarchy, only : cluster_dendrogram, dendrogram, &
        is_isomorphic, from_mlab_linkage, to_mlab_linkage, cluster_success
    implicit none
    real(dp), allocatable :: z(:,:), matlab(:,:), again(:,:)
    integer, allocatable :: t1(:), t2(:)
    type(cluster_dendrogram) :: result
    integer :: n, mode, io, i, status, code, code2, iso, roundtrip

    do
        read(*,*,iostat=io) n, mode
        if (io < 0) exit
        if (io /= 0 .or. n < 2) error stop 'reference driver invalid header'
        allocate(z(n-1,4), t1(n), t2(n))
        do i = 1,n-1
            read(*,*,iostat=io) z(i,:)
            if (io /= 0) error stop 'reference driver invalid linkage row'
        end do
        read(*,*,iostat=io) t1
        if (io /= 0) error stop 'reference driver invalid T1'
        read(*,*,iostat=io) t2
        if (io /= 0) error stop 'reference driver invalid T2'

        select case(mode)
        case(0)
            result = dendrogram(z, status)
        case(1)
            result = dendrogram(z, status, count_sort='ascending')
        case(2)
            result = dendrogram(z, status, count_sort='descending')
        case(3)
            result = dendrogram(z, status, distance_sort='ascending')
        case(4)
            result = dendrogram(z, status, distance_sort='descending')
        case default
            error stop 'reference driver invalid mode'
        end select
        matlab = to_mlab_linkage(z,code)
        if (code /= cluster_success) error stop 'MATLAB export failure'
        again = from_mlab_linkage(matlab,code2)
        roundtrip = 0
        if (code2 == cluster_success) then
            if (all(shape(z) == shape(again))) then
                if (all(again == z)) roundtrip = 1
            end if
        end if
        iso = merge(1,0,is_isomorphic(t1,t2))
        write(*,'(3(i0,1x))') status,iso,roundtrip
        write(*,'(*(i0,1x))') result%leaves
        do i = 1,n-1
            write(*,'(8(es25.16,1x))') result%icoord(i,:),result%dcoord(i,:)
        end do
        deallocate(z,t1,t2)
    end do
end program cluster_hierarchy_layout_reference_driver
