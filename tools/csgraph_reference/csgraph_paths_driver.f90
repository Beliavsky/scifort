! SPDX-License-Identifier: MIT
! Validate SciFort Bellman-Ford and breadth-first traversals using test graphs.
program csgraph_paths_driver
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csgraph_from_dense, &
        breadth_first_order, bellman_ford, csgraph_success, csgraph_negative_cycle
    implicit none
    type(csr_graph) :: graph
    real(dp), allocatable :: a(:,:), distances(:)
    integer, allocatable :: order(:), predecessors(:)
    integer :: cases, case_id, n, i, status, k

    read(*,*) cases
    do case_id=1,cases
        read(*,*) n
        allocate(a(n,n))
        do i=1,n
            read(*,*) a(i,:)
        end do
        graph = csgraph_from_dense(a,status)
        if (status /= csgraph_success) error stop 'invalid input'
        write(*,'(i0)') n
        do i=1,n
            call bellman_ford(graph,i,distances,predecessors,status)
            select case (status)
            case (csgraph_success)
                write(*,'(*(es26.17,1x))') distances
            case (csgraph_negative_cycle)
                write(*,'(a)') 'NEGATIVE_CYCLE'
            case default
                error stop 'unexpected Bellman-Ford status'
            end select
        end do
        do k=1,2
            call breadth_first_order(graph,1,order,predecessors,status,directed=k==1)
            if (status /= csgraph_success) error stop 'unexpected BFS status'
            write(*,'(i0)') size(order)
            write(*,'(*(i0,1x))') order
            write(*,'(*(i0,1x))') predecessors
        end do
        deallocate(a)
    end do
end program csgraph_paths_driver
