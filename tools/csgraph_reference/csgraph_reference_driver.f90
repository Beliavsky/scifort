! SPDX-License-Identifier: MIT
! Read dense graph cases from stdin and emit SciFort CSR algorithm results.
program csgraph_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csgraph_from_dense, &
        connected_components, dijkstra, minimum_spanning_tree, csgraph_success
    implicit none
    type(csr_graph) :: graph, forest
    real(dp), allocatable :: a(:,:), distance(:)
    integer, allocatable :: labels(:), weak(:), strong(:), predecessors(:)
    integer :: n_cases, case_id, n, i, status, nw, ns

    read(*,*) n_cases
    do case_id = 1, n_cases
        read(*,*) n
        allocate(a(n,n))
        do i = 1, n
            read(*,*) a(i,:)
        end do
        graph = csgraph_from_dense(a, status)
        if (status /= csgraph_success) error stop 'graph constructor failed'
        call connected_components(graph,nw,weak,status,directed=.true.,connection='weak')
        if (status /= csgraph_success) error stop 'weak components failed'
        call connected_components(graph,ns,strong,status,directed=.true.,connection='strong')
        if (status /= csgraph_success) error stop 'strong components failed'
        forest = minimum_spanning_tree(graph,status)
        if (status /= csgraph_success) error stop 'minimum spanning forest failed'
        write(*,'(3(i0,1x))') n,nw,ns
        write(*,'(*(i0,1x))') weak
        write(*,'(*(i0,1x))') strong
        write(*,'(i0,1x,es26.17)') size(forest%data),sum(forest%data)
        if (any(a < 0.0_dp)) then
            write(*,'(a)') 'NEGATIVE'
        else
            write(*,'(a)') 'POSITIVE'
            do i = 1, n
                call dijkstra(graph,i,distance,predecessors,status)
                if (status /= csgraph_success) error stop 'directed dijkstra failed'
                write(*,'(*(es26.17,1x))') distance
            end do
            do i = 1, n
                call dijkstra(graph,i,distance,predecessors,status,directed=.false.)
                if (status /= csgraph_success) error stop 'undirected dijkstra failed'
                write(*,'(*(es26.17,1x))') distance
            end do
        end if
        deallocate(a)
    end do
end program csgraph_reference_driver
