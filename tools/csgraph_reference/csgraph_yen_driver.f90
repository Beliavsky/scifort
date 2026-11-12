! SPDX-License-Identifier: MIT
! Protocol: trial count; per trial, n m source sink K directed unweighted;
! followed by m lines of one-based source target finite weight.
program csgraph_yen_driver
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csgraph_from_edges, yen, csgraph_success
    implicit none
    type(csr_graph) :: graph
    integer, allocatable :: source_edges(:), sink_edges(:), predecessors(:,:)
    real(dp), allocatable :: weights(:), costs(:)
    integer :: trials, t, n, m, source, sink, k, directed, unweighted, code, status, i, j

    read(*,*) trials
    do t = 1, trials
        read(*,*) n,m,source,sink,k,directed,unweighted
        allocate(source_edges(m),sink_edges(m),weights(m))
        do i = 1, m
            read(*,*) source_edges(i),sink_edges(i),weights(i)
        end do
        graph = csgraph_from_edges(n,source_edges,sink_edges,weights,code)
        if (code /= csgraph_success) error stop 'invalid graph in reference driver'
        call yen(graph,source,sink,k,costs,predecessors,status, &
            directed=directed/=0,unweighted=unweighted/=0)
        write(*,'(i0,1x,i0)') status,size(costs)
        do i = 1, size(costs)
            write(*,'(es25.17,1x,*(i0,1x))') costs(i),(predecessors(i,j),j=1,n)
        end do
        deallocate(source_edges,sink_edges,weights)
    end do
end program csgraph_yen_driver
