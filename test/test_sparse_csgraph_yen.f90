! SPDX-License-Identifier: MIT
! Native examples and edge cases for K shortest simple paths (SciPy csgraph.yen).
program test_sparse_csgraph_yen
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csgraph_from_edges, yen, &
        path_from_predecessors, csgraph_success, csgraph_invalid_input, csgraph_negative_cycle
    implicit none
    type(csr_graph) :: graph
    real(dp), allocatable :: distances(:)
    integer, allocatable :: predecessors(:,:), path(:)
    integer :: status, i, j
    real(dp), parameter :: tol = 1.0e-12_dp

    ! Two equal-cost alternatives coexist in the Yen candidate heap.
    graph = csgraph_from_edges(4,[1,1,2,2,3],[2,3,4,3,4], &
        [1.0_dp,2.0_dp,1.0_dp,1.0_dp,3.0_dp])
    call yen(graph,1,4,8,distances,predecessors,status)
    call require(status == csgraph_success,'directed success')
    call require(size(distances)==3,'exact number of simple routes')
    call require(all(shape(predecessors)==[3,4]),'predecessor shape')
    call require(all(abs(distances-[2.0_dp,5.0_dp,5.0_dp])<tol),'directed ordered costs')
    call path_from_predecessors(1,4,predecessors(1,:),path,status)
    call require(status==csgraph_success .and. all(path==[1,2,4]),'first path vertices')
    call path_from_predecessors(1,4,predecessors(2,:),path,status)
    call require(status==csgraph_success .and. all(path==[1,2,3,4]),'second path vertices')
    call path_from_predecessors(1,4,predecessors(3,:),path,status)
    call require(status==csgraph_success .and. all(path==[1,3,4]),'third path vertices')
    call yen(graph,1,4,2,distances,predecessors,status)
    call require(status==csgraph_success .and. size(distances)==2,'K smaller than path count')

    ! Johnson reweighting allows negative edges in the absence of cycles.
    graph = csgraph_from_edges(4,[1,1,1,2,2,3],[2,3,4,3,4,4], &
        [-3.0_dp,2.0_dp,6.0_dp,1.0_dp,4.0_dp,0.0_dp])
    call yen(graph,1,4,7,distances,predecessors,status)
    call require(status==csgraph_success .and. size(distances)==4,'negative DAG paths')
    call require(all(abs(distances-[-2.0_dp,1.0_dp,2.0_dp,6.0_dp])<tol), &
        'negative DAG distances')
    call yen(graph,1,4,7,distances,predecessors,status,unweighted=.true.)
    call require(status==csgraph_success .and. size(distances)==4,'unweighted signed graph')
    call require(all(abs(distances-[1.0_dp,2.0_dp,2.0_dp,3.0_dp])<tol), &
        'unit cost regardless of signed stored weights')

    graph = csgraph_from_edges(4,[1,2,3],[2,3,4], &
        [1.0_dp,2.0_dp,3.0_dp])
    call yen(graph,1,4,4,distances,predecessors,status,directed=.false.)
    call require(status==csgraph_success .and. size(distances)==1,'undirected chain')
    call require(abs(distances(1)-6.0_dp)<tol,'undirected cost')
    call yen(graph,4,1,4,distances,predecessors,status,directed=.false.)
    call require(status==csgraph_success .and. size(distances)==1,'reverse undirected path')
    call require(abs(distances(1)-6.0_dp)<tol,'reverse undirected cost')
    call yen(graph,4,1,4,distances,predecessors,status)
    call require(status==csgraph_success .and. size(distances)==0,'directed unreachable')
    call require(all(shape(predecessors)==[0,4]),'unreachable predecessor shape')

    ! Explicit zeros and parallel edges represent one vertex sequence.
    graph = csgraph_from_edges(3,[1,1,1,2,2,3],[2,2,1,3,3,3], &
        [2.0_dp,0.0_dp,17.0_dp,0.0_dp,2.0_dp,8.0_dp])
    call yen(graph,1,3,7,distances,predecessors,status)
    call require(status==csgraph_success .and. size(distances)==1,'deduplicated parallel paths')
    call require(distances(1)==0.0_dp,'zero-weight arcs remain traversable')
    call yen(graph,1,1,5,distances,predecessors,status)
    call require(status==csgraph_success .and. size(distances)==1,'same source and sink')
    call require(distances(1)==0.0_dp .and. all(predecessors==0), &
        'zero-length source/sink path')
    call yen(graph,1,3,0,distances,predecessors,status)
    call require(status==csgraph_success .and. size(distances)==0,'K zero')
    call require(all(shape(predecessors)==[0,3]),'K zero predecessor shape')

    ! Negative edges in the undirected interpretation create negative cycles.
    graph = csgraph_from_edges(2,[1],[2],[-1.0_dp])
    call yen(graph,1,2,3,distances,predecessors,status,directed=.false.)
    call require(status==csgraph_negative_cycle,'undirected negative edge cycle')
    call yen(graph,1,2,3,distances,predecessors,status,directed=.false.,unweighted=.true.)
    call require(status==csgraph_success .and. size(distances)==1,'unweighted ignores negative cycle')

    ! Check cycles even when they are disconnected from both endpoints.
    graph = csgraph_from_edges(5,[1,4,5],[2,5,4], &
        [3.0_dp,-2.0_dp,1.0_dp])
    call yen(graph,1,2,3,distances,predecessors,status)
    call require(status==csgraph_negative_cycle,'disconnected negative cycle')
    call yen(graph,1,2,0,distances,predecessors,status)
    call require(status==csgraph_success,'K zero bypasses cycle validation')

    graph = csgraph_from_edges(4,[1],[2],[1.0_dp])
    call yen(graph,0,2,2,distances,predecessors,status)
    call require(status==csgraph_invalid_input,'invalid source')
    call yen(graph,1,5,2,distances,predecessors,status)
    call require(status==csgraph_invalid_input,'invalid sink')
    call yen(graph,1,2,-1,distances,predecessors,status)
    call require(status==csgraph_invalid_input,'negative K')
    call yen(graph,3,4,huge(i),distances,predecessors,status)
    call require(status==csgraph_success .and. size(distances)==0, &
        'large K with unreachable target does not allocate K paths')

    ! An equal-cost layered DAG has eight distinct simple routes; selecting
    ! five exercises candidate-heap tie-breaking and dynamic path storage.
    graph = csgraph_from_edges(8, &
        [1,1,2,2,3,3,4,4,5,5,6,7], &
        [2,3,4,5,4,5,6,7,6,7,8,8], &
        [1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
         1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp])
    call yen(graph,1,8,5,distances,predecessors,status)
    call require(status==csgraph_success .and. size(distances)==5,'layered graph K=5')
    call require(all(distances==4.0_dp),'layered paths all same cost')
    do i = 1, size(distances)
        do j = i+1, size(distances)
            call require(any(predecessors(i,:)/=predecessors(j,:)), &
                'distinct equally weighted vertex routes')
        end do
    end do

    graph = csr_graph()
    call yen(graph,1,2,2,distances,predecessors,status)
    call require(status==csgraph_invalid_input,'unallocated graph')
    graph = csgraph_from_edges(0,[integer::],[integer::],[real(dp)::])
    call yen(graph,1,1,2,distances,predecessors,status)
    call require(status==csgraph_invalid_input,'empty graph invalid source')

    print '(a)', 'test_sparse_csgraph_yen: PASS'
contains
    subroutine require(condition,description)
        logical, intent(in) :: condition !! Expected invariant.
        character(len=*), intent(in) :: description !! Failure diagnostic.
        if (.not. condition) then
            print '(a)', 'FAILED: '//description
            error stop 1
        end if
    end subroutine require
end program test_sparse_csgraph_yen
