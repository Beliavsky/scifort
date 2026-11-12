! SPDX-License-Identifier: MIT
program test_sparse_csgraph
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csgraph_from_dense, &
        csgraph_from_edges, csgraph_to_dense, valid_csgraph, connected_components, breadth_first_order, &
        dijkstra, bellman_ford, minimum_spanning_tree, &
        csgraph_success, csgraph_invalid_input, csgraph_negative_cycle
    implicit none
    type(csr_graph) :: graph, forest, other
    real(dp) :: matrix(6,6), a(3,3)
    real(dp), allocatable :: d(:), dense(:,:)
    integer, allocatable :: labels(:), predecessors(:), visit_order(:)
    integer :: status, n_components

    matrix = 0.0_dp
    matrix(1,2) = 2.0_dp
    matrix(2,3) = 1.0_dp
    matrix(1,3) = 10.0_dp
    matrix(3,2) = 4.0_dp
    matrix(4,5) = 3.0_dp
    matrix(5,4) = 2.0_dp
    graph = csgraph_from_dense(matrix, status)
    call require(status == csgraph_success, 'construct dense graph')
    call require(valid_csgraph(graph), 'valid CSR')
    dense = csgraph_to_dense(graph,status)
    call require(status == csgraph_success .and. all(dense == matrix), 'dense round trip')
    call connected_components(graph,n_components,labels,status,directed=.true.,connection='weak')
    call require(status == csgraph_success .and. n_components == 3, 'weak components')
    call require(all(labels == [1,1,1,2,2,3]), 'weak labels')
    call connected_components(graph,n_components,labels,status,directed=.true.,connection='strong')
    call require(status == csgraph_success .and. n_components == 4, 'strong components')
    call require(labels(2) == labels(3) .and. labels(4) == labels(5), 'strong groups')
    call require(labels(1) /= labels(2) .and. labels(6) /= labels(1), 'strong singleton')
    call connected_components(graph,n_components,labels,status,directed=.false.,connection='strong')
    call require(status == csgraph_success .and. n_components == 3, 'undirected components')

    call dijkstra(graph,1,d,predecessors,status)
    call require(status == csgraph_success, 'directed dijkstra')
    call require(all(abs(d(:3)-[0.0_dp,2.0_dp,3.0_dp]) < 1.0e-12_dp), 'shortest path')
    call require(all(predecessors(:3) == [0,1,2]), 'shortest path predecessor')
    call require(all(.not. ieee_is_finite(d(4:))), 'unreachable vertices')
    call dijkstra(graph,3,d,predecessors,status,directed=.false.)
    call require(status == csgraph_success, 'undirected dijkstra')
    call require(all(abs(d(:3)-[3.0_dp,1.0_dp,0.0_dp]) < 1.0e-12_dp), 'reverse paths')
    call dijkstra(graph,1,d,predecessors,status,limit=2.5_dp)
    call require(status == csgraph_success .and. .not. ieee_is_finite(d(3)), 'distance limit')
    call dijkstra(graph,1,d,predecessors,status,limit=-1.0_dp)
    call require(status == csgraph_invalid_input .and. size(d) == 0, 'reject limit')

    call breadth_first_order(graph,1,visit_order,predecessors,status)
    call require(status == csgraph_success, 'directed BFS')
    call require(all(visit_order == [1,2,3]), 'BFS order')
    call require(all(predecessors == [0,1,1,0,0,0]), 'BFS predecessor')
    call breadth_first_order(graph,3,visit_order,predecessors,status,directed=.false.)
    call require(status == csgraph_success, 'undirected BFS')
    call require(all(visit_order == [3,2,1]), 'undirected BFS order')
    call bellman_ford(graph,1,d,predecessors,status)
    call require(status == csgraph_success, 'Bellman-Ford positive graph')
    call require(all(abs(d(:3)-[0.0_dp,2.0_dp,3.0_dp]) < 1.0e-12_dp), &
        'Bellman-Ford matches Dijkstra')

    forest = minimum_spanning_tree(graph,status)
    call require(status == csgraph_success, 'minimum spanning forest')
    call require(size(forest%data) == 3, 'forest edge count')
    call require(abs(sum(forest%data)-5.0_dp) < 1.0e-12_dp, 'forest weight')
    call connected_components(forest,n_components,labels,status,directed=.false.)
    call require(status == csgraph_success .and. n_components == 3, 'forest connectivity')

    other = csgraph_from_edges(3,[1,2],[2,3],[0.0_dp,3.0_dp],status)
    call require(status == csgraph_success, 'explicit zero arc')
    call connected_components(other,n_components,labels,status)
    call require(status == csgraph_success .and. n_components == 1, 'zero arc is edge')
    call dijkstra(other,1,d,predecessors,status)
    call require(status == csgraph_success, 'zero-weight shortest path')
    call require(all(abs(d-[0.0_dp,0.0_dp,3.0_dp]) < 1.0e-12_dp), 'zero-weight distance')

    a = 0.0_dp
    a(1,2) = -3.0_dp
    a(2,3) = 2.0_dp
    graph = csgraph_from_dense(a,status)
    forest = minimum_spanning_tree(graph,status)
    call require(status == csgraph_success .and. size(forest%data) == 2, 'negative MST')
    call require(abs(sum(forest%data)+1.0_dp) < 1.0e-12_dp, 'negative MST weight')
    call dijkstra(graph,1,d,predecessors,status)
    call require(status == csgraph_invalid_input, 'negative dijkstra rejected')
    call bellman_ford(graph,1,d,predecessors,status)
    call require(status == csgraph_success, 'negative-edge Bellman-Ford')
    call require(all(abs(d-[0.0_dp,-3.0_dp,-1.0_dp]) < 1.0e-12_dp), 'negative shortest paths')
    a(3,1) = -2.0_dp
    graph = csgraph_from_dense(a,status)
    call bellman_ford(graph,1,d,predecessors,status)
    call require(status == csgraph_negative_cycle .and. size(d) == 0, 'negative cycle detected')
    a(3,1) = 0.0_dp

    graph = csgraph_from_edges(0,[integer::],[integer::],[real(dp)::],status)
    call require(status == csgraph_success .and. valid_csgraph(graph), 'zero-node graph')
    call connected_components(graph,n_components,labels,status,connection='strong')
    call require(status == csgraph_success .and. n_components == 0, 'zero-node components')
    forest = minimum_spanning_tree(graph,status)
    call require(status == csgraph_success .and. forest%n == 0, 'empty forest')
    graph = csgraph_from_edges(3,[1],[4],[1.0_dp],status)
    call require(status == csgraph_invalid_input .and. .not. valid_csgraph(graph), 'bad endpoint')
    graph = csgraph_from_dense(a,status)
    graph%indptr(1) = 2
    call require(.not. valid_csgraph(graph), 'bad offsets rejected')
    call connected_components(graph,n_components,labels,status)
    call require(status == csgraph_invalid_input .and. size(labels) == 0, 'bad graph rejected')

    print *, 'sparse csgraph regression PASS'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok !! Condition to verify.
        character(len=*), intent(in) :: message !! Diagnostic on failure.
        if (.not. ok) then
            print *, 'FAILED: ',message
            error stop 1
        end if
    end subroutine require
end program test_sparse_csgraph
