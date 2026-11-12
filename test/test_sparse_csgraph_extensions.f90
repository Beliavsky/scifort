! SPDX-License-Identifier: MIT
program test_sparse_csgraph_extensions
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csgraph_from_dense, csgraph_from_edges, csgraph_to_dense, &
        floyd_warshall, johnson, shortest_path, laplacian, depth_first_order, &
        breadth_first_tree, depth_first_tree, reconstruct_path, path_from_predecessors, &
        csgraph_success, csgraph_invalid_input, csgraph_negative_cycle
    implicit none
    type(csr_graph) :: graph, tree
    real(dp) :: a(5,5), b(3,3), expected(3,3), tol
    real(dp), allocatable :: distances(:,:), lap(:,:), diagonal(:), dense(:,:)
    integer, allocatable :: predecessors(:,:), order(:), parents(:), path(:)
    integer :: status

    tol = 5.0e-13_dp
    a = 0.0_dp
    a(1,2) = 2.0_dp
    a(1,3) = 10.0_dp
    a(2,3) = -3.0_dp
    a(3,4) = 4.0_dp
    a(4,1) = 1.0_dp
    graph = csgraph_from_dense(a,status)
    call require(status == csgraph_success,'construct directed weighted graph')
    call floyd_warshall(graph,distances,predecessors,status)
    call require(status == csgraph_success,'Floyd-Warshall negative edges')
    call require(all(abs(distances(1,1:4)-[0.0_dp,2.0_dp,-1.0_dp,3.0_dp]) < tol), &
        'Floyd-Warshall distances')
    call require(.not. ieee_is_finite(distances(1,5)), 'Floyd-Warshall isolated node')
    call require(all(predecessors(1,1:4)==[0,1,2,3]),'Floyd-Warshall predecessors')
    call path_from_predecessors(1,4,predecessors(1,:),path,status)
    call require(status == csgraph_success .and. all(path==[1,2,3,4]),'Floyd-Warshall path reconstruction')
    call path_from_predecessors(1,5,predecessors(1,:),path,status)
    call require(status == csgraph_success .and. size(path)==0,'unreachable path')
    call path_from_predecessors(1,1,predecessors(1,:),path,status)
    call require(status == csgraph_success .and. all(path==[1]),'zero-length path')
    call johnson(graph,distances,predecessors,status)
    call require(status == csgraph_success,'Johnson negative-edge all pairs')
    call require(all(abs(distances(1,1:4)-[0.0_dp,2.0_dp,-1.0_dp,3.0_dp]) < tol),'Johnson distances')
    call require(.not. ieee_is_finite(distances(5,1)),'Johnson isolated row')
    call require(all(predecessors(1,1:4)==[0,1,2,3]),'Johnson predecessors')
    call johnson(graph,distances,predecessors,status,indices=[4,1,4])
    call require(status == csgraph_success,'Johnson selected source rows')
    call require(all(shape(distances)==[3,5]),'Johnson selected shape')
    call require(all(abs(distances(2,1:4)-[0.0_dp,2.0_dp,-1.0_dp,3.0_dp])<tol),'selected source values')
    call require(all(distances(1,:)==distances(3,:)),'repeated source rows')
    call johnson(graph,distances,predecessors,status,indices=[0])
    call require(status==csgraph_invalid_input .and. size(distances)==0,'invalid selected source')
    call johnson(graph,distances,predecessors,status,indices=[integer::])
    call require(status==csgraph_success .and. all(shape(distances)==[0,5]),'empty selected sources')

    call floyd_warshall(graph,distances,predecessors,status,unweighted=.true.)
    call require(status==csgraph_success .and. distances(1,3)==1.0_dp,'Floyd-Warshall unweighted')
    call johnson(graph,distances,predecessors,status,unweighted=.true.)
    call require(status==csgraph_success .and. distances(1,4)==2.0_dp,'Johnson unweighted')
    call floyd_warshall(graph,distances,predecessors,status,directed=.false.,unweighted=.true.)
    call require(status==csgraph_success .and. distances(3,1)==1.0_dp,'Floyd-Warshall undirected')
    call johnson(graph,distances,predecessors,status,directed=.false.,unweighted=.true.)
    call require(status==csgraph_success .and. distances(3,1)==1.0_dp,'Johnson undirected')
    ! An undirected negative edge gives a negative round trip.
    call johnson(graph,distances,predecessors,status,directed=.false.)
    call require(status==csgraph_negative_cycle .and. size(distances)==0,'Johnson undirected negative cycle')
    call floyd_warshall(graph,distances,predecessors,status,directed=.false.)
    call require(status==csgraph_negative_cycle .and. size(distances)==0,'Floyd-Warshall undirected negative cycle')

    ! Disconnected negative cycle must be detected by all-pairs algorithms.
    graph = csgraph_from_edges(5,[1,4,5],[2,5,4],[1.0_dp,-2.0_dp,1.0_dp],status)
    call johnson(graph,distances,predecessors,status)
    call require(status==csgraph_negative_cycle,'Johnson disconnected negative cycle')
    call floyd_warshall(graph,distances,predecessors,status)
    call require(status==csgraph_negative_cycle,'Floyd-Warshall disconnected negative cycle')
    call johnson(graph,distances,predecessors,status,unweighted=.true.)
    call require(status==csgraph_success,'Johnson unweighted ignores signed edges')
    call floyd_warshall(graph,distances,predecessors,status,unweighted=.true.)
    call require(status==csgraph_success,'Floyd unweighted ignores signed edges')

    ! Parallel arcs: paths choose the minimum edge, while Laplacians sum weights.
    graph = csgraph_from_edges(3,[1,1,2,3,1],[2,2,3,1,1], &
        [3.0_dp,0.0_dp,4.0_dp,2.0_dp,99.0_dp],status)
    call floyd_warshall(graph,distances,predecessors,status)
    call require(status==csgraph_success .and. distances(1,3)==4.0_dp,'Floyd stored zero parallel arc')
    call johnson(graph,distances,predecessors,status)
    call require(status==csgraph_success .and. distances(1,3)==4.0_dp,'Johnson stored zero parallel arc')
    call laplacian(graph,lap,diagonal,status,use_out_degree=.true.)
    call require(status==csgraph_success,'out-degree graph Laplacian')
    expected = reshape([3.0_dp,0.0_dp,-2.0_dp,-3.0_dp,4.0_dp,0.0_dp, &
        0.0_dp,-4.0_dp,2.0_dp],[3,3])
    call require(all(abs(lap-expected)<tol),'Laplacian sums parallel arcs, ignores loops')
    call require(all(abs(diagonal-[3.0_dp,4.0_dp,2.0_dp])<tol),'out-degree vector')
    call laplacian(graph,lap,diagonal,status)
    call require(status==csgraph_success,'in-degree Laplacian')
    call require(all(abs(diagonal-[2.0_dp,3.0_dp,4.0_dp])<tol),'in-degree vector')
    call laplacian(graph,lap,diagonal,status,symmetrized=.true.)
    call require(status==csgraph_success,'symmetrized Laplacian')
    call require(all(abs(lap-transpose(lap))<tol),'symmetric Laplacian')
    call laplacian(graph,lap,diagonal,status,symmetrized=.true.,normed=.true.)
    call require(status==csgraph_success,'symmetrically normalized Laplacian')
    call require(all(abs(diag_of(lap)-1.0_dp)<tol),'normalized diagonal ones')
    call require(all(abs(lap-transpose(lap))<tol),'normalized symmetric Laplacian')

    b = 0.0_dp
    b(1,1) = 2.0_dp
    b(2,3) = 4.0_dp
    graph = csgraph_from_dense(b,status)
    call laplacian(graph,lap,diagonal,status,normed=.true.,use_out_degree=.true.)
    call require(status==csgraph_success,'normalized isolated vertices')
    call require(all(abs(diagonal-[1.0_dp,2.0_dp,1.0_dp])<tol),'isolate scales')
    call require(all(abs(diag_of(lap)-[0.0_dp,1.0_dp,0.0_dp])<tol),'normalized isolate diag')
    b(2,3)=-4.0_dp
    graph = csgraph_from_dense(b,status)
    call laplacian(graph,lap,diagonal,status,normed=.true.,use_out_degree=.true.)
    call require(status==csgraph_invalid_input .and. size(lap)==0, 'reject negative normalized degree')
    call laplacian(graph,lap,diagonal,status,use_out_degree=.true.)
    call require(status==csgraph_success .and. lap(2,2)==-4.0_dp,'unnormalized signed Laplacian')

    graph = csgraph_from_edges(5,[1,1,2,2,3,4],[2,3,4,5,5,5], &
        [4.0_dp,5.0_dp,6.0_dp,7.0_dp,8.0_dp,9.0_dp],status)
    call depth_first_order(graph,1,order,parents,status)
    call require(status==csgraph_success .and. all(order==[1,2,4,5,3]),'DFS discovery ordering')
    call require(all(parents==[0,1,1,2,4]),'DFS parent pointers')
    call path_from_predecessors(1,5,parents,path,status)
    call require(status==csgraph_success .and. all(path==[1,2,4,5]),'DFS path')
    call depth_first_order(graph,5,order,parents,status,directed=.false.)
    call require(status==csgraph_success .and. size(order)==5,'undirected DFS reachability')
    call depth_first_order(graph,6,order,parents,status)
    call require(status==csgraph_invalid_input .and. size(order)==0,'invalid DFS root')
    call path_from_predecessors(1,3,[0,3,2],path,status)
    call require(status==csgraph_invalid_input .and. size(path)==0,'cyclic invalid predecessor chain')


    ! The unified shortest-path dispatcher covers four named algorithms,
    ! automatic selection, selected source rows, and unit-edge traversal.
    graph = csgraph_from_dense(a,status)
    call shortest_path(graph,distances,predecessors,status)
    call require(status==csgraph_success .and. distances(1,4)==3.0_dp,'auto shortest paths')
    call shortest_path(graph,distances,predecessors,status,method='FW',indices=[4,1])
    call require(status==csgraph_success .and. all(shape(distances)==[2,5]),'FW selected rows')
    call require(distances(2,4)==3.0_dp,'FW selected shortest path')
    call shortest_path(graph,distances,predecessors,status,method='J',indices=[1])
    call require(status==csgraph_success .and. distances(1,3)==-1.0_dp,'J selected source')
    call shortest_path(graph,distances,predecessors,status,method='BF',indices=[1])
    call require(status==csgraph_success .and. distances(1,3)==-1.0_dp,'BF selected source')
    call shortest_path(graph,distances,predecessors,status,method='D')
    call require(status==csgraph_invalid_input .and. size(distances)==0,'D rejects signed arcs')
    call shortest_path(graph,distances,predecessors,status,method='D',unweighted=.true.)
    call require(status==csgraph_success .and. distances(1,4)==2.0_dp,'D unit-weight mode')
    call shortest_path(graph,distances,predecessors,status,method='nonsense')
    call require(status==csgraph_invalid_input .and. size(distances)==0,'invalid shortest path method')
    call shortest_path(graph,distances,predecessors,status,indices=[9])
    call require(status==csgraph_invalid_input,'invalid shortest path indices')

    graph = csgraph_from_edges(1,[1],[1],[-1.0_dp],status)
    call floyd_warshall(graph,distances,predecessors,status)
    call require(status==csgraph_negative_cycle,'strict negative self-loop detection')
    call johnson(graph,distances,predecessors,status)
    call require(status==csgraph_negative_cycle,'Johnson negative self-loop detection')

    graph = csgraph_from_edges(5,[1,1,2,2,3,4],[2,3,4,5,5,5], &
        [4.0_dp,5.0_dp,6.0_dp,7.0_dp,8.0_dp,9.0_dp],status)
    tree = breadth_first_tree(graph,1,status)
    call require(status==csgraph_success .and. size(tree%data)==4,'BFS spanning tree count')
    dense = csgraph_to_dense(tree,status)
    call require(dense(1,2)==4.0_dp .and. dense(2,5)==7.0_dp,'BFS tree weights')
    call require(dense(4,5)==0.0_dp .and. dense(1,3)==5.0_dp,'BFS tree level order')
    tree = depth_first_tree(graph,1,status)
    call require(status==csgraph_success .and. size(tree%data)==4,'DFS spanning tree count')
    dense = csgraph_to_dense(tree,status)
    call require(dense(4,5)==9.0_dp .and. dense(1,3)==5.0_dp,'DFS tree weights')
    call require(dense(2,5)==0.0_dp,'DFS tree depth order')
    tree = breadth_first_tree(graph,3,status,directed=.false.)
    call require(status==csgraph_success .and. tree%n==5,'undirected traversal tree')
    dense = csgraph_to_dense(tree,status)
    call require(dense(3,1)==5.0_dp,'undirected tree reversed edge retains weight')
    tree = reconstruct_path(graph,[0,1,1,2,4],status)
    call require(status==csgraph_success .and. size(tree%data)==4,'reconstruct predecessor tree')
    dense = csgraph_to_dense(tree,status)
    call require(dense(1,2)==4.0_dp .and. dense(4,5)==9.0_dp,'reconstructed original weights')
    tree = reconstruct_path(graph,[0,1,3,2,4],status)
    call require(status==csgraph_invalid_input,'reject cyclic predecessor tree')
    tree = reconstruct_path(graph,[0,1,1,1,4],status)
    call require(status==csgraph_invalid_input,'reject predecessor absent from graph')
    tree = reconstruct_path(graph,[3,1,0,2,4],status,directed=.false.)
    call require(status==csgraph_success,'undirected reverse predecessor edge')
    dense = csgraph_to_dense(tree,status)
    call require(dense(3,1)==5.0_dp .and. dense(1,2)==4.0_dp,'undirected tree original weights')
    tree = depth_first_tree(graph,6,status)
    call require(status==csgraph_invalid_input,'invalid tree root')

    graph = csgraph_from_edges(0,[integer::],[integer::],[real(dp)::],status)
    call johnson(graph,distances,predecessors,status)
    call require(status==csgraph_success .and. all(shape(distances)==[0,0]),'empty Johnson')
    call floyd_warshall(graph,distances,predecessors,status)
    call require(status==csgraph_success .and. all(shape(distances)==[0,0]),'empty Floyd')
    call laplacian(graph,lap,diagonal,status)
    call require(status==csgraph_success .and. all(shape(lap)==[0,0]),'empty Laplacian')
    print *, 'sparse csgraph extensions regression PASS'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! Condition to verify.
        character(len=*), intent(in) :: message !! Diagnostic on failure.
        if (.not. ok) then
            print *, 'FAILED: ',message
            error stop 1
        end if
    end subroutine require

    function diag_of(a) result(diagonal)
        real(dp), intent(in) :: a(:,:) !! Square matrix.
        real(dp) :: diagonal(size(a,1))
        integer :: j
        do j=1,size(a,1)
            diagonal(j)=a(j,j)
        end do
    end function diag_of
end program test_sparse_csgraph_extensions
