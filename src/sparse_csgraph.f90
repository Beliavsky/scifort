! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent compressed-sparse-row graph algorithms corresponding to a
! subset of scipy.sparse.csgraph. All vertex and CSR offsets are one-based.
module scifort_sparse_csgraph
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: csgraph_success = 0
    integer, parameter, public :: csgraph_invalid_input = 1
    integer, parameter, public :: csgraph_negative_cycle = 2
    integer, parameter, public :: csgraph_no_full_matching = 3
    public :: csr_graph, csgraph_from_dense, csgraph_from_edges
    public :: csgraph_to_dense, valid_csgraph
    public :: connected_components, breadth_first_order
    public :: dijkstra, bellman_ford, minimum_spanning_tree
    public :: floyd_warshall, johnson, laplacian, depth_first_order, reconstruct_path, path_from_predecessors
    public :: shortest_path, breadth_first_tree, depth_first_tree
    public :: reverse_cuthill_mckee, maximum_flow
    public :: csr_bipartite_graph, bipartite_from_edges, bipartite_from_dense, &
        bipartite_from_csgraph, valid_bipartite_graph
    public :: maximum_bipartite_matching, structural_rank
    public :: csr_weighted_bipartite_graph, weighted_bipartite_from_edges, &
        weighted_bipartite_from_dense, weighted_bipartite_from_csgraph, &
        valid_weighted_bipartite_graph, minimum_weight_full_bipartite_matching
    public :: min_weight_full_bipartite_matching
    public :: laplacian_sparse, yen

    ! Explicit stored zero-weight arcs remain graph edges for traversal and
    ! shortest paths; dense input uses zero to denote absence of an edge.
    type :: csr_graph
        integer :: n = 0
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
        real(dp), allocatable :: data(:)
    end type csr_graph

    ! Structural rectangular CSR for bipartite matching. There are no values:
    ! explicitly stored zero-weight connections are preserved by the edge constructor.
    type :: csr_bipartite_graph
        integer :: n_rows = 0
        integer :: n_cols = 0
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
    end type csr_bipartite_graph

    ! Rectangular weighted biadjacency; the structural type above remains
    ! unchanged for callers that need unweighted maximum-cardinality matches.
    type :: csr_weighted_bipartite_graph
        integer :: n_rows = 0
        integer :: n_cols = 0
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
        real(dp), allocatable :: data(:)
    end type csr_weighted_bipartite_graph

    ! One simple path for Yen's algorithm.  Keep selected arc costs so that
    ! root prefixes can be reused without ambiguities from parallel arcs.
    type :: yen_path
        integer, allocatable :: vertices(:)
        real(dp), allocatable :: edge_weights(:)
        real(dp) :: cost = 0.0_dp
    end type yen_path

    type :: graph_edge
        integer :: source = 0
        integer :: target = 0
        real(dp) :: weight = 0.0_dp
    end type graph_edge

contains

    function csgraph_from_dense(a, status) result(graph)
        real(dp), intent(in) :: a(:,:) !! Square finite adjacency matrix; zero means no arc.
        integer, intent(out), optional :: status !! Zero on success, invalid input otherwise.
        type(csr_graph) :: graph
        integer :: n, m, i, j, k

        if (present(status)) status = csgraph_invalid_input
        n = size(a,1)
        if (size(a,2) /= n .or. .not. all(ieee_is_finite(a))) return
        if (count(a /= 0.0_dp, kind=int64) > int(huge(m)-1,int64)) return
        m = count(a /= 0.0_dp)
        graph%n = n
        allocate(graph%indptr(n+1), graph%indices(m), graph%data(m))
        k = 0
        do i = 1, n
            graph%indptr(i) = k+1
            do j = 1, n
                if (a(i,j) == 0.0_dp) cycle
                k = k+1
                graph%indices(k) = j
                graph%data(k) = a(i,j)
            end do
        end do
        graph%indptr(n+1) = m+1
        if (present(status)) status = csgraph_success
    end function csgraph_from_dense

    function csgraph_from_edges(n, source, target, weight, status) result(graph)
        integer, intent(in) :: n !! Number of vertices, may be zero.
        integer, intent(in) :: source(:) !! One-based source vertices.
        integer, intent(in) :: target(:) !! One-based destination vertices.
        real(dp), intent(in) :: weight(:) !! Finite edge weights; zeros allowed.
        integer, intent(out), optional :: status !! Zero on success, invalid input otherwise.
        type(csr_graph) :: graph
        integer, allocatable :: cursor(:)
        integer :: m, k, i, p

        if (present(status)) status = csgraph_invalid_input
        if (n < 0) return
        m = size(source)
        if (size(target) /= m .or. size(weight) /= m) return
        if (m == huge(m)) return
        if (any(source < 1) .or. any(source > n) .or. &
            any(target < 1) .or. any(target > n)) return
        if (.not. all(ieee_is_finite(weight))) return
        graph%n = n
        allocate(graph%indptr(n+1), graph%indices(m), graph%data(m), cursor(n))
        graph%indptr = 0
        graph%indptr(1) = 1
        do k = 1, m
            graph%indptr(source(k)+1) = graph%indptr(source(k)+1)+1
        end do
        do i = 1, n
            graph%indptr(i+1) = graph%indptr(i+1)+graph%indptr(i)
        end do
        cursor = graph%indptr(:n)
        do k = 1, m
            p = cursor(source(k))
            graph%indices(p) = target(k)
            graph%data(p) = weight(k)
            cursor(source(k)) = p+1
        end do
        if (present(status)) status = csgraph_success
    end function csgraph_from_edges

    pure logical function valid_csgraph(graph) result(ok)
        type(csr_graph), intent(in) :: graph !! CSR with one-based indices and offsets.
        integer :: n, m

        ok = .false.
        n = graph%n
        if (n < 0) return
        if (.not. allocated(graph%indptr) .or. .not. allocated(graph%indices) .or. &
            .not. allocated(graph%data)) return
        m = size(graph%indices)
        if (m == huge(m)) return
        if (size(graph%indptr) /= n+1 .or. size(graph%data) /= m) return
        if (graph%indptr(1) /= 1 .or. graph%indptr(n+1) /= m+1) return
        if (any(graph%indptr < 1) .or. any(graph%indptr > m+1)) return
        if (any(graph%indptr(2:) < graph%indptr(:n))) return
        if (any(graph%indices < 1) .or. any(graph%indices > n)) return
        if (.not. all(ieee_is_finite(graph%data))) return
        ok = .true.
    end function valid_csgraph

    function csgraph_to_dense(graph, status) result(a)
        type(csr_graph), intent(in) :: graph !! Valid compressed-sparse-row graph.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: a(:,:)
        integer :: i, k

        if (present(status)) status = csgraph_invalid_input
        allocate(a(0,0))
        if (.not. valid_csgraph(graph)) return
        deallocate(a)
        allocate(a(graph%n,graph%n))
        a = 0.0_dp
        do i = 1, graph%n
            do k = graph%indptr(i), graph%indptr(i+1)-1
                a(i,graph%indices(k)) = graph%data(k)
            end do
        end do
        if (present(status)) status = csgraph_success
    end function csgraph_to_dense

    ! Construct the arc-reversed CSR in O(n + m) time and storage.
    subroutine reverse_graph(graph, reverse)
        type(csr_graph), intent(in) :: graph !! Valid directed CSR graph.
        type(csr_graph), intent(out) :: reverse !! Edge-reversed graph.
        integer, allocatable :: cursor(:)
        integer :: n, m, i, k, j, p

        n = graph%n
        m = size(graph%indices)
        reverse%n = n
        allocate(reverse%indptr(n+1), reverse%indices(m), reverse%data(m), cursor(n))
        reverse%indptr = 0
        reverse%indptr(1) = 1
        do k = 1, m
            j = graph%indices(k)
            reverse%indptr(j+1) = reverse%indptr(j+1)+1
        end do
        do i = 1, n
            reverse%indptr(i+1) = reverse%indptr(i+1)+reverse%indptr(i)
        end do
        cursor = reverse%indptr(:n)
        do i = 1, n
            do k = graph%indptr(i), graph%indptr(i+1)-1
                j = graph%indices(k)
                p = cursor(j)
                reverse%indices(p) = i
                reverse%data(p) = graph%data(k)
                cursor(j) = p+1
            end do
        end do
    end subroutine reverse_graph

    subroutine connected_components(graph, n_components, labels, status, directed, connection)
        type(csr_graph), intent(in) :: graph !! Valid directed or undirected CSR graph.
        integer, intent(out) :: n_components !! Number of weak/strong components.
        integer, allocatable, intent(out) :: labels(:) !! One-based component labels.
        integer, intent(out) :: status !! Zero on success.
        logical, intent(in), optional :: directed !! Default true.
        character(len=*), intent(in), optional :: connection !! 'weak' (default) or 'strong'.
        type(csr_graph) :: reverse
        integer, allocatable :: queue(:), path(:), cursor(:), finish(:)
        logical, allocatable :: visited(:)
        integer :: n, i, top, k, v, u, head, tail, n_finished, c, which
        logical :: strong
        character(len=:), allocatable :: rule

        n_components = 0
        status = csgraph_invalid_input
        allocate(labels(0))
        if (.not. valid_csgraph(graph)) return
        rule = 'weak'
        if (present(connection)) rule = trim(connection)
        if (rule /= 'weak' .and. rule /= 'strong') return
        strong = rule == 'strong'
        if (present(directed)) then
            if (.not. directed) strong = .false.
        end if
        n = graph%n
        call reverse_graph(graph, reverse)
        deallocate(labels)
        allocate(labels(n), queue(n))
        labels = 0
        if (.not. strong) then
            ! A weak component is a BFS on outgoing and incoming arcs.
            do i = 1, n
                if (labels(i) /= 0) cycle
                n_components = n_components+1
                head = 1
                tail = 1
                queue(1) = i
                labels(i) = n_components
                do while (head <= tail)
                    u = queue(head)
                    head = head+1
                    do which = 1, 2
                        if (which == 1) then
                            do k = graph%indptr(u), graph%indptr(u+1)-1
                                v = graph%indices(k)
                                if (labels(v) /= 0) cycle
                                tail = tail+1
                                queue(tail) = v
                                labels(v) = n_components
                            end do
                        else
                            do k = reverse%indptr(u), reverse%indptr(u+1)-1
                                v = reverse%indices(k)
                                if (labels(v) /= 0) cycle
                                tail = tail+1
                                queue(tail) = v
                                labels(v) = n_components
                            end do
                        end if
                    end do
                end do
            end do
        else
            ! Iterative Kosaraju: finish-order DFS, then reverse-arc flood.
            allocate(visited(n), path(n), cursor(n), finish(n))
            visited = .false.
            n_finished = 0
            do i = 1, n
                if (visited(i)) cycle
                top = 1
                path(top) = i
                cursor(top) = graph%indptr(i)
                visited(i) = .true.
                do while (top > 0)
                    u = path(top)
                    if (cursor(top) < graph%indptr(u+1)) then
                        v = graph%indices(cursor(top))
                        cursor(top) = cursor(top)+1
                        if (visited(v)) cycle
                        top = top+1
                        path(top) = v
                        cursor(top) = graph%indptr(v)
                        visited(v) = .true.
                    else
                        n_finished = n_finished+1
                        finish(n_finished) = u
                        top = top-1
                    end if
                end do
            end do
            do i = n, 1, -1
                v = finish(i)
                if (labels(v) /= 0) cycle
                n_components = n_components+1
                head = 1
                tail = 1
                queue(1) = v
                labels(v) = n_components
                do while (head <= tail)
                    u = queue(head)
                    head = head+1
                    do k = reverse%indptr(u), reverse%indptr(u+1)-1
                        c = reverse%indices(k)
                        if (labels(c) /= 0) cycle
                        tail = tail+1
                        queue(tail) = c
                        labels(c) = n_components
                    end do
                end do
            end do
        end if
        status = csgraph_success
    end subroutine connected_components

    subroutine breadth_first_order(graph, start, order, predecessors, status, directed)
        type(csr_graph), intent(in) :: graph !! Valid graph in one-based CSR format.
        integer, intent(in) :: start !! One-based start node.
        integer, allocatable, intent(out) :: order(:) !! Visited nodes in BFS order.
        integer, allocatable, intent(out) :: predecessors(:) !! One-based predecessor; zero for root/unvisited.
        integer, intent(out) :: status !! Zero on success.
        logical, intent(in), optional :: directed !! Default true; false follows reverse arcs too.
        type(csr_graph) :: reverse
        logical, allocatable :: visited(:)
        integer, allocatable :: queue(:)
        integer :: n, head, tail, u, k, v, which
        logical :: is_directed

        status = csgraph_invalid_input
        allocate(order(0),predecessors(0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        if (start < 1 .or. start > n) return
        is_directed = .true.
        if (present(directed)) is_directed = directed
        if (.not. is_directed) call reverse_graph(graph,reverse)
        allocate(visited(n),queue(n))
        visited = .false.
        deallocate(predecessors)
        allocate(predecessors(n))
        predecessors = 0
        head = 1
        tail = 1
        queue(1) = start
        visited(start) = .true.
        do while (head <= tail)
            u = queue(head)
            head = head+1
            do which = 1, 2
                if (which == 2 .and. is_directed) exit
                if (which == 1) then
                    do k = graph%indptr(u), graph%indptr(u+1)-1
                        v = graph%indices(k)
                        if (visited(v)) cycle
                        visited(v) = .true.
                        predecessors(v) = u
                        tail = tail+1
                        queue(tail) = v
                    end do
                else
                    do k = reverse%indptr(u), reverse%indptr(u+1)-1
                        v = reverse%indices(k)
                        if (visited(v)) cycle
                        visited(v) = .true.
                        predecessors(v) = u
                        tail = tail+1
                        queue(tail) = v
                    end do
                end if
            end do
        end do
        order = queue(:tail)
        status = csgraph_success
    end subroutine breadth_first_order

    ! Binary min-heap. Duplicate records are deliberately allowed, with
    ! stale records ignored after extraction. Ties choose the lower vertex.
    pure logical function heap_less(key_a, vertex_a, key_b, vertex_b) result(ok)
        real(dp), intent(in) :: key_a, key_b !! Distances.
        integer, intent(in) :: vertex_a, vertex_b !! One-based vertices.
        ok = key_a < key_b
        if (key_a == key_b) ok = vertex_a < vertex_b
    end function heap_less

    pure subroutine heap_push(keys, vertices, heap_size, key, vertex)
        real(dp), intent(inout) :: keys(:) !! Heap distance keys.
        integer, intent(inout) :: vertices(:) !! Heap vertex IDs.
        integer, intent(inout) :: heap_size !! Current number of heap entries.
        real(dp), intent(in) :: key !! New tentative distance.
        integer, intent(in) :: vertex !! One-based vertex to enqueue.
        integer :: i, parent

        heap_size = heap_size+1
        i = heap_size
        do while (i > 1)
            parent = i/2
            if (.not. heap_less(key,vertex,keys(parent),vertices(parent))) exit
            keys(i) = keys(parent)
            vertices(i) = vertices(parent)
            i = parent
        end do
        keys(i) = key
        vertices(i) = vertex
    end subroutine heap_push

    pure subroutine heap_pop(keys, vertices, heap_size, key, vertex)
        real(dp), intent(inout) :: keys(:) !! Heap distance keys.
        integer, intent(inout) :: vertices(:) !! Heap vertex IDs.
        integer, intent(inout) :: heap_size !! Current number of heap entries.
        real(dp), intent(out) :: key !! Minimum tentative distance.
        integer, intent(out) :: vertex !! Minimum vertex.
        real(dp) :: last_key
        integer :: last_vertex, i, child

        key = keys(1)
        vertex = vertices(1)
        last_key = keys(heap_size)
        last_vertex = vertices(heap_size)
        heap_size = heap_size-1
        if (heap_size == 0) return
        i = 1
        do while (i <= heap_size/2)
            child = 2*i
            if (child < heap_size) then
                if (heap_less(keys(child+1),vertices(child+1), &
                    keys(child),vertices(child))) child = child+1
            end if
            if (.not. heap_less(keys(child),vertices(child),last_key,last_vertex)) exit
            keys(i) = keys(child)
            vertices(i) = vertices(child)
            i = child
        end do
        keys(i) = last_key
        vertices(i) = last_vertex
    end subroutine heap_pop

    subroutine relax_row(graph, u, current, bound, distances, predecessors, &
        heap_keys, heap_vertices, heap_size)
        type(csr_graph), intent(in) :: graph !! CSR graph in traversal orientation.
        integer, intent(in) :: u !! Current source node.
        real(dp), intent(in) :: current !! Finalized source distance.
        real(dp), intent(in) :: bound !! Maximum permitted path cost.
        real(dp), intent(inout) :: distances(:) !! Tentative best distances.
        integer, intent(inout) :: predecessors(:) !! Previous node on best path.
        real(dp), intent(inout) :: heap_keys(:) !! Min-heap distance storage.
        integer, intent(inout) :: heap_vertices(:) !! Min-heap vertex storage.
        integer, intent(inout) :: heap_size !! Min-heap entry count.
        integer :: k, v
        real(dp) :: trial

        do k = graph%indptr(u), graph%indptr(u+1)-1
            v = graph%indices(k)
            trial = current+graph%data(k)
            if (.not. ieee_is_finite(trial)) cycle
            if (trial >= distances(v) .or. trial > bound) cycle
            distances(v) = trial
            predecessors(v) = u
            call heap_push(heap_keys,heap_vertices,heap_size,trial,v)
        end do
    end subroutine relax_row

    subroutine dijkstra(graph, source, distances, predecessors, status, directed, limit)
        type(csr_graph), intent(in) :: graph !! Nonnegative-weight CSR graph.
        integer, intent(in) :: source !! One-based source vertex.
        real(dp), allocatable, intent(out) :: distances(:) !! Infinity for unreachable vertices.
        integer, allocatable, intent(out) :: predecessors(:) !! One-based; zero for source/unreachable.
        integer, intent(out) :: status !! Zero on success.
        logical, intent(in), optional :: directed !! Default true; false uses arcs both directions.
        real(dp), intent(in), optional :: limit !! Maximum path length, default infinity.
        type(csr_graph) :: reverse
        real(dp), allocatable :: heap_keys(:)
        integer, allocatable :: heap_vertices(:)
        real(dp) :: inf, bound, current
        integer :: n, m, capacity, heap_size, u
        logical :: is_directed

        status = csgraph_invalid_input
        allocate(distances(0), predecessors(0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        if (source < 1 .or. source > n) return
        if (any(graph%data < 0.0_dp)) return
        inf = ieee_value(0.0_dp,ieee_positive_inf)
        bound = inf
        if (present(limit)) then
            if (ieee_is_nan(limit)) return
            if (limit < 0.0_dp) return
            bound = limit
        end if
        is_directed = .true.
        if (present(directed)) is_directed = directed
        if (.not. is_directed) call reverse_graph(graph,reverse)
        m = size(graph%data)
        ! Each arc is relaxed at most once per traversal orientation.
        if (2_int64*int(m,int64)+1_int64 > int(huge(m),int64)) return
        capacity = 2*m+1
        allocate(heap_keys(capacity), heap_vertices(capacity))
        deallocate(distances,predecessors)
        allocate(distances(n),predecessors(n))
        distances = inf
        predecessors = 0
        distances(source) = 0.0_dp
        heap_size = 0
        call heap_push(heap_keys,heap_vertices,heap_size,0.0_dp,source)
        do while (heap_size > 0)
            call heap_pop(heap_keys,heap_vertices,heap_size,current,u)
            if (current > distances(u)) cycle
            if (current > bound) exit
            call relax_row(graph,u,current,bound,distances,predecessors, &
                heap_keys,heap_vertices,heap_size)
            if (.not. is_directed) then
                call relax_row(reverse,u,current,bound,distances,predecessors, &
                    heap_keys,heap_vertices,heap_size)
            end if
        end do
        status = csgraph_success
    end subroutine dijkstra

    subroutine bellman_ford(graph, source, distances, predecessors, status, directed)
        type(csr_graph), intent(in) :: graph !! Finite-weight CSR graph; negative edges allowed.
        integer, intent(in) :: source !! One-based starting vertex.
        real(dp), allocatable, intent(out) :: distances(:) !! Infinity for unreachable vertices.
        integer, allocatable, intent(out) :: predecessors(:) !! One-based; zero for root/unreachable.
        integer, intent(out) :: status !! Success, invalid input, or reachable negative cycle.
        logical, intent(in), optional :: directed !! Default true; false follows arcs in both directions.
        type(csr_graph) :: reverse
        real(dp) :: inf
        integer :: n, iteration
        logical :: is_directed, changed

        status = csgraph_invalid_input
        allocate(distances(0),predecessors(0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        if (source < 1 .or. source > n) return
        is_directed = .true.
        if (present(directed)) is_directed = directed
        if (.not. is_directed) call reverse_graph(graph,reverse)
        inf = ieee_value(0.0_dp,ieee_positive_inf)
        deallocate(distances,predecessors)
        allocate(distances(n),predecessors(n))
        distances = inf
        predecessors = 0
        distances(source) = 0.0_dp
        do iteration = 1, n-1
            changed = .false.
            call relax_bellman_ford(graph,distances,predecessors,changed,.true.)
            if (.not. is_directed) call relax_bellman_ford(reverse,distances,predecessors,changed,.true.)
            if (.not. changed) exit
        end do
        changed = .false.
        call relax_bellman_ford(graph,distances,predecessors,changed,.false.)
        if (.not. is_directed) call relax_bellman_ford(reverse,distances,predecessors,changed,.false.)
        if (changed) then
            status = csgraph_negative_cycle
            deallocate(distances,predecessors)
            allocate(distances(0),predecessors(0))
            return
        end if
        status = csgraph_success
    end subroutine bellman_ford

    subroutine relax_bellman_ford(graph, distances, predecessors, changed, update)
        type(csr_graph), intent(in) :: graph !! Directed CSR adjacency for one relaxation pass.
        real(dp), intent(inout) :: distances(:) !! Current best distances.
        integer, intent(inout) :: predecessors(:) !! Current best predecessors.
        logical, intent(inout) :: changed !! True if any arc improves a shortest path.
        logical, intent(in) :: update !! If false, only check for negative cycles.
        integer :: u, v, k
        real(dp) :: candidate

        do u = 1, graph%n
            if (.not. ieee_is_finite(distances(u))) cycle
            do k = graph%indptr(u), graph%indptr(u+1)-1
                v = graph%indices(k)
                candidate = distances(u)+graph%data(k)
                if (.not. ieee_is_finite(candidate)) cycle
                if (candidate >= distances(v)) cycle
                changed = .true.
                if (.not. update) return
                distances(v) = candidate
                predecessors(v) = u
            end do
        end do
    end subroutine relax_bellman_ford

    pure logical function edge_less(a,b) result(ok)
        type(graph_edge), intent(in) :: a, b !! Edges compared by weight and end points.
        ok = a%weight < b%weight
        if (a%weight == b%weight) then
            ok = a%source < b%source
            if (a%source == b%source) ok = a%target < b%target
        end if
    end function edge_less

    subroutine sort_edges(edges)
        type(graph_edge), intent(inout) :: edges(:) !! Edges sorted ascending by weight.
        type(graph_edge), allocatable :: temp(:)
        integer :: n, width, left, middle, right, i, j, k

        n = size(edges)
        allocate(temp(n))
        width = 1
        do while (width < n)
            left = 1
            do while (left <= n)
                middle = min(n,left+width-1)
                right = min(n,middle+width)
                i = left
                j = middle+1
                do k = left, right
                    if (i > middle) then
                        temp(k) = edges(j)
                        j = j+1
                    else if (j > right) then
                        temp(k) = edges(i)
                        i = i+1
                    else if (edge_less(edges(j),edges(i))) then
                        temp(k) = edges(j)
                        j = j+1
                    else
                        temp(k) = edges(i)
                        i = i+1
                    end if
                end do
                left = right+1
            end do
            edges = temp
            if (width > n/2) exit
            width = width*2
        end do
    end subroutine sort_edges

    pure integer function disjoint_root(parent, v) result(root)
        integer, intent(in) :: parent(:) !! Union-find parent links.
        integer, intent(in) :: v !! Vertex to query.

        root = v
        do while (parent(root) /= root)
            root = parent(root)
        end do
    end function disjoint_root

    function minimum_spanning_tree(graph, status) result(forest)
        type(csr_graph), intent(in) :: graph !! Finite-weight directed or undirected graph.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_graph) :: forest
        type(graph_edge), allocatable :: edges(:)
        integer, allocatable :: parent(:), rank(:), sources(:), targets(:)
        real(dp), allocatable :: weights(:)
        integer :: n, i, j, k, m, a, b, num, code

        if (present(status)) status = csgraph_invalid_input
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        m = count(graph%data /= 0.0_dp)
        allocate(edges(m),parent(n),rank(n),sources(max(0,n-1)), &
            targets(max(0,n-1)),weights(max(0,n-1)))
        k = 0
        do i = 1, n
            do j = graph%indptr(i), graph%indptr(i+1)-1
                if (graph%data(j) == 0.0_dp) cycle
                k = k+1
                edges(k)%source = i
                edges(k)%target = graph%indices(j)
                edges(k)%weight = graph%data(j)
            end do
        end do
        call sort_edges(edges)
        do i = 1, n
            parent(i) = i
        end do
        rank = 0
        num = 0
        do k = 1, m
            a = disjoint_root(parent,edges(k)%source)
            b = disjoint_root(parent,edges(k)%target)
            if (a == b) cycle
            num = num+1
            sources(num) = edges(k)%source
            targets(num) = edges(k)%target
            weights(num) = edges(k)%weight
            if (rank(a) < rank(b)) then
                parent(a) = b
            else if (rank(a) > rank(b)) then
                parent(b) = a
            else
                parent(b) = a
                rank(a) = rank(a)+1
            end if
        end do
        forest = csgraph_from_edges(n,sources(:num),targets(:num),weights(:num),code)
        if (present(status)) status = code
    end function minimum_spanning_tree


    ! Floyd-Warshall computes all-pairs paths in O(n**3) time and O(n**2)
    ! storage. Stored zero-weight arcs are respected (unlike dense zero gaps).
    subroutine floyd_warshall(graph, distances, predecessors, status, directed, unweighted)
        type(csr_graph), intent(in) :: graph !! Valid one-based CSR graph.
        real(dp), allocatable, intent(out) :: distances(:,:) !! n by n shortest-path distances.
        integer, allocatable, intent(out) :: predecessors(:,:) !! One-based previous nodes; zero if absent.
        integer, intent(out) :: status !! Success, invalid input, or negative cycle.
        logical, intent(in), optional :: directed !! Default true; false allows both arc directions.
        logical, intent(in), optional :: unweighted !! Default false; true gives every arc unit cost.
        real(dp) :: inf, weight, candidate
        integer :: n, i, j, k, p
        logical :: is_directed, is_unweighted

        status = csgraph_invalid_input
        allocate(distances(0,0), predecessors(0,0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        is_directed = .true.
        if (present(directed)) is_directed = directed
        is_unweighted = .false.
        if (present(unweighted)) is_unweighted = unweighted
        inf = ieee_value(0.0_dp, ieee_positive_inf)
        deallocate(distances, predecessors)
        allocate(distances(n,n), predecessors(n,n))
        distances = inf
        predecessors = 0
        do i = 1, n
            distances(i,i) = 0.0_dp
        end do
        do i = 1, n
            do p = graph%indptr(i), graph%indptr(i+1)-1
                j = graph%indices(p)
                weight = graph%data(p)
                if (is_unweighted) weight = 1.0_dp
                if (weight < distances(i,j)) then
                    distances(i,j) = weight
                    if (i /= j) predecessors(i,j) = i
                end if
                if (.not. is_directed) then
                    if (weight < distances(j,i)) then
                        distances(j,i) = weight
                        if (i /= j) predecessors(j,i) = j
                    end if
                end if
            end do
        end do
        do k = 1, n
            do i = 1, n
                if (.not. ieee_is_finite(distances(i,k))) cycle
                do j = 1, n
                    if (.not. ieee_is_finite(distances(k,j))) cycle
                    candidate = distances(i,k) + distances(k,j)
                    if (.not. ieee_is_finite(candidate)) cycle
                    if (candidate >= distances(i,j)) cycle
                    distances(i,j) = candidate
                    predecessors(i,j) = predecessors(k,j)
                end do
            end do
        end do
        do i = 1, n
            if (distances(i,i) < 0.0_dp) then
                status = csgraph_negative_cycle
                deallocate(distances, predecessors)
                allocate(distances(0,0), predecessors(0,0))
                return
            end if
        end do
        status = csgraph_success
    end subroutine floyd_warshall

    ! Potential relaxation for Johnson: all vertices initially have potential
    ! zero, as if they were connected to a synthetic source by zero arcs.
    subroutine relax_potentials(graph, potentials, changed, update)
        type(csr_graph), intent(in) :: graph !! Directed CSR edge list.
        real(dp), intent(inout) :: potentials(:) !! Current Bellman-Ford potentials.
        logical, intent(inout) :: changed !! Whether any edge can be relaxed.
        logical, intent(in) :: update !! True to update the potentials.
        integer :: u, v, k
        real(dp) :: candidate

        do u = 1, graph%n
            do k = graph%indptr(u), graph%indptr(u+1)-1
                v = graph%indices(k)
                candidate = potentials(u) + graph%data(k)
                if (.not. ieee_is_finite(candidate)) cycle
                if (candidate >= potentials(v)) cycle
                changed = .true.
                if (.not. update) return
                potentials(v) = candidate
            end do
        end do
    end subroutine relax_potentials

    ! Relax the edges in the reweighted graph for one Dijkstra queue item.
    subroutine relax_johnson_row(graph, u, current, potentials, unweighted, &
        distances, predecessors, heap_keys, heap_vertices, heap_size)
        type(csr_graph), intent(in) :: graph !! CSR edges in traversal direction.
        integer, intent(in) :: u !! Node being expanded.
        real(dp), intent(in) :: current !! Reweighted distance to u.
        real(dp), intent(in) :: potentials(:) !! Johnson vertex potentials.
        logical, intent(in) :: unweighted !! All edges cost one when true.
        real(dp), intent(inout) :: distances(:) !! Tentative reweighted distances.
        integer, intent(inout) :: predecessors(:) !! Previous vertex on each path.
        real(dp), intent(inout) :: heap_keys(:) !! Binary min-heap keys.
        integer, intent(inout) :: heap_vertices(:) !! Binary min-heap vertices.
        integer, intent(inout) :: heap_size !! Number of heap entries.
        real(dp) :: weight, trial
        integer :: k, v

        do k = graph%indptr(u), graph%indptr(u+1)-1
            v = graph%indices(k)
            weight = graph%data(k)
            if (unweighted) then
                weight = 1.0_dp
            else
                ! In exact arithmetic every reduced weight is nonnegative.
                weight = max(0.0_dp, (weight + potentials(u)) - potentials(v))
            end if
            trial = current + weight
            if (.not. ieee_is_finite(trial)) cycle
            if (trial >= distances(v)) cycle
            distances(v) = trial
            predecessors(v) = u
            call heap_push(heap_keys, heap_vertices, heap_size, trial, v)
        end do
    end subroutine relax_johnson_row

    ! Johnson's algorithm supports negative arcs but not negative cycles.
    ! Optional source indices select an arbitrary ordered subset of rows.
    subroutine johnson(graph, distances, predecessors, status, directed, unweighted, indices)
        type(csr_graph), intent(in) :: graph !! Valid one-based CSR graph.
        real(dp), allocatable, intent(out) :: distances(:,:) !! Selected-source by n path costs.
        integer, allocatable, intent(out) :: predecessors(:,:) !! One-based parents; zero if absent.
        integer, intent(out) :: status !! Success, invalid input, or negative cycle.
        logical, intent(in), optional :: directed !! Default true; false uses arcs in either direction.
        logical, intent(in), optional :: unweighted !! Default false; true makes each arc unit weight.
        integer, intent(in), optional :: indices(:) !! Selected one-based source nodes, default all.
        type(csr_graph) :: reverse
        real(dp), allocatable :: potential(:), reweighted(:), heap_keys(:)
        integer, allocatable :: sources(:), parent(:), heap_vertices(:)
        real(dp) :: inf, current
        integer :: n, m, capacity, row, source, u, heap_size, iteration, j
        logical :: is_directed, is_unweighted, changed

        status = csgraph_invalid_input
        allocate(distances(0,0), predecessors(0,0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        m = size(graph%data)
        if (2_int64*int(m,int64)+1_int64 > int(huge(m),int64)) return
        if (present(indices)) then
            if (any(indices < 1) .or. any(indices > n)) return
            sources = indices
        else
            allocate(sources(n))
            do j = 1, n
                sources(j) = j
            end do
        end if
        is_directed = .true.
        if (present(directed)) is_directed = directed
        is_unweighted = .false.
        if (present(unweighted)) is_unweighted = unweighted
        if (.not. is_directed) call reverse_graph(graph,reverse)
        allocate(potential(n))
        potential = 0.0_dp
        if (.not. is_unweighted) then
            do iteration = 1, n-1
                changed = .false.
                call relax_potentials(graph,potential,changed,.true.)
                if (.not. is_directed) call relax_potentials(reverse,potential,changed,.true.)
                if (.not. changed) exit
            end do
            changed = .false.
            call relax_potentials(graph,potential,changed,.false.)
            if (.not. is_directed) call relax_potentials(reverse,potential,changed,.false.)
            if (changed) then
                status = csgraph_negative_cycle
                return
            end if
        end if
        inf = ieee_value(0.0_dp,ieee_positive_inf)
        capacity = 2*m+1
        allocate(heap_keys(capacity),heap_vertices(capacity),reweighted(n),parent(n))
        deallocate(distances,predecessors)
        allocate(distances(size(sources),n),predecessors(size(sources),n))
        distances = inf
        predecessors = 0
        do row = 1, size(sources)
            source = sources(row)
            reweighted = inf
            parent = 0
            reweighted(source) = 0.0_dp
            heap_size = 0
            call heap_push(heap_keys,heap_vertices,heap_size,0.0_dp,source)
            do while (heap_size > 0)
                call heap_pop(heap_keys,heap_vertices,heap_size,current,u)
                if (current > reweighted(u)) cycle
                call relax_johnson_row(graph,u,current,potential,is_unweighted, &
                    reweighted,parent,heap_keys,heap_vertices,heap_size)
                if (.not. is_directed) then
                    call relax_johnson_row(reverse,u,current,potential,is_unweighted, &
                        reweighted,parent,heap_keys,heap_vertices,heap_size)
                end if
            end do
            predecessors(row,:) = parent
            do j = 1, n
                if (.not. ieee_is_finite(reweighted(j))) cycle
                distances(row,j) = reweighted(j) - potential(source) + potential(j)
            end do
            distances(row,source) = 0.0_dp
            predecessors(row,source) = 0
        end do
        status = csgraph_success
    end subroutine johnson

    ! Return a dense real(dp) graph Laplacian, including isolated vertices.
    ! In-degree (column sum) matches scipy.sparse.csgraph.laplacian's default.
    subroutine laplacian(graph, lap, diagonal, status, normed, use_out_degree, symmetrized)
        type(csr_graph), intent(in) :: graph !! Finite CSR edge weights.
        real(dp), allocatable, intent(out) :: lap(:,:) !! Dense n-by-n Laplacian.
        real(dp), allocatable, intent(out) :: diagonal(:) !! Degrees or normalized scaling factors.
        integer, intent(out) :: status !! Invalid on negative normalized degree or overflow.
        logical, intent(in), optional :: normed !! Default false; symmetric normalization if true.
        logical, intent(in), optional :: use_out_degree !! Default false uses incoming degrees.
        logical, intent(in), optional :: symmetrized !! Default false; true adds A and transpose(A).
        real(dp), allocatable :: degree(:)
        integer :: n, i, j, k
        logical :: normalize, outgoing, symmetrize

        status = csgraph_invalid_input
        allocate(lap(0,0), diagonal(0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        normalize = .false.
        if (present(normed)) normalize = normed
        outgoing = .false.
        if (present(use_out_degree)) outgoing = use_out_degree
        symmetrize = .false.
        if (present(symmetrized)) symmetrize = symmetrized
        deallocate(lap,diagonal)
        allocate(lap(n,n),diagonal(n),degree(n))
        lap = 0.0_dp
        do i = 1, n
            do k = graph%indptr(i), graph%indptr(i+1)-1
                j = graph%indices(k)
                if (i == j) cycle
                lap(i,j) = lap(i,j) - graph%data(k)
                if (symmetrize) lap(j,i) = lap(j,i) - graph%data(k)
            end do
        end do
        if (outgoing) then
            degree = -sum(lap,dim=2)
        else
            degree = -sum(lap,dim=1)
        end if
        if (.not. all(ieee_is_finite(lap)) .or. .not. all(ieee_is_finite(degree))) then
            deallocate(lap,diagonal)
            allocate(lap(0,0),diagonal(0))
            return
        end if
        if (normalize) then
            if (any(degree < 0.0_dp)) then
                deallocate(lap,diagonal)
                allocate(lap(0,0),diagonal(0))
                return
            end if
            diagonal = 1.0_dp
            where (degree > 0.0_dp) diagonal = sqrt(degree)
            do i = 1, n
                do j = 1, n
                    if (i == j) cycle
                    lap(i,j) = (lap(i,j)/diagonal(i))/diagonal(j)
                end do
                lap(i,i) = 0.0_dp
                if (degree(i) /= 0.0_dp) lap(i,i) = 1.0_dp
            end do
        else
            diagonal = degree
            do i = 1, n
                lap(i,i) = degree(i)
            end do
        end if
        if (.not. all(ieee_is_finite(lap))) then
            deallocate(lap,diagonal)
            allocate(lap(0,0),diagonal(0))
            return
        end if
        status = csgraph_success
    end subroutine laplacian

    ! Iterative depth-first traversal: no recursion or stack-size dependency.
    ! In undirected mode, outgoing arcs are explored before incoming arcs.
    subroutine depth_first_order(graph, start, order, predecessors, status, directed)
        type(csr_graph), intent(in) :: graph !! Valid one-based CSR graph.
        integer, intent(in) :: start !! One-based initial vertex.
        integer, allocatable, intent(out) :: order(:) !! Visited nodes, DFS discovery order.
        integer, allocatable, intent(out) :: predecessors(:) !! Parent IDs, zero for root/unvisited.
        integer, intent(out) :: status !! Zero on success.
        logical, intent(in), optional :: directed !! Default true; false follows incoming arcs too.
        type(csr_graph) :: reverse
        integer, allocatable :: stack(:), next_out(:), next_in(:), discovery(:)
        logical, allocatable :: seen(:)
        integer :: n, top, found, u, k, v, count_seen
        logical :: is_directed

        status = csgraph_invalid_input
        allocate(order(0),predecessors(0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        if (start < 1 .or. start > n) return
        is_directed = .true.
        if (present(directed)) is_directed = directed
        if (.not. is_directed) call reverse_graph(graph,reverse)
        allocate(stack(n),next_out(n),next_in(n),seen(n),discovery(n))
        next_out = graph%indptr(:n)
        if (.not. is_directed) next_in = reverse%indptr(:n)
        deallocate(predecessors)
        allocate(predecessors(n))
        predecessors = 0
        seen = .false.
        top = 1
        stack(1) = start
        seen(start) = .true.
        count_seen = 1
        discovery(1) = start
        do while (top > 0)
            u = stack(top)
            found = 0
            do while (next_out(u) < graph%indptr(u+1))
                k = next_out(u)
                next_out(u) = k+1
                v = graph%indices(k)
                if (seen(v)) cycle
                found = v
                exit
            end do
            if (found == 0 .and. .not. is_directed) then
                do while (next_in(u) < reverse%indptr(u+1))
                    k = next_in(u)
                    next_in(u) = k+1
                    v = reverse%indices(k)
                    if (seen(v)) cycle
                    found = v
                    exit
                end do
            end if
            if (found == 0) then
                top = top-1
                cycle
            end if
            seen(found) = .true.
            predecessors(found) = u
            count_seen = count_seen+1
            discovery(count_seen) = found
            top = top+1
            stack(top) = found
        end do
        order = discovery(:count_seen)
        status = csgraph_success
    end subroutine depth_first_order

    ! Return the one-based vertex sequence from a source to a target using
    ! either a single-source predecessor vector or one selected APSP row.
    subroutine path_from_predecessors(source, target, predecessors, path, status)
        integer, intent(in) :: source !! One-based source vertex.
        integer, intent(in) :: target !! One-based target vertex.
        integer, intent(in) :: predecessors(:) !! One-based predecessor vector; zero for absent.
        integer, allocatable, intent(out) :: path(:) !! Source-to-target nodes, empty if unreachable.
        integer, intent(out) :: status !! Zero on a valid vector and indices.
        integer, allocatable :: reverse_order(:)
        logical, allocatable :: seen(:)
        integer :: n, v, count_nodes, i

        status = csgraph_invalid_input
        allocate(path(0))
        n = size(predecessors)
        if (source < 1 .or. source > n .or. target < 1 .or. target > n) return
        if (any(predecessors < 0) .or. any(predecessors > n)) return
        allocate(reverse_order(n),seen(n))
        seen = .false.
        v = target
        count_nodes = 0
        do while (v /= 0)
            if (seen(v)) return
            seen(v) = .true.
            count_nodes = count_nodes+1
            reverse_order(count_nodes) = v
            if (v == source) then
                deallocate(path)
                allocate(path(count_nodes))
                do i = 1, count_nodes
                    path(i) = reverse_order(count_nodes+1-i)
                end do
                status = csgraph_success
                return
            end if
            v = predecessors(v)
        end do
        status = csgraph_success
    end subroutine path_from_predecessors


    ! All-pairs or selected-source shortest path interface. This exposes the
    ! same algorithms as scipy.sparse.csgraph.shortest_path; source IDs and
    ! predecessor IDs are one-based, with zero for unavailable parents.
    subroutine shortest_path(graph, distances, predecessors, status, method, directed, unweighted, indices)
        type(csr_graph), intent(in) :: graph !! Valid one-based CSR graph.
        real(dp), allocatable, intent(out) :: distances(:,:) !! Selected source by n shortest distances.
        integer, allocatable, intent(out) :: predecessors(:,:) !! One-based predecessor matrix.
        integer, intent(out) :: status !! Success, invalid input, or negative cycle.
        character(len=*), intent(in), optional :: method !! 'auto', 'FW', 'D', 'BF', 'J'.
        logical, intent(in), optional :: directed !! Default true; false allows reverse arcs.
        logical, intent(in), optional :: unweighted !! Default false; each arc costs one if true.
        integer, intent(in), optional :: indices(:) !! Ordered one-based source IDs, default all.
        type(csr_graph) :: working
        real(dp), allocatable :: all_distances(:,:), vector_distances(:)
        integer, allocatable :: all_predecessors(:,:), vector_predecessors(:), sources(:)
        integer :: n, m, i, row
        logical :: is_directed, is_unweighted
        character(len=:), allocatable :: algorithm
        real(dp) :: inf

        status = csgraph_invalid_input
        allocate(distances(0,0),predecessors(0,0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        m = size(graph%data)
        if (present(indices)) then
            if (any(indices < 1) .or. any(indices > n)) return
            sources = indices
        else
            allocate(sources(n))
            do i = 1, n
                sources(i) = i
            end do
        end if
        algorithm = 'auto'
        if (present(method)) algorithm = trim(method)
        select case (algorithm)
        case ('auto', 'FW', 'D', 'BF', 'J')
        case default
            return
        end select
        is_directed = .true.
        if (present(directed)) is_directed = directed
        is_unweighted = .false.
        if (present(unweighted)) is_unweighted = unweighted
        if (algorithm == 'auto') then
            algorithm = 'J'
            if (real(m,dp) >= 0.25_dp*real(n,dp)**2 .and. n <= 256) algorithm = 'FW'
            if (size(sources) < n) algorithm = 'J'
        end if
        select case (algorithm)
        case ('FW')
            call floyd_warshall(graph,all_distances,all_predecessors,status, &
                directed=is_directed,unweighted=is_unweighted)
            if (status /= csgraph_success) return
            distances = all_distances(sources,:)
            predecessors = all_predecessors(sources,:)
        case ('J')
            call johnson(graph,all_distances,all_predecessors,status, &
                directed=is_directed,unweighted=is_unweighted,indices=sources)
            if (status /= csgraph_success) return
            distances = all_distances
            predecessors = all_predecessors
        case ('D', 'BF')
            working = graph
            if (is_unweighted) working%data = 1.0_dp
            if (algorithm == 'D' .and. any(working%data < 0.0_dp)) return
            inf = ieee_value(0.0_dp,ieee_positive_inf)
            deallocate(distances,predecessors)
            allocate(distances(size(sources),n),predecessors(size(sources),n))
            distances = inf
            predecessors = 0
            do row = 1, size(sources)
                if (algorithm == 'D') then
                    call dijkstra(working,sources(row),vector_distances,vector_predecessors, &
                        status,directed=is_directed)
                else
                    call bellman_ford(working,sources(row),vector_distances,vector_predecessors, &
                        status,directed=is_directed)
                end if
                if (status /= csgraph_success) then
                    deallocate(distances,predecessors)
                    allocate(distances(0,0),predecessors(0,0))
                    return
                end if
                distances(row,:) = vector_distances
                predecessors(row,:) = vector_predecessors
            end do
        end select
        status = csgraph_success
    end subroutine shortest_path

    ! Extract a spanning traversal tree, keeping original edge weights.
    ! In undirected mode, a reverse edge supplies the weight if needed.
    ! Reconstruct the oriented CSR tree represented by a predecessor vector,
    ! matching scipy.sparse.csgraph.reconstruct_path (with one-based indices).
    function reconstruct_path(graph, predecessors, status, directed) result(tree)
        type(csr_graph), intent(in) :: graph !! Valid original weighted CSR graph.
        integer, intent(in) :: predecessors(:) !! One-based predecessor per vertex; zero for roots.
        integer, intent(out) :: status !! Zero on success, invalid if path or CSR is malformed.
        logical, intent(in), optional :: directed !! Default true; false allows reverse arcs.
        type(csr_graph) :: tree
        integer, allocatable :: sources(:), targets(:), visited_from(:)
        real(dp), allocatable :: weights(:)
        integer :: n, m, i, j, k, u, v, count_edges
        logical :: found, is_directed

        status = csgraph_invalid_input
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        if (size(predecessors) /= n) return
        if (any(predecessors < 0) .or. any(predecessors > n)) return
        is_directed = .true.
        if (present(directed)) is_directed = directed
        allocate(visited_from(n))
        visited_from = 0
        ! A valid predecessor forest must not contain parent cycles.
        do i = 1, n
            v = i
            do while (v /= 0)
                if (visited_from(v) == i) return
                if (visited_from(v) /= 0) exit
                visited_from(v) = i
                v = predecessors(v)
            end do
        end do
        m = count(predecessors /= 0)
        allocate(sources(m),targets(m),weights(m))
        count_edges = 0
        do v = 1, n
            u = predecessors(v)
            if (u == 0) cycle
            found = .false.
            do k = graph%indptr(u),graph%indptr(u+1)-1
                if (graph%indices(k) /= v) cycle
                count_edges = count_edges+1
                sources(count_edges) = u
                targets(count_edges) = v
                weights(count_edges) = graph%data(k)
                found = .true.
                exit
            end do
            if (.not. found .and. .not. is_directed) then
                do k = graph%indptr(v),graph%indptr(v+1)-1
                    if (graph%indices(k) /= u) cycle
                    count_edges = count_edges+1
                    sources(count_edges) = u
                    targets(count_edges) = v
                    weights(count_edges) = graph%data(k)
                    found = .true.
                    exit
                end do
            end if
            if (.not. found) return
        end do
        tree = csgraph_from_edges(n,sources,targets,weights,status)
    end function reconstruct_path

    function traversal_tree(graph, start, status, depth_first, directed) result(tree)
        type(csr_graph), intent(in) :: graph !! Valid CSR graph.
        integer, intent(in) :: start !! One-based root vertex.
        integer, intent(out) :: status !! Zero on success.
        logical, intent(in) :: depth_first !! DFS if true, BFS otherwise.
        logical, intent(in), optional :: directed !! Default true.
        type(csr_graph) :: tree
        integer, allocatable :: order(:), parents(:)

        if (depth_first) then
            call depth_first_order(graph,start,order,parents,status,directed)
        else
            call breadth_first_order(graph,start,order,parents,status,directed)
        end if
        if (status /= csgraph_success) return
        tree = reconstruct_path(graph,parents,status,directed)
    end function traversal_tree

    function breadth_first_tree(graph, start, status, directed) result(tree)
        type(csr_graph), intent(in) :: graph !! Valid one-based CSR graph.
        integer, intent(in) :: start !! One-based root.
        integer, intent(out) :: status !! Zero on success.
        logical, intent(in), optional :: directed !! Default true.
        type(csr_graph) :: tree

        tree = traversal_tree(graph,start,status,.false.,directed)
    end function breadth_first_tree

    function depth_first_tree(graph, start, status, directed) result(tree)
        type(csr_graph), intent(in) :: graph !! Valid one-based CSR graph.
        integer, intent(in) :: start !! One-based root.
        integer, intent(out) :: status !! Zero on success.
        logical, intent(in), optional :: directed !! Default true.
        type(csr_graph) :: tree

        tree = traversal_tree(graph,start,status,.true.,directed)
    end function depth_first_tree

    ! ---------------------------------------------------------------
    ! Sparse graph orderings, flow, and rectangular structural matching
    ! ---------------------------------------------------------------

    subroutine reverse_cuthill_mckee(graph, permutation, status, symmetric_mode)
        type(csr_graph), intent(in) :: graph !! Valid square CSR; only structural entries count.
        integer, allocatable, intent(out) :: permutation(:) !! One-based new-to-old vertex order.
        integer, intent(out) :: status !! csgraph_success or csgraph_invalid_input.
        logical, intent(in), optional :: symmetric_mode !! True: input already structurally symmetric.
        type(csr_graph) :: reversed
        integer, allocatable :: offsets(:), neighbors(:), degree(:), mark(:), queue(:), &
            forward(:), heap_vertices(:)
        real(dp), allocatable :: heap_keys(:)
        logical, allocatable :: assigned(:)
        integer :: n, m, used, i, k, j, t, which, u, v, root, candidate, next_root, &
            height, next_height, nordered, head, tail, heap_size
        logical :: sym
        real(dp) :: heap_key
        status = csgraph_invalid_input
        allocate(permutation(0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        m = size(graph%indices)
        if (m > huge(m)/2) return
        sym = .false.
        if (present(symmetric_mode)) sym = symmetric_mode
        if (.not. sym) call reverse_graph(graph,reversed)
        allocate(offsets(n+1), neighbors(2*m), degree(n), mark(n))
        degree = 0
        mark = 0
        used = 0
        do i = 1, n
            offsets(i) = used+1
            do which = 1, 2
                if (sym .and. which == 2) cycle
                if (which == 1) then
                    do k = graph%indptr(i), graph%indptr(i+1)-1
                        j = graph%indices(k)
                        if (j == i .or. mark(j) == i) cycle
                        mark(j) = i
                        used = used+1
                        neighbors(used) = j
                    end do
                else
                    do k = reversed%indptr(i), reversed%indptr(i+1)-1
                        j = reversed%indices(k)
                        if (j == i .or. mark(j) == i) cycle
                        mark(j) = i
                        used = used+1
                        neighbors(used) = j
                    end do
                end if
            end do
            degree(i) = used-offsets(i)+1
        end do
        offsets(n+1) = used+1
        allocate(assigned(n), forward(n), queue(n), heap_vertices(n), heap_keys(n))
        assigned = .false.
        nordered = 0
        do while (nordered < n)
            root = 0
            do i = 1, n
                if (assigned(i)) cycle
                if (root == 0) then
                    root = i
                else if (degree(i) < degree(root)) then
                    root = i
                end if
            end do
            ! Pseudo-peripheral root refinement: do not assume a particular
            ! SciPy tie-breaking permutation, only the structural objective.
            call rcm_farthest(root, offsets, neighbors, degree, assigned, height, candidate)
            do t = 1, n
                call rcm_farthest(candidate, offsets, neighbors, degree, assigned, next_height, next_root)
                if (next_height <= height) exit
                root = candidate
                height = next_height
                candidate = next_root
            end do
            root = candidate
            head = 1
            tail = 1
            queue(1) = root
            assigned(root) = .true.
            do while (head <= tail)
                u = queue(head)
                head = head+1
                nordered = nordered+1
                forward(nordered) = u
                heap_size = 0
                do k = offsets(u), offsets(u+1)-1
                    v = neighbors(k)
                    if (assigned(v)) cycle
                    assigned(v) = .true.
                    call heap_push(heap_keys,heap_vertices,heap_size,real(degree(v),dp),v)
                end do
                do while (heap_size > 0)
                    call heap_pop(heap_keys,heap_vertices,heap_size,heap_key,v)
                    tail = tail+1
                    queue(tail) = v
                end do
            end do
        end do
        deallocate(permutation)
        allocate(permutation(n))
        permutation = forward(n:1:-1)
        status = csgraph_success
    end subroutine reverse_cuthill_mckee

    subroutine rcm_farthest(root, offsets, neighbors, degree, assigned, maxlevel, candidate)
        integer, intent(in) :: root !! BFS starting vertex.
        integer, intent(in) :: offsets(:) !! CSR adjacency offsets.
        integer, intent(in) :: neighbors(:) !! CSR adjacency vertices.
        integer, intent(in) :: degree(:) !! Distinct-neighbor structural degree.
        logical, intent(in) :: assigned(:) !! Already ordered vertices.
        integer, intent(out) :: maxlevel !! Last BFS level.
        integer, intent(out) :: candidate !! Minimum-degree last-level vertex.
        integer, allocatable :: queue(:), dist(:)
        integer :: n, head, tail, u, v, k
        n = size(assigned)
        allocate(queue(n), dist(n))
        dist = -1
        dist(root) = 0
        head = 1
        tail = 1
        queue(1) = root
        candidate = root
        maxlevel = 0
        do while (head <= tail)
            u = queue(head)
            head = head+1
            if (dist(u) > maxlevel) then
                maxlevel = dist(u)
                candidate = u
            else if (dist(u) == maxlevel) then
                if (degree(u) < degree(candidate) .or. &
                    (degree(u) == degree(candidate) .and. u < candidate)) candidate = u
            end if
            do k = offsets(u), offsets(u+1)-1
                v = neighbors(k)
                if (assigned(v) .or. dist(v) >= 0) cycle
                dist(v) = dist(u)+1
                tail = tail+1
                queue(tail) = v
            end do
        end do
    end subroutine rcm_farthest

    function bipartite_from_edges(n_rows, n_cols, rows, columns, status) result(graph)
        integer, intent(in) :: n_rows !! Number of left vertices; >= 0.
        integer, intent(in) :: n_cols !! Number of right vertices; >= 0.
        integer, intent(in) :: rows(:) !! One-based left endpoints (duplicates allowed).
        integer, intent(in) :: columns(:) !! One-based right endpoints.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_bipartite_graph) :: graph
        integer, allocatable :: cursor(:)
        integer :: m, i, k, p
        if (present(status)) status = csgraph_invalid_input
        if (n_rows < 0 .or. n_cols < 0) return
        if (n_rows == huge(n_rows) .or. n_cols == huge(n_cols)) return
        m = size(rows)
        if (m == huge(m) .or. size(columns) /= m) return
        if (any(rows < 1) .or. any(rows > n_rows) .or. &
            any(columns < 1) .or. any(columns > n_cols)) return
        graph%n_rows = n_rows
        graph%n_cols = n_cols
        allocate(graph%indptr(n_rows+1), graph%indices(m), cursor(n_rows))
        graph%indptr = 0
        graph%indptr(1) = 1
        do k = 1, m
            graph%indptr(rows(k)+1) = graph%indptr(rows(k)+1)+1
        end do
        do i = 1, n_rows
            graph%indptr(i+1) = graph%indptr(i+1)+graph%indptr(i)
        end do
        cursor = graph%indptr(:n_rows)
        do k = 1, m
            p = cursor(rows(k))
            graph%indices(p) = columns(k)
            cursor(rows(k)) = p+1
        end do
        if (present(status)) status = csgraph_success
    end function bipartite_from_edges

    function bipartite_from_dense(a, status) result(graph)
        real(dp), intent(in) :: a(:,:) !! Rectangular matrix; nonzero entries become edges.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_bipartite_graph) :: graph
        integer, allocatable :: rows(:), cols(:)
        integer :: m, i, j, k
        if (present(status)) status = csgraph_invalid_input
        if (.not. all(ieee_is_finite(a))) return
        if (count(a /= 0.0_dp,kind=int64) >= int(huge(m),int64)) return
        m = count(a /= 0.0_dp)
        allocate(rows(m),cols(m))
        k = 0
        do i = 1, size(a,1)
            do j = 1, size(a,2)
                if (a(i,j) == 0.0_dp) cycle
                k = k+1
                rows(k) = i
                cols(k) = j
            end do
        end do
        graph = bipartite_from_edges(size(a,1),size(a,2),rows,cols,status)
    end function bipartite_from_dense

    function bipartite_from_csgraph(square, status) result(graph)
        type(csr_graph), intent(in) :: square !! Square CSR graph; stored zeros remain edges.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_bipartite_graph) :: graph
        if (present(status)) status = csgraph_invalid_input
        if (.not. valid_csgraph(square)) return
        graph%n_rows = square%n
        graph%n_cols = square%n
        graph%indptr = square%indptr
        graph%indices = square%indices
        if (present(status)) status = csgraph_success
    end function bipartite_from_csgraph

    pure logical function valid_bipartite_graph(graph) result(ok)
        type(csr_bipartite_graph), intent(in) :: graph !! Rectangular one-based CSR pattern.
        integer :: n, m
        ok = .false.
        n = graph%n_rows
        if (n < 0 .or. graph%n_cols < 0) return
        if (n == huge(n)) return
        if (.not. allocated(graph%indptr) .or. .not. allocated(graph%indices)) return
        m = size(graph%indices)
        if (m == huge(m)) return
        if (size(graph%indptr) /= n+1) return
        if (graph%indptr(1) /= 1 .or. graph%indptr(n+1) /= m+1) return
        if (any(graph%indptr < 1) .or. any(graph%indptr > m+1)) return
        if (any(graph%indptr(2:) < graph%indptr(:n))) return
        if (any(graph%indices < 1) .or. any(graph%indices > graph%n_cols)) return
        ok = .true.
    end function valid_bipartite_graph

    subroutine maximum_bipartite_matching(graph, matching, status, perm_type)
        type(csr_bipartite_graph), intent(in) :: graph !! Bipartite CSR pattern, rows vs columns.
        integer, allocatable, intent(out) :: matching(:) !! 1-based partner, zero for unmatched.
        integer, intent(out) :: status !! Zero on success.
        character(len=*), intent(in), optional :: perm_type !! 'row' (default) or 'column'.
        integer, allocatable :: row_match(:), col_match(:), level(:), queue(:), cursor(:), &
            stack_row(:), stack_col(:)
        integer :: nr, nc, head, tail, i, u, v, mate, k, shortest
        character(len=:), allocatable :: rule
        status = csgraph_invalid_input
        allocate(matching(0))
        if (.not. valid_bipartite_graph(graph)) return
        rule = 'row'
        if (present(perm_type)) rule = trim(perm_type)
        if (rule /= 'row' .and. rule /= 'column') return
        nr = graph%n_rows
        nc = graph%n_cols
        allocate(row_match(nr),col_match(nc),level(nr),queue(nr),cursor(nr), &
            stack_row(nr),stack_col(nr))
        row_match = 0
        col_match = 0
        do
            level = -1
            head = 1
            tail = 0
            do i = 1, nr
                if (row_match(i) /= 0) cycle
                level(i) = 0
                tail = tail+1
                queue(tail) = i
            end do
            shortest = huge(shortest)
            do while (head <= tail)
                u = queue(head)
                head = head+1
                if (level(u) >= shortest) cycle
                do k = graph%indptr(u),graph%indptr(u+1)-1
                    v = graph%indices(k)
                    mate = col_match(v)
                    if (mate == 0) then
                        shortest = level(u)+1
                    else if (level(mate) < 0) then
                        level(mate) = level(u)+1
                        tail = tail+1
                        queue(tail) = mate
                    end if
                end do
            end do
            if (shortest == huge(shortest)) exit
            cursor = graph%indptr(:nr)
            do i = 1, nr
                if (row_match(i) /= 0) cycle
                if (level(i) /= 0) cycle
                if (hk_augment(i,graph,row_match,col_match,level,cursor, &
                    shortest,stack_row,stack_col)) continue
            end do
        end do
        deallocate(matching)
        if (rule == 'row') then
            matching = col_match
        else
            matching = row_match
        end if
        status = csgraph_success
    end subroutine maximum_bipartite_matching

    logical function hk_augment(start, graph, row_match, col_match, level, cursor, &
        shortest, stack_row, stack_col) result(found)
        integer, intent(in) :: start !! Free starting row.
        type(csr_bipartite_graph), intent(in) :: graph !! Valid pattern.
        integer, intent(inout) :: row_match(:) !! Row-to-column matches.
        integer, intent(inout) :: col_match(:) !! Column-to-row matches.
        integer, intent(inout) :: level(:) !! Hopcroft--Karp BFS levels (pruned in DFS).
        integer, intent(inout) :: cursor(:) !! Current outgoing edge for each row.
        integer, intent(in) :: shortest !! Shortest augmenting path depth.
        integer, intent(inout) :: stack_row(:) !! Preallocated DFS stack.
        integer, intent(inout) :: stack_col(:) !! Chosen alternating edges.
        integer :: top, u, v, mate, k, t
        logical :: descend
        found = .false.
        top = 1
        stack_row(1) = start
        do while (top > 0)
            u = stack_row(top)
            descend = .false.
            do while (cursor(u) < graph%indptr(u+1))
                k = cursor(u)
                cursor(u) = cursor(u)+1
                v = graph%indices(k)
                mate = col_match(v)
                if (mate == 0) then
                    if (level(u)+1 /= shortest) cycle
                    row_match(u) = v
                    col_match(v) = u
                    do t = top-1, 1, -1
                        u = stack_row(t)
                        v = stack_col(t)
                        row_match(u) = v
                        col_match(v) = u
                    end do
                    found = .true.
                    return
                end if
                if (level(mate) /= level(u)+1) cycle
                stack_col(top) = v
                top = top+1
                stack_row(top) = mate
                descend = .true.
                exit
            end do
            if (descend) cycle
            level(u) = -1
            top = top-1
        end do
    end function hk_augment

    subroutine structural_rank(graph, rank, status)
        type(csr_bipartite_graph), intent(in) :: graph !! Rectangular structural CSR pattern.
        integer, intent(out) :: rank !! Maximum transversal cardinality.
        integer, intent(out) :: status !! Zero on success.
        integer, allocatable :: matching(:)
        rank = 0
        call maximum_bipartite_matching(graph, matching, status, perm_type='column')
        if (status == csgraph_success) rank = count(matching /= 0)
    end subroutine structural_rank

    subroutine maximum_flow(graph, source, sink, flow_value, flow, status)
        type(csr_graph), intent(in) :: graph !! Directed CSR with nonnegative integral capacities.
        integer, intent(in) :: source !! One-based source, distinct from sink.
        integer, intent(in) :: sink !! One-based sink.
        integer(int64), intent(out) :: flow_value !! Maximum integral source-to-sink flow.
        type(csr_graph), intent(out) :: flow !! Sparse antisymmetric *net* flow matrix.
        integer, intent(out) :: status !! Zero on success, invalid input otherwise.
        integer, allocatable :: head(:), nxt(:), dest(:), edge_src(:), level(:), queue(:), &
            cursor(:), path_vertex(:), path_edge(:), rows(:), cols(:), marks(:), touched(:), &
            out_rows(:), out_cols(:)
        integer(int64), allocatable :: residual(:), capacities(:), totals(:)
        real(dp), allocatable :: values(:), out_values(:)
        type(csr_graph) :: raw
        integer :: n, m, k, u, v, e, depth, qhead, qtail, used, ntouched, count_out, t
        integer(int64) :: delta, weight, net
        real(dp), parameter :: max_exact = 9007199254740991.0_dp
        flow_value = 0_int64
        status = csgraph_invalid_input
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        if (source < 1 .or. source > n .or. sink < 1 .or. sink > n .or. source == sink) return
        m = size(graph%indices)
        if (m > huge(m)/2) return
        if (any(graph%data < 0.0_dp) .or. any(graph%data > max_exact)) return
        if (any(graph%data /= anint(graph%data))) return
        allocate(head(n),nxt(2*m),dest(2*m),edge_src(m),level(n),queue(n), &
            cursor(n),path_vertex(n),path_edge(n),residual(2*m),capacities(m))
        head = 0
        do u = 1, n
            do k = graph%indptr(u),graph%indptr(u+1)-1
                edge_src(k) = u
                v = graph%indices(k)
                capacities(k) = int(graph%data(k),int64)
                e = 2*k-1
                dest(e) = v
                residual(e) = capacities(k)
                nxt(e) = head(u)
                head(u) = e
                dest(e+1) = u
                residual(e+1) = 0_int64
                nxt(e+1) = head(v)
                head(v) = e+1
            end do
        end do
        do
            level = -1
            level(source) = 0
            qhead = 1
            qtail = 1
            queue(1) = source
            do while (qhead <= qtail)
                u = queue(qhead)
                qhead = qhead+1
                e = head(u)
                do while (e /= 0)
                    v = dest(e)
                    if (residual(e) > 0_int64 .and. level(v) < 0) then
                        level(v) = level(u)+1
                        qtail = qtail+1
                        queue(qtail) = v
                    end if
                    e = nxt(e)
                end do
            end do
            if (level(sink) < 0) exit
            cursor = head
            u = source
            depth = 0
            do
                if (u == sink) then
                    delta = huge(delta)
                    do k = 1, depth
                        delta = min(delta,residual(path_edge(k)))
                    end do
                    if (delta <= 0_int64) exit ! Defensive: should be unreachable.
                    if (flow_value > huge(flow_value)-delta) return
                    do k = 1, depth
                        e = path_edge(k)
                        residual(e) = residual(e)-delta
                        ! Reverse capacity plus delta cannot exceed the initial
                        ! forward capacity for each original stored arc.
                        if (mod(e,2) == 1) then
                            residual(e+1) = residual(e+1)+delta
                        else
                            residual(e-1) = residual(e-1)+delta
                        end if
                    end do
                    flow_value = flow_value+delta
                    u = source
                    depth = 0
                    cycle
                end if
                e = cursor(u)
                do while (e /= 0)
                    if (residual(e) > 0_int64 .and. level(dest(e)) == level(u)+1) exit
                    e = nxt(e)
                end do
                cursor(u) = e
                if (e /= 0) then
                    depth = depth+1
                    path_vertex(depth) = u
                    path_edge(depth) = e
                    u = dest(e)
                else
                    level(u) = -1
                    if (depth == 0) exit
                    u = path_vertex(depth)
                    cursor(u) = nxt(path_edge(depth))
                    depth = depth-1
                end if
            end do
        end do
        ! The signed flow associated with each stored forward arc is exactly
        ! capacity - residual(forward). Combine opposite arcs so the returned
        ! matrix is antisymmetric, even for antiparallel or parallel inputs.
        allocate(rows(2*m),cols(2*m),values(2*m))
        used = 0
        do k = 1, m
            weight = capacities(k)-residual(2*k-1)
            if (weight == 0_int64) cycle
            used = used+1
            rows(used) = edge_src(k)
            cols(used) = graph%indices(k)
            values(used) = real(weight,dp)
            used = used+1
            rows(used) = graph%indices(k)
            cols(used) = edge_src(k)
            values(used) = -real(weight,dp)
        end do
        raw = csgraph_from_edges(n,rows(:used),cols(:used),values(:used))
        allocate(marks(n),touched(n),totals(n),out_rows(used),out_cols(used),out_values(used))
        marks = 0
        totals = 0_int64
        count_out = 0
        do u = 1, n
            ntouched = 0
            do k = raw%indptr(u),raw%indptr(u+1)-1
                v = raw%indices(k)
                if (marks(v) /= u) then
                    marks(v) = u
                    ntouched = ntouched+1
                    touched(ntouched) = v
                    totals(v) = 0_int64
                end if
                weight = int(raw%data(k),int64)
                if (weight > 0_int64) then
                    if (totals(v) > huge(weight)-weight) return
                else if (weight < 0_int64) then
                    if (totals(v) < -huge(weight)-weight) return
                end if
                totals(v) = totals(v)+weight
            end do
            do t = 1, ntouched
                v = touched(t)
                net = totals(v)
                if (net == 0_int64) cycle
                if (abs(real(net,dp)) > max_exact) return
                count_out = count_out+1
                out_rows(count_out) = u
                out_cols(count_out) = v
                out_values(count_out) = real(net,dp)
            end do
        end do
        flow = csgraph_from_edges(n,out_rows(:count_out),out_cols(:count_out), &
            out_values(:count_out),status)
    end subroutine maximum_flow


    ! Sparse Laplacian with O(n + m) extra space, including normalization and
    ! A + transpose(A). Every row contains distinct column indices; repeated
    ! edges are summed, and off-diagonal cancellations are omitted.
    subroutine laplacian_sparse(graph, lap, diagonal, status, normed, use_out_degree, symmetrized)
        type(csr_graph), intent(in) :: graph !! Valid one-based weighted square CSR graph.
        type(csr_graph), intent(out) :: lap !! Sparse Laplacian, with explicit nonzero diagonal.
        real(dp), allocatable, intent(out) :: diagonal(:) !! Degrees or normalization scales.
        integer, intent(out) :: status !! Success or invalid input/overflow.
        logical, intent(in), optional :: normed !! If true, normalize by degree square roots.
        logical, intent(in), optional :: use_out_degree !! If true, sum row degrees (default columns).
        logical, intent(in), optional :: symmetrized !! If true, use A + transpose(A).
        type(csr_graph) :: reverse
        real(dp), allocatable :: degree(:), accum(:), scale(:)
        integer, allocatable :: mark(:), touched(:)
        integer :: n, m, cap, i, j, k, t, ntouched, used, side
        real(dp) :: value
        logical :: normalize, outgoing, symmetrize

        status = csgraph_invalid_input
        allocate(diagonal(0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        m = size(graph%indices)
        if (n == huge(n)) return
        normalize = .false.
        if (present(normed)) normalize = normed
        outgoing = .false.
        if (present(use_out_degree)) outgoing = use_out_degree
        symmetrize = .false.
        if (present(symmetrized)) symmetrize = symmetrized
        if (symmetrize .and. m > (huge(m)-n)/2) return
        if (.not. symmetrize .and. m > huge(m)-n) return
        cap = n+m
        if (symmetrize) cap = n+2*m
        allocate(degree(n),scale(n),accum(n),mark(n),touched(n))
        degree = 0.0_dp
        do i = 1, n
            do k = graph%indptr(i), graph%indptr(i+1)-1
                j = graph%indices(k)
                if (i == j) cycle
                value = graph%data(k)
                if (symmetrize) then
                    degree(i) = degree(i)+value
                    degree(j) = degree(j)+value
                else if (outgoing) then
                    degree(i) = degree(i)+value
                else
                    degree(j) = degree(j)+value
                end if
            end do
        end do
        if (.not. all(ieee_is_finite(degree))) return
        if (normalize .and. any(degree < 0.0_dp)) return
        deallocate(diagonal)
        allocate(diagonal(n))
        if (normalize) then
            diagonal = 1.0_dp
            where (degree > 0.0_dp) diagonal = sqrt(degree)
        else
            diagonal = degree
        end if
        scale = diagonal
        if (.not. normalize) scale = 1.0_dp
        if (symmetrize) call reverse_graph(graph,reverse)
        lap%n = n
        allocate(lap%indptr(n+1),lap%indices(cap),lap%data(cap))
        used = 0
        mark = 0
        accum = 0.0_dp
        do i = 1, n
            lap%indptr(i) = used+1
            ntouched = 0
            do side = 1, 1+merge(1,0,symmetrize)
                if (side == 1) then
                    do k = graph%indptr(i), graph%indptr(i+1)-1
                        j = graph%indices(k)
                        if (j == i) cycle
                        if (mark(j) /= i) then
                            mark(j) = i
                            ntouched = ntouched+1
                            touched(ntouched) = j
                            accum(j) = 0.0_dp
                        end if
                        accum(j) = accum(j)-graph%data(k)
                    end do
                else
                    do k = reverse%indptr(i), reverse%indptr(i+1)-1
                        j = reverse%indices(k)
                        if (j == i) cycle
                        if (mark(j) /= i) then
                            mark(j) = i
                            ntouched = ntouched+1
                            touched(ntouched) = j
                            accum(j) = 0.0_dp
                        end if
                        accum(j) = accum(j)-reverse%data(k)
                    end do
                end if
            end do
            do t = 1, ntouched
                j = touched(t)
                value = accum(j)
                if (normalize) value = (value/scale(i))/scale(j)
                if (.not. ieee_is_finite(value)) return
                if (value == 0.0_dp) cycle
                used = used+1
                lap%indices(used) = j
                lap%data(used) = value
            end do
            if (normalize) then
                value = 0.0_dp
                if (degree(i) /= 0.0_dp) value = 1.0_dp
            else
                value = degree(i)
            end if
            if (value /= 0.0_dp) then
                used = used+1
                lap%indices(used) = i
                lap%data(used) = value
            end if
        end do
        lap%indptr(n+1) = used+1
        lap%indices = lap%indices(:used)
        lap%data = lap%data(:used)
        status = csgraph_success
    end subroutine laplacian_sparse

    ! Construct weighted rectangular biadjacency without losing explicit
    ! zero-weight edges. Duplicate edges are retained; the assignment solver
    ! naturally chooses the least expensive duplicate.
    function weighted_bipartite_from_edges(n_rows, n_cols, rows, columns, weights, status) result(graph)
        integer, intent(in) :: n_rows !! Number of left vertices (>= 0).
        integer, intent(in) :: n_cols !! Number of right vertices (>= 0).
        integer, intent(in) :: rows(:) !! One-based row indices.
        integer, intent(in) :: columns(:) !! One-based column indices.
        real(dp), intent(in) :: weights(:) !! Finite costs; zeros and negative values allowed.
        integer, intent(out), optional :: status !! Success or invalid input.
        type(csr_weighted_bipartite_graph) :: graph
        integer, allocatable :: cursor(:)
        integer :: m, i, k, pos

        if (present(status)) status = csgraph_invalid_input
        if (n_rows < 0 .or. n_cols < 0) return
        if (n_rows == huge(n_rows) .or. n_cols == huge(n_cols)) return
        m = size(rows)
        if (m == huge(m) .or. size(columns) /= m .or. size(weights) /= m) return
        if (any(rows < 1) .or. any(rows > n_rows) .or. &
            any(columns < 1) .or. any(columns > n_cols)) return
        if (.not. all(ieee_is_finite(weights))) return
        graph%n_rows = n_rows
        graph%n_cols = n_cols
        allocate(graph%indptr(n_rows+1),graph%indices(m),graph%data(m),cursor(n_rows))
        graph%indptr = 0
        graph%indptr(1) = 1
        do k = 1, m
            graph%indptr(rows(k)+1) = graph%indptr(rows(k)+1)+1
        end do
        do i = 1, n_rows
            graph%indptr(i+1) = graph%indptr(i+1)+graph%indptr(i)
        end do
        cursor = graph%indptr(:n_rows)
        do k = 1, m
            pos = cursor(rows(k))
            graph%indices(pos) = columns(k)
            graph%data(pos) = weights(k)
            cursor(rows(k)) = pos+1
        end do
        if (present(status)) status = csgraph_success
    end function weighted_bipartite_from_edges

    function weighted_bipartite_from_dense(a, status) result(graph)
        real(dp), intent(in) :: a(:,:) !! Finite rectangular costs; zero denotes no edge.
        integer, intent(out), optional :: status !! Success or invalid input.
        type(csr_weighted_bipartite_graph) :: graph
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: weights(:)
        integer :: i, j, k, m
        if (present(status)) status = csgraph_invalid_input
        if (.not. all(ieee_is_finite(a))) return
        if (count(a /= 0.0_dp,kind=int64) >= int(huge(m),int64)) return
        m = count(a /= 0.0_dp)
        allocate(rows(m),cols(m),weights(m))
        k = 0
        do i = 1, size(a,1)
            do j = 1, size(a,2)
                if (a(i,j) == 0.0_dp) cycle
                k = k+1
                rows(k) = i
                cols(k) = j
                weights(k) = a(i,j)
            end do
        end do
        graph = weighted_bipartite_from_edges(size(a,1),size(a,2),rows,cols,weights,status)
    end function weighted_bipartite_from_dense

    function weighted_bipartite_from_csgraph(square, status) result(graph)
        type(csr_graph), intent(in) :: square !! Valid CSR graph, including explicit zero arcs.
        integer, intent(out), optional :: status !! Success or invalid input.
        type(csr_weighted_bipartite_graph) :: graph
        if (present(status)) status = csgraph_invalid_input
        if (.not. valid_csgraph(square)) return
        graph%n_rows = square%n
        graph%n_cols = square%n
        graph%indptr = square%indptr
        graph%indices = square%indices
        graph%data = square%data
        if (present(status)) status = csgraph_success
    end function weighted_bipartite_from_csgraph

    pure logical function valid_weighted_bipartite_graph(graph) result(ok)
        type(csr_weighted_bipartite_graph), intent(in) :: graph !! Weighted rectangular CSR graph.
        integer :: n, m
        ok = .false.
        n = graph%n_rows
        if (n < 0 .or. graph%n_cols < 0 .or. n == huge(n)) return
        if (.not. allocated(graph%indptr) .or. .not. allocated(graph%indices) .or. &
            .not. allocated(graph%data)) return
        m = size(graph%indices)
        if (m == huge(m)) return
        if (size(graph%indptr) /= n+1 .or. size(graph%data) /= m) return
        if (graph%indptr(1) /= 1 .or. graph%indptr(n+1) /= m+1) return
        if (any(graph%indptr < 1) .or. any(graph%indptr > m+1)) return
        if (any(graph%indptr(2:) < graph%indptr(:n))) return
        if (any(graph%indices < 1) .or. any(graph%indices > graph%n_cols)) return
        if (.not. all(ieee_is_finite(graph%data))) return
        ok = .true.
    end function valid_weighted_bipartite_graph

    ! Transpose to run the augmenting-path solver with its smaller partition
    ! always on the left. This is an O(n + m) CSR conversion.
    subroutine transpose_weighted_bipartite(original, transposed)
        type(csr_weighted_bipartite_graph), intent(in) :: original !! Valid weighted CSR.
        type(csr_weighted_bipartite_graph), intent(out) :: transposed !! CSR of its transpose.
        integer, allocatable :: cursor(:)
        integer :: i, j, k, pos, m
        m = size(original%indices)
        transposed%n_rows = original%n_cols
        transposed%n_cols = original%n_rows
        allocate(transposed%indptr(transposed%n_rows+1),transposed%indices(m), &
            transposed%data(m),cursor(transposed%n_rows))
        transposed%indptr = 0
        transposed%indptr(1) = 1
        do k = 1, m
            j = original%indices(k)
            transposed%indptr(j+1) = transposed%indptr(j+1)+1
        end do
        do i = 1, transposed%n_rows
            transposed%indptr(i+1) = transposed%indptr(i+1)+transposed%indptr(i)
        end do
        cursor = transposed%indptr(:transposed%n_rows)
        do i = 1, original%n_rows
            do k = original%indptr(i),original%indptr(i+1)-1
                j = original%indices(k)
                pos = cursor(j)
                transposed%indices(pos) = i
                transposed%data(pos) = original%data(k)
                cursor(j) = pos+1
            end do
        end do
    end subroutine transpose_weighted_bipartite

    ! SciPy-compatible name using the existing native sparse solver. Keep the
    ! original entry point and argument ordering for backward compatibility.
    subroutine min_weight_full_bipartite_matching(graph, row_ind, col_ind, status, maximize, total_weight)
        implicit none
        type(csr_weighted_bipartite_graph), intent(in) :: graph !! Finite weighted rectangular CSR.
        integer, allocatable, intent(out) :: row_ind(:) !! Ascending one-based matched rows.
        integer, allocatable, intent(out) :: col_ind(:) !! Corresponding one-based matched columns.
        integer, intent(out) :: status !! Success (0), invalid (1), or no full matching (3).
        logical, intent(in), optional :: maximize !! Maximize total weight; default false.
        real(dp), intent(out), optional :: total_weight !! Sum of selected original costs; 0 on failure.
        real(dp) :: total

        call minimum_weight_full_bipartite_matching(graph,row_ind,col_ind,total,status,maximize)
        if (present(total_weight)) total_weight = total
    end subroutine min_weight_full_bipartite_matching

    ! Successive shortest augmenting paths with feasible dual potentials and a
    ! binary heap. Sparse memory O(n_rows + n_cols + m), rather than building
    ! an n_rows-by-n_cols dense cost array. Ties may differ from SciPy LAPJVsp.
    subroutine minimum_weight_full_bipartite_matching(graph, row_ind, col_ind, total_weight, status, maximize)
        type(csr_weighted_bipartite_graph), intent(in) :: graph !! Finite weighted rectangular CSR.
        integer, allocatable, intent(out) :: row_ind(:) !! Ascending one-based original row IDs.
        integer, allocatable, intent(out) :: col_ind(:) !! Corresponding one-based original columns.
        real(dp), intent(out) :: total_weight !! Sum of original selected edge costs.
        integer, intent(out) :: status !! Zero, invalid, or no full matching.
        logical, intent(in), optional :: maximize !! If true, find maximum-weight full matching.
        type(csr_weighted_bipartite_graph) :: working
        real(dp), allocatable :: row_dual(:), col_dual(:), row_dist(:), col_dist(:), &
            heap_keys(:), chosen_cost(:)
        integer, allocatable :: row_match(:), col_match(:), prev_row(:), prev_edge(:), &
            heap_vertex(:), inverse(:)
        logical, allocatable :: expanded(:), settled(:)
        integer :: nl, nr, m, start, i, j, k, heap_size, terminal, old, count_out
        real(dp) :: infinity, current, delta, sign_weight
        logical :: transposed, maximum, arithmetic_failure

        status = csgraph_invalid_input
        total_weight = 0.0_dp
        allocate(row_ind(0),col_ind(0))
        if (.not. valid_weighted_bipartite_graph(graph)) return
        nl = graph%n_rows
        nr = graph%n_cols
        if (min(nl,nr) == 0) then
            status = csgraph_success
            return
        end if
        transposed = nl > nr
        if (transposed) then
            call transpose_weighted_bipartite(graph,working)
        else
            working = graph
        end if
        nl = working%n_rows
        nr = working%n_cols
        m = size(working%indices)
        maximum = .false.
        if (present(maximize)) maximum = maximize
        sign_weight = 1.0_dp
        if (maximum) sign_weight = -1.0_dp
        infinity = ieee_value(0.0_dp,ieee_positive_inf)
        allocate(row_dual(nl),col_dual(nr),row_dist(nl),col_dist(nr), &
            row_match(nl),col_match(nr),prev_row(nr),prev_edge(nr), &
            expanded(nl),settled(nr),heap_keys(max(1,m)), &
            heap_vertex(max(1,m)),chosen_cost(nl))
        col_dual = 0.0_dp
        row_match = 0
        col_match = 0
        chosen_cost = 0.0_dp
        do i = 1, nl
            if (working%indptr(i) == working%indptr(i+1)) then
                status = csgraph_no_full_matching
                return
            end if
            row_dual(i) = infinity
            do k = working%indptr(i),working%indptr(i+1)-1
                row_dual(i) = min(row_dual(i),sign_weight*working%data(k))
            end do
        end do
        do start = 1, nl
            col_dist = infinity
            row_dist = infinity
            expanded = .false.
            settled = .false.
            prev_row = 0
            prev_edge = 0
            heap_size = 0
            arithmetic_failure = .false.
            expanded(start) = .true.
            row_dist(start) = 0.0_dp
            call relax_assignment_row(start,0.0_dp)
            if (arithmetic_failure) return
            terminal = 0
            delta = 0.0_dp
            do while (heap_size > 0)
                call heap_pop(heap_keys,heap_vertex,heap_size,current,j)
                if (settled(j) .or. current /= col_dist(j)) cycle
                settled(j) = .true.
                if (col_match(j) == 0) then
                    terminal = j
                    delta = current
                    exit
                end if
                i = col_match(j)
                expanded(i) = .true.
                row_dist(i) = current
                call relax_assignment_row(i,current)
                if (arithmetic_failure) return
            end do
            if (terminal == 0) then
                status = csgraph_no_full_matching
                return
            end if
            ! Potential update makes all edges along the discovered shortest
            ! augmenting path tight and retains nonnegative reduced costs.
            do i = 1, nl
                if (expanded(i)) row_dual(i) = row_dual(i)+(delta-row_dist(i))
            end do
            do j = 1, nr
                if (settled(j)) col_dual(j) = col_dual(j)-(delta-col_dist(j))
            end do
            if (.not. all(ieee_is_finite(row_dual)) .or. &
                .not. all(ieee_is_finite(col_dual))) return
            j = terminal
            do
                i = prev_row(j)
                old = row_match(i)
                row_match(i) = j
                col_match(j) = i
                chosen_cost(i) = working%data(prev_edge(j))
                if (old == 0) exit
                j = old
            end do
        end do
        total_weight = sum(chosen_cost)
        if (.not. ieee_is_finite(total_weight)) then
            total_weight = 0.0_dp
            return
        end if
        deallocate(row_ind,col_ind)
        allocate(row_ind(nl),col_ind(nl))
        if (transposed) then
            ! Transposed row i denotes original column i. Enumerate the
            ! selected original rows to satisfy SciPy's sorted-row contract.
            allocate(inverse(graph%n_rows))
            inverse = 0
            do i = 1, nl
                inverse(row_match(i)) = i
            end do
            count_out = 0
            do j = 1, graph%n_rows
                if (inverse(j) == 0) cycle
                count_out = count_out+1
                row_ind(count_out) = j
                col_ind(count_out) = inverse(j)
            end do
        else
            do i = 1, nl
                row_ind(i) = i
                col_ind(i) = row_match(i)
            end do
        end if
        status = csgraph_success
    contains
        subroutine relax_assignment_row(row, distance)
            integer, intent(in) :: row !! Expanded one-based row.
            real(dp), intent(in) :: distance !! Path length to this row.
            integer :: edge, column
            real(dp) :: red, trial
            do edge = working%indptr(row),working%indptr(row+1)-1
                column = working%indices(edge)
                if (settled(column)) cycle
                red = (sign_weight*working%data(edge)-row_dual(row))-col_dual(column)
                if (.not. ieee_is_finite(red)) then
                    arithmetic_failure = .true.
                    return
                end if
                red = max(0.0_dp,red)
                trial = distance+red
                if (.not. ieee_is_finite(trial)) then
                    arithmetic_failure = .true.
                    return
                end if
                if (trial >= col_dist(column)) cycle
                col_dist(column) = trial
                prev_row(column) = row
                prev_edge(column) = edge
                call heap_push(heap_keys,heap_vertex,heap_size,trial,column)
            end do
        end subroutine relax_assignment_row
    end subroutine minimum_weight_full_bipartite_matching


    ! Yen's K shortest *simple* source-to-target paths.  Rows of predecessors
    ! describe only vertices on each path; all other entries are zero.
    ! The output may contain fewer than k paths when routes are exhausted.
    subroutine yen(graph, source, sink, k, distances, predecessors, status, directed, unweighted)
        type(csr_graph), intent(in) :: graph !! Finite, one-based CSR graph.
        integer, intent(in) :: source !! One-based starting vertex.
        integer, intent(in) :: sink !! One-based ending vertex.
        integer, intent(in) :: k !! Maximum number of simple paths; may be zero.
        real(dp), allocatable, intent(out) :: distances(:) !! Ascending path costs (size <= k).
        integer, allocatable, intent(out) :: predecessors(:,:) !! (number of paths, graph%n), zero if absent.
        integer, intent(out) :: status !! 0 success, 1 invalid input, 2 negative cycle.
        logical, intent(in), optional :: directed !! Default true; false follows either arc orientation.
        logical, intent(in), optional :: unweighted !! Default false; true assigns every edge unit cost.
        type(csr_graph) :: reverse
        type(yen_path), allocatable :: accepted(:), candidates(:), grown(:)
        type(yen_path) :: spur_path, candidate
        real(dp), allocatable :: potentials(:)
        logical, allocatable :: blocked(:), banned_next(:)
        logical :: is_directed, is_unweighted, changed, duplicate
        integer :: n, m, iteration, npaths, heap_count, spur, path_length, i, j, v, code

        status = csgraph_invalid_input
        allocate(distances(0), predecessors(0,0))
        if (.not. valid_csgraph(graph)) return
        n = graph%n
        m = size(graph%indices)
        if (source < 1 .or. source > n .or. sink < 1 .or. sink > n .or. k < 0) return
        ! Each stored arc can improve a Dijkstra tentative distance at most once.
        if (2_int64 * int(m,int64) + 1_int64 > int(huge(m),int64)) return
        is_directed = .true.
        if (present(directed)) is_directed = directed
        is_unweighted = .false.
        if (present(unweighted)) is_unweighted = unweighted
        deallocate(predecessors)
        allocate(predecessors(0,n))
        if (k == 0) then
            status = csgraph_success
            return
        end if
        if (.not. is_directed) call reverse_graph(graph,reverse)

        ! Bellman-Ford with a super-source supplies globally feasible Johnson
        ! potentials.  This also detects disconnected negative cycles.  With
        ! unweighted edges no reweighting or cycle check is required.
        allocate(potentials(n))
        potentials = 0.0_dp
        if (.not. is_unweighted .and. any(graph%data < 0.0_dp)) then
            do iteration = 1, n-1
                changed = .false.
                call relax_potentials(graph,potentials,changed,.true.)
                if (.not. is_directed) call relax_potentials(reverse,potentials,changed,.true.)
                if (.not. changed) exit
            end do
            changed = .false.
            call relax_potentials(graph,potentials,changed,.false.)
            if (.not. is_directed) call relax_potentials(reverse,potentials,changed,.false.)
            if (changed) then
                status = csgraph_negative_cycle
                return
            end if
        end if

        ! Grow accepted routes only when found: K can be very large even
        ! for a graph with no source-to-target route.
        allocate(accepted(1), candidates(0), blocked(n), banned_next(n))
        blocked = .false.
        banned_next = .false.
        call yen_spur_path(graph,reverse,is_directed,is_unweighted,potentials, &
            source,sink,blocked,banned_next,accepted(1),code)
        if (code /= csgraph_success) then
            status = code
            return
        end if
        npaths = 0
        if (allocated(accepted(1)%vertices)) npaths = 1
        heap_count = 0

        do while (npaths > 0 .and. npaths < k)
            path_length = size(accepted(npaths)%vertices)
            do spur = 1, path_length-1
                ! Vertices before the spur in the previous best route are
                ! excluded to guarantee a simple concatenated route.
                blocked = .false.
                blocked(accepted(npaths)%vertices(:spur-1)) = .true.
                banned_next = .false.
                do i = 1, npaths
                    if (size(accepted(i)%vertices) <= spur) cycle
                    if (.not. all(accepted(i)%vertices(:spur) == &
                        accepted(npaths)%vertices(:spur))) cycle
                    banned_next(accepted(i)%vertices(spur+1)) = .true.
                end do
                call yen_spur_path(graph,reverse,is_directed,is_unweighted,potentials, &
                    accepted(npaths)%vertices(spur),sink,blocked,banned_next,spur_path,code)
                if (code /= csgraph_success) then
                    status = code
                    return
                end if
                if (.not. allocated(spur_path%vertices)) cycle
                candidate%vertices = [accepted(npaths)%vertices(:spur-1),spur_path%vertices]
                candidate%edge_weights = [accepted(npaths)%edge_weights(:spur-1),spur_path%edge_weights]
                candidate%cost = sum(candidate%edge_weights)
                if (.not. ieee_is_finite(candidate%cost)) then
                    status = csgraph_invalid_input
                    return
                end if
                duplicate = .false.
                do j = 1, npaths
                    if (yen_same_path(candidate,accepted(j))) then
                        duplicate = .true.
                        exit
                    end if
                end do
                if (.not. duplicate) then
                    do j = 1, heap_count
                        if (yen_same_path(candidate,candidates(j))) then
                            duplicate = .true.
                            exit
                        end if
                    end do
                end if
                if (.not. duplicate) call yen_candidate_push(candidates,heap_count,candidate)
            end do
            if (heap_count == 0) exit
            if (npaths == size(accepted)) then
                allocate(grown(int(min(int(k,int64),2_int64*int(size(accepted),int64)))))
                grown(:npaths) = accepted(:npaths)
                call move_alloc(grown,accepted)
            end if
            npaths = npaths+1
            call yen_candidate_pop(candidates,heap_count,accepted(npaths))
        end do

        deallocate(distances,predecessors)
        allocate(distances(npaths),predecessors(npaths,n))
        predecessors = 0
        do i = 1, npaths
            distances(i) = accepted(i)%cost
            do j = 2, size(accepted(i)%vertices)
                v = accepted(i)%vertices(j)
                predecessors(i,v) = accepted(i)%vertices(j-1)
            end do
        end do
        status = csgraph_success
    end subroutine yen

    ! A shortest allowed spur, using a binary heap with reweighted arcs.
    ! The supplied banned_next mask applies only to outgoing arcs of source.
    subroutine yen_spur_path(graph,reverse,is_directed,is_unweighted,potentials, &
        source,sink,blocked,banned_next,path,status)
        type(csr_graph), target, intent(in) :: graph !! Forward adjacency.
        type(csr_graph), target, intent(in) :: reverse !! Reversed adjacency if needed.
        logical, intent(in) :: is_directed !! Whether to ignore reversed arcs.
        logical, intent(in) :: is_unweighted !! Treat each traversed arc as cost one.
        real(dp), intent(in) :: potentials(:) !! Johnson feasible potentials.
        integer, intent(in) :: source !! Spur vertex.
        integer, intent(in) :: sink !! Destination vertex.
        logical, intent(in) :: blocked(:) !! Excluded root-prefix vertices.
        logical, intent(in) :: banned_next(:) !! Excluded next vertices at spur.
        type(yen_path), intent(out) :: path !! Empty when there is no route.
        integer, intent(out) :: status !! Success or invalid numerical input.
        type(csr_graph), pointer :: adjacency
        real(dp), allocatable :: dist(:), prev_weight(:), heap_keys(:)
        integer, allocatable :: parent(:), heap_vertices(:), reversed_path(:)
        logical, allocatable :: finalized(:)
        real(dp) :: current, weight, reduced, trial, inf
        integer :: n, m, heap_count, u, v, edge, direction, n_directions, length, j

        status = csgraph_success
        n = graph%n
        m = size(graph%data)
        if (blocked(source) .or. blocked(sink)) return
        if (source == sink) then
            path%vertices = [source]
            allocate(path%edge_weights(0))
            path%cost = 0.0_dp
            return
        end if
        allocate(dist(n),prev_weight(n),parent(n),finalized(n), &
            heap_keys(2*m+1),heap_vertices(2*m+1))
        inf = ieee_value(0.0_dp,ieee_positive_inf)
        dist = inf
        parent = 0
        prev_weight = 0.0_dp
        finalized = .false.
        dist(source) = 0.0_dp
        heap_count = 0
        call heap_push(heap_keys,heap_vertices,heap_count,0.0_dp,source)
        n_directions = 1
        if (.not. is_directed) n_directions = 2
        do while (heap_count > 0)
            call heap_pop(heap_keys,heap_vertices,heap_count,current,u)
            if (finalized(u) .or. current > dist(u)) cycle
            finalized(u) = .true.
            if (u == sink) exit
            do direction = 1, n_directions
                if (direction == 1) then
                    adjacency => graph
                else
                    adjacency => reverse
                end if
                do edge = adjacency%indptr(u), adjacency%indptr(u+1)-1
                    v = adjacency%indices(edge)
                    if (blocked(v) .or. finalized(v)) cycle
                    if (u == source .and. banned_next(v)) cycle
                    weight = adjacency%data(edge)
                    if (is_unweighted) weight = 1.0_dp
                    reduced = weight
                    if (.not. is_unweighted) reduced = (weight+potentials(u))-potentials(v)
                    if (.not. ieee_is_finite(reduced)) then
                        status = csgraph_invalid_input
                        return
                    end if
                    ! Permit tiny negative reduced costs from subtraction
                    ! rounding, but do not hide numerically infeasible edges.
                    if (reduced < -64.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
                        abs(weight),abs(potentials(u)),abs(potentials(v)))) then
                        status = csgraph_invalid_input
                        return
                    end if
                    reduced = max(0.0_dp,reduced)
                    trial = current+reduced
                    if (.not. ieee_is_finite(trial)) then
                        status = csgraph_invalid_input
                        return
                    end if
                    if (trial >= dist(v)) cycle
                    dist(v) = trial
                    parent(v) = u
                    prev_weight(v) = weight
                    call heap_push(heap_keys,heap_vertices,heap_count,trial,v)
                end do
            end do
        end do
        if (.not. finalized(sink)) return
        allocate(reversed_path(n))
        length = 1
        reversed_path(1) = sink
        u = sink
        do while (u /= source)
            v = parent(u)
            if (v <= 0 .or. length == n) then
                status = csgraph_invalid_input
                return
            end if
            length = length+1
            reversed_path(length) = v
            u = v
        end do
        allocate(path%vertices(length),path%edge_weights(length-1))
        path%vertices = reversed_path(length:1:-1)
        do j = 1, length-1
            path%edge_weights(j) = prev_weight(path%vertices(j+1))
        end do
        path%cost = sum(path%edge_weights)
        if (.not. ieee_is_finite(path%cost)) status = csgraph_invalid_input
    end subroutine yen_spur_path

    pure logical function yen_same_path(a,b) result(same)
        type(yen_path), intent(in) :: a !! Candidate path.
        type(yen_path), intent(in) :: b !! Candidate or accepted path.
        same = .false.
        if (.not. allocated(a%vertices) .or. .not. allocated(b%vertices)) return
        if (size(a%vertices) /= size(b%vertices)) return
        same = all(a%vertices == b%vertices)
    end function yen_same_path

    ! Deterministic lexicographic vertex ordering for equal-cost alternatives.
    pure logical function yen_path_less(a,b) result(less)
        type(yen_path), intent(in) :: a !! First simple path.
        type(yen_path), intent(in) :: b !! Second simple path.
        integer :: i
        less = a%cost < b%cost
        if (a%cost < b%cost .or. a%cost > b%cost) return
        do i = 1, min(size(a%vertices),size(b%vertices))
            if (a%vertices(i) == b%vertices(i)) cycle
            less = a%vertices(i) < b%vertices(i)
            return
        end do
        less = size(a%vertices) < size(b%vertices)
    end function yen_path_less

    subroutine yen_candidate_push(heap,count,item)
        type(yen_path), allocatable, intent(inout) :: heap(:) !! Candidate binary min heap.
        integer, intent(inout) :: count !! Number of active candidate paths.
        type(yen_path), intent(in) :: item !! Newly discovered path.
        type(yen_path), allocatable :: grown(:)
        integer :: i, parent, capacity
        if (count == size(heap)) then
            capacity = int(min(int(huge(count),int64), &
                max(8_int64,2_int64*int(size(heap),int64))))
            allocate(grown(capacity))
            grown(:count) = heap(:count)
            call move_alloc(grown,heap)
        end if
        count = count+1
        i = count
        do while (i > 1)
            parent = i/2
            if (.not. yen_path_less(item,heap(parent))) exit
            heap(i) = heap(parent)
            i = parent
        end do
        heap(i) = item
    end subroutine yen_candidate_push

    subroutine yen_candidate_pop(heap,count,item)
        type(yen_path), intent(inout) :: heap(:) !! Candidate binary min heap.
        integer, intent(inout) :: count !! Number of active candidate paths.
        type(yen_path), intent(out) :: item !! Cheapest candidate.
        type(yen_path) :: last
        integer :: i, child
        item = heap(1)
        last = heap(count)
        count = count-1
        if (count == 0) return
        i = 1
        do while (i <= count/2)
            child = 2*i
            if (child < count) then
                if (yen_path_less(heap(child+1),heap(child))) child = child+1
            end if
            if (.not. yen_path_less(heap(child),last)) exit
            heap(i) = heap(child)
            i = child
        end do
        heap(i) = last
    end subroutine yen_candidate_pop

end module scifort_sparse_csgraph
