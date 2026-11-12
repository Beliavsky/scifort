! SPDX-License-Identifier: MIT
program test_sparse_csgraph_advanced
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csr_bipartite_graph, &
        csgraph_from_edges, csgraph_from_dense, csgraph_to_dense, &
        csgraph_success, csgraph_invalid_input, &
        reverse_cuthill_mckee, maximum_flow, maximum_bipartite_matching, &
        structural_rank, valid_bipartite_graph, bipartite_from_edges, &
        bipartite_from_dense, bipartite_from_csgraph
    implicit none
    type(csr_graph) :: graph, flow, invalid_graph
    type(csr_bipartite_graph) :: bi, invalid_bi
    real(dp), allocatable :: net(:,:)
    real(dp) :: a(4,4)
    integer, allocatable :: perm(:), matched(:)
    integer :: status, rank
    integer(int64) :: value

    graph = csgraph_from_edges(6, [3,6,1,5,2], [6,1,5,2,4], &
        [1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp],status)
    call require(status == csgraph_success,'path input constructed')
    call reverse_cuthill_mckee(graph,perm,status)
    call require(status == csgraph_success,'reverse Cuthill--McKee on directed pattern')
    call check_permutation(perm,6)
    call require(bandwidth(graph,perm) == 1,'reordered path has optimal bandwidth 1')
    call reverse_cuthill_mckee(graph,perm,status,symmetric_mode=.true.)
    call require(status == csgraph_success,'directed symmetric-mode call')
    call check_permutation(perm,6)
    graph = csgraph_from_edges(3, [1,1,2,3], [2,2,3,3], &
        [0.0_dp,2.0_dp,1.0_dp,1.0_dp],status)
    call reverse_cuthill_mckee(graph,perm,status)
    call check_permutation(perm,3)
    call require(bandwidth(graph,perm)==1,'stored zeros and repeated neighbors')
    graph = csgraph_from_edges(0,[integer::],[integer::],[real(dp)::],status)
    call reverse_cuthill_mckee(graph,perm,status)
    call require(status==csgraph_success .and. size(perm)==0,'empty ordering')
    call reverse_cuthill_mckee(invalid_graph,perm,status)
    call require(status==csgraph_invalid_input .and. size(perm)==0,'invalid graph ordering')

    bi = bipartite_from_edges(2,3,[1,2,2],[3,1,2],status)
    call require(status==csgraph_success .and. valid_bipartite_graph(bi),'rectangular input')
    call maximum_bipartite_matching(bi,matched,status,perm_type='row')
    call require(status==csgraph_success .and. size(matched)==3,'matching rows-to-columns')
    call require(count(matched/=0)==2,'matching cardinality two')
    call check_matching(bi,matched,'row')
    call maximum_bipartite_matching(bi,matched,status,perm_type='column')
    call require(status==csgraph_success .and. size(matched)==2,'matching columns-to-rows')
    call check_matching(bi,matched,'column')
    call structural_rank(bi,rank,status)
    call require(status==csgraph_success .and. rank==2,'rectangular structural rank')
    bi = bipartite_from_edges(4,3,[1,2,3,4,4,4],[1,1,2,3,3,3],status)
    call maximum_bipartite_matching(bi,matched,status,'column')
    call check_matching(bi,matched,'column')
    call require(count(matched/=0)==3,'duplicates do not alter rank')
    bi = bipartite_from_edges(0,4,[integer::],[integer::],status)
    call maximum_bipartite_matching(bi,matched,status,'row')
    call require(size(matched)==4 .and. all(matched==0),'empty left partition')
    call maximum_bipartite_matching(bi,matched,status,'column')
    call require(size(matched)==0,'empty left matching column mode')
    bi = bipartite_from_edges(3,0,[integer::],[integer::],status)
    call maximum_bipartite_matching(bi,matched,status,'column')
    call require(size(matched)==3 .and. all(matched==0),'empty right partition')
    bi = bipartite_from_edges(2,2,[1,1,2],[1,2,1],status)
    call maximum_bipartite_matching(bi,matched,status,'row')
    call require(count(matched/=0)==2,'augmenting path finds perfect matching')
    call maximum_bipartite_matching(bi,matched,status,'invalid')
    call require(status==csgraph_invalid_input .and. size(matched)==0,'invalid permutation type')
    call maximum_bipartite_matching(invalid_bi,matched,status)
    call require(status==csgraph_invalid_input .and. size(matched)==0,'invalid structural CSR')
    bi = bipartite_from_edges(3,2,[0],[1],status)
    call require(status==csgraph_invalid_input,'invalid left vertex')
    a = 0.0_dp
    a(1,2) = 7.0_dp
    a(3,4) = -1.0_dp
    bi = bipartite_from_dense(a,status)
    call structural_rank(bi,rank,status)
    call require(rank == 2 .and. status==csgraph_success,'dense bipartite conversion')
    graph = csgraph_from_dense(a,status)
    bi = bipartite_from_csgraph(graph,status)
    call structural_rank(bi,rank,status)
    call require(rank==2,'square CSR structural conversion')

    graph = csgraph_from_edges(6, &
        [1,1,2,3,2,4,3,5,4,5], &
        [2,3,3,2,4,3,5,4,6,6], &
        [16.0_dp,13.0_dp,10.0_dp,4.0_dp,12.0_dp,9.0_dp,14.0_dp,7.0_dp,20.0_dp,4.0_dp],status)
    call maximum_flow(graph,1,6,value,flow,status)
    call require(status==csgraph_success .and. value==23_int64,'classic maximum flow 23')
    call check_flow(graph,1,6,value,flow)
    graph = csgraph_from_edges(3,[1,1,2,1],[2,2,3,3], &
        [3.0_dp,4.0_dp,5.0_dp,1.0_dp],status)
    call maximum_flow(graph,1,3,value,flow,status)
    call require(status==csgraph_success .and. value==6_int64,'parallel-arc maximum flow')
    call check_flow(graph,1,3,value,flow)
    net = csgraph_to_dense(flow,status)
    call require(net(1,2)==5.0_dp .and. net(2,1)==-5.0_dp,'net flow combines parallel arcs')
    graph = csgraph_from_edges(4,[1,2,2,3,3,2,4],[2,1,3,2,4,4,4], &
        [4.0_dp,5.0_dp,3.0_dp,7.0_dp,2.0_dp,0.0_dp,99.0_dp],status)
    call maximum_flow(graph,1,4,value,flow,status)
    call require(status==csgraph_success .and. value==2_int64,'opposite arcs and loops')
    call check_flow(graph,1,4,value,flow)
    graph = csgraph_from_edges(3,[1],[2],[0.0_dp],status)
    call maximum_flow(graph,1,3,value,flow,status)
    call require(status==csgraph_success .and. value==0_int64,'unreachable zero-flow sink')
    call check_flow(graph,1,3,value,flow)
    graph = csgraph_from_edges(2,[1],[2],[-1.0_dp],status)
    call maximum_flow(graph,1,2,value,flow,status)
    call require(status==csgraph_invalid_input,'negative capacity rejected')
    graph = csgraph_from_edges(2,[1],[2],[1.5_dp],status)
    call maximum_flow(graph,1,2,value,flow,status)
    call require(status==csgraph_invalid_input,'fractional capacity rejected')
    graph = csgraph_from_edges(2,[1],[2],[2.0_dp],status)
    call maximum_flow(graph,1,1,value,flow,status)
    call require(status==csgraph_invalid_input,'identical source and sink rejected')
    call maximum_flow(graph,0,2,value,flow,status)
    call require(status==csgraph_invalid_input,'invalid source rejected')
    call maximum_flow(invalid_graph,1,2,value,flow,status)
    call require(status==csgraph_invalid_input,'invalid CSR maximum flow')
    graph = csgraph_from_edges(2,[1],[2],[9007199254740991.0_dp],status)
    call maximum_flow(graph,1,2,value,flow,status)
    call require(status==csgraph_success .and. &
        value==9007199254740991_int64,'largest exactly representable capacity')
    call check_flow(graph,1,2,value,flow)
    graph = csgraph_from_edges(2,[1,1],[2,2], &
        [9007199254740991.0_dp,9007199254740991.0_dp],status)
    call maximum_flow(graph,1,2,value,flow,status)
    call require(status==csgraph_invalid_input,'unrepresentable merged net flow rejected')
    print '(a)', 'test_sparse_csgraph_advanced: PASS'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print '(a)', 'FAILED: '//trim(message)
            error stop 1
        end if
    end subroutine require

    subroutine check_permutation(p,n)
        integer, intent(in) :: p(:), n
        integer :: v
        call require(size(p)==n,'permutation length')
        do v=1,n
            call require(count(p==v)==1,'permutation covers all vertices')
        end do
    end subroutine check_permutation

    integer function bandwidth(g,p) result(width)
        type(csr_graph), intent(in) :: g
        integer, intent(in) :: p(:)
        integer, allocatable :: inverse(:)
        integer :: x, y, e
        allocate(inverse(g%n))
        do x=1,g%n
            inverse(p(x))=x
        end do
        width=0
        do x=1,g%n
            do e=g%indptr(x),g%indptr(x+1)-1
                y=g%indices(e)
                width=max(width,abs(inverse(x)-inverse(y)))
            end do
        end do
    end function bandwidth

    subroutine check_matching(g,match,rule)
        type(csr_bipartite_graph), intent(in) :: g
        integer, intent(in) :: match(:)
        character(len=*), intent(in) :: rule
        integer :: r, c, k
        if (rule=='row') then
            do c=1,g%n_cols
                r=match(c)
                if (r==0) cycle
                call require(r>=1 .and. r<=g%n_rows,'row match range')
                call require(count(match==r)==1,'matching rows unique')
                call require(any(g%indices(g%indptr(r):g%indptr(r+1)-1)==c), &
                    'matched edge exists')
            end do
        else
            do r=1,g%n_rows
                c=match(r)
                if (c==0) cycle
                call require(c>=1 .and. c<=g%n_cols,'column match range')
                call require(count(match==c)==1,'matching columns unique')
                call require(any(g%indices(g%indptr(r):g%indptr(r+1)-1)==c), &
                    'matched edge exists')
            end do
        end if
    end subroutine check_matching

    subroutine check_flow(g,s,t,f,actual)
        type(csr_graph), intent(in) :: g,actual
        integer, intent(in) :: s,t
        integer(int64), intent(in) :: f
        real(dp), allocatable :: dense(:,:),caps(:,:), netflow(:)
        integer :: u,k
        dense=csgraph_to_dense(actual)
        allocate(caps(g%n,g%n),netflow(g%n))
        caps=0.0_dp
        do u=1,g%n
            do k=g%indptr(u),g%indptr(u+1)-1
                caps(u,g%indices(k))=caps(u,g%indices(k))+g%data(k)
            end do
        end do
        call require(all(abs(dense+transpose(dense))<1.0e-12_dp),'flow is antisymmetric')
        call require(all(dense<=caps+1.0e-12_dp),'flow respects directed capacities')
        netflow=sum(dense,dim=2)
        call require(abs(netflow(s)-real(f,dp))<1.0e-12_dp,'source flow balance')
        call require(abs(netflow(t)+real(f,dp))<1.0e-12_dp,'sink flow balance')
        do u=1,g%n
            if (u==s .or. u==t) cycle
            call require(abs(netflow(u))<1.0e-12_dp,'interior conservation')
        end do
    end subroutine check_flow
end program test_sparse_csgraph_advanced
