! SPDX-License-Identifier: MIT
program test_sparse_csgraph_weighted
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csr_weighted_bipartite_graph, &
        weighted_bipartite_from_edges, weighted_bipartite_from_dense, &
        weighted_bipartite_from_csgraph, valid_weighted_bipartite_graph, &
        minimum_weight_full_bipartite_matching, laplacian, laplacian_sparse, &
        csgraph_from_edges, csgraph_to_dense, csgraph_success, &
        csgraph_invalid_input, csgraph_no_full_matching
    implicit none
    type(csr_graph) :: graph, sparse_lap
    type(csr_weighted_bipartite_graph) :: bipartite, invalid_bi
    real(dp), allocatable :: degrees(:)
    integer, allocatable :: ri(:), ci(:)
    integer :: status, n, i
    real(dp) :: total

    bipartite = weighted_bipartite_from_dense(reshape([ &
        3.0_dp,4.0_dp,10.0_dp,3.0_dp,3.0_dp,1.0_dp, &
        6.0_dp,5.0_dp,8.0_dp],[3,3]),status)
    call require(status == csgraph_success .and. valid_weighted_bipartite_graph(bipartite), &
        'dense weighted constructor')
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status)
    call require(status == csgraph_success .and. all(ri == [1,2,3]),'square weighted rows')
    call require(all(ci == [1,3,2]) .and. total == 9.0_dp,'square minimum cost = 9')
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status,maximize=.true.)
    call require(status == csgraph_success .and. total == 19.0_dp,'square maximum weight = 19')

    bipartite = weighted_bipartite_from_edges(2,3,[1,1,2,2],[1,2,1,3], &
        [0.0_dp,4.0_dp,3.0_dp,-1.0_dp],status)
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status)
    call require(status == csgraph_success .and. total == -1.0_dp,'stored zero and negative weights')
    call require(all(ri == [1,2]) .and. all(ci == [1,3]),'rectangular rows <= cols')
    bipartite = weighted_bipartite_from_edges(3,2,[1,1,2,2,3,3], &
        [1,2,1,2,1,2],[4.0_dp,1.0_dp,8.0_dp,7.0_dp,2.0_dp,6.0_dp],status)
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status)
    call require(status == csgraph_success .and. total == 3.0_dp,'transpose minimum matching')
    call require(all(ri == [1,3]) .and. all(ci == [2,1]),'sorted original row indices')
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status,maximize=.true.)
    call require(status == csgraph_success .and. total == 14.0_dp,'transpose maximum matching')

    bipartite = weighted_bipartite_from_edges(2,2,[1,1,2,2],[1,1,1,2], &
        [20.0_dp,1.0_dp,2.0_dp,3.0_dp],status)
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status)
    call require(status == csgraph_success .and. total == 4.0_dp,'parallel edges choose cheaper')
    bipartite = weighted_bipartite_from_edges(2,2,[1,2],[1,1],[1.0_dp,2.0_dp],status)
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status)
    call require(status == csgraph_no_full_matching .and. size(ri) == 0, 'infeasible assignment status')
    bipartite = weighted_bipartite_from_edges(0,3,[integer::],[integer::],[real(dp)::],status)
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status)
    call require(status == csgraph_success .and. size(ri) == 0, 'empty left partition')
    bipartite = weighted_bipartite_from_edges(3,0,[integer::],[integer::],[real(dp)::],status)
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status)
    call require(status == csgraph_success .and. size(ri) == 0, 'empty right partition')
    bipartite = weighted_bipartite_from_edges(2,3,[0],[2],[1.0_dp],status)
    call require(status == csgraph_invalid_input, 'invalid weighted row')
    bipartite = weighted_bipartite_from_edges(2,3,[1],[2], &
        [ieee_value(0.0_dp,ieee_quiet_nan)],status)
    call require(status == csgraph_invalid_input, 'NaN weighted cost rejected')
    call minimum_weight_full_bipartite_matching(invalid_bi,ri,ci,total,status)
    call require(status == csgraph_invalid_input, 'unallocated weighted CSR rejected')
    n = 2000
    bipartite = weighted_bipartite_from_edges(n,n,[(i,i=1,n)],[(i,i=1,n)], &
        [(1.0_dp,i=1,n)],status)
    call minimum_weight_full_bipartite_matching(bipartite,ri,ci,total,status)
    call require(status == csgraph_success .and. total == real(n,dp), &
        'large sparse full weighted matching')
    call require(all(ri == ci) .and. size(ri) == n,'large diagonal weighted matching')
    graph = csgraph_from_edges(2,[1,2],[2,1],[0.0_dp,3.0_dp],status)
    bipartite = weighted_bipartite_from_csgraph(graph,status)
    call require(status == csgraph_success .and. all(bipartite%data == [0.0_dp,3.0_dp]), &
        'square weighted conversion preserves zeros')

    graph = csgraph_from_edges(3,[1,1,2,3,1,2],[2,2,3,1,1,2], &
        [3.0_dp,2.0_dp,4.0_dp,-1.0_dp,30.0_dp,20.0_dp],status)
    call check_laplacian(graph)
    graph = csgraph_from_edges(3,[1,1,1,2],[2,2,3,1], &
        [1.0_dp,-1.0_dp,2.0_dp,5.0_dp],status)
    call check_laplacian(graph)
    graph = csgraph_from_edges(0,[integer::],[integer::],[real(dp)::],status)
    call check_laplacian(graph)
    graph = csgraph_from_edges(5,[1,2],[2,3],[2.0_dp,4.0_dp],status)
    call check_laplacian(graph)

    ! This would require ~800 MB if the implementation were accidentally dense.
    n = 12000
    graph = csgraph_from_edges(n,[(i,i=1,n-1)],[(i,i=2,n)], &
        [(1.0_dp,i=1,n-1)])
    call laplacian_sparse(graph,sparse_lap,degrees,status,use_out_degree=.true.)
    call require(status == csgraph_success .and. sparse_lap%n == n, 'large sparse Laplacian')
    call require(size(sparse_lap%data) == 2*(n-1), 'large sparse storage O(m+n)')
    call require(degrees(1) == 1.0_dp .and. degrees(n) == 0.0_dp, 'large sparse degree')
    print '(a)', 'test_sparse_csgraph_weighted: PASS'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! Condition to check.
        character(len=*), intent(in) :: message !! Failure diagnostic.
        if (.not. ok) then
            print '(a)', 'FAILED: '//trim(message)
            error stop 1
        end if
    end subroutine require

    subroutine check_laplacian(input)
        type(csr_graph), intent(in) :: input !! Valid square weighted CSR.
        integer :: a, b, c, result_code
        real(dp), allocatable :: dense_reference(:,:), expected_degree(:), got_degree(:), sparse_dense(:,:)
        type(csr_graph) :: result
        do a = 0, 1
            do b = 0, 1
                do c = 0, 1
                    call laplacian(input,dense_reference,expected_degree,result_code, &
                        normed=a==1,use_out_degree=b==1,symmetrized=c==1)
                    call laplacian_sparse(input,result,got_degree,status, &
                        normed=a==1,use_out_degree=b==1,symmetrized=c==1)
                    call require(status == result_code, 'sparse Laplacian status parity')
                    if (status /= csgraph_success) cycle
                    call require(all(abs(got_degree-expected_degree) < 1.0e-11_dp), &
                        'sparse Laplacian degree parity')
                    sparse_dense = csgraph_to_dense(result,result_code)
                    call require(result_code == csgraph_success, 'valid sparse output CSR')
                    call require(all(abs(sparse_dense-dense_reference) < 1.0e-11_dp), &
                        'sparse and dense Laplacian numerical parity')
                    call require(size(result%indices) <= 2*size(input%indices)+input%n, &
                        'sparse output bounded storage')
                end do
            end do
        end do
    end subroutine check_laplacian
end program test_sparse_csgraph_weighted
