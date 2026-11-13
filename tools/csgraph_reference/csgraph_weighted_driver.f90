! SPDX-License-Identifier: MIT
! Line-oriented reference protocol: W nr nc m maximize | L n m normed outgoing symmetrized.
program csgraph_weighted_driver
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csr_weighted_bipartite_graph, &
        weighted_bipartite_from_edges, min_weight_full_bipartite_matching, &
        csgraph_from_edges, laplacian_sparse, csgraph_to_dense, csgraph_success
    implicit none
    type(csr_graph) :: graph, sparse_lap
    type(csr_weighted_bipartite_graph) :: weighted
    integer, allocatable :: rows(:), cols(:), ri(:), ci(:)
    real(dp), allocatable :: values(:), diag(:), dense(:,:)
    real(dp) :: cost
    integer :: trials, t, nr, nc, m, i, j, status, a, b, c
    character(len=1) :: kind
    read(*,*) trials
    do t = 1, trials
        read(*,*) kind, nr, nc, m, a
        allocate(rows(m),cols(m),values(m))
        do i = 1, m
            read(*,*) rows(i),cols(i),values(i)
        end do
        if (kind == 'W') then
            weighted = weighted_bipartite_from_edges(nr,nc,rows,cols,values,status)
            if (status == csgraph_success) then
                call min_weight_full_bipartite_matching(weighted,ri,ci,status,a==1,cost)
            end if
            write(*,'(i0,1x,es25.16,1x,i0)') status,cost,merge(size(ri),0,status==csgraph_success)
            if (status == csgraph_success) then
                do i = 1, size(ri)
                    write(*,'(i0,1x,i0)') ri(i),ci(i)
                end do
            end if
        else if (kind == 'L') then
            read(*,*) b, c
            graph = csgraph_from_edges(nr,rows,cols,values,status)
            if (status == csgraph_success) then
                call laplacian_sparse(graph,sparse_lap,diag,status,normed=a==1, &
                    use_out_degree=b==1,symmetrized=c==1)
            end if
            write(*,'(i0,1x,i0)') status,nr
            if (status == csgraph_success) then
                dense = csgraph_to_dense(sparse_lap,status)
                do i = 1, nr
                    write(*,'(es25.16)') diag(i)
                end do
                do i = 1, nr
                    write(*,'(*(es25.16,1x))') dense(i,:)
                end do
            end if
        else
            error stop 'unknown reference task'
        end if
        deallocate(rows,cols,values)
    end do
end program csgraph_weighted_driver
