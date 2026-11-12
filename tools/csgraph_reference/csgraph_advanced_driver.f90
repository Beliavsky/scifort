! SPDX-License-Identifier: MIT
! Text protocol for deterministic randomized SciPy comparisons.
program csgraph_advanced_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csr_bipartite_graph, &
        csgraph_from_edges, csgraph_to_dense, bipartite_from_edges, &
        reverse_cuthill_mckee, maximum_flow, maximum_bipartite_matching, &
        structural_rank, csgraph_success
    implicit none
    type(csr_graph) :: graph, flow
    type(csr_bipartite_graph) :: bipartite
    real(dp), allocatable :: weights(:), dense(:,:)
    integer, allocatable :: rows(:), cols(:), matching(:), matching2(:), perm(:)
    integer :: trials, trial, n, m, a1, a2, i, j, status, status2, rank
    integer(int64) :: value
    character(len=1) :: task
    read(*,*) trials
    do trial = 1, trials
        read(*,*) task, n, m, a1, a2
        allocate(rows(m),cols(m),weights(m))
        do i = 1, m
            read(*,*) rows(i), cols(i), weights(i)
        end do
        select case (task)
        case ('F')
            graph = csgraph_from_edges(n,rows,cols,weights,status)
            call maximum_flow(graph,a1,a2,value,flow,status)
            write(*,'(i0,1x,i0)') status,value
            if (status == csgraph_success) then
                dense = csgraph_to_dense(flow,status2)
                do i = 1, n
                    write(*,'(*(i0,1x))') (nint(dense(i,j),kind=int64),j=1,n)
                end do
                deallocate(dense)
            end if
        case ('R')
            graph = csgraph_from_edges(n,rows,cols,weights,status)
            call reverse_cuthill_mckee(graph,perm,status,symmetric_mode=(a1==1))
            write(*,'(*(i0,1x))') status,perm
        case ('M')
            bipartite = bipartite_from_edges(n,a1,rows,cols,status)
            call maximum_bipartite_matching(bipartite,matching,status,'row')
            write(*,'(*(i0,1x))') status,matching
            call maximum_bipartite_matching(bipartite,matching2,status,'column')
            write(*,'(*(i0,1x))') status,matching2
            call structural_rank(bipartite,rank,status)
            write(*,'(*(i0,1x))') status,rank
        case default
            error stop 'unknown reference case'
        end select
        deallocate(rows,cols,weights)
    end do
end program csgraph_advanced_driver
