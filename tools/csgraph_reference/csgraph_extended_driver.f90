! SPDX-License-Identifier: MIT
! Data-only bridge for cross-checking the native graph algorithms with SciPy.
program csgraph_extended_driver
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csgraph_from_dense, &
        floyd_warshall, johnson, laplacian, depth_first_order, &
        csgraph_success
    implicit none
    type(csr_graph) :: graph
    real(dp), allocatable :: a(:,:), distances(:,:), l(:,:), diagonal(:)
    integer, allocatable :: predecessors(:,:), order(:), parents(:)
    integer :: cases, case_id, n, i, direction, unweight, sym, out, norm, status

    read(*,*) cases
    do case_id = 1, cases
        read(*,*) n
        allocate(a(n,n))
        do i = 1, n
            read(*,*) a(i,:)
        end do
        graph = csgraph_from_dense(a,status)
        if (status /= csgraph_success) error stop 'invalid reference graph'
        write(*,'(i0)') n
        do direction = 0, 1
            do unweight = 0, 1
                call floyd_warshall(graph,distances,predecessors,status, &
                    directed=direction==1,unweighted=unweight==1)
                write(*,'(i0)') status
                if (status == csgraph_success) then
                    do i = 1, n
                        write(*,'(*(es26.17,1x))') distances(i,:)
                    end do
                    do i = 1, n
                        write(*,'(*(i0,1x))') predecessors(i,:)
                    end do
                end if
                call johnson(graph,distances,predecessors,status, &
                    directed=direction==1,unweighted=unweight==1)
                write(*,'(i0)') status
                if (status == csgraph_success) then
                    do i = 1, n
                        write(*,'(*(es26.17,1x))') distances(i,:)
                    end do
                    do i = 1, n
                        write(*,'(*(i0,1x))') predecessors(i,:)
                    end do
                end if
            end do
            call depth_first_order(graph,1,order,parents,status,directed=direction==1)
            if (status /= csgraph_success) error stop 'invalid DFS graph'
            write(*,'(i0)') size(order)
            write(*,'(*(i0,1x))') order
            write(*,'(*(i0,1x))') parents
        end do
        do sym = 0, 1
            do out = 0, 1
                do norm = 0, 1
                    call laplacian(graph,l,diagonal,status, &
                        normed=norm==1,use_out_degree=out==1,symmetrized=sym==1)
                    write(*,'(i0)') status
                    if (status == csgraph_success) then
                        write(*,'(*(es26.17,1x))') diagonal
                        do i = 1, n
                            write(*,'(*(es26.17,1x))') l(i,:)
                        end do
                    end if
                end do
            end do
        end do
        deallocate(a)
    end do
end program csgraph_extended_driver
