! SPDX-License-Identifier: MIT
program csgraph_dispatch_driver
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_graph, csgraph_from_dense, csgraph_to_dense, &
        shortest_path, breadth_first_tree, depth_first_tree, csgraph_success
    implicit none
    type(csr_graph) :: graph, tree
    real(dp), allocatable :: a(:,:), d(:,:), tree_dense(:,:)
    integer, allocatable :: pred(:,:)
    integer :: n, case_id, cases, i, j, which, direct, uw, code
    character(len=4), parameter :: algorithms(5) = [character(len=4)::'auto','FW','J','BF','D']

    read(*,*) cases
    do case_id=1,cases
        read(*,*) n
        allocate(a(n,n))
        do i=1,n
            read(*,*) a(i,:)
        end do
        graph = csgraph_from_dense(a,code)
        if (code /= csgraph_success) error stop 'invalid graph'
        write(*,'(i0)') n
        do direct=0,1
            do uw=0,1
                do which=1,5
                    call shortest_path(graph,d,pred,code,method=trim(algorithms(which)), &
                        directed=direct==1,unweighted=uw==1,indices=[1,n])
                    write(*,'(i0)') code
                    if (code /= csgraph_success) cycle
                    do i=1,2
                        write(*,'(*(es26.17,1x))') d(i,:)
                    end do
                    do i=1,2
                        write(*,'(*(i0,1x))') pred(i,:)
                    end do
                end do
            end do
            do which=1,2
                if (which==1) then
                    tree = breadth_first_tree(graph,1,code,directed=direct==1)
                else
                    tree = depth_first_tree(graph,1,code,directed=direct==1)
                end if
                write(*,'(i0)') code
                if (code /= csgraph_success) cycle
                tree_dense = csgraph_to_dense(tree,code)
                if (code /= csgraph_success) error stop 'invalid tree'
                do j=1,n
                    write(*,'(*(es26.17,1x))') tree_dense(j,:)
                end do
            end do
        end do
        deallocate(a)
    end do
end program csgraph_dispatch_driver
