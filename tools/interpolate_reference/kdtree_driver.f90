! SPDX-License-Identifier: MIT
! Line-oriented driver to compare independent SciFort KD-tree with scipy.spatial.cKDTree.
program scifort_kdtree_driver
    use scifort_kinds, only: dp
    use scifort_interpolate_kdtree, only: kd_tree, kd_tree_type
    implicit none
    type(kd_tree_type) :: tree
    real(dp), allocatable :: x(:,:), q(:,:), d2(:)
    integer, allocatable :: ids(:), ball(:)
    integer :: n, dim, m, k, i, j, code
    real(dp) :: radius
    read(*,*) n,dim,m,k,radius
    allocate(x(n,dim),q(m,dim))
    do i=1,n
        read(*,*) x(i,:)
    end do
    do i=1,m
        read(*,*) q(i,:)
    end do
    tree=kd_tree(x,status=code)
    if (code/=0) error stop 'invalid test index'
    do i=1,m
        call tree%query(q(i,:),ids,distance2=d2,k=k,status=code)
        if (code/=0) error stop 'invalid kNN query'
        do j=1,k
            write(*,'(i0,1x,es26.17)') ids(j),d2(j)
        end do
        call tree%query_ball_point(q(i,:),radius,ball,status=code)
        if (code/=0) error stop 'invalid radius query'
        write(*,'(i0)') size(ball)
        do j=1,size(ball)
            write(*,'(i0)') ball(j)
        end do
    end do
end program scifort_kdtree_driver
