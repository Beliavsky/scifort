! SPDX-License-Identifier: MIT
! Test oracle driver: emits each Fortran Hausdorff result and witness indices.
program hausdorff_driver
    use scifort_kinds, only : dp
    use scifort_spatial, only : directed_hausdorff, hausdorff_distance
    implicit none
    integer :: nc, n, m, ndim, k, i, j, code, forward_flag
    real(dp), allocatable :: u(:,:), v(:,:)
    real(dp) :: d
    logical :: forward

    read(*,*) nc
    do k=1,nc
        read(*,*) n,m,ndim
        allocate(u(n,ndim),v(m,ndim))
        do i=1,n
            read(*,*) u(i,:)
        end do
        do j=1,m
            read(*,*) v(j,:)
        end do
        d = directed_hausdorff(u,v,index_u=i,index_v=j,status=code)
        write(*,'(es28.17e3,1x,3(i0,1x))') d,i,j,code
        d = directed_hausdorff(v,u,index_u=i,index_v=j,status=code)
        write(*,'(es28.17e3,1x,3(i0,1x))') d,i,j,code
        d = hausdorff_distance(u,v,index_u=i,index_v=j,from_u=forward,status=code)
        forward_flag=0
        if (forward) forward_flag=1
        write(*,'(es28.17e3,1x,4(i0,1x))') d,i,j,forward_flag,code
        d = directed_hausdorff(u,v,index_u=i,index_v=j,status=code,early_break=.false.)
        write(*,'(es28.17e3,1x,3(i0,1x))') d,i,j,code
        deallocate(u,v)
    end do
end program hausdorff_driver
