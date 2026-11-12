! SPDX-License-Identifier: MIT
! Test oracle driver: evaluates great-circle interpolation at scalar/vector t.
program slerp_driver
    use scifort_kinds, only : dp
    use scifort_spatial, only : geometric_slerp
    implicit none
    integer :: nc, ndim, nt, k, i, code_array, code_scalar
    real(dp), allocatable :: start(:), finish(:), times(:), point(:), points(:,:)

    read(*,*) nc
    do k=1,nc
        read(*,*) ndim,nt
        allocate(start(ndim),finish(ndim),times(nt))
        read(*,*) start
        read(*,*) finish
        if (nt>0) read(*,*) times
        points = geometric_slerp(start,finish,times,status=code_array)
        if (nt>0) then
            point = geometric_slerp(start,finish,times(1),status=code_scalar)
        else
            code_scalar=0
            allocate(point(0))
        end if
        write(*,'(2(i0,1x))') code_array,code_scalar
        do i=1,nt
            write(*,'(*(es28.17e3,1x))') points(i,:)
        end do
        if (nt>0) write(*,'(*(es28.17e3,1x))') point
        deallocate(start,finish,times,point,points)
    end do
end program slerp_driver
