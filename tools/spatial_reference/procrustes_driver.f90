! SPDX-License-Identifier: MIT
! Input/output bridge for independent scipy.spatial.procrustes comparisons.
program procrustes_driver
    use scifort_kinds, only : dp
    use scifort_spatial, only : procrustes, orthogonal_procrustes
    implicit none
    integer :: cases, scenario, n, d, i, code
    real(dp), allocatable :: a(:,:), b(:,:), x(:,:), y(:,:), rotation(:,:)
    real(dp) :: disparity, scale

    read(*,*) cases
    do scenario=1,cases
        read(*,*) n,d
        allocate(a(n,d),b(n,d))
        do i=1,n
            read(*,*) a(i,:)
        end do
        do i=1,n
            read(*,*) b(i,:)
        end do
        call procrustes(a,b,x,y,disparity,status=code)
        write(*,'(i0,1x,es26.17e3)') code,disparity
        if (code == 0) then
            do i=1,n
                write(*,'(*(es26.17e3,1x))') x(i,:)
            end do
            do i=1,n
                write(*,'(*(es26.17e3,1x))') y(i,:)
            end do
        end if
        call orthogonal_procrustes(a,b,rotation,scale,status=code)
        write(*,'(i0,1x,es26.17e3)') code,scale
        if (code == 0) then
            do i=1,d
                write(*,'(*(es26.17e3,1x))') rotation(i,:)
            end do
        end if
        deallocate(a,b,x,y,rotation)
    end do
end program procrustes_driver
