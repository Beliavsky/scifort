! SPDX-License-Identifier: MIT
program halfspace_driver
    use scifort_kinds, only : dp
    use scifort_spatial, only : halfspace_intersection, halfspace_intersection_type
    implicit none
    type(halfspace_intersection_type) :: poly
    real(dp), allocatable :: hs(:,:), interior(:)
    integer :: ncases, n, d, i, k, code

    read(*,*) ncases
    do k=1,ncases
        read(*,*) n,d
        allocate(hs(n,d+1),interior(d))
        read(*,*) interior
        do i=1,n
            read(*,*) hs(i,:)
        end do
        poly = halfspace_intersection(hs,interior,status=code)
        if (code == 0) then
            write(*,'(3(i0,1x))') code,size(poly%intersections,1),size(poly%dual_facets,2)
            write(*,'(2(es26.17e3,1x))') poly%dual_area,poly%dual_volume
            do i=1,size(poly%intersections,1)
                write(*,'(*(es26.17e3,1x))') poly%intersections(i,:)
            end do
            do i=1,n
                write(*,'(*(es26.17e3,1x))') poly%dual_points(i,:)
            end do
        else
            write(*,'(3(i0,1x))') code,0,0
        end if
        deallocate(hs,interior)
    end do
end program halfspace_driver
