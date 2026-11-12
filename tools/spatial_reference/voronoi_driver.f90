! SPDX-License-Identifier: MIT
! Console bridge for independent scipy.spatial.SphericalVoronoi comparisons.
program voronoi_driver
    use scifort_kinds, only : dp
    use scifort_spatial, only : spherical_voronoi, spherical_voronoi_type
    implicit none
    type(spherical_voronoi_type) :: diagram
    real(dp), allocatable :: points(:,:), center(:), areas(:)
    integer, allocatable :: indices(:)
    integer :: cases, scenario, n, d, i, code, nvertices
    real(dp) :: radius, threshold

    read(*,*) cases
    do scenario=1,cases
        read(*,*) n,d,radius,threshold
        allocate(points(n,d),center(d))
        read(*,*) center
        do i=1,n
            read(*,*) points(i,:)
        end do
        diagram = spherical_voronoi(points,radius=radius,center=center,threshold=threshold,status=code)
        nvertices = 0
        if (code==0) nvertices=size(diagram%vertices,1)
        write(*,'(i0,1x,i0)') code,nvertices
        if (code == 0) then
            do i=1,nvertices
                write(*,'(*(es26.17e3,1x))') diagram%vertices(i,:)
            end do
            areas = diagram%calculate_areas()
            do i=1,n
                write(*,'(es26.17e3)') areas(i)
            end do
            do i=1,n
                indices = diagram%region(i)
                write(*,'(*(i0,1x))') size(indices),indices
            end do
        end if
        deallocate(points,center)
    end do
end program voronoi_driver
