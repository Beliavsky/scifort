! SPDX-License-Identifier: MIT
! Stream protocol for checking 3D tetrahedralization against SciPy/Qhull.
program delaunay_3d_driver
    use scifort_kinds, only: dp
    use scifort_interpolate_scattered, only: delaunay_3d
    implicit none
    integer :: n,i,status
    real(dp), allocatable :: points(:,:)
    integer, allocatable :: tetrahedra(:,:)
    read(*,*) n
    allocate(points(n,3))
    do i=1,n
        read(*,*) points(i,:)
    end do
    call delaunay_3d(points,tetrahedra,status)
    print '(i0)',status
    if (status/=0) stop
    print '(i0)',size(tetrahedra,2)
    do i=1,size(tetrahedra,2)
        print '(4(i0,1x))',tetrahedra(:,i)
    end do
end program delaunay_3d_driver
