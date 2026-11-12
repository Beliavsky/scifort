! SPDX-License-Identifier: MIT
! Batch driver for independent SciPy ConvexHull numerical comparisons.
program convex_hull_driver
    use scifort_kinds, only : dp
    use scifort_spatial, only : convex_hull, convex_hull_type
    implicit none
    type(convex_hull_type) :: hull
    real(dp), allocatable :: points(:,:)
    integer :: nt, caseid, n, d, i, j, status

    read(*,*) nt
    do caseid=1,nt
        read(*,*) n,d
        allocate(points(n,d))
        do i=1,n
            read(*,*) points(i,:)
        end do
        hull = convex_hull(points,status=status)
        if (hull%valid()) then
            print '(3(i0,1x))', status,size(hull%vertices),size(hull%simplices,2)
            print '(2(es26.17,1x))', hull%area,hull%volume
            print '(*(i0,1x))', hull%vertices
            do j=1,size(hull%simplices,2)
                print '(*(i0,1x))', hull%simplices(:,j)
                print '(*(es26.17,1x))', hull%equations(:,j)
                print '(*(i0,1x))', hull%neighbors(:,j)
            end do
        else
            print '(3(i0,1x))', status,0,0
        end if
        deallocate(points)
    end do
end program convex_hull_driver
