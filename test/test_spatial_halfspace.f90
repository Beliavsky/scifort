! SPDX-License-Identifier: MIT
program test_spatial_halfspace
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial, only : halfspace_intersection, halfspace_intersection_type, &
        convex_hull, convex_hull_type, spatial_distance_success, spatial_distance_invalid_input
    implicit none
    type(halfspace_intersection_type) :: poly
    type(convex_hull_type) :: hull
    real(dp) :: square(6,3), cuboid(6,4), octahedron(8,4), simplex(4,4), &
                expected_square(4,2), expected_cube(8,3), expected_octa(6,3), &
                expected_tetra(4,3), center(3), moved(6,4), short(3,3)
    integer :: i, j, k, m, code

    square(1,:) = [1.0_dp,0.0_dp,-1.0_dp]
    square(2,:) = [-1.0_dp,0.0_dp,-1.0_dp]
    square(3,:) = [0.0_dp,1.0_dp,-2.0_dp]
    square(4,:) = [0.0_dp,-1.0_dp,-2.0_dp]
    square(5,:) = [1.0_dp,0.0_dp,-3.0_dp] ! redundant
    square(6,:) = square(1,:)              ! repeated face
    poly = halfspace_intersection(square,[0.0_dp,0.0_dp],status=code)
    if (code /= spatial_distance_success .or. .not. poly%valid()) error stop 'rectangle status'
    expected_square(1,:) = [-1.0_dp,-2.0_dp]
    expected_square(2,:) = [-1.0_dp,2.0_dp]
    expected_square(3,:) = [1.0_dp,-2.0_dp]
    expected_square(4,:) = [1.0_dp,2.0_dp]
    call check_vertices(poly%intersections,expected_square,1.0e-11_dp)
    if (size(poly%dual_facets,2) /= 4) error stop 'rectangle polar facet count'
    if (maxval(abs(poly%dual_points(1,:)-[1.0_dp,0.0_dp])) > 1.0e-12_dp) &
        error stop 'rectangle polar point'

    cuboid(1,:) = [1.0_dp,0.0_dp,0.0_dp,-2.0_dp]
    cuboid(2,:) = [-1.0_dp,0.0_dp,0.0_dp,-2.0_dp]
    cuboid(3,:) = [0.0_dp,1.0_dp,0.0_dp,-3.0_dp]
    cuboid(4,:) = [0.0_dp,-1.0_dp,0.0_dp,-3.0_dp]
    cuboid(5,:) = [0.0_dp,0.0_dp,1.0_dp,-4.0_dp]
    cuboid(6,:) = [0.0_dp,0.0_dp,-1.0_dp,-4.0_dp]
    m = 0
    do i=-1,1,2
        do j=-1,1,2
            do k=-1,1,2
                m = m+1
                expected_cube(m,:) = [2.0_dp*real(i,dp),3.0_dp*real(j,dp),4.0_dp*real(k,dp)]
            end do
        end do
    end do
    poly = halfspace_intersection(cuboid,[0.0_dp,0.0_dp,0.0_dp],status=code)
    if (code /= 0 .or. .not. poly%valid()) error stop 'cuboid status'
    call check_vertices(poly%intersections,expected_cube,1.0e-10_dp)
    hull = convex_hull(poly%intersections)
    if (.not. hull%valid()) error stop 'cuboid output hull'
    if (abs(hull%volume-192.0_dp) > 1.0e-9_dp) error stop 'cuboid volume'
    ! Every dual facet equation supports its vertices and contains all
    ! polar points on its inner side, regardless of triangulation order.
    do i=1,size(poly%dual_facets,2)
        if (any(poly%dual_facets(:,i)<1) .or. &
            any(poly%dual_facets(:,i)>size(cuboid,1))) error stop 'dual indices'
        if (abs(norm2(poly%dual_equations(1:3,i))-1.0_dp)>1.0e-11_dp) &
            error stop 'dual equation normal'
        do j=1,size(poly%dual_points,1)
            if (dot_product(poly%dual_equations(1:3,i),poly%dual_points(j,:)) + &
                poly%dual_equations(4,i)>1.0e-10_dp) error stop 'dual containment'
        end do
        do j=1,3
            k = poly%dual_facets(j,i)
            if (abs(dot_product(poly%dual_equations(1:3,i),poly%dual_points(k,:)) + &
                poly%dual_equations(4,i)) > 1.0e-10_dp) error stop 'dual support'
        end do
    end do

    ! The polar of an octahedron has coplanar square facets; each is
    ! triangulated by the hull builder but corresponds to ONE primal vertex.
    m = 0
    do i=-1,1,2
        do j=-1,1,2
            do k=-1,1,2
                m = m+1
                octahedron(m,1:3) = real([i,j,k],dp)
            end do
        end do
    end do
    octahedron(:,4) = -1.0_dp
    expected_octa(1,:) = [1.0_dp,0.0_dp,0.0_dp]
    expected_octa(2,:) = [-1.0_dp,0.0_dp,0.0_dp]
    expected_octa(3,:) = [0.0_dp,1.0_dp,0.0_dp]
    expected_octa(4,:) = [0.0_dp,-1.0_dp,0.0_dp]
    expected_octa(5,:) = [0.0_dp,0.0_dp,1.0_dp]
    expected_octa(6,:) = [0.0_dp,0.0_dp,-1.0_dp]
    poly = halfspace_intersection(octahedron,[0.0_dp,0.0_dp,0.0_dp],status=code)
    if (code /= 0) error stop 'octahedron status'
    call check_vertices(poly%intersections,expected_octa,1.0e-10_dp)
    if (size(poly%dual_facets,2) /= 12) error stop 'cube polar triangulation'

    simplex(1,:) = [-1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    simplex(2,:) = [0.0_dp,-1.0_dp,0.0_dp,0.0_dp]
    simplex(3,:) = [0.0_dp,0.0_dp,-1.0_dp,0.0_dp]
    simplex(4,:) = [1.0_dp,1.0_dp,1.0_dp,-1.0_dp]
    expected_tetra(1,:) = [0.0_dp,0.0_dp,0.0_dp]
    expected_tetra(2,:) = [1.0_dp,0.0_dp,0.0_dp]
    expected_tetra(3,:) = [0.0_dp,1.0_dp,0.0_dp]
    expected_tetra(4,:) = [0.0_dp,0.0_dp,1.0_dp]
    poly = halfspace_intersection(simplex,[0.2_dp,0.2_dp,0.2_dp],status=code)
    if (code /= 0) error stop 'simplex status'
    call check_vertices(poly%intersections,expected_tetra,1.0e-10_dp)

    center = [120.0_dp,-75.0_dp,500.0_dp]
    moved = cuboid
    do i=1,6
        moved(i,4) = cuboid(i,4)*1.0e-4_dp - dot_product(cuboid(i,1:3),center)
    end do
    poly = halfspace_intersection(moved,center,status=code)
    if (code /= 0) error stop 'translated small box status'
    call check_vertices(poly%intersections,spread(center,dim=1,ncopies=8)+ &
        expected_cube*1.0e-4_dp,5.0e-10_dp)

    ! Strict interior, unboundedness, invalid data, and wrong shapes.
    poly = halfspace_intersection(square,[1.0_dp,0.0_dp],status=code)
    if (code /= spatial_distance_invalid_input .or. poly%valid()) error stop 'boundary point'
    poly = halfspace_intersection(square(1:3,:),[0.0_dp,0.0_dp],status=code)
    if (code /= spatial_distance_invalid_input) error stop 'unbounded input'
    poly = halfspace_intersection(cuboid(1:5,:),[0.0_dp,0.0_dp,0.0_dp],status=code)
    if (code /= spatial_distance_invalid_input) error stop '3d unbounded input'
    poly = halfspace_intersection(square(1:2,:),[0.0_dp,0.0_dp],status=code)
    if (code /= spatial_distance_invalid_input) error stop 'few constraints'
    short = square(1:3,:)
    short(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    poly = halfspace_intersection(short,[0.0_dp,0.0_dp],status=code)
    if (code /= spatial_distance_invalid_input) error stop 'nonfinite input'
    poly = halfspace_intersection(cuboid,[0.0_dp,0.0_dp],status=code)
    if (code /= spatial_distance_invalid_input) error stop 'dimension mismatch'
    print '(a)', 'test_spatial_halfspace PASS'
contains
    subroutine check_vertices(actual,expected,tol)
        real(dp), intent(in) :: actual(:,:) !! actual vertices (n,d)
        real(dp), intent(in) :: expected(:,:) !! reference vertices (n,d)
        real(dp), intent(in) :: tol !! maximum Euclidean discrepancy
        integer :: i,j
        if (any(shape(actual)/=shape(expected))) then
            print *, 'actual/expected vertex counts', shape(actual), shape(expected)
            error stop 'vertex array shape'
        end if
        do i=1,size(actual,1)
            if (minval([(norm2(actual(i,:)-expected(j,:)),j=1,size(expected,1))])>tol) &
                error stop 'unmatched vertex'
        end do
        do i=1,size(expected,1)
            if (minval([(norm2(expected(i,:)-actual(j,:)),j=1,size(actual,1))])>tol) &
                error stop 'missing vertex'
        end do
    end subroutine check_vertices
end program test_spatial_halfspace
