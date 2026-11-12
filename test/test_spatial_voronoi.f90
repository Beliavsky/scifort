! SPDX-License-Identifier: MIT
program test_spatial_voronoi
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial, only : spherical_voronoi, spherical_voronoi_type, &
        spatial_distance_success, spatial_distance_invalid_input
    implicit none
    real(dp), parameter :: pi=acos(-1.0_dp)
    type(spherical_voronoi_type) :: vor
    real(dp) :: oct(6,3), square(4,2), circle(5,2), cube(8,3), center(3), t(6,3), regular(4,3)
    real(dp) :: angles(5), gaps(5), e_areas(5)
    real(dp), allocatable :: areas(:)
    integer, allocatable :: region(:)
    integer :: i, j, k, m, code

    oct = 0.0_dp
    oct(1,:)=[1.0_dp,0.0_dp,0.0_dp]
    oct(2,:)=[-1.0_dp,0.0_dp,0.0_dp]
    oct(3,:)=[0.0_dp,1.0_dp,0.0_dp]
    oct(4,:)=[0.0_dp,-1.0_dp,0.0_dp]
    oct(5,:)=[0.0_dp,0.0_dp,1.0_dp]
    oct(6,:)=[0.0_dp,0.0_dp,-1.0_dp]
    vor = spherical_voronoi(oct, status=code)
    if (code /= spatial_distance_success .or. .not. vor%valid()) error stop 'octahedron status'
    if (any(shape(vor%vertices) /= [8,3])) error stop 'octahedron vertices'
    if (any(shape(vor%simplices) /= [3,8])) error stop 'octahedron hull faces'
    do i=1,6
        region = vor%region(i)
        if (size(region)/=4) error stop 'octahedron region degree'
        if (any(region<1) .or. any(region>8)) error stop 'octahedron region indices'
    end do
    if (size(vor%region(0)) /= 0 .or. size(vor%region(7))/=0) error stop 'invalid generator'
    areas = vor%calculate_areas()
    if (maxval(abs(areas-2.0_dp*pi/3.0_dp)) > 5.0e-13_dp) then
        print *, areas
        error stop 'octahedron area'
    end if
    if (abs(sum(areas)-4.0_dp*pi)>1.0e-12_dp) error stop 'total octahedron area'

    ! Planar quadrilateral facets in a cube are triangulated, with
    ! duplicate dual Voronoi vertices across coplanar hull triangles.
    m = 0
    do i=-1,1,2
        do j=-1,1,2
            do k=-1,1,2
                m = m+1
                cube(m,:) = real([i,j,k],dp)/sqrt(3.0_dp)
            end do
        end do
    end do
    vor = spherical_voronoi(cube,status=code)
    if (code /= 0 .or. any(shape(vor%vertices)/=[12,3])) &
        error stop 'cube coplanar facet triangulation'
    areas = vor%calculate_areas()
    if (maxval(abs(areas-0.5_dp*pi))>1.0e-12_dp) error stop 'cube equal areas'

    regular(1,:)=[1.0_dp,1.0_dp,1.0_dp]/sqrt(3.0_dp)
    regular(2,:)=[1.0_dp,-1.0_dp,-1.0_dp]/sqrt(3.0_dp)
    regular(3,:)=[-1.0_dp,1.0_dp,-1.0_dp]/sqrt(3.0_dp)
    regular(4,:)=[-1.0_dp,-1.0_dp,1.0_dp]/sqrt(3.0_dp)
    vor = spherical_voronoi(regular, status=code)
    if (code /= 0) error stop 'tetrahedron status'
    areas = vor%calculate_areas()
    if (maxval(abs(areas-pi))>1.0e-12_dp) error stop 'tetrahedron equal areas'

    center = [10.0_dp,-20.0_dp,30.0_dp]
    do i=1,6
        t(i,:) = center + 7.0_dp*oct(i,:)
    end do
    vor = spherical_voronoi(t,radius=7.0_dp,center=center,status=code)
    if (code /= 0) error stop 'translated sphere status'
    areas = vor%calculate_areas()
    if (maxval(abs(areas-49.0_dp*2.0_dp*pi/3.0_dp))>1.0e-11_dp) &
        error stop 'translated sphere areas'
    if (maxval(abs(vor%vertices(:,1)-center(1)))>7.01_dp) &
        error stop 'translated sphere vertices'

    square(1,:) = [1.0_dp,0.0_dp]
    square(2,:) = [0.0_dp,1.0_dp]
    square(3,:) = [-1.0_dp,0.0_dp]
    square(4,:) = [0.0_dp,-1.0_dp]
    vor = spherical_voronoi(square, status=code)
    if (code/=0 .or. any(shape(vor%vertices)/=[4,2])) error stop 'circle status'
    areas = vor%calculate_areas()
    if (maxval(abs(areas-0.5_dp*pi))>1.0e-12_dp) error stop 'circle equal arcs'
    vor = spherical_voronoi(square*3.0_dp,radius=3.0_dp,status=code)
    areas = vor%calculate_areas()
    if (code/=0 .or. maxval(abs(areas-1.5_dp*pi))>1.0e-12_dp) &
        error stop 'scaled circle'

    ! Unequally spaced circle sites: SciPy hull-edge arcs differ from the
    ! actual per-generator Voronoi arc lengths, except for regular spacing.
    angles = [0.0_dp, 0.2_dp, 1.6_dp, 3.4_dp, 5.1_dp]
    circle(:,1) = cos(angles)
    circle(:,2) = sin(angles)
    vor = spherical_voronoi(circle,status=code)
    if (code/=0) error stop 'irregular circle status'
    gaps(1:4) = angles(2:5)-angles(1:4)
    gaps(5) = 2.0_dp*pi-angles(5)+angles(1)
    areas = vor%calculate_areas()
    if (maxval(abs(areas-gaps))>2.0e-13_dp) error stop 'SciPy circle hull-edge arcs'
    areas = vor%region_measures()
    do i=1,5
        e_areas(i) = 0.5_dp*(gaps(modulo(i+3,5)+1)+gaps(i))
    end do
    if (maxval(abs(areas-e_areas))>2.0e-13_dp) then
        print *, areas,e_areas
        error stop 'circular per-generator arc lengths'
    end if

    ! Coplanar points on a great circle do not define a 3D spherical hull.
    t = 0.0_dp
    do i=1,6
        t(i,1) = cos(2.0_dp*pi*real(i-1,dp)/6.0_dp)
        t(i,2) = sin(2.0_dp*pi*real(i-1,dp)/6.0_dp)
    end do
    vor = spherical_voronoi(t,status=code)
    if (code/=spatial_distance_invalid_input) error stop 'coplanar 3D points'

    vor = spherical_voronoi(oct(:3,:),status=code)
    if (code/=spatial_distance_invalid_input .or. vor%valid()) &
        error stop 'too few points'
    if (size(vor%calculate_areas()) /= 0) error stop 'invalid areas'
    t = oct
    t(2,:) = t(1,:)
    vor = spherical_voronoi(t,status=code)
    if (code/=spatial_distance_invalid_input) error stop 'duplicate generator'
    t = oct
    t(2,:) = 2.0_dp*t(2,:)
    vor = spherical_voronoi(t,status=code)
    if (code/=spatial_distance_invalid_input) error stop 'off-sphere generator'
    t = oct
    t(2,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    vor = spherical_voronoi(t,status=code)
    if (code/=spatial_distance_invalid_input) error stop 'nonfinite generator'
    vor = spherical_voronoi(oct,center=[1.0_dp,2.0_dp],status=code)
    if (code/=spatial_distance_invalid_input) error stop 'wrong center shape'
    vor = spherical_voronoi(oct,radius=-1.0_dp,status=code)
    if (code/=spatial_distance_invalid_input) error stop 'bad radius'
    vor = spherical_voronoi(oct,threshold=-1.0_dp,status=code)
    if (code/=spatial_distance_invalid_input) error stop 'bad threshold'
    print '(a)', 'test_spatial_voronoi PASS'
end program test_spatial_voronoi
