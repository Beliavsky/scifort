! SPDX-License-Identifier: MIT
program test_spatial_convex_hull
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial, only : convex_hull, convex_hull_type, &
        spatial_distance_success, spatial_distance_invalid_input
    implicit none
    type(convex_hull_type) :: hull
    real(dp) :: square(7,2), cube(10,3), shifted(10,3), flat(5,3), eps
    real(dp) :: n(3), dot, tetra(4,3), expected_volume
    integer :: i, j, k, code, idx, s, f, v

    eps = 4.0e-12_dp
    square = 0.0_dp
    square(1,:) = [0.0_dp,0.0_dp]
    square(2,:) = [1.0_dp,0.0_dp]
    square(3,:) = [1.0_dp,1.0_dp]
    square(4,:) = [0.0_dp,1.0_dp]
    square(5,:) = [0.4_dp,0.2_dp]
    square(6,:) = square(3,:)
    square(7,:) = [0.5_dp,0.0_dp]
    hull = convex_hull(square,status=code)
    if (code /= spatial_distance_success .or. .not. hull%valid()) &
        error stop 'square status'
    if (size(hull%vertices) /= 4 .or. any(shape(hull%simplices) /= [2,4])) &
        error stop 'square vertices/simplices'
    if (maxval(abs([hull%area-4.0_dp,hull%volume-1.0_dp]))>eps) &
        error stop 'square area/volume'
    if (any(hull%vertices<1) .or. any(hull%vertices>4)) error stop 'square vertices'
    if (any(hull%neighbors<1) .or. any(hull%neighbors>4)) error stop '2d adjacency'
    do i=1,4
        if (hull%simplices(1,i) /= hull%vertices(i)) error stop '2d simplex order'
        if (hull%simplices(2,i) /= hull%vertices(mod(i,4)+1)) error stop '2d edge order'
        if (hull%neighbors(1,i) /= mod(i,4)+1) error stop '2d next edge'
        if (hull%neighbors(2,i) /= mod(i+2,4)+1) error stop '2d previous edge'
        do j=1,7
            if (dot_product(hull%equations(1:2,i),square(j,:)) + &
                hull%equations(3,i) > eps) error stop '2d outward equation'
        end do
    end do

    s = 0
    do i=-1,1,2
        do j=-1,1,2
            do k=-1,1,2
                s = s+1
                cube(s,:) = real([i,j,k],dp)
            end do
        end do
    end do
    cube(9,:) = [0.0_dp,0.0_dp,0.0_dp]
    cube(10,:) = cube(5,:)
    hull = convex_hull(cube,status=code)
    if (code/=0 .or. .not. hull%valid()) error stop 'cube status'
    if (size(hull%vertices)/=8 .or. any(shape(hull%simplices)/=[3,12])) &
        error stop 'cube topology'
    if (abs(hull%area-24.0_dp)>eps .or. abs(hull%volume-8.0_dp)>eps) &
        error stop 'cube area/volume'
    if (any(hull%vertices /= [(i,i=1,8)])) error stop 'cube vertex index order'
    do f=1,size(hull%simplices,2)
        do v=1,3
            idx = hull%neighbors(v,f)
            if (idx<1 .or. idx>size(hull%simplices,2)) error stop '3d closed facets'
            if (count([(any(hull%simplices(:,idx)==hull%simplices(k,f)), k=1,3)]) /= 2) &
                error stop '3d shared facet edge'
        end do
        n = hull%equations(1:3,f)
        if (abs(norm2(n)-1.0_dp) > eps) error stop 'face unit normal'
        do j=1,size(cube,1)
            dot = sum(n*cube(j,:)) + hull%equations(4,f)
            if (dot > eps) error stop '3d outward equation'
        end do
    end do

    ! Irregular tetrahedron: a closed four-triangle hull with known volume.
    tetra(1,:) = [0.0_dp,0.0_dp,0.0_dp]
    tetra(2,:) = [2.0_dp,0.0_dp,0.0_dp]
    tetra(3,:) = [0.0_dp,3.0_dp,0.0_dp]
    tetra(4,:) = [0.0_dp,0.0_dp,4.0_dp]
    expected_volume = 4.0_dp
    hull = convex_hull(tetra,status=code)
    if (code/=0 .or. .not. hull%valid()) error stop 'tetrahedron status'
    if (size(hull%vertices)/=4 .or. size(hull%simplices,2)/=4) &
        error stop 'tetrahedron topology'
    if (abs(hull%volume-expected_volume)>eps) error stop 'tetrahedron volume'

    shifted = 5.0_dp*cube
    shifted(:,1) = shifted(:,1) + 12.0_dp
    shifted(:,2) = shifted(:,2) - 40.0_dp
    shifted(:,3) = shifted(:,3) + 110.0_dp
    hull = convex_hull(shifted,status=code)
    if (code /= 0) error stop 'translated hull status'
    if (abs(hull%area-600.0_dp)>eps .or. abs(hull%volume-1000.0_dp)>eps) &
        error stop 'translated volume/surface area'

    hull = convex_hull(square(1:2,:),status=code)
    if (code /= spatial_distance_invalid_input .or. hull%valid()) error stop '2d few points'
    hull = convex_hull(cube(1:3,:),status=code)
    if (code /= spatial_distance_invalid_input) error stop '3d few points'
    flat = 0.0_dp
    do i=1,5
        flat(i,:) = [real(i,dp),real(i*i,dp),0.0_dp]
    end do
    hull = convex_hull(flat,status=code)
    if (code /= spatial_distance_invalid_input) error stop 'flat 3d points'
    square(:,2) = 0.0_dp
    hull = convex_hull(square,status=code)
    if (code /= spatial_distance_invalid_input) error stop 'collinear 2d points'
    cube(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    hull = convex_hull(cube,status=code)
    if (code /= spatial_distance_invalid_input) error stop 'nonfinite coordinates'
    print '(a)', 'test_spatial_convex_hull PASS'
end program test_spatial_convex_hull
