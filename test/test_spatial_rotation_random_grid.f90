! SPDX-License-Identifier: MIT
program test_spatial_rotation_random_grid
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, rotation_from_rotvec, rotation_mean
    use scifort_spatial_rotation_stack, only : rotation_rng_type, rotation_rng, &
        rotation_random, rotation_stack_random, rotation_stack_type
    use scifort_spatial_rotation_grid, only : rotation_grid_type, &
        rotation_grid_from_quat, rotation_grid_from_rotvec, &
        rotation_grid_from_matrix, rotation_grid_from_mrp, rotation_grid_from_euler, rotation_grid_random
    implicit none
    type(rotation_rng_type) :: a, b, uninitialized
    type(rotation_stack_type) :: stack, averages, bad_stack
    type(rotation_grid_type) :: grid, singleton, expanded, other, bad
    type(rotation_type) :: r, r0, r1, r2, rotations(2,3)
    real(dp), allocatable :: q(:,:,:), q0(:,:), result(:,:,:), matrices(:,:,:,:), rv(:,:,:)
    real(dp), allocatable :: scratch(:,:,:), mrps(:,:,:), eulers(:,:,:)
    real(dp) :: vectors(2,1,3), unit(3), moments(3), second(3), x(3)
    real(dp) :: weights(2), tolerance, pi, nan
    integer :: i, j, status

    tolerance = 2.0e-13_dp
    pi = acos(-1.0_dp)
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    a = rotation_rng(987654321_int64)
    b = rotation_rng(987654321_int64)
    if (.not. a%valid() .or. .not. b%valid()) error stop 'rng construction'
    stack = rotation_stack_random(20, a, status)
    if (status /= 0 .or. .not. stack%valid() .or. stack%size() /= 20) error stop 'stack random'
    q0 = stack%as_quat()
    do i = 1, 20
        r = rotation_random(b, status)
        if (status /= 0 .or. .not. r%valid()) error stop 'one random'
        if (maxval(abs(r%as_quat()-q0(i,:))) > tolerance) error stop 'random stream splitting'
    end do
    r = rotation_random(a)
    r0 = rotation_random(b)
    if (maxval(abs(r%as_quat()-r0%as_quat())) > tolerance) error stop 'random stream continuation'
    bad_stack = rotation_stack_random(-1, status=status)
    if (status /= 1 .or. bad_stack%valid()) error stop 'negative random count'
    bad_stack = rotation_stack_random(2, uninitialized, status)
    if (status /= 1 .or. bad_stack%valid()) error stop 'uninitialized random generator'
    r = rotation_random(uninitialized, status)
    if (status /= 1 .or. r%valid()) error stop 'invalid random scalar'
    stack = rotation_stack_random(0, b, status)
    if (status /= 0 .or. .not. stack%valid() .or. stack%size() /= 0) error stop 'empty random stack'

    ! Check that the rotated z-axis is isotropic, not just that q is unit length.
    a = rotation_rng(123456789_int64)
    moments = 0.0_dp
    second = 0.0_dp
    unit = [0.0_dp, 0.0_dp, 1.0_dp]
    do i = 1, 12000
        r = rotation_random(a)
        x = r%apply(unit)
        moments = moments + x
        second = second + x*x
    end do
    moments = moments/12000.0_dp
    second = second/12000.0_dp
    if (maxval(abs(moments)) > 0.035_dp) error stop 'random isotropy means'
    if (maxval(abs(second-1.0_dp/3.0_dp)) > 0.035_dp) error stop 'random isotropy moments'

    ! The random grid fills the first Fortran axis fastest and shares the RNG.
    a = rotation_rng(-943_int64)
    b = rotation_rng(-943_int64)
    grid = rotation_grid_random([2,3],a,status)
    if (status /= 0 .or. .not. grid%valid()) error stop 'random grid'
    if (any(grid%shape() /= [2,3])) error stop 'grid shape'
    q = grid%as_quat()
    do j = 1, 3
        do i = 1, 2
            r = rotation_random(b)
            if (maxval(abs(r%as_quat()-q(i,j,:))) > tolerance) error stop 'grid stream order'
        end do
    end do
    other = rotation_grid_from_quat(q,status=status)
    if (status /= 0 .or. any(other%shape()/=[2,3])) error stop 'quaternion grid'
    if (maxval(abs(q-other%as_quat())) > tolerance) error stop 'quat grid roundtrip'
    rv = grid%as_rotvec()
    other = rotation_grid_from_rotvec(rv,status=status)
    if (status /= 0) error stop 'rotvec grid'
    do j = 1, 3
        do i = 1, 2
            r = grid%at(i,j)
            r0 = other%at(i,j)
            if (maxval(abs(r%as_matrix()-r0%as_matrix())) > tolerance) error stop 'rotvec grid roundtrip'
        end do
    end do
    matrices = grid%as_matrix()
    if (any(shape(matrices) /= [2,3,3,3])) error stop 'grid matrix shape'
    other = rotation_grid_from_matrix(matrices,status)
    if (status /= 0) error stop 'from matrix grid status'
    if (maxval(abs(other%as_matrix()-matrices)) > tolerance) error stop 'matrix grid roundtrip'
    mrps = grid%as_mrp()
    other = rotation_grid_from_mrp(mrps,status)
    if (status /= 0) error stop 'from mrp grid status'
    if (maxval(abs(other%as_matrix()-matrices)) > tolerance) error stop 'mrp grid roundtrip'
    eulers = grid%as_euler('ZYX', status=status)
    if (status /= 0) error stop 'as euler grid status'
    other = rotation_grid_from_euler('ZYX',eulers,status=status)
    if (status /= 0) error stop 'from euler grid status'
    if (maxval(abs(other%as_matrix()-matrices)) > tolerance) error stop 'euler grid roundtrip'
    do j = 1, 3
        do i = 1, 2
            r = grid%at(i,j,status)
            if (status /= 0) error stop 'at valid'
            if (maxval(abs(matrices(i,j,:,:)-r%as_matrix())) > tolerance) error stop 'matrix values'
        end do
    end do
    r = grid%at(3,2,status)
    if (status /= 1 .or. r%valid()) error stop 'invalid grid index'
    expanded = grid%compose(grid%inverse(),status)
    if (status /= 0 .or. maxval(abs(expanded%magnitude())) > tolerance) error stop 'inverse composition'
    unit = [1.0_dp,0.0_dp,0.0_dp]
    result = grid%apply(unit,status=status)
    if (status /= 0 .or. any(shape(result) /= [2,3,3])) error stop 'vector grid apply'
    vectors(1,1,:) = [1.0_dp,2.0_dp,3.0_dp]
    vectors(2,1,:) = [3.0_dp,2.0_dp,1.0_dp]
    result = grid%apply(vectors,status=status)
    if (status /= 0 .or. any(shape(result) /= [2,3,3])) error stop 'broadcast apply'
    do j = 1, 3
        do i = 1, 2
            r = grid%at(i,j)
            if (maxval(abs(result(i,j,:)-r%apply(vectors(i,1,:)))) > tolerance) error stop 'broadcast values'
        end do
    end do
    singleton = rotation_grid_from_rotvec(reshape([0.0_dp,0.0_dp,pi/4.0_dp, &
        0.0_dp,pi/5.0_dp,0.0_dp,pi/6.0_dp,0.0_dp,0.0_dp], [1,3,3]))
    expanded = grid%compose(singleton,status)
    if (status /= 0 .or. any(expanded%shape() /= [2,3])) error stop 'compose broadcast'
    do j = 1, 3
        do i = 1, 2
            r0 = grid%at(i,j)
            r1 = singleton%at(1,j)
            r = r0%compose(r1)
            r2 = expanded%at(i,j)
            if (maxval(abs(r2%as_matrix()-r%as_matrix())) > tolerance) &
                error stop 'compose broadcasting values'
        end do
    end do
    expanded = grid%compose(rotation_identity(),status)
    if (status /= 0 .or. maxval(abs(expanded%as_quat()-q)) > tolerance) error stop 'scalar compose'
    other = grid%inverse()
    result = other%apply(grid%apply(vectors), status=status)
    if (status /= 0) error stop 'inverse vector apply'
    do j = 1, 3
        do i = 1, 2
            if (maxval(abs(result(i,j,:)-vectors(i,1,:))) > tolerance) error stop 'inverse apply values'
        end do
    end do

    weights = [1.0_dp, 2.0_dp]
    averages = grid%mean(1,weights,status)
    if (status /= 0 .or. averages%size() /= 3) error stop 'weighted axis-one mean'
    do j = 1, 3
        rotations(:,j) = [grid%at(1,j),grid%at(2,j)]
        r = rotation_mean(rotations(:,j),weights)
        r0 = averages%at(j)
        if (maxval(abs(r0%as_matrix()-r%as_matrix())) > tolerance) error stop 'mean axis one'
    end do
    averages = grid%mean(2,status=status)
    if (status /= 0 .or. averages%size() /= 2) error stop 'axis two mean'
    do i = 1, 2
        r = rotation_mean(rotations(i,:))
        r0 = averages%at(i)
        if (maxval(abs(r0%as_matrix()-r%as_matrix())) > tolerance) error stop 'mean axis two'
    end do

    ! No partially valid grid may escape an error path.
    bad = rotation_grid_from_rotvec(reshape([0.0_dp,1.0_dp],[1,1,2]),status=status)
    if (status /= 1 .or. bad%valid()) error stop 'invalid rotvec width'
    q(2,3,4) = nan
    bad = rotation_grid_from_quat(q,status=status)
    if (status /= 1 .or. bad%valid()) error stop 'invalid quaternion in grid'
    bad = rotation_grid_from_euler('ZZZ',eulers,status=status)
    if (status /= 1 .or. bad%valid()) error stop 'invalid euler grid sequence'
    bad = rotation_grid_from_matrix(matrices(:,:,1:2,:),status)
    if (status /= 1 .or. bad%valid()) error stop 'invalid matrix grid shape'
    bad = rotation_grid_from_mrp(mrps(:,:,1:2),status)
    if (status /= 1 .or. bad%valid()) error stop 'invalid mrp grid shape'
    bad = rotation_grid_random([0,3],status=status)
    if (status /= 1 .or. bad%valid()) error stop 'invalid random grid dims'
    bad = rotation_grid_random([2,3],uninitialized,status)
    if (status /= 1 .or. bad%valid()) error stop 'invalid grid rng'
    scratch = reshape([(real(i,dp),i=1,24)], [2,4,3])
    result = grid%apply(scratch,status=status)
    if (status /= 1 .or. size(result) /= 0) error stop 'invalid apply shape'
    bad = rotation_grid_from_rotvec(reshape([0.0_dp,0.0_dp,0.0_dp, &
        0.0_dp,0.0_dp,0.0_dp],[2,1,3]))
    expanded = grid%compose(bad,status=status)
    if (status /= 0 .or. any(expanded%shape()/=[2,3])) error stop 'other-axis broadcast'
    scratch = reshape([(real(i,dp),i=1,18)], [3,2,3])
    result = grid%apply(scratch,status=status)
    if (status /= 1 .or. size(result) /= 0) error stop 'incompatible axes'
    averages = grid%mean(3,status=status)
    if (status /= 1 .or. averages%valid()) error stop 'invalid mean axis'
    averages = grid%mean(1, weights=[0.0_dp,0.0_dp], status=status)
    if (status /= 1 .or. averages%valid()) error stop 'all-zero weights'
    print '(a)', 'test_spatial_rotation_random_grid PASS'
end program test_spatial_rotation_random_grid
