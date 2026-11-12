! SPDX-License-Identifier: MIT
program test_spatial_rotation_stack
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_quat, rotation_identity
    use scifort_spatial_rotation_stack, only : rotation_stack_type, &
        rotation_stack_identity, rotation_stack_from_rotations, &
        rotation_stack_from_quat, rotation_stack_from_rotvec, &
        rotation_stack_from_matrix, rotation_stack_from_euler, &
        rotation_stack_from_mrp, rotation_stack_concatenate
    implicit none
    type(rotation_stack_type) :: stack, singles, values, subset, combined, bad, temp
    type(rotation_type) :: r, invalid
    real(dp) :: pi, rv(3,3), v(3), vectors(3,3), w(3), mrps(3,3)
    real(dp), allocatable :: result(:,:), quats(:,:), m3(:,:,:), angles(:,:), reals(:)
    logical, allocatable :: locks(:), equal(:)
    integer :: status

    pi = acos(-1.0_dp)
    rv(1,:) = [0.0_dp,0.0_dp,0.0_dp]
    rv(2,:) = [0.0_dp,0.0_dp,pi/2.0_dp]
    rv(3,:) = [pi,0.0_dp,0.0_dp]
    stack = rotation_stack_from_rotvec(rv,status=status)
    if (status /= 0 .or. .not. stack%valid() .or. stack%size() /= 3) error stop 'make stack'
    quats = stack%as_quat()
    if (any(shape(quats) /= [3,4])) error stop 'quaternion shape'
    values = rotation_stack_from_quat(quats,status=status)
    if (status /= 0) error stop 'from_quat'
    equal = stack%approx_equal(values)
    if (size(equal) /= 3 .or. .not. all(equal)) error stop 'quat roundtrip'
    m3 = stack%as_matrix()
    if (any(shape(m3) /= [3,3,3])) error stop 'matrix shape'
    values = rotation_stack_from_matrix(m3,status=status)
    if (status /= 0 .or. .not. all(stack%approx_equal(values))) error stop 'matrix roundtrip'
    mrps = stack%as_mrp()
    values = rotation_stack_from_mrp(mrps,status=status)
    if (status /= 0 .or. .not. all(stack%approx_equal(values))) error stop 'mrp roundtrip'

    v = [1.0_dp,0.0_dp,0.0_dp]
    result = stack%apply(v,status=status)
    if (status /= 0 .or. any(shape(result) /= [3,3])) error stop 'apply shape'
    call close(result(1,:),v,1.0e-14_dp)
    call close(result(2,:),[0.0_dp,1.0_dp,0.0_dp],1.0e-14_dp)
    call close(result(3,:),v,1.0e-14_dp)
    vectors(1,:) = [1.0_dp,0.0_dp,0.0_dp]
    vectors(2,:) = [0.0_dp,1.0_dp,0.0_dp]
    vectors(3,:) = [0.0_dp,0.0_dp,1.0_dp]
    result = stack%apply(vectors,status=status)
    if (status /= 0 .or. any(shape(result) /= [3,3])) error stop 'many vectors'
    call close(result(2,:),[-1.0_dp,0.0_dp,0.0_dp],1.0e-14_dp)
    temp = rotation_stack_identity(1)
    result = temp%apply(vectors,status=status)
    if (status /= 0 .or. any(shape(result) /= [3,3])) error stop 'singleton broadcast'
    if (maxval(abs(result-vectors)) > 1.0e-14_dp) error stop 'singleton result'
    result = stack%apply(vectors(1:1,:),status=status)
    if (status /= 0 .or. any(shape(result) /= [3,3])) error stop 'one vector row'
    if (maxval(abs(result-stack%apply(v))) > 1.0e-14_dp) error stop 'broadcast value'
    ! The inverse must undo forward rotation on a per-row basis.
    result = stack%apply(vectors)
    temp = stack%inverse(status=status)
    if (status /= 0) error stop 'inverse status'
    result = temp%apply(result)
    if (maxval(abs(result-vectors)) > 3.0e-14_dp) error stop 'inverse undo'

    values = stack%compose(stack%inverse(),status=status)
    if (status /= 0) error stop 'compose inverse'
    if (maxval(abs(values%magnitude())) > 3.0e-14_dp) error stop 'identity rotation'
    singles = rotation_stack_from_rotvec(rv(:1,:))
    combined = stack%compose(singles,status=status)
    if (status /= 0 .or. .not. all(combined%approx_equal(stack))) error stop 'right broadcast'
    combined = singles%compose(stack,status=status)
    if (status /= 0 .or. .not. all(combined%approx_equal(stack))) error stop 'left broadcast'
    r = rotation_identity()
    combined = stack%compose(r,status=status)
    if (status /= 0 .or. .not. all(combined%approx_equal(stack))) error stop 'scalar compose'

    subset = stack%slice([3,1,3],status=status)
    if (status /= 0 .or. subset%size() /= 3) error stop 'slice size'
    if (.not. all(subset%approx_equal(rotation_stack_from_rotvec(rv([3,1,3],:))))) &
        error stop 'slice order'
    r = stack%at(2,status=status)
    if (status /= 0) error stop 'at'
    call close(r%as_rotvec(),rv(2,:),1.0e-14_dp)
    r = stack%at(4,status=status)
    if (status /= 1 .or. r%valid()) error stop 'out-of-range at'
    combined = rotation_stack_concatenate([subset,stack],status=status)
    if (status /= 0 .or. combined%size() /= 6) error stop 'concatenation'
    temp = combined%slice([4,5,6])
    if (.not. all(temp%approx_equal(stack))) error stop 'concat suffix'

    values = stack%power(2.0_dp,status=status)
    if (status /= 0) error stop 'power status'
    combined = stack%compose(stack)
    if (.not. all(values%approx_equal(combined,atol=1.0e-6_dp))) error stop 'power squared'
    values = stack%power(0.5_dp,status=status)
    if (status /= 0) error stop 'root status'
    combined = values%compose(values)
    if (.not. all(combined%approx_equal(stack,atol=1.0e-6_dp))) error stop 'root squared'

    values = rotation_stack_from_euler('zyx',reshape([0.0_dp,0.0_dp,0.0_dp, &
        20.0_dp,10.0_dp,30.0_dp, 90.0_dp,0.0_dp,0.0_dp],[3,3],order=[2,1]),degrees=.true.,status=status)
    if (status /= 0 .or. values%size() /= 3) error stop 'euler creation'
    angles = values%as_euler('zyx',degrees=.true.,gimbal_lock=locks,status=status)
    if (status /= 0 .or. any(shape(angles) /= [3,3])) error stop 'euler extraction'
    if (any(locks)) error stop 'unexpected gimbal lock'
    temp = rotation_stack_from_euler('zyx',angles,degrees=.true.)
    if (.not. all(temp%approx_equal(values))) error stop 'euler roundtrip'
    reals = stack%magnitude()
    if (size(reals) /= 3) error stop 'magnitudes'
    w = [2.0_dp,3.0_dp,4.0_dp]
    r = stack%mean(weights=w,status=status)
    if (status /= 0 .or. .not. r%valid()) error stop 'mean'

    ! Invalid inputs fail without creating partially valid arrays.
    bad = rotation_stack_from_quat(reshape([1.0_dp,0.0_dp],[2,1]),status=status)
    if (status /= 1 .or. bad%valid()) error stop 'wrong quaternion shape'
    bad = rotation_stack_from_rotvec(rv(:,1:2),status=status)
    if (status /= 1 .or. bad%valid()) error stop 'wrong rotvec shape'
    quats = stack%as_quat()
    quats(2,:) = 0.0_dp
    bad = rotation_stack_from_quat(quats,status=status)
    if (status /= 1 .or. bad%valid()) error stop 'invalid middle quaternion'
    quats = stack%as_quat(scalar_first=.true.)
    values = rotation_stack_from_quat(quats,scalar_first=.true.,status=status)
    if (status /= 0 .or. .not. all(values%approx_equal(stack))) error stop 'scalar-first quaternion'
    equal = stack%approx_equal(singles,atol=-1.0_dp,status=status)
    if (status /= 1 .or. size(equal) /= 0) error stop 'negative equality threshold'
    bad = stack%slice([1,4],status=status)
    if (status /= 1 .or. bad%valid()) error stop 'invalid slice'
    bad = rotation_stack_from_rotvec(rv(:2,:))
    values = stack%compose(bad,status=status)
    if (status /= 1 .or. values%valid()) error stop 'invalid broadcast composition'
    result = stack%apply(vectors(:2,:),status=status)
    if (status /= 1 .or. size(result,1) /= 0) error stop 'invalid broadcast apply'
    result = stack%apply([1.0_dp,2.0_dp],status=status)
    if (status /= 1 .or. size(result,1) /= 0) error stop 'invalid vector size'
    values = stack%power(ieee_value(0.0_dp,ieee_quiet_nan),status=status)
    if (status /= 1 .or. values%valid()) error stop 'invalid power'
    invalid = rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,0.0_dp])
    values = rotation_stack_from_rotations([rotation_identity(),invalid],status=status)
    if (status /= 1 .or. values%valid()) error stop 'invalid component'
    values = rotation_stack_identity(0,status=status)
    if (status /= 0 .or. .not. values%valid() .or. values%size() /= 0) error stop 'empty identity'
    values = rotation_stack_identity(-1,status=status)
    if (status /= 1 .or. values%valid()) error stop 'negative identity'

    print '(a)', 'test_spatial_rotation_stack PASS'
contains
    subroutine close(actual,expected,tol)
        real(dp), intent(in) :: actual(:) !! vector to check
        real(dp), intent(in) :: expected(:) !! expected vector
        real(dp), intent(in) :: tol !! max absolute error
        if (size(actual) /= size(expected)) error stop 'close shape'
        if (.not. all(abs(actual-expected) <= tol)) then
            print *, 'Mismatch:',actual,expected
            error stop 'close values'
        end if
    end subroutine close
end program test_spatial_rotation_stack
