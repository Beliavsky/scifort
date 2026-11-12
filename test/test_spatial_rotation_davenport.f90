! SPDX-License-Identifier: MIT
! Rotation.from_davenport/as_davenport, including batched stacks and grids.
program test_spatial_rotation_davenport
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_davenport, &
        rotation_from_euler, spatial_rotation_invalid_input
    use scifort_spatial_rotation_stack, only : rotation_stack_type, &
        rotation_stack_from_davenport
    use scifort_spatial_rotation_grid, only : rotation_grid_type, &
        rotation_grid_from_davenport
    implicit none
    type(rotation_type) :: r, back, ref
    type(rotation_stack_type) :: stack, bad_stack
    type(rotation_grid_type) :: grid, bad_grid
    real(dp) :: axes(3,3), scaled(3,3), angles(3), recovered(3)
    real(dp) :: samples(3,3), lattice(2,2,3)
    real(dp), allocatable :: batch(:,:), all_angles(:,:,:)
    logical :: locked
    logical, allocatable :: flags(:), grid_flags(:,:)
    integer :: status, i, j, k
    axes(1,:) = [1.0_dp,0.0_dp,0.0_dp]
    axes(2,:) = [0.0_dp,1.0_dp,0.0_dp]
    axes(3,:) = [0.5_dp,0.0_dp,sqrt(3.0_dp)/2.0_dp]
    angles = [0.3_dp,0.5_dp,-0.8_dp]
    do k = 1, 2
        if (k == 1) then
            r = rotation_from_davenport(axes,'extrinsic',angles,status=status)
            recovered = r%as_davenport(axes,'extrinsic',gimbal_lock=locked,status=status)
            back = rotation_from_davenport(axes,'extrinsic',recovered,status=status)
        else
            r = rotation_from_davenport(axes,'intrinsic',angles,status=status)
            recovered = r%as_davenport(axes,'intrinsic',gimbal_lock=locked,status=status)
            back = rotation_from_davenport(axes,'intrinsic',recovered,status=status)
        end if
        if (status /= 0 .or. locked) error stop 'regular Davenport construction/extraction'
        call same_rotation(r,back,4.0e-14_dp)
        call close(recovered,angles,8.0e-14_dp)
    end do
    ! Equivalent to a standard Euler sequence for orthogonal global axes.
    axes(3,:) = [0.0_dp,0.0_dp,1.0_dp]
    r = rotation_from_davenport(axes,'e',angles)
    ref = rotation_from_euler('xyz',angles)
    call same_rotation(r,ref,3.0e-15_dp)
    r = rotation_from_davenport(axes,'i',angles)
    ref = rotation_from_euler('XYZ',angles)
    call same_rotation(r,ref,3.0e-15_dp)
    r = rotation_from_davenport([0.0_dp,0.0_dp,3.0_dp],'e',90.0_dp,degrees=.true.,status=status)
    if (status /= 0) error stop 'single axis constructor'
    call close(r%apply([1.0_dp,0.0_dp,0.0_dp]),[0.0_dp,1.0_dp,0.0_dp],4.0e-15_dp)
    r = rotation_from_davenport(axes(1:2,:),'e',angles(1:2),status=status)
    ref = rotation_from_euler('xy',angles(1:2))
    call same_rotation(r,ref,3.0e-15_dp)
    ! Scaled axes, including dimensions far outside the square-safe range.
    scaled(1,:) = axes(1,:)*1.0e200_dp
    scaled(2,:) = axes(2,:)*1.0e-200_dp
    scaled(3,:) = axes(3,:)*(-3.0_dp)
    r = rotation_from_davenport(scaled,'e',angles,status=status)
    if (status /= 0) error stop 'scaled axes invalid'
    ref = rotation_from_davenport(axes,'e',[angles(1),angles(2),-angles(3)])
    call same_rotation(r,ref,4.0e-15_dp)
    ! Gimbal lock, third angle set to zero, while reconstructing SO(3).
    r = rotation_from_davenport(axes,'e',[0.6_dp,pi_over_2(),-1.2_dp])
    recovered = r%as_davenport(axes,'e',gimbal_lock=locked,status=status)
    if (status /= 0 .or. .not. locked .or. abs(recovered(3)) > 5.0e-14_dp) &
        error stop 'Davenport gimbal-lock convention'
    back = rotation_from_davenport(axes,'e',recovered)
    call same_rotation(r,back,4.0e-14_dp)
    ! Vectorized stack and grid accessors.
    samples(1,:) = [0.2_dp, 0.4_dp, -0.3_dp]
    samples(2,:) = [-0.8_dp, 0.1_dp, 0.9_dp]
    samples(3,:) = [0.5_dp, -0.2_dp, -0.7_dp]
    stack = rotation_stack_from_davenport(axes,'i',samples,status=status)
    if (status /= 0 .or. .not. stack%valid() .or. stack%size() /= 3) error stop 'Davenport stack invalid'
    batch = stack%as_davenport(axes,'i',gimbal_lock=flags,status=status)
    if (status /= 0 .or. any(flags)) error stop 'Davenport stack accessor'
    do i = 1, 3
        ref = rotation_from_davenport(axes,'i',samples(i,:))
        back = rotation_from_davenport(axes,'i',batch(i,:))
        call same_rotation(ref,back,5.0e-14_dp)
    end do
    lattice(1,1,:) = samples(1,:)
    lattice(2,1,:) = samples(2,:)
    lattice(1,2,:) = samples(3,:)
    lattice(2,2,:) = 0.0_dp
    grid = rotation_grid_from_davenport(axes,'e',lattice,status=status)
    if (status /= 0 .or. .not. grid%valid()) error stop 'Davenport grid invalid'
    all_angles = grid%as_davenport(axes,'e',gimbal_lock=grid_flags,status=status)
    if (status /= 0) error stop 'Davenport grid accessor'
    do j=1,2
        do i=1,2
            ref = rotation_from_davenport(axes,'e',lattice(i,j,:))
            back = rotation_from_davenport(axes,'e',all_angles(i,j,:))
            call same_rotation(ref,back,6.0e-14_dp)
        end do
    end do
    ! Reject invalid sequence order, zero/nonfinite axes, and nonorthogonal axes.
    r = rotation_from_davenport(axes,'other',angles,status=status)
    if (status /= spatial_rotation_invalid_input .or. r%valid()) error stop 'invalid order accepted'
    scaled = axes
    scaled(2,:) = [1.0_dp,1.0_dp,0.0_dp]
    r = rotation_from_davenport(scaled,'e',angles,status=status)
    if (status /= 1 .or. r%valid()) error stop 'nonperpendicular axes accepted'
    r = rotation_from_davenport(axes,'e',[0.1_dp,0.2_dp],status=status)
    if (status /= 1 .or. r%valid()) error stop 'wrong angle count accepted'
    scaled = axes
    scaled(3,:) = 0.0_dp
    r = rotation_from_davenport(scaled,'e',angles,status=status)
    if (status /= 1 .or. r%valid()) error stop 'zero axis accepted'
    scaled = axes
    scaled(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    r = rotation_from_davenport(scaled,'e',angles,status=status)
    if (status /= 1 .or. r%valid()) error stop 'NaN axis accepted'
    r = rotation_from_davenport(axes,'e',angles)
    recovered = r%as_davenport(axes(1:2,:),'e',status=status)
    if (status /= 1) error stop 'invalid extraction axes accepted'
    bad_stack = rotation_stack_from_davenport(scaled,'e',samples,status=status)
    if (status /= 1 .or. bad_stack%valid()) error stop 'invalid stack accepted'
    bad_grid = rotation_grid_from_davenport(scaled,'e',lattice,status=status)
    if (status /= 1 .or. bad_grid%valid()) error stop 'invalid grid accepted'
    print '(a)', 'test_spatial_rotation_davenport PASS'
contains
    pure function pi_over_2() result(value)
        real(dp) :: value
        value = acos(-1.0_dp)/2.0_dp
    end function pi_over_2
    subroutine close(a,b,tol)
        real(dp), intent(in) :: a(:) !! actual values
        real(dp), intent(in) :: b(:) !! expected values
        real(dp), intent(in) :: tol !! comparison tolerance
        if (size(a) /= size(b)) error stop 'Davenport shape mismatch'
        if (.not. all(abs(a-b) <= tol)) then
            print *, a, b
            error stop 'Davenport angle mismatch'
        end if
    end subroutine close
    subroutine same_rotation(a,b,tol)
        type(rotation_type), intent(in) :: a !! reference rotation
        type(rotation_type), intent(in) :: b !! compared rotation
        real(dp), intent(in) :: tol !! quaternion coordinate tolerance
        real(dp) :: x(4), y(4)
        if (.not. a%valid() .or. .not. b%valid()) error stop 'invalid rotation'
        x = a%as_quat()
        y = b%as_quat()
        if (min(maxval(abs(x-y)),maxval(abs(x+y))) > tol) then
            print *, x, y
            error stop 'Davenport rotation mismatch'
        end if
    end subroutine same_rotation
end program test_spatial_rotation_davenport
