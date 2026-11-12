! SPDX-License-Identifier: MIT
program test_spatial_rigid_transform
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, rotation_from_euler, rotation_from_quat
    use scifort_spatial_rigid_transform
    implicit none
    type(rotation_type) :: r, rotations(3)
    type(rigid_transform_type) :: tf, tf2, tf3, mean_tf, tmp_tf, unit_tf
    type(rigid_transform_stack_type) :: stack, singleton, broadcasted
    real(dp) :: m(4,4), points(4,3), trans(3,3), exp6(6), dq(8), out(3)
    real(dp), allocatable :: arr(:,:), mats(:,:,:)
    integer :: i, status

    unit_tf = rigid_identity()
    r = rotation_from_euler('z',[90.0_dp],degrees=.true.)
    tf = rigid_from_components([2.0_dp,-1.0_dp,3.0_dp],r,status)
    if (status/=0 .or. .not. tf%valid()) error stop 'construction'
    call check3(tf%apply([1.0_dp,0.0_dp,0.0_dp]),[2.0_dp,0.0_dp,3.0_dp],1.0e-13_dp)
    call check3(tf%apply(tf%apply([1.0_dp,2.0_dp,3.0_dp]),inverse=.true.), &
        [1.0_dp,2.0_dp,3.0_dp],1.0e-13_dp)
    tf2 = tf%inverse()
    tf3 = tf%compose(tf2)
    call check_matrix(tf3%as_matrix(),unit_tf%as_matrix(),1.0e-13_dp)
    tf3 = tf%compose(rigid_from_translation([0.2_dp,0.4_dp,0.6_dp]))
    call check3(tf3%apply([1.0_dp,2.0_dp,3.0_dp]), &
        tf%apply([1.2_dp,2.4_dp,3.6_dp]),1.0e-12_dp)

    m = tf%as_matrix()
    tf2 = rigid_from_matrix(m,status)
    if (status/=0) error stop 'matrix input'
    call check_matrix(tf2%as_matrix(),m,1.0e-13_dp)
    if (maxval(abs(m(4,:)-[0.0_dp,0.0_dp,0.0_dp,1.0_dp])) > 0.0_dp) error stop 'matrix bottom row'
    exp6 = tf%as_exp_coords()
    tf2 = rigid_from_exp_coords(exp6,status)
    if (status/=0) error stop 'exp coords'
    call check_matrix(tf2%as_matrix(),m,1.0e-12_dp)
    tf2 = rigid_from_exp_coords([0.0_dp,0.0_dp,0.0_dp,1.0_dp,-2.0_dp,3.0_dp])
    call check3(tf2%translation(),[1.0_dp,-2.0_dp,3.0_dp],1.0e-14_dp)
    tf2 = rigid_from_exp_coords([1.0e-12_dp,2.0e-12_dp,-3.0e-12_dp,0.5_dp,0.4_dp,-0.3_dp])
    tf3 = rigid_from_exp_coords(tf2%as_exp_coords())
    call check_matrix(tf2%as_matrix(),tf3%as_matrix(),1.0e-13_dp)
    tf2 = rigid_from_exp_coords([3.141592653589793_dp-1.0e-7_dp,0.0_dp,0.0_dp, &
        1.0_dp,2.0_dp,3.0_dp])
    tf3 = rigid_from_exp_coords(tf2%as_exp_coords())
    call check_matrix(tf2%as_matrix(),tf3%as_matrix(),5.0e-12_dp)

    dq = tf%as_dual_quat()
    tf2 = rigid_from_dual_quat(dq,status=status)
    if (status/=0) error stop 'dual-quaternion input'
    call check_matrix(tf2%as_matrix(),m,1.0e-13_dp)
    dq = tf%as_dual_quat(scalar_first=.true.)
    tf2 = rigid_from_dual_quat(dq,scalar_first=.true.)
    call check_matrix(tf2%as_matrix(),m,1.0e-13_dp)
    tf2 = rigid_from_dual_quat(5.0_dp*dq,scalar_first=.true.)
    call check_matrix(tf2%as_matrix(),m,1.0e-13_dp)
    tf2 = tf%power(0.5_dp)
    tf3 = tf2%compose(tf2)
    call check_matrix(tf3%as_matrix(),m,1.0e-12_dp)
    tf2 = tf%power(-1.0_dp)
    tf3 = tf%inverse()
    call check_matrix(tf2%as_matrix(),tf3%as_matrix(),1.0e-13_dp)
    tf2 = tf%power(0.0_dp)
    call check_matrix(tf2%as_matrix(),unit_tf%as_matrix(),1.0e-13_dp)

    points(1,:) = [0.0_dp,0.0_dp,0.0_dp]
    points(2,:) = [1.0_dp,0.0_dp,0.0_dp]
    points(3,:) = [2.0_dp,3.0_dp,4.0_dp]
    points(4,:) = [-1.0_dp,0.5_dp,0.7_dp]
    arr = tf%apply(points)
    do i=1,4
        call check3(arr(i,:),tf%apply(points(i,:)),1.0e-14_dp)
    end do

    trans(1,:) = [0.0_dp,0.0_dp,0.0_dp]
    trans(2,:) = [1.0_dp,2.0_dp,3.0_dp]
    trans(3,:) = [-1.0_dp,0.0_dp,1.0_dp]
    rotations(1) = rotation_identity()
    rotations(2) = r
    rotations(3) = rotation_from_euler('x',[45.0_dp],degrees=.true.)
    stack = rigid_stack_from_components(trans,rotations,status)
    if (status/=0 .or. .not. stack%valid() .or. stack%size()/=3) error stop 'stack construct'
    mats = stack%as_matrix()
    singleton = rigid_stack_from_components(trans(1:1,:),rotations(1:1))
    if (.not. singleton%valid()) error stop 'singleton'
    broadcasted = stack%compose(singleton)
    if (.not. broadcasted%valid()) error stop 'broadcast compose'
    broadcasted = singleton%compose(stack)
    if (.not. broadcasted%valid()) error stop 'broadcast inverse compose'
    broadcasted = stack%inverse()
    if (.not. broadcasted%valid()) error stop 'inverse stack'
    do i=1,3
        tf2 = stack%at(i)
        call check_matrix(tf2%as_matrix(),mats(i,:,:),1.0e-14_dp)
        tmp_tf = broadcasted%at(i)
        tf3 = tmp_tf%compose(tf2)
        call check_matrix(tf3%as_matrix(),unit_tf%as_matrix(),1.0e-13_dp)
    end do
    broadcasted = rigid_stack_from_matrix(mats,status)
    if (status/=0) error stop 'matrix stack'
    arr = broadcasted%as_exp_coords()
    broadcasted = rigid_stack_from_exp_coords(arr,status)
    if (status/=0) error stop 'exp stack'
    tmp_tf = broadcasted%at(2)
    call check_matrix(tmp_tf%as_matrix(),mats(2,:,:),1.0e-12_dp)
    arr = broadcasted%as_dual_quat()
    broadcasted = rigid_stack_from_dual_quat(arr,status=status)
    if (status/=0) error stop 'dual stack'
    tmp_tf = broadcasted%at(3)
    call check_matrix(tmp_tf%as_matrix(),mats(3,:,:),1.0e-12_dp)
    arr = stack%apply([1.0_dp,0.0_dp,0.0_dp])
    do i=1,3
        tmp_tf = stack%at(i)
        call check3(arr(i,:),tmp_tf%apply([1.0_dp,0.0_dp,0.0_dp]),1.0e-13_dp)
    end do
    arr = stack%apply(trans)
    do i=1,3
        tmp_tf = stack%at(i)
        call check3(arr(i,:),tmp_tf%apply(trans(i,:)),1.0e-13_dp)
    end do
    mean_tf = stack%mean(weights=[0.0_dp,1.0_dp,0.0_dp],status=status)
    if (status/=0) error stop 'mean status'
    call check_matrix(mean_tf%as_matrix(),mats(2,:,:),1.0e-12_dp)
    mean_tf = stack%mean(weights=[0.25_dp,0.5_dp,0.25_dp])
    if (.not. mean_tf%valid()) error stop 'weighted mean'
    call check3(mean_tf%translation(),[0.25_dp,1.0_dp,1.75_dp],1.0e-13_dp)

    tf2 = rigid_from_translation([0.0_dp,1.0_dp],status)
    if (status==0 .or. tf2%valid()) error stop 'short translation'
    tf2 = rigid_from_exp_coords([1.0_dp,2.0_dp,3.0_dp],status)
    if (status==0 .or. tf2%valid()) error stop 'short exp'
    dq=0.0_dp
    tf2 = rigid_from_dual_quat(dq,status=status)
    if (status==0 .or. tf2%valid()) error stop 'zero dual quat'
    m=unit_tf%as_matrix()
    m(4,1)=1.0_dp
    tf2=rigid_from_matrix(m,status)
    if (status==0 .or. tf2%valid()) error stop 'invalid bottom row'
    m=unit_tf%as_matrix()
    m(1,1)=2.0_dp
    tf2=rigid_from_matrix(m,status)
    if (status==0 .or. tf2%valid()) error stop 'nonrotation input'
    tf2=rigid_from_translation([ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,0.0_dp],status)
    if (status==0 .or. tf2%valid()) error stop 'nan translation'
    tf2=tf%power(ieee_value(0.0_dp,ieee_quiet_nan),status)
    if (status==0 .or. tf2%valid()) error stop 'nan exponent'
    stack=rigid_stack_from_components(trans,rotations(1:2),status)
    if (status==0 .or. stack%valid()) error stop 'bad broadcast'
    stack=rigid_stack_from_components(trans,rotations)
    mean_tf=stack%mean(weights=[0.0_dp,0.0_dp,0.0_dp],status=status)
    if (status==0 .or. mean_tf%valid()) error stop 'zero weights'
    broadcasted = stack%compose(rigid_stack_from_components(trans(1:2,:),rotations(1:2)))
    if (broadcasted%valid()) error stop 'invalid compose broadcast'
    tmp_tf = stack%at(0)
    if (tmp_tf%valid()) error stop 'invalid index'
    print '(A)', 'PASS: rigid transform SE(3), exp/log, dual quaternions, stacks, errors'
contains
    subroutine check3(a,b,tol)
        real(dp), intent(in) :: a(:) !! calculated vector
        real(dp), intent(in) :: b(:) !! expected vector
        real(dp), intent(in) :: tol !! absolute tolerance
        if (size(a)/=3 .or. size(b)/=3) error stop 'check3 shape'
        if (.not. all(abs(a-b) <= tol)) then
            print *, 'actual=',a, ' expected=',b
            error stop 'check3 values'
        end if
    end subroutine check3
    subroutine check_matrix(a,b,tol)
        real(dp), intent(in) :: a(:,:) !! calculated matrix
        real(dp), intent(in) :: b(:,:) !! expected matrix
        real(dp), intent(in) :: tol !! absolute tolerance
        if (any(shape(a)/=[4,4]) .or. any(shape(b)/=[4,4])) error stop 'matrix shape'
        if (.not. all(abs(a-b) <= tol)) then
            print *, 'max error=',maxval(abs(a-b))
            error stop 'matrix values'
        end if
    end subroutine check_matrix
end program test_spatial_rigid_transform
