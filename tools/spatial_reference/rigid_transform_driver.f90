! SPDX-License-Identifier: MIT
! Streaming numerical reference driver for SciPy RigidTransform comparisons.
program rigid_transform_driver
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_quat
    use scifort_spatial_rigid_transform
    implicit none
    type(rotation_type) :: ra,rb,rot(2)
    type(rigid_transform_type) :: a,b,c,avg
    type(rigid_transform_stack_type) :: stack,wide
    real(dp) :: q1(4), q2(4), t1(3), t2(3), point(3), exponent
    real(dp) :: translations(2,3), weights(2), du(8)
    real(dp) :: m(4,4)
    integer :: i,n,code
    read(*,*) n
    do i=1,n
        read(*,*) q1,t1,q2,t2,point,exponent
        ra = rotation_from_quat(q1)
        rb = rotation_from_quat(q2)
        a = rigid_from_components(t1,ra)
        b = rigid_from_components(t2,rb)
        if (.not. a%valid() .or. .not. b%valid()) error stop 'invalid input'
        c = a%compose(b)
        write(*, '(*(ES25.16E3,1X))') a%as_matrix()
        write(*, '(*(ES25.16E3,1X))') a%as_exp_coords()
        write(*, '(*(ES25.16E3,1X))') a%as_dual_quat()
        write(*, '(*(ES25.16E3,1X))') c%as_matrix()
        c = a%inverse()
        write(*, '(*(ES25.16E3,1X))') c%as_matrix()
        c = a%power(exponent,code)
        if (code/=0) error stop 'power'
        write(*, '(*(ES25.16E3,1X))') c%as_matrix()
        write(*, '(*(ES25.16E3,1X))') a%apply(point),a%apply(point,inverse=.true.)
        c = rigid_from_exp_coords(a%as_exp_coords(),code)
        if (code/=0) error stop 'exp input'
        write(*, '(*(ES25.16E3,1X))') c%as_matrix()
        du = a%as_dual_quat()
        c = rigid_from_dual_quat(du,status=code)
        if (code/=0) error stop 'dual input'
        write(*, '(*(ES25.16E3,1X))') c%as_matrix()
        m = a%as_matrix()
        c = rigid_from_matrix(m,status=code)
        if (code/=0) error stop 'matrix input'
        write(*, '(*(ES25.16E3,1X))') c%as_matrix()
        translations(1,:) = t1
        translations(2,:) = t2
        rot = [ra,rb]
        stack = rigid_stack_from_components(translations,rot,status=code)
        if (code/=0) error stop 'stack input'
        weights = [1.0_dp,2.0_dp]
        avg = stack%mean(weights=weights,status=code)
        if (code/=0) error stop 'stack mean'
        write(*, '(*(ES25.16E3,1X))') avg%as_matrix()
        write(*, '(*(ES25.16E3,1X))') stack%apply(point)
        wide = stack%compose(stack%inverse())
        if (.not. wide%valid()) error stop 'stack compose'
        write(*, '(*(ES25.16E3,1X))') wide%as_matrix()
        wide = rigid_stack_from_exp_coords(stack%as_exp_coords())
        if (.not. wide%valid()) error stop 'stack exp'
        write(*, '(*(ES25.16E3,1X))') wide%as_matrix()
    end do
end program rigid_transform_driver
