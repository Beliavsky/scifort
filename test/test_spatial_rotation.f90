! SPDX-License-Identifier: MIT
program test_spatial_rotation
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, &
        rotation_from_quat, rotation_from_rotvec, rotation_from_matrix, &
        spatial_rotation_success, spatial_rotation_invalid_input
    implicit none
    type(rotation_type) :: r, s, combined, unit, invalid
    real(dp) :: m(3,3), v(3), q(4), theta, big, tiny
    integer :: status, k

    unit = rotation_identity()
    call close(unit%apply([2.0_dp,3.0_dp,4.0_dp]), &
        [2.0_dp,3.0_dp,4.0_dp],0.0_dp)
    call close(unit%as_rotvec(),[0.0_dp,0.0_dp,0.0_dp],0.0_dp)
    theta = acos(-1.0_dp)/2.0_dp
    r = rotation_from_rotvec([0.0_dp,0.0_dp,theta],status=status)
    if (status /= spatial_rotation_success .or. .not. r%valid()) &
        error stop 'construct from rotvec'
    call close(r%apply([1.0_dp,0.0_dp,0.0_dp]), &
        [0.0_dp,1.0_dp,0.0_dp],2.0e-15_dp)
    call close(r%as_rotvec(),[0.0_dp,0.0_dp,theta],2.0e-15_dp)
    call close(r%as_rotvec(degrees=.true.),[0.0_dp,0.0_dp,90.0_dp],2.0e-13_dp)
    call close(r%apply([3.0_dp,4.0_dp,5.0_dp],inverse=.true.), &
        [4.0_dp,-3.0_dp,5.0_dp],3.0e-15_dp)
    s = r%inverse()
    call close(s%as_quat(),[0.0_dp,0.0_dp, &
        -sin(theta/2.0_dp),cos(theta/2.0_dp)],2.0e-15_dp)

    m = r%as_matrix()
    call close(reshape(m,[9]),reshape(transpose(transpose(m)),[9]),0.0_dp)
    s = rotation_from_matrix(m,status=status)
    if (status /= 0) error stop 'matrix conversion status'
    call close(s%apply([1.0_dp,0.0_dp,0.0_dp]), &
        [0.0_dp,1.0_dp,0.0_dp],2.0e-15_dp)
    call close(reshape(matmul(transpose(m),m),[9]), &
        reshape(identity(),[9]),2.0e-15_dp)
    if (abs(s%magnitude()-theta)>2.0e-15_dp) error stop 'magnitude'

    s = rotation_from_rotvec([90.0_dp,0.0_dp,0.0_dp],degrees=.true.)
    combined = r%compose(s)
    v = [0.5_dp,1.0_dp,-0.2_dp]
    call close(combined%apply(v),r%apply(s%apply(v)),2.0e-15_dp)
    s = r%inverse()
    call close(s%apply(r%apply(v)),v,2.0e-15_dp)
    combined = combined%compose(combined%inverse())
    call close(combined%as_quat(canonical=.true.), &
        [0.0_dp,0.0_dp,0.0_dp,1.0_dp],2.0e-15_dp)

    q = [1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    s = rotation_from_quat(q, status=status)
    if (status /= 0) error stop 'pi rotation quaternion'
    call close(s%as_rotvec(),[acos(-1.0_dp),0.0_dp,0.0_dp],3.0e-15_dp)
    s = rotation_from_matrix(s%as_matrix(),status=status)
    if (status /= 0) error stop 'matrix pi conversion'
    call close(s%apply([0.0_dp,3.0_dp,0.0_dp]), &
        [0.0_dp,-3.0_dp,0.0_dp],2.0e-15_dp)

    s = rotation_from_quat([-2.0_dp,0.0_dp,0.0_dp,0.0_dp])
    call close(s%as_quat(canonical=.true.),[1.0_dp,0.0_dp,0.0_dp,0.0_dp],0.0_dp)
    s = rotation_from_quat([3.0_dp,1.0_dp,0.0_dp,0.0_dp],scalar_first=.true.)
    call close(s%as_quat(scalar_first=.true.), &
        [3.0_dp,1.0_dp,0.0_dp,0.0_dp]/sqrt(10.0_dp),1.0e-15_dp)

    ! Quaternion magnitude normalization is stable for arbitrary scaling.
    big = 1.0e300_dp
    tiny = 1.0e-300_dp
    s = rotation_from_quat(big*r%as_quat(),status=status)
    if (status /= 0) error stop 'large quaternion'
    call close(s%apply(v),r%apply(v),2.0e-15_dp)
    s = rotation_from_quat(tiny*r%as_quat(),status=status)
    if (status /= 0) error stop 'small quaternion'
    call close(s%apply(v),r%apply(v),2.0e-15_dp)

    do k=1,3
        v = 0.0_dp
        v(k) = 1.0e-12_dp
        s = rotation_from_rotvec(v,status=status)
        if (status /= 0) error stop 'small rotvec'
        call close(s%as_rotvec(),v,1.0e-26_dp)
    end do

    invalid = rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,0.0_dp],status=status)
    if (status /= spatial_rotation_invalid_input .or. invalid%valid()) &
        error stop 'zero quaternion invalid'
    if (.not. all(ieee_is_nan(invalid%as_quat()))) error stop 'invalid quaternion output'
    invalid = rotation_from_quat([1.0_dp,0.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'wrong quaternion shape'
    invalid = rotation_from_rotvec([1.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'wrong rotvec shape'
    invalid = rotation_from_rotvec([1.0_dp, &
        ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'nonfinite rotvec'
    m = identity()
    m(1,1) = -1.0_dp
    invalid = rotation_from_matrix(m,status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'reflection invalid'
    m = 2.0_dp*identity()
    invalid = rotation_from_matrix(m,status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'nonrotation invalid'
    m = identity()
    invalid = rotation_from_matrix(m(:2,:),status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'wrong matrix shape'

    print '(a)', 'test_spatial_rotation PASS'
contains
    subroutine close(a,b,tol)
        real(dp), intent(in) :: a(:) !! actual values
        real(dp), intent(in) :: b(:) !! expected values
        real(dp), intent(in) :: tol !! absolute tolerance
        if (size(a) /= size(b)) error stop 'comparison shape'
        if (.not. all(abs(a-b) <= tol)) then
            print *, a,b
            error stop 'spatial rotation value discrepancy'
        end if
    end subroutine close

    function identity() result(m)
        real(dp) :: m(3,3)
        integer :: i
        m = 0.0_dp
        do i=1,3
            m(i,i) = 1.0_dp
        end do
    end function identity
end program test_spatial_rotation
