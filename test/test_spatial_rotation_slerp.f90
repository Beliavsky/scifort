! SPDX-License-Identifier: MIT
program test_spatial_rotation_slerp
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, &
        rotation_from_quat, rotation_from_rotvec, rotation_slerp, rotation_slerp_path
    implicit none
    type(rotation_type) :: r1, r2, r3, bad, knots(4), temp
    type(rotation_type), allocatable :: result(:)
    real(dp) :: q(4), t(4), qtime(7), pi
    integer :: code, j
    pi = acos(-1.0_dp)

    r1 = rotation_identity()
    r2 = rotation_from_rotvec([0.0_dp,0.0_dp,pi])
    r3 = rotation_slerp(r1,r2,0.5_dp,status=code)
    if (code /= 0 .or. .not. r3%valid()) error stop 'slerp status'
    call close_vector(r3%apply([1.0_dp,0.0_dp,0.0_dp]), &
        [0.0_dp,1.0_dp,0.0_dp],3.0e-15_dp)
    r3 = rotation_slerp(r1,r2,0.0_dp)
    call close(r3%as_matrix(),r1%as_matrix(),0.0_dp)
    r3 = rotation_slerp(r1,r2,1.0_dp)
    call close(r3%as_matrix(),r2%as_matrix(),2.0e-15_dp)
    r3 = rotation_slerp(r1,r2,-0.5_dp)
    call close_vector(r3%apply([1.0_dp,0.0_dp,0.0_dp]), &
        [0.0_dp,-1.0_dp,0.0_dp],3.0e-15_dp)

    ! Reversing the quaternion sign must not change the physical path.
    q = -r2%as_quat()
    temp = rotation_from_quat(q)
    r3 = rotation_slerp(r1,temp,0.5_dp)
    call close_vector(r3%apply([1.0_dp,0.0_dp,0.0_dp]), &
        [0.0_dp,1.0_dp,0.0_dp],3.0e-15_dp)

    ! Rotation sequences use nonuniform knots; queries need not be sorted.
    t = [-2.0_dp,0.0_dp,3.0_dp,10.0_dp]
    knots(1) = rotation_from_rotvec([0.0_dp,0.0_dp,0.0_dp])
    knots(2) = rotation_from_rotvec([0.0_dp,0.0_dp,pi/3.0_dp])
    knots(3) = rotation_from_rotvec([0.0_dp,0.0_dp,2.0_dp*pi/3.0_dp])
    knots(4) = rotation_from_rotvec([0.0_dp,0.0_dp,pi])
    qtime = [10.0_dp,-2.0_dp,1.5_dp,3.0_dp,-1.0_dp,7.0_dp,0.0_dp]
    call rotation_slerp_path(t,knots,qtime,result,code)
    if (code /= 0 .or. size(result) /= 7) error stop 'slerp path shape'
    do j=1,7
        r3 = rotation_slerp(knots(find_segment(qtime(j),t)), &
                            knots(find_segment(qtime(j),t)+1), &
             (qtime(j)-t(find_segment(qtime(j),t))) / &
             (t(find_segment(qtime(j),t)+1)-t(find_segment(qtime(j),t))))
        call close(result(j)%as_matrix(),r3%as_matrix(),3.0e-15_dp)
    end do

    ! Invalid paths fail atomically and return an empty result.
    call rotation_slerp_path([0.0_dp,0.0_dp],knots(:2),[0.5_dp],result,code)
    if (code /= 1 .or. size(result) /= 0) error stop 'duplicate time'
    call rotation_slerp_path(t,knots,[10.1_dp],result,code)
    if (code /= 1 .or. size(result) /= 0) error stop 'time out of range'
    call rotation_slerp_path(t,knots,[-3.0_dp],result,code)
    if (code /= 1 .or. size(result) /= 0) error stop 'time below range'
    call rotation_slerp_path(t,knots(:3),qtime,result,code)
    if (code /= 1 .or. size(result) /= 0) error stop 'mismatched knots'
    bad = rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,0.0_dp])
    r3 = rotation_slerp(r1,bad,0.4_dp,status=code)
    if (code /= 1 .or. r3%valid()) error stop 'invalid endpoint'
    r3 = rotation_slerp(r1,r2,ieee_value(0.0_dp,ieee_quiet_nan),status=code)
    if (code /= 1 .or. r3%valid()) error stop 'nonfinite fraction'

    print '(a)', 'test_spatial_rotation_slerp PASS'
contains
    integer function find_segment(x,times) result(k)
        real(dp), intent(in) :: x !! query time in [times(1), times(n)]
        real(dp), intent(in) :: times(:) !! ascending knots
        k = 1
        do while (k < size(times)-1)
            if (x < times(k+1)) exit
            k = k+1
        end do
    end function find_segment

    subroutine close_vector(a,b,tol)
        real(dp), intent(in) :: a(:) !! actual vector
        real(dp), intent(in) :: b(:) !! expected vector
        real(dp), intent(in) :: tol !! absolute tolerance
        if (size(a) /= size(b)) error stop 'vector shape'
        if (.not. all(abs(a-b)<=tol)) error stop 'slerp vector discrepancy'
    end subroutine close_vector

    subroutine close(a,b,tol)
        real(dp), intent(in) :: a(:,:) !! actual matrix
        real(dp), intent(in) :: b(:,:) !! expected matrix
        real(dp), intent(in) :: tol !! absolute tolerance
        if (any(shape(a) /= shape(b))) error stop 'comparison shape'
        if (.not. all(abs(a-b) <= tol)) then
            print *, a,b
            error stop 'slerp value discrepancy'
        end if
    end subroutine close
end program test_spatial_rotation_slerp
