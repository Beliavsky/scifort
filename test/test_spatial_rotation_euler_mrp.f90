! SPDX-License-Identifier: MIT
program test_spatial_rotation_euler_mrp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_euler, &
        rotation_from_mrp, rotation_from_quat, spatial_rotation_success, spatial_rotation_invalid_input
    implicit none
    type(rotation_type) :: r, recovered, invalid
    real(dp) :: angles(3), mrp(3), pi
    logical :: gimbal
    integer :: status, i
    character(len=3), parameter :: axes(24) = [character(len=3) :: &
        'xyz','xzy','yxz','yzx','zxy','zyx', &
        'xyx','xzx','yxy','yzy','zxz','zyz', &
        'XYZ','XZY','YXZ','YZX','ZXY','ZYX', &
        'XYX','XZX','YXY','YZY','ZXZ','ZYZ']

    pi = acos(-1.0_dp)
    r = rotation_from_euler('xyz', [90.0_dp,0.0_dp,0.0_dp], &
                            degrees=.true.,status=status)
    if (status /= spatial_rotation_success) error stop 'Euler constructor'
    call close(r%apply([0.0_dp,1.0_dp,0.0_dp]), &
               [0.0_dp,0.0_dp,1.0_dp],3.0e-15_dp)
    call close(r%as_euler('xyz',degrees=.true.), &
               [90.0_dp,0.0_dp,0.0_dp],2.0e-13_dp)

    do i=1,24
        angles = [0.21_dp,-0.63_dp,1.28_dp]
        r = rotation_from_euler(axes(i),angles,status=status)
        if (status /= 0) error stop 'Euler seq constructor'
        angles = r%as_euler(axes(i),gimbal_lock=gimbal,status=status)
        if (status /= 0 .or. gimbal) error stop 'Euler extraction'
        recovered = rotation_from_euler(axes(i),angles,status=status)
        if (status /= 0) error stop 'Euler reconstructed'
        call same_rotation(r,recovered,4.0e-15_dp)
        if (axes(i)(1:1) == axes(i)(3:3)) then
            angles = [0.44_dp,pi,0.71_dp]
        else
            angles = [0.44_dp,0.5_dp*pi,0.71_dp]
        end if
        r = rotation_from_euler(axes(i),angles,status=status)
        angles = r%as_euler(axes(i),gimbal_lock=gimbal,status=status)
        if (status /= 0 .or. .not. gimbal) error stop 'Euler gimbal lock'
        if (abs(angles(3)) > 1.0e-12_dp) error stop 'Euler gimbal last angle'
        recovered = rotation_from_euler(axes(i),angles)
        call same_rotation(r,recovered,3.0e-14_dp)
    end do
    r = rotation_from_euler('z',[pi/2.0_dp],status=status)
    if (status /= 0) error stop 'one Euler axis'
    call close(r%apply([1.0_dp,0.0_dp,0.0_dp]), &
               [0.0_dp,1.0_dp,0.0_dp],2.0e-15_dp)
    r = rotation_from_euler('ZX',[0.0_dp,pi/2.0_dp],status=status)
    if (status /= 0) error stop 'two Euler axes'
    call close(r%apply([0.0_dp,1.0_dp,0.0_dp]), &
               [0.0_dp,0.0_dp,1.0_dp],2.0e-15_dp)
    invalid = rotation_from_euler('xYz',[0.0_dp,0.0_dp,0.0_dp],status=status)
    if (status /= spatial_rotation_invalid_input .or. invalid%valid()) error stop 'mixed case accepted'
    invalid = rotation_from_euler('xxy',[0.0_dp,0.0_dp,0.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'adjacent repeated axis accepted'
    invalid = rotation_from_euler('xyz',[0.0_dp,0.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'Euler wrong angles size'
    invalid = rotation_from_euler('xyz',[0.0_dp, ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp], &
                                  status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'Euler NaN accepted'
    angles = r%as_euler('xxy',status=status)
    if (status /= 1 .or. .not. all(ieee_is_nan(angles))) error stop 'Euler invalid seq extraction'
    angles = r%as_euler('xy',status=status)
    if (status /= 1 .or. .not. all(ieee_is_nan(angles))) error stop 'Euler wrong extraction seq'

    mrp = [0.0_dp,0.0_dp,0.4142135623730950488_dp]
    r = rotation_from_mrp(mrp,status=status)
    if (status /= 0 .or. .not. r%valid()) error stop 'MRP constructor'
    call close(r%apply([1.0_dp,0.0_dp,0.0_dp]), &
               [0.0_dp,1.0_dp,0.0_dp],5.0e-15_dp)
    call close(r%as_mrp(),mrp,2.0e-15_dp)
    do i=1,9
        mrp = [(2.0_dp**(i-6)), -0.5_dp*(2.0_dp**(i-6)), 0.3_dp*(2.0_dp**(i-6))]
        r = rotation_from_mrp(mrp,status=status)
        if (status /= 0) error stop 'MRP scaled constructor'
        recovered = rotation_from_mrp(r%as_mrp(),status=status)
        call same_rotation(r,recovered,5.0e-15_dp)
        if (sqrt(sum(r%as_mrp()**2)) > 1.0_dp+1.0e-14_dp) &
            error stop 'MRP not in principal shadow set'
    end do
    mrp = [1.0e200_dp,-2.0e200_dp,3.0e200_dp]
    r = rotation_from_mrp(mrp,status=status)
    if (status /= 0) error stop 'very large MRP invalid'
    call same_rotation(r,rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,-1.0_dp]),1.0e-185_dp)
    invalid = rotation_from_mrp([1.0_dp,2.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'bad MRP size'
    invalid = rotation_from_mrp([1.0_dp, ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'nonfinite MRP'
    if (.not. all(ieee_is_nan(invalid%as_mrp()))) error stop 'invalid MRP accessor'
    print '(a)', 'test_spatial_rotation_euler_mrp PASS'
contains
    subroutine close(a,b,tol)
        real(dp), intent(in) :: a(:) !! actual values
        real(dp), intent(in) :: b(:) !! target values
        real(dp), intent(in) :: tol !! absolute tolerance
        if (size(a) /= size(b)) error stop 'shape'
        if (.not. all(abs(a-b) < tol)) then
            print *, a,b
            error stop 'Euler/MRP mismatch'
        end if
    end subroutine close
    subroutine same_rotation(a,b,tol)
        type(rotation_type), intent(in) :: a !! first rotation
        type(rotation_type), intent(in) :: b !! second rotation
        real(dp), intent(in) :: tol !! error tolerance
        real(dp) :: x(4),y(4)
        x = a%as_quat()
        y = b%as_quat()
        if (.not. a%valid() .or. .not. b%valid()) error stop 'invalid rotation'
        if (min(maxval(abs(x-y)),maxval(abs(x+y))) > tol) then
            print *, x,y
            error stop 'reconstructed quaternion mismatch'
        end if
    end subroutine same_rotation
end program test_spatial_rotation_euler_mrp
