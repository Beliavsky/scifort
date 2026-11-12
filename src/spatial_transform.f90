! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent 3D active rotation utilities and piecewise quaternion Slerp.
! Quaternion convention: scalar-last (x,y,z,w), right handed.
! See docs/spatial_transform.md for mathematical definitions and limitations.
! Davenport-angle quaternion half-angle method follows Shuster-Markley and
! Bernardes-Viollet, cross-checked/adapted from SciPy 1.17 (BSD-3-Clause).
! See LICENSES/scipy-spatial-transform-BSD-3-Clause.txt for attribution.
module scifort_spatial_transform
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: spatial_rotation_success = 0
    integer, parameter, public :: spatial_rotation_invalid_input = 1
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: rotation_type
        private
        real(dp) :: q(4) = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
        logical :: defined = .true.
    contains
        procedure :: valid => rotation_valid
        procedure :: as_quat => rotation_as_quat
        procedure :: as_matrix => rotation_as_matrix
        procedure :: as_rotvec => rotation_as_rotvec
        procedure :: as_euler => rotation_as_euler
        procedure :: as_mrp => rotation_as_mrp
        procedure :: as_davenport => rotation_as_davenport
        procedure :: apply => rotation_apply
        procedure :: inverse => rotation_inverse
        procedure :: compose => rotation_compose
        procedure :: magnitude => rotation_magnitude
    end type rotation_type

    public :: rotation_identity, rotation_from_quat, rotation_from_rotvec
    public :: rotation_from_matrix, rotation_from_euler, rotation_from_mrp
    public :: rotation_from_davenport

    interface rotation_from_davenport
        module procedure rotation_from_davenport_axes
        module procedure rotation_from_davenport_axis
    end interface rotation_from_davenport
    public :: rotation_slerp, rotation_slerp_path
    public :: rotation_align_vectors, rotation_mean

contains

    pure function invalid_rotation() result(r)
        type(rotation_type) :: r
        r%q = ieee_value(0.0_dp, ieee_quiet_nan)
        r%defined = .false.
    end function invalid_rotation

    pure subroutine set_status(status, code)
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        integer, intent(in) :: code !! status value to return
        if (present(status)) status = code
    end subroutine set_status

    pure logical function rotation_valid(self) result(ok)
        class(rotation_type), intent(in) :: self !! rotation to check
        ok = self%defined .and. all(ieee_is_finite(self%q))
    end function rotation_valid

    pure function rotation_identity() result(r)
        type(rotation_type) :: r
        r%q = [0.0_dp,0.0_dp,0.0_dp,1.0_dp]
        r%defined = .true.
    end function rotation_identity

    ! Scale before taking the norm: otherwise very small/large input
    ! quaternions can underflow or overflow even when their direction is valid.
    pure subroutine normalize_quaternion(q, success)
        real(dp), intent(inout) :: q(4) !! quaternion, normalized in place
        logical, intent(out) :: success !! true if finite and nonzero
        real(dp) :: scale, mag
        success = .false.
        if (.not. all(ieee_is_finite(q))) return
        scale = maxval(abs(q))
        if (scale <= 0.0_dp) return
        q = q / scale
        mag = sqrt(sum(q*q))
        q = q / mag
        success = .true.
    end subroutine normalize_quaternion

    function rotation_from_quat(quat, scalar_first, status) result(r)
        real(dp), intent(in) :: quat(:) !! four quaternion components; nonzero, finite
        logical, intent(in), optional :: scalar_first !! w,x,y,z order when true; default false
        integer, intent(out), optional :: status !! zero on success; one on invalid input
        type(rotation_type) :: r
        logical :: scalar_head, ok
        r = invalid_rotation()
        call set_status(status, spatial_rotation_invalid_input)
        if (size(quat) /= 4) return
        scalar_head = .false.
        if (present(scalar_first)) scalar_head = scalar_first
        if (scalar_head) then
            r%q = [quat(2:4), quat(1)]
        else
            r%q = quat
        end if
        call normalize_quaternion(r%q, ok)
        if (.not. ok) then
            r = invalid_rotation()
            return
        end if
        r%defined = .true.
        call set_status(status, spatial_rotation_success)
    end function rotation_from_quat

    function rotation_from_rotvec(rotvec, degrees, status) result(r)
        real(dp), intent(in) :: rotvec(:) !! axis multiplied by angle, length three
        logical, intent(in), optional :: degrees !! interpret angle in degrees; default radians
        integer, intent(out), optional :: status !! zero on success; one on invalid input
        type(rotation_type) :: r
        real(dp) :: v(3), scale, angle, sine_factor, x2
        logical :: ok
        r = invalid_rotation()
        call set_status(status, spatial_rotation_invalid_input)
        if (size(rotvec) /= 3) return
        if (.not. all(ieee_is_finite(rotvec))) return
        v = rotvec
        if (present(degrees)) then
            if (degrees) v = v*(pi/180.0_dp)
        end if
        if (.not. all(ieee_is_finite(v))) return
        scale = maxval(abs(v))
        if (scale == 0.0_dp) then
            r = rotation_identity()
            call set_status(status, spatial_rotation_success)
            return
        end if
        angle = scale*sqrt(sum((v/scale)**2))
        if (.not. ieee_is_finite(angle)) return
        if (angle < 1.0e-4_dp) then
            x2 = angle*angle
            sine_factor = 0.5_dp - x2/48.0_dp + x2*x2/3840.0_dp
        else
            sine_factor = sin(0.5_dp*angle)/angle
        end if
        r%q(1:3) = v*sine_factor
        r%q(4) = cos(0.5_dp*angle)
        call normalize_quaternion(r%q, ok)
        if (.not. ok) then
            r = invalid_rotation()
            return
        end if
        r%defined = .true.
        call set_status(status, spatial_rotation_success)
    end function rotation_from_rotvec

    function rotation_from_matrix(matrix, status) result(r)
        real(dp), intent(in) :: matrix(:,:) !! 3x3 proper orthogonal active rotation
        integer, intent(out), optional :: status !! zero on success; one on invalid matrix
        type(rotation_type) :: r
        real(dp) :: m(3,3), e(3,3), trace, s, det
        logical :: ok
        integer :: i
        r = invalid_rotation()
        call set_status(status, spatial_rotation_invalid_input)
        if (size(matrix,1) /= 3 .or. size(matrix,2) /= 3) return
        if (.not. all(ieee_is_finite(matrix))) return
        m = matrix
        e = matmul(transpose(m),m)
        do i=1,3
            e(i,i) = e(i,i)-1.0_dp
        end do
        det = dot_product(m(:,1),cross3(m(:,2),m(:,3)))
        ! A first version accepts rotations, but does not project arbitrary
        ! near-rotation matrices onto SO(3) as SciPy does.
        if (maxval(abs(e)) > 1.0e-8_dp .or. abs(det-1.0_dp) > 1.0e-8_dp) return
        trace = m(1,1)+m(2,2)+m(3,3)
        if (trace > 0.0_dp) then
            s = 2.0_dp*sqrt(max(0.0_dp,1.0_dp+trace))
            r%q(4) = 0.25_dp*s
            r%q(1) = (m(3,2)-m(2,3))/s
            r%q(2) = (m(1,3)-m(3,1))/s
            r%q(3) = (m(2,1)-m(1,2))/s
        else if (m(1,1) >= m(2,2) .and. m(1,1) >= m(3,3)) then
            s = 2.0_dp*sqrt(max(0.0_dp,1.0_dp+m(1,1)-m(2,2)-m(3,3)))
            r%q(4) = (m(3,2)-m(2,3))/s
            r%q(1) = 0.25_dp*s
            r%q(2) = (m(1,2)+m(2,1))/s
            r%q(3) = (m(1,3)+m(3,1))/s
        else if (m(2,2) >= m(3,3)) then
            s = 2.0_dp*sqrt(max(0.0_dp,1.0_dp+m(2,2)-m(1,1)-m(3,3)))
            r%q(4) = (m(1,3)-m(3,1))/s
            r%q(1) = (m(1,2)+m(2,1))/s
            r%q(2) = 0.25_dp*s
            r%q(3) = (m(2,3)+m(3,2))/s
        else
            s = 2.0_dp*sqrt(max(0.0_dp,1.0_dp+m(3,3)-m(1,1)-m(2,2)))
            r%q(4) = (m(2,1)-m(1,2))/s
            r%q(1) = (m(1,3)+m(3,1))/s
            r%q(2) = (m(2,3)+m(3,2))/s
            r%q(3) = 0.25_dp*s
        end if
        call normalize_quaternion(r%q,ok)
        if (.not. ok) then
            r = invalid_rotation()
            return
        end if
        r%defined = .true.
        call set_status(status, spatial_rotation_success)
    end function rotation_from_matrix

    pure function cross3(a,b) result(c)
        real(dp), intent(in) :: a(3) !! first vector
        real(dp), intent(in) :: b(3) !! second vector
        real(dp) :: c(3)
        c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), &
             a(1)*b(2)-a(2)*b(1)]
    end function cross3

    pure function rotation_as_quat(self, canonical, scalar_first) result(quat)
        class(rotation_type), intent(in) :: self !! input rotation
        logical, intent(in), optional :: canonical !! choose deterministic q/-q sign; default false
        logical, intent(in), optional :: scalar_first !! put scalar component first; default false
        real(dp) :: quat(4)
        logical :: fix_sign, scalar_head
        integer :: k
        quat = self%q
        if (.not. self%valid()) then
            quat = ieee_value(0.0_dp,ieee_quiet_nan)
            return
        end if
        fix_sign = .false.
        if (present(canonical)) fix_sign = canonical
        if (fix_sign) then
            if (quat(4) < 0.0_dp) then
                quat = -quat
            else if (quat(4) == 0.0_dp) then
                do k=1,3
                    if (quat(k) == 0.0_dp) cycle
                    if (quat(k) < 0.0_dp) quat = -quat
                    exit
                end do
            end if
        end if
        scalar_head = .false.
        if (present(scalar_first)) scalar_head = scalar_first
        if (scalar_head) quat = [quat(4),quat(1:3)]
    end function rotation_as_quat

    pure function rotation_as_matrix(self) result(m)
        class(rotation_type), intent(in) :: self !! input rotation
        real(dp) :: m(3,3)
        real(dp) :: x, y, z, w
        if (.not. self%valid()) then
            m = ieee_value(0.0_dp,ieee_quiet_nan)
            return
        end if
        x = self%q(1)
        y = self%q(2)
        z = self%q(3)
        w = self%q(4)
        m(1,:) = [1.0_dp-2.0_dp*(y*y+z*z), 2.0_dp*(x*y-z*w), 2.0_dp*(x*z+y*w)]
        m(2,:) = [2.0_dp*(x*y+z*w), 1.0_dp-2.0_dp*(x*x+z*z), 2.0_dp*(y*z-x*w)]
        m(3,:) = [2.0_dp*(x*z-y*w), 2.0_dp*(y*z+x*w), 1.0_dp-2.0_dp*(x*x+y*y)]
    end function rotation_as_matrix

    pure function rotation_as_rotvec(self, degrees) result(v)
        class(rotation_type), intent(in) :: self !! input rotation
        logical, intent(in), optional :: degrees !! return degrees, default radians
        real(dp) :: v(3), q(4), s, angle, factor, s2
        if (.not. self%valid()) then
            v = ieee_value(0.0_dp,ieee_quiet_nan)
            return
        end if
        q = self%q
        if (q(4) < 0.0_dp) q = -q
        s = sqrt(sum(q(1:3)**2))
        angle = 2.0_dp*atan2(s,q(4))
        if (s < 1.0e-4_dp) then
            s2 = s*s
            factor = 2.0_dp + s2/3.0_dp + 3.0_dp*s2*s2/20.0_dp
        else
            factor = angle/s
        end if
        v = q(1:3)*factor
        if (present(degrees)) then
            if (degrees) v = v*(180.0_dp/pi)
        end if
    end function rotation_as_rotvec

    pure function rotation_apply(self, vector, inverse) result(v)
        class(rotation_type), intent(in) :: self !! active 3D rotation
        real(dp), intent(in) :: vector(:) !! three Cartesian vector components
        logical, intent(in), optional :: inverse !! apply inverse rotation; default false
        real(dp) :: v(3), qv(3), t(3), w
        logical :: backward
        if (size(vector) /= 3 .or. .not. self%valid()) then
            v = ieee_value(0.0_dp,ieee_quiet_nan)
            return
        end if
        if (.not. all(ieee_is_finite(vector))) then
            v = ieee_value(0.0_dp,ieee_quiet_nan)
            return
        end if
        backward = .false.
        if (present(inverse)) backward = inverse
        qv = self%q(1:3)
        w = self%q(4)
        if (backward) qv = -qv
        t = 2.0_dp*cross3(qv,vector)
        v = vector + w*t + cross3(qv,t)
    end function rotation_apply

    pure function rotation_inverse(self) result(r)
        class(rotation_type), intent(in) :: self !! input rotation
        type(rotation_type) :: r
        if (.not. self%valid()) then
            r = invalid_rotation()
            return
        end if
        r%q(1:3) = -self%q(1:3)
        r%q(4) = self%q(4)
    end function rotation_inverse

    pure function rotation_compose(self, other) result(r)
        class(rotation_type), intent(in) :: self !! rotation applied second
        type(rotation_type), intent(in) :: other !! rotation applied first
        type(rotation_type) :: r
        real(dp) :: a(4), b(4)
        logical :: ok
        if (.not. self%valid() .or. .not. other%valid()) then
            r = invalid_rotation()
            return
        end if
        a = self%q
        b = other%q
        r%q(1:3) = a(4)*b(1:3) + b(4)*a(1:3) + cross3(a(1:3),b(1:3))
        r%q(4) = a(4)*b(4)-dot_product(a(1:3),b(1:3))
        call normalize_quaternion(r%q, ok)
        if (.not. ok) r = invalid_rotation()
    end function rotation_compose

    pure function rotation_magnitude(self) result(angle)
        class(rotation_type), intent(in) :: self !! input rotation
        real(dp) :: angle
        if (.not. self%valid()) then
            angle = ieee_value(0.0_dp,ieee_quiet_nan)
        else
            angle = 2.0_dp*atan2(sqrt(sum(self%q(1:3)**2)),abs(self%q(4)))
        end if
    end function rotation_magnitude


    ! Validate axis strings without changing SciPy's case-sensitive distinction:
    ! lowercase = extrinsic (fixed axes), uppercase = intrinsic (moving axes).
    pure logical function valid_euler_sequence(seq, full_length) result(ok)
        character(len=*), intent(in) :: seq !! axis string, case-sensitive
        logical, intent(in) :: full_length !! require exactly 3 axes when true
        logical :: lowercase, uppercase
        integer :: i, v
        ok = .false.
        if (len(seq) < 1 .or. len(seq) > 3) return
        if (full_length .and. len(seq) /= 3) return
        lowercase = verify(seq, 'xyz') == 0
        uppercase = verify(seq, 'XYZ') == 0
        if (.not. (lowercase .or. uppercase)) return
        do i = 2, len(seq)
            v = iachar(seq(i:i)) - iachar(seq(i-1:i-1))
            if (v == 0) return
        end do
        ok = .true.
    end function valid_euler_sequence

    pure integer function axis_number(ch) result(axis)
        character(len=1), intent(in) :: ch !! x, y, z, X, Y, or Z
        select case(ch)
        case('x', 'X')
            axis = 1
        case('y', 'Y')
            axis = 2
        case('z', 'Z')
            axis = 3
        case default
            axis = 0
        end select
    end function axis_number

    ! Elementary active rotation, implemented directly as a quaternion.
    pure function axis_rotation(axis, angle) result(r)
        integer, intent(in) :: axis !! Cartesian axis (1:3)
        real(dp), intent(in) :: angle !! rotation angle in radians
        type(rotation_type) :: r
        r = rotation_identity()
        r%q = 0.0_dp
        r%q(axis) = sin(0.5_dp*angle)
        r%q(4) = cos(0.5_dp*angle)
    end function axis_rotation

    function rotation_from_euler(seq, angles, degrees, status) result(r)
        character(len=*), intent(in) :: seq !! 1-3 lowercase extrinsic or uppercase intrinsic axes
        real(dp), intent(in) :: angles(:) !! one angle for each axis in seq
        logical, intent(in), optional :: degrees !! angles in degrees when true; default false
        integer, intent(out), optional :: status !! success or invalid-input code
        type(rotation_type) :: r, elementary
        real(dp) :: angle
        logical :: intrinsic, deg
        integer :: i, axis
        r = invalid_rotation()
        call set_status(status, spatial_rotation_invalid_input)
        if (.not. valid_euler_sequence(seq, .false.)) return
        if (size(angles) /= len(seq)) return
        if (.not. all(ieee_is_finite(angles))) return
        intrinsic = verify(seq, 'XYZ') == 0
        deg = .false.
        if (present(degrees)) deg = degrees
        r = rotation_identity()
        do i = 1, len(seq)
            axis = axis_number(seq(i:i))
            angle = angles(i)
            if (deg) angle = angle*(pi/180.0_dp)
            if (.not. ieee_is_finite(angle)) then
                r = invalid_rotation()
                return
            end if
            elementary = axis_rotation(axis, angle)
            if (intrinsic) then
                r = r%compose(elementary)
            else
                r = elementary%compose(r)
            end if
        end do
        call set_status(status, spatial_rotation_success)
    end function rotation_from_euler

    ! Recover a single active elementary rotation angle from its matrix.
    ! Used at intrinsic gimbal lock to force the *last* returned Euler angle
    ! to zero, matching SciPy's convention rather than only reconstructing SO(3).
    pure function matrix_axis_angle(m, axis) result(angle)
        real(dp), intent(in) :: m(3,3) !! near-elementary rotation matrix
        integer, intent(in) :: axis !! one-based rotation axis
        real(dp) :: angle
        select case(axis)
        case(1)
            angle = atan2(m(3,2),m(2,2))
        case(2)
            angle = atan2(m(1,3),m(1,1))
        case default
            angle = atan2(m(2,1),m(1,1))
        end select
    end function matrix_axis_angle

    ! For any of the 12 extrinsic conventions, recover a,b,c from the
    ! active-rotation matrix. Intrinsic X-Y-Z is reverse extrinsic z-y-x.
    ! This uses direct rotation-matrix identities, with a separate case for
    ! repeated axes (proper Euler angles) and for Tait-Bryan angles.
    function rotation_as_euler(self, seq, degrees, gimbal_lock, status) result(angles)
        class(rotation_type), intent(in) :: self !! rotation to represent
        character(len=*), intent(in) :: seq !! three axes; uppercase intrinsic, lowercase extrinsic
        logical, intent(in), optional :: degrees !! output degrees when true; default radians
        logical, intent(out), optional :: gimbal_lock !! true when first and third angles are not unique
        integer, intent(out), optional :: status !! success or invalid-input code
        real(dp) :: angles(3)
        real(dp) :: m(3,3), a, b, c, hyp, rm(3,3), rmid(3,3)
        integer :: axes(3), i, j, k, parity, d, nex(4)
        logical :: intrinsic, repeated, locked
        character(len=3) :: fixed_axes
        type(rotation_type) :: middle
        angles = ieee_value(0.0_dp, ieee_quiet_nan)
        call set_status(status, spatial_rotation_invalid_input)
        if (present(gimbal_lock)) gimbal_lock = .false.
        if (.not. self%valid()) return
        if (.not. valid_euler_sequence(seq, .true.)) return
        intrinsic = verify(seq, 'XYZ') == 0
        if (intrinsic) then
            fixed_axes = seq(3:3)//seq(2:2)//seq(1:1)
        else
            fixed_axes = seq
        end if
        do d = 1, 3
            axes(d) = axis_number(fixed_axes(d:d))
        end do
        i = axes(1)
        parity = 0
        if (axes(2) /= mod(i,3)+1) parity = 1
        nex = [2,3,1,2]
        j = nex(i+parity)
        k = nex(i-parity+1)
        repeated = axes(1) == axes(3)
        m = self%as_matrix()
        ! Identical mathematical orientations may yield tiny nonzero terms
        ! at gimbal lock because of floating-point matrix construction.
        locked = .false.
        if (repeated) then
            hyp = sqrt(m(i,j)**2+m(i,k)**2)
            locked = hyp <= 1.0e-7_dp
            b = atan2(hyp,m(i,i))
            if (.not. locked) then
                a = atan2(m(i,j),m(i,k))
                c = atan2(m(j,i),-m(k,i))
            else
                a = atan2(-m(j,k),m(j,j))
                c = 0.0_dp
            end if
        else
            hyp = sqrt(m(i,i)**2+m(j,i)**2)
            locked = hyp <= 1.0e-7_dp
            b = atan2(-m(k,i),hyp)
            if (.not. locked) then
                a = atan2(m(k,j),m(k,k))
                c = atan2(m(j,i),m(i,i))
            else
                a = atan2(-m(j,k),m(j,j))
                c = 0.0_dp
            end if
        end if
        if (parity == 1) then
            a = -a
            b = -b
            c = -c
        end if
        if (repeated .and. b < 0.0_dp) then
            ! Proper Euler angles use middle angles in [0,pi]. This
            ! flips the middle axis and shifts the outer two by pi.
            b = -b
            if (.not. locked) then
                a = a + pi
                c = c - pi
            else
                ! With outer angle c=0, recover a anew at the chosen
                ! positive middle angle (the same singularity gauge).
                middle = axis_rotation(axes(2), b)
                rmid = middle%as_matrix()
                rm = matmul(transpose(rmid),m)
                a = matrix_axis_angle(rm,axes(1))
                c = 0.0_dp
            end if
        end if
        if (intrinsic) then
            if (locked) then
                ! At a singularity, the first and third Euler angles can
                ! be combined. Set extrinsic first=0, so intrinsic third=0.
                middle = axis_rotation(axes(2), b)
                rmid = middle%as_matrix()
                rm = matmul(m,transpose(rmid))
                c = matrix_axis_angle(rm,axes(3))
                a = 0.0_dp
            end if
            angles = [c,b,a]
        else
            angles = [a,b,c]
        end if
        ! Wrap both outer angles to SciPy's principal [-pi,pi] interval.
        angles(1) = atan2(sin(angles(1)), cos(angles(1)))
        angles(3) = atan2(sin(angles(3)), cos(angles(3)))
        if (present(degrees)) then
            if (degrees) angles = angles*(180.0_dp/pi)
        end if
        if (present(gimbal_lock)) gimbal_lock = locked
        call set_status(status, spatial_rotation_success)
    end function rotation_as_euler


    ! For a Davenport sequence the adjacent axes must be orthogonal, but
    ! the first and third axes need not be either parallel or orthogonal.
    ! Axes are supplied as (number_of_axes, 3) in sequence order.
    pure subroutine validate_davenport_axes(axes, order, unit, intrinsic, ok)
        real(dp), intent(in) :: axes(:,:) !! 1-3 Cartesian axes, one per row
        character(len=*), intent(in) :: order !! 'e', 'extrinsic', 'i', or 'intrinsic'
        real(dp), intent(out) :: unit(3,3) !! unit vectors; only first n rows used
        logical, intent(out) :: intrinsic !! true for intrinsic/body-fixed rotations
        logical, intent(out) :: ok !! whether input is admissible
        real(dp) :: scale, v(3), mag
        integer :: i, n
        unit = 0.0_dp
        intrinsic = .false.
        ok = .false.
        n = size(axes,1)
        if (size(axes,2) /= 3 .or. n < 1 .or. n > 3) return
        select case (order)
        case ('e', 'extrinsic')
            intrinsic = .false.
        case ('i', 'intrinsic')
            intrinsic = .true.
        case default
            return
        end select
        if (.not. all(ieee_is_finite(axes))) return
        do i = 1, n
            scale = maxval(abs(axes(i,:)))
            if (scale <= 0.0_dp) return
            v = axes(i,:) / scale
            mag = sqrt(sum(v*v))
            unit(i,:) = v/mag
        end do
        do i = 2, n
            if (abs(dot_product(unit(i-1,:),unit(i,:))) > 1.0e-8_dp) return
        end do
        ok = .true.
    end subroutine validate_davenport_axes

    ! Generalized Euler (Davenport) rotations. Composition uses the same
    ! active-rotation and intrinsic/extrinsic conventions as rotation_from_euler.
    function rotation_from_davenport_axes(axes, order, angles, degrees, status) result(r)
        real(dp), intent(in) :: axes(:,:) !! axes in order, shape (n,3), n=1..3
        character(len=*), intent(in) :: order !! 'e'/'extrinsic' or 'i'/'intrinsic'
        real(dp), intent(in) :: angles(:) !! one radian angle per input axis
        logical, intent(in), optional :: degrees !! input angles in degrees, default false
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rotation_type) :: r, elementary
        real(dp) :: u(3,3), v(3), angle
        logical :: intrinsic, ok, deg
        integer :: i, n, code
        r = invalid_rotation()
        call set_status(status, spatial_rotation_invalid_input)
        call validate_davenport_axes(axes, order, u, intrinsic, ok)
        if (.not. ok) return
        n = size(axes,1)
        if (size(angles) /= n) return
        if (.not. all(ieee_is_finite(angles))) return
        deg = .false.
        if (present(degrees)) deg = degrees
        r = rotation_identity()
        do i = 1, n
            angle = angles(i)
            if (deg) angle = angle*(pi/180.0_dp)
            if (.not. ieee_is_finite(angle)) then
                r = invalid_rotation()
                return
            end if
            v = u(i,:)*angle
            elementary = rotation_from_rotvec(v, status=code)
            if (code /= spatial_rotation_success) then
                r = invalid_rotation()
                return
            end if
            if (intrinsic) then
                r = r%compose(elementary)
            else
                r = elementary%compose(r)
            end if
        end do
        call set_status(status, spatial_rotation_success)
    end function rotation_from_davenport_axes

    ! Single Cartesian axis and a scalar angle, the rank-one SciPy form.
    function rotation_from_davenport_axis(axis, order, angle, degrees, status) result(r)
        real(dp), intent(in) :: axis(:) !! one Cartesian axis, length three
        character(len=*), intent(in) :: order !! 'e'/'extrinsic' or 'i'/'intrinsic'
        real(dp), intent(in) :: angle !! one rotation angle
        logical, intent(in), optional :: degrees !! true for degrees
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_type) :: r
        r = invalid_rotation()
        call set_status(status, spatial_rotation_invalid_input)
        if (size(axis) /= 3) return
        r = rotation_from_davenport_axes(reshape(axis,[1,3]),order,[angle],degrees,status)
    end function rotation_from_davenport_axis

    ! Extract generalized Euler angles using the quaternion factorization of
    ! Shuster and Markley's Davenport sequence. A rotation about the middle
    ! axis first maps the oblique third axis into an equivalent proper-Euler
    ! reference frame. The remaining half-angle sums and differences follow
    ! from the vector/scalar projections of the transformed quaternion.
    ! The singular gauge uses a zero THIRD angle, as SciPy does.
    function rotation_as_davenport(self, axes, order, degrees, gimbal_lock, status) result(angles)
        class(rotation_type), intent(in) :: self !! rotation to decompose
        real(dp), intent(in) :: axes(:,:) !! 3 axis vectors, shape (3,3)
        character(len=*), intent(in) :: order !! 'e'/'extrinsic' or 'i'/'intrinsic'
        logical, intent(in), optional :: degrees !! output degrees, default radians
        logical, intent(out), optional :: gimbal_lock !! true for degenerate decomposition
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        real(dp) :: angles(3), u(3,3), n1(3), n2(3), n3(3), cross12(3)
        real(dp) :: lambda, q_lam(4), transformed(4), av, bv, cv, dv
        real(dp) :: half_sum, half_difference, middle, first, third
        logical :: intrinsic, ok, reverse_middle, locked
        type(rotation_type) :: middle_rotation, combined
        angles = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(gimbal_lock)) gimbal_lock = .false.
        call set_status(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (size(axes,1) /= 3 .or. size(axes,2) /= 3) return
        call validate_davenport_axes(axes,order,u,intrinsic,ok)
        if (.not. ok) return
        n1 = u(1,:)
        n2 = u(2,:)
        n3 = u(3,:)
        if (intrinsic) then
            n1 = u(3,:)
            n3 = u(1,:)
        end if
        cross12 = cross3(n1,n2)
        lambda = atan2(dot_product(n3,cross12),dot_product(n3,n1))
        reverse_middle = lambda < 0.0_dp
        if (reverse_middle) then
            n2 = -n2
            cross12 = -cross12
            lambda = -lambda
        end if
        q_lam(1:3) = sin(0.5_dp*lambda)*n2
        q_lam(4) = cos(0.5_dp*lambda)
        middle_rotation = rotation_from_quat(q_lam)
        combined = middle_rotation%compose(self)
        transformed = combined%as_quat()
        av = transformed(4)
        bv = dot_product(transformed(1:3),n1)
        cv = dot_product(transformed(1:3),n2)
        dv = dot_product(transformed(1:3),cross12)
        middle = 2.0_dp*atan2(hypot(cv,dv),hypot(av,bv))
        half_sum = atan2(bv,av)
        half_difference = atan2(dv,cv)
        locked = middle <= 1.0e-7_dp .or. pi-middle <= 1.0e-7_dp
        if (locked) then
            if (middle <= 1.0e-7_dp) then
                first = 2.0_dp*half_sum
            else if (intrinsic) then
                first = 2.0_dp*half_difference
            else
                first = -2.0_dp*half_difference
            end if
            third = 0.0_dp
        else if (intrinsic) then
            first = half_sum + half_difference
            third = half_sum - half_difference
        else
            first = half_sum - half_difference
            third = half_sum + half_difference
        end if
        middle = middle-lambda
        if (reverse_middle) middle = -middle
        angles = [atan2(sin(first),cos(first)), &
                  atan2(sin(middle),cos(middle)), &
                  atan2(sin(third),cos(third))]
        if (present(degrees)) then
            if (degrees) angles = angles*(180.0_dp/pi)
        end if
        if (present(gimbal_lock)) gimbal_lock = locked
        call set_status(status, spatial_rotation_success)
    end function rotation_as_davenport

    ! MRPs sigma = axis*tan(angle/4); q = [2*sigma, 1-|sigma|^2]
    ! divided by (1+|sigma|^2). The scale split avoids intermediate overflow.
    function rotation_from_mrp(mrp, status) result(r)
        real(dp), intent(in) :: mrp(:) !! three modified Rodrigues parameters
        integer, intent(out), optional :: status !! success or invalid-input code
        type(rotation_type) :: r
        real(dp) :: scale, ss, invscale, den, u(3)
        logical :: ok
        r = invalid_rotation()
        call set_status(status, spatial_rotation_invalid_input)
        if (size(mrp) /= 3) return
        if (.not. all(ieee_is_finite(mrp))) return
        scale = maxval(abs(mrp))
        if (scale <= 1.0_dp) then
            ss = sum(mrp*mrp)
            den = 1.0_dp + ss
            r%q(1:3) = (2.0_dp/den)*mrp
            r%q(4) = (1.0_dp-ss)/den
        else
            invscale = 1.0_dp/scale
            u = mrp/scale
            ss = sum(u*u)
            den = ss + invscale*invscale
            r%q(1:3) = (2.0_dp*invscale/den)*u
            r%q(4) = (invscale*invscale-ss)/den
        end if
        call normalize_quaternion(r%q,ok)
        if (.not. ok) then
            r = invalid_rotation()
            return
        end if
        r%defined = .true.
        call set_status(status, spatial_rotation_success)
    end function rotation_from_mrp

    pure function rotation_as_mrp(self) result(mrp)
        class(rotation_type), intent(in) :: self !! rotation to represent
        real(dp) :: mrp(3)
        real(dp) :: q(4)
        mrp = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. self%valid()) return
        q = self%q
        if (q(4) < 0.0_dp) q = -q
        mrp = q(1:3)/(1.0_dp+q(4))
    end function rotation_as_mrp

    function rotation_slerp(start, finish, fraction, status) result(r)
        type(rotation_type), intent(in) :: start !! starting 3D rotation
        type(rotation_type), intent(in) :: finish !! ending 3D rotation
        real(dp), intent(in) :: fraction !! interpolation fraction; extrapolation permitted
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rotation_type) :: r
        real(dp) :: a(4), b(4), dot, omega, factor, sn
        logical :: ok
        r = invalid_rotation()
        call set_status(status,spatial_rotation_invalid_input)
        if (.not. start%valid() .or. .not. finish%valid()) return
        if (.not. ieee_is_finite(fraction)) return
        a = start%q
        b = finish%q
        dot = sum(a*b)
        ! Quaternions q and -q describe the same physical rotation.
        ! Following the hemisphere of a chooses the shortest SO(3) arc.
        if (dot < 0.0_dp) then
            b = -b
            dot = -dot
        end if
        dot = min(1.0_dp,max(0.0_dp,dot))
        if (dot > 1.0_dp-1.0e-10_dp) then
            r%q = (1.0_dp-fraction)*a + fraction*b
        else
            omega = acos(dot)
            sn = sin(omega)
            factor = sin((1.0_dp-fraction)*omega)/sn
            r%q = factor*a + (sin(fraction*omega)/sn)*b
        end if
        call normalize_quaternion(r%q,ok)
        if (.not. ok) then
            r = invalid_rotation()
            return
        end if
        r%defined = .true.
        call set_status(status,spatial_rotation_success)
    end function rotation_slerp

    subroutine rotation_slerp_path(times, rotations, query_times, values, status)
        real(dp), intent(in) :: times(:) !! strictly increasing keyframe times, at least 2
        type(rotation_type), intent(in) :: rotations(:) !! one rotation per keyframe time
        real(dp), intent(in) :: query_times(:) !! interpolation times in the keyframe range
        type(rotation_type), allocatable, intent(out) :: values(:) !! interpolated rotations
        integer, intent(out), optional :: status !! zero on success; one if any input is invalid
        integer :: n, j, lo, hi, mid, k, code
        real(dp) :: time, weight
        call set_status(status,spatial_rotation_invalid_input)
        allocate(values(0))
        n = size(times)
        if (n < 2 .or. size(rotations) /= n) return
        if (.not. all(ieee_is_finite(times))) return
        if (.not. all(ieee_is_finite(query_times))) return
        if (any(times(2:n) <= times(1:n-1))) return
        do j=1,n
            if (.not. rotations(j)%valid()) return
        end do
        if (any(query_times < times(1)) .or. any(query_times > times(n))) return
        deallocate(values)
        allocate(values(size(query_times)))
        do j=1,size(query_times)
            time = query_times(j)
            lo = 1
            hi = n
            do while (hi-lo > 1)
                mid = (lo+hi)/2
                if (times(mid) <= time) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            k = min(lo,n-1)
            weight = (time-times(k))/(times(k+1)-times(k))
            values(j) = rotation_slerp(rotations(k),rotations(k+1),weight,code)
            if (code /= spatial_rotation_success) then
                deallocate(values)
                allocate(values(0))
                return
            end if
        end do
        call set_status(status,spatial_rotation_success)
    end subroutine rotation_slerp_path

    ! Symmetric four-by-four Jacobi eigensolver. The largest eigenvector is
    ! a unit quaternion for both Markley's mean and Davenport's Wahba fit.
    ! Uses only bounded-size local storage and does not square a condition
    ! number by forming a normal matrix of an already ill-conditioned matrix.
    pure subroutine largest_eigenvector4(matrix, vector, success)
        real(dp), intent(in) :: matrix(4,4) !! real symmetric input matrix
        real(dp), intent(out) :: vector(4) !! normalized leading eigenvector
        logical, intent(out) :: success !! Jacobi convergence flag
        real(dp) :: a(4,4), eigenvectors(4,4), apq, tau, t, c, s, oldp, oldq, val
        real(dp) :: scale
        integer :: p, q, k, sweep, best
        logical :: changed

        success = .false.
        vector = 0.0_dp
        if (.not. all(ieee_is_finite(matrix))) return
        scale = maxval(abs(matrix))
        if (scale <= 0.0_dp) return
        a = matrix/scale
        eigenvectors = 0.0_dp
        do k=1,4
            eigenvectors(k,k) = 1.0_dp
        end do
        do sweep=1,80
            changed = .false.
            do p=1,3
                do q=p+1,4
                    apq = a(p,q)
                    if (abs(apq) <= 16.0_dp*epsilon(1.0_dp)) cycle
                    tau = (a(q,q)-a(p,p))/(2.0_dp*apq)
                    if (tau >= 0.0_dp) then
                        t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
                    else
                        t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
                    end if
                    c = 1.0_dp/sqrt(1.0_dp+t*t)
                    s = t*c
                    a(p,p) = a(p,p)-t*apq
                    a(q,q) = a(q,q)+t*apq
                    a(p,q) = 0.0_dp
                    a(q,p) = 0.0_dp
                    do k=1,4
                        if (k /= p .and. k /= q) then
                            oldp = a(k,p)
                            oldq = a(k,q)
                            a(k,p) = c*oldp-s*oldq
                            a(p,k) = a(k,p)
                            a(k,q) = s*oldp+c*oldq
                            a(q,k) = a(k,q)
                        end if
                        oldp = eigenvectors(k,p)
                        oldq = eigenvectors(k,q)
                        eigenvectors(k,p) = c*oldp-s*oldq
                        eigenvectors(k,q) = s*oldp+c*oldq
                    end do
                    changed = .true.
                end do
            end do
            if (.not. changed) exit
        end do
        if (sweep > 80) return
        best = maxloc([a(1,1),a(2,2),a(3,3),a(4,4)],dim=1)
        vector = eigenvectors(:,best)
        val = sqrt(sum(vector*vector))
        if (val <= 0.0_dp .or. .not. ieee_is_finite(val)) return
        vector = vector/val
        if (vector(4) < 0.0_dp) vector = -vector
        success = .true.
    end subroutine largest_eigenvector4

    ! Markley chordal L2 mean: maximize the weighted sum of squared
    ! quaternion inner products. Quaternion antipodes contribute equally.
    function rotation_mean(rotations, weights, status) result(r)
        type(rotation_type), intent(in) :: rotations(:) !! source rotations, at least one
        real(dp), intent(in), optional :: weights(:) !! nonnegative weights; default unit weights
        integer, intent(out), optional :: status !! zero success, one invalid inputs
        type(rotation_type) :: r
        real(dp) :: matrix(4,4), q(4), weight, scale
        logical :: ok
        integer :: i, j, k, n

        r = invalid_rotation()
        call set_status(status,spatial_rotation_invalid_input)
        n = size(rotations)
        if (n < 1) return
        scale = 1.0_dp
        if (present(weights)) then
            if (size(weights) /= n) return
            if (.not. all(ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            scale = maxval(weights)
            if (scale <= 0.0_dp) return
        end if
        matrix = 0.0_dp
        do i=1,n
            if (.not. rotations(i)%valid()) return
            weight = 1.0_dp
            if (present(weights)) weight = weights(i)/scale
            q = rotations(i)%q
            do j=1,4
                do k=j,4
                    matrix(j,k) = matrix(j,k)+weight*q(j)*q(k)
                end do
            end do
        end do
        do j=2,4
            do k=1,j-1
                matrix(j,k) = matrix(k,j)
            end do
        end do
        call largest_eigenvector4(matrix,q,ok)
        if (.not. ok) return
        r%q = q
        r%defined = .true.
        call set_status(status,spatial_rotation_success)
    end function rotation_mean

    ! Stable length and normalization for a three-vector, including inputs
    ! of very small/large nonzero magnitude.
    pure subroutine unit_vector3(vector, unit, success)
        real(dp), intent(in) :: vector(3) !! finite input vector
        real(dp), intent(out) :: unit(3) !! unit direction
        logical, intent(out) :: success !! vector has positive magnitude
        real(dp) :: scale
        unit = 0.0_dp
        success = .false.
        if (.not. all(ieee_is_finite(vector))) return
        scale = maxval(abs(vector))
        if (scale <= 0.0_dp) return
        unit = vector/scale
        unit = unit/sqrt(sum(unit*unit))
        success = .true.
    end subroutine unit_vector3

    pure function shortest_vector_rotation(a, b) result(r)
        real(dp), intent(in) :: a(3) !! destination vector
        real(dp), intent(in) :: b(3) !! source vector
        type(rotation_type) :: r
        real(dp) :: u(3), v(3), axis(3), dot
        logical :: oka, okb, okq
        integer :: k
        r = invalid_rotation()
        call unit_vector3(a,u,oka)
        call unit_vector3(b,v,okb)
        if (.not. oka .or. .not. okb) return
        dot = max(-1.0_dp,min(1.0_dp,dot_product(u,v)))
        axis = cross3(v,u)
        if (dot < 0.0_dp) then
            ! 1+dot suffers catastrophic cancellation near antipodes.
            ! The half-angle formula based on ||u+v|| retains the tiny
            ! departure from a 180-degree rotation.
            if (maxval(abs(axis)) <= 8.0_dp*epsilon(1.0_dp)) then
                k = minloc(abs(v),dim=1)
                axis = 0.0_dp
                axis(k) = 1.0_dp
                axis = cross3(v,axis)
            end if
            axis = axis/sqrt(sum(axis*axis))
            r%q = [axis*(0.5_dp*sqrt(sum((u-v)**2))), &
                   0.5_dp*sqrt(sum((u+v)**2))]
        else
            r%q = [axis,1.0_dp+dot]
        end if
        call normalize_quaternion(r%q,okq)
        if (.not. okq) then
            r = invalid_rotation()
            return
        end if
        r%defined = .true.
    end function shortest_vector_rotation

    ! Compute the inverse of a 3x3 positive-definite information matrix.
    ! The scale-free matrix makes the inversion insensitive to units.
    pure subroutine invert_information3(matrix, inverse, success)
        real(dp), intent(in) :: matrix(3,3) !! symmetric positive-definite matrix
        real(dp), intent(out) :: inverse(3,3) !! inverse on success
        logical, intent(out) :: success !! matrix is invertible
        real(dp) :: work(3,6), pivot, scale, temp(6)
        integer :: i, k, j
        inverse = ieee_value(0.0_dp,ieee_quiet_nan)
        success = .false.
        scale = maxval(abs(matrix))
        if (scale <= 0.0_dp .or. .not. ieee_is_finite(scale)) return
        work = 0.0_dp
        work(:,1:3) = matrix/scale
        do i=1,3
            work(i,i+3) = 1.0_dp
        end do
        do j=1,3
            k = j-1+maxloc(abs(work(j:3,j)),dim=1)
            if (abs(work(k,j)) <= 64.0_dp*epsilon(1.0_dp)) return
            if (k /= j) then
                temp = work(j,:)
                work(j,:) = work(k,:)
                work(k,:) = temp
            end if
            pivot = work(j,j)
            work(j,:) = work(j,:)/pivot
            do i=1,3
                if (i == j) cycle
                pivot = work(i,j)
                work(i,:) = work(i,:)-pivot*work(j,:)
            end do
        end do
        inverse = work(:,4:6)/scale
        success = all(ieee_is_finite(inverse))
    end subroutine invert_information3

    ! Weighted Wahba/Kabsch rotation mapping source rows b onto target rows a.
    ! Uses the leading eigenvector of Davenport's 4x4 quaternion matrix,
    ! without requiring LAPACK. The infinite-weight case is handled as an
    ! exact primary-vector constraint followed by rotation about that axis.
    function rotation_align_vectors(a, b, weights, rssd, sensitivity, status) result(r)
        real(dp), intent(in) :: a(:,:) !! target vectors, shape (n,3)
        real(dp), intent(in) :: b(:,:) !! source vectors, shape (n,3)
        real(dp), intent(in), optional :: weights(:) !! strictly positive weights; at most one +Inf
        real(dp), intent(out), optional :: rssd !! root weighted sum of squared residual distances
        real(dp), intent(out), optional :: sensitivity(3,3) !! sensitivity of rotation-vector fit; finite multi-vector case only
        integer, intent(out), optional :: status !! zero success; one invalid or non-identifiable inputs
        type(rotation_type) :: r, primary, axial
        real(dp) :: h(3,3), kappa(4,4), q(4), wscale, ascale, bscale, c, s
        real(dp) :: x(3), y(3), z(3), u(3), info(3,3), inv_info(3,3)
        real(dp) :: w(size(a,1)), loss, harmonic, residual(3), axis_angle(3)
        integer :: n, i, j, l, primary_idx
        logical :: ok, infinite_weights(size(a,1))

        r = invalid_rotation()
        call set_status(status,spatial_rotation_invalid_input)
        if (present(rssd)) rssd = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(sensitivity)) sensitivity = ieee_value(0.0_dp,ieee_quiet_nan)
        n = size(a,1)
        if (n < 1 .or. size(a,2) /= 3) return
        if (size(b,1) /= n .or. size(b,2) /= 3) return
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        w = 1.0_dp
        if (present(weights)) then
            if (size(weights) /= n) return
            if (any(ieee_is_nan(weights))) return
            if (any(weights <= 0.0_dp)) return
            w = weights
        end if
        infinite_weights = .not. ieee_is_finite(w)
        if (count(infinite_weights) > 1) return
        if (present(sensitivity)) then
            if (n < 2 .or. any(infinite_weights)) return
        end if
        if (any(infinite_weights)) then
            primary_idx = findloc(infinite_weights,.true.,dim=1)
            primary = shortest_vector_rotation(a(primary_idx,:),b(primary_idx,:))
            if (.not. primary%valid()) return
            call unit_vector3(a(primary_idx,:),u,ok)
            if (.not. ok) return
            ! Choose the optimal twist about the already aligned primary axis.
            c = 0.0_dp
            s = 0.0_dp
            wscale = 0.0_dp
            ascale = maxval(abs(a))
            bscale = maxval(abs(b))
            do i=1,n
                if (.not. infinite_weights(i)) wscale = max(wscale,w(i))
            end do
            if (wscale > 0.0_dp) then
                do i=1,n
                    if (infinite_weights(i)) cycle
                    x = a(i,:)/ascale
                    y = primary%apply(b(i,:)/bscale)
                    x = x-dot_product(u,x)*u
                    y = y-dot_product(u,y)*u
                    c = c+(w(i)/wscale)*dot_product(x,y)
                    s = s+(w(i)/wscale)*dot_product(x,cross3(u,y))
                end do
            end if
            axis_angle = u*atan2(s,c)
            axial = rotation_from_rotvec(axis_angle)
            r = axial%compose(primary)
        else if (n == 1) then
            r = shortest_vector_rotation(a(1,:),b(1,:))
        else
            ascale = maxval(abs(a))
            bscale = maxval(abs(b))
            if (ascale <= 0.0_dp .or. bscale <= 0.0_dp) return
            wscale = maxval(w)
            h = 0.0_dp
            do i=1,n
                x = a(i,:)/ascale
                y = b(i,:)/bscale
                do j=1,3
                    do l=1,3
                        h(j,l) = h(j,l)+(w(i)/wscale)*x(j)*y(l)
                    end do
                end do
            end do
            kappa = 0.0_dp
            do j=1,3
                do l=1,3
                    kappa(j,l) = h(j,l)+h(l,j)
                end do
                kappa(j,j) = kappa(j,j)-h(1,1)-h(2,2)-h(3,3)
            end do
            kappa(1:3,4) = [h(3,2)-h(2,3),h(1,3)-h(3,1),h(2,1)-h(1,2)]
            kappa(4,1:3) = kappa(1:3,4)
            kappa(4,4) = h(1,1)+h(2,2)+h(3,3)
            call largest_eigenvector4(kappa,q,ok)
            if (.not. ok) return
            r%q = q
            r%defined = .true.
        end if
        if (.not. r%valid()) return

        if (present(rssd)) then
            ascale = max(maxval(abs(a)),maxval(abs(b)))
            wscale = maxval(w,mask=.not. infinite_weights)
            if (wscale <= 0.0_dp) then
                rssd = 0.0_dp
            else
                loss = 0.0_dp
                do i=1,n
                    if (infinite_weights(i)) cycle
                    residual = a(i,:)/ascale-r%apply(b(i,:)/ascale)
                    loss = loss+(w(i)/wscale)*sum(residual*residual)
                end do
                rssd = (ascale*sqrt(wscale))*sqrt(max(0.0_dp,loss))
            end if
        end if
        if (present(sensitivity)) then
            ascale = maxval(abs(a))
            bscale = maxval(abs(b))
            if (ascale <= 0.0_dp .or. bscale <= 0.0_dp) then
                r = invalid_rotation()
                return
            end if
            wscale = maxval(w)
            info = 0.0_dp
            do i=1,n
                x = a(i,:)/ascale
                y = r%apply(b(i,:)/bscale)
                info = info+(w(i)/wscale)*dot_product(x,y)*identity3()
                do j=1,3
                    do l=1,3
                        info(j,l) = info(j,l)-0.5_dp*(w(i)/wscale)*(x(j)*y(l)+y(j)*x(l))
                    end do
                end do
            end do
            call invert_information3(info,inv_info,ok)
            if (.not. ok) then
                r = invalid_rotation()
                return
            end if
            ! SciPy normalizes inverse information by the arithmetic mean
            ! of supplied weights (scale-independent under common scaling).
            harmonic = sum(w/wscale)/real(n,dp)
            sensitivity = (harmonic/ascale)*(inv_info/bscale)
        end if
        call set_status(status,spatial_rotation_success)
    end function rotation_align_vectors

    pure function identity3() result(m)
        real(dp) :: m(3,3)
        m = 0.0_dp
        m(1,1) = 1.0_dp
        m(2,2) = 1.0_dp
        m(3,3) = 1.0_dp
    end function identity3

end module scifort_spatial_transform
