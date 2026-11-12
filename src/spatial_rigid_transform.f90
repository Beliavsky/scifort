! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Native SE(3) rigid transforms. Conventions match SciPy RigidTransform:
! x' = R*x + t, active rotations, quaternion order (x,y,z,w), and
! exponential coordinates (rx,ry,rz,vx,vy,vz).
module scifort_spatial_rigid_transform
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, &
        rotation_from_matrix, rotation_from_rotvec, rotation_from_quat, rotation_mean
    implicit none
    private

    integer, parameter, public :: spatial_rigid_success = 0
    integer, parameter, public :: spatial_rigid_invalid_input = 1
    real(dp), parameter :: small_angle = 1.0e-4_dp

    type, public :: rigid_transform_type
        private
        type(rotation_type) :: r
        real(dp) :: t(3) = 0.0_dp
        logical :: defined = .true.
    contains
        procedure :: valid => rigid_valid
        procedure :: rotation => rigid_rotation
        procedure :: translation => rigid_translation
        procedure :: as_matrix => rigid_as_matrix
        procedure :: as_exp_coords => rigid_as_exp_coords
        procedure :: as_dual_quat => rigid_as_dual_quat
        procedure, private :: rigid_apply_vector
        procedure, private :: rigid_apply_points
        generic :: apply => rigid_apply_vector, rigid_apply_points
        procedure :: inverse => rigid_inverse
        procedure :: compose => rigid_compose
        procedure :: power => rigid_power
    end type rigid_transform_type

    ! A one-dimensional transform stack, with size-one broadcasting.
    type, public :: rigid_transform_stack_type
        private
        type(rigid_transform_type), allocatable :: values(:)
        logical :: defined = .false.
    contains
        procedure :: valid => rigid_stack_valid
        procedure :: size => rigid_stack_size
        procedure :: at => rigid_stack_at
        procedure :: as_matrix => rigid_stack_as_matrix
        procedure :: as_exp_coords => rigid_stack_as_exp_coords
        procedure :: as_dual_quat => rigid_stack_as_dual_quat
        procedure :: translation => rigid_stack_translation
        procedure, private :: rigid_stack_apply_vector
        procedure, private :: rigid_stack_apply_points
        generic :: apply => rigid_stack_apply_vector, rigid_stack_apply_points
        procedure :: inverse => rigid_stack_inverse
        procedure, private :: rigid_stack_compose_stack
        procedure, private :: rigid_stack_compose_single
        generic :: compose => rigid_stack_compose_stack, rigid_stack_compose_single
        procedure :: power => rigid_stack_power
        procedure :: mean => rigid_stack_mean
    end type rigid_transform_stack_type

    public :: rigid_identity, rigid_from_components, rigid_from_translation, rigid_from_rotation
    public :: rigid_from_matrix, rigid_from_exp_coords, rigid_from_dual_quat
    public :: rigid_stack_from_components, rigid_stack_from_matrix, rigid_stack_from_exp_coords
    public :: rigid_stack_from_dual_quat

contains

    pure subroutine set_status(status, code)
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        integer, intent(in) :: code !! status to return
        if (present(status)) status = code
    end subroutine set_status

    pure function nan3() result(v)
        real(dp) :: v(3)
        v = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan3

    pure function invalid_rigid() result(tf)
        type(rigid_transform_type) :: tf
        tf%defined = .false.
        tf%t = nan3()
    end function invalid_rigid

    pure logical function rigid_valid(self) result(ok)
        class(rigid_transform_type), intent(in) :: self !! SE(3) transform
        ok = self%defined .and. self%r%valid() .and. all(ieee_is_finite(self%t))
    end function rigid_valid

    pure function rigid_identity() result(tf)
        type(rigid_transform_type) :: tf
        tf%r = rotation_identity()
        tf%t = 0.0_dp
        tf%defined = .true.
    end function rigid_identity

    function rigid_from_components(translation, rotation, status) result(tf)
        real(dp), intent(in) :: translation(:) !! translation components; length 3
        type(rotation_type), intent(in) :: rotation !! active rotation, applied before translation
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_type) :: tf
        tf = invalid_rigid()
        call set_status(status, spatial_rigid_invalid_input)
        if (size(translation) /= 3) return
        if (.not. all(ieee_is_finite(translation))) return
        if (.not. rotation%valid()) return
        tf%t = translation
        tf%r = rotation
        tf%defined = .true.
        call set_status(status, spatial_rigid_success)
    end function rigid_from_components

    function rigid_from_translation(translation, status) result(tf)
        real(dp), intent(in) :: translation(:) !! Cartesian translation, length 3
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_type) :: tf
        tf = rigid_from_components(translation, rotation_identity(), status)
    end function rigid_from_translation

    function rigid_from_rotation(rotation, status) result(tf)
        type(rotation_type), intent(in) :: rotation !! active rotation
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_type) :: tf
        tf = rigid_from_components([0.0_dp,0.0_dp,0.0_dp], rotation, status)
    end function rigid_from_rotation

    function rigid_from_matrix(matrix, status) result(tf)
        real(dp), intent(in) :: matrix(:,:) !! 4x4 homogeneous rigid transform
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_type) :: tf
        type(rotation_type) :: r
        integer :: code
        tf = invalid_rigid()
        call set_status(status, spatial_rigid_invalid_input)
        if (size(matrix,1) /= 4 .or. size(matrix,2) /= 4) return
        if (.not. all(ieee_is_finite(matrix))) return
        if (maxval(abs(matrix(4,:) - [0.0_dp,0.0_dp,0.0_dp,1.0_dp])) > 1.0e-10_dp) return
        r = rotation_from_matrix(matrix(1:3,1:3), status=code)
        if (code /= 0) return
        tf = rigid_from_components(matrix(1:3,4), r, status)
    end function rigid_from_matrix

    function rigid_rotation(self) result(r)
        class(rigid_transform_type), intent(in) :: self !! SE(3) transform
        type(rotation_type) :: r
        r = self%r
        if (.not. self%valid()) r = rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,0.0_dp])
    end function rigid_rotation

    pure function rigid_translation(self) result(t)
        class(rigid_transform_type), intent(in) :: self !! SE(3) transform
        real(dp) :: t(3)
        t = self%t
        if (.not. self%valid()) t = nan3()
    end function rigid_translation

    pure function rigid_as_matrix(self) result(matrix)
        class(rigid_transform_type), intent(in) :: self !! SE(3) transform
        real(dp) :: matrix(4,4)
        matrix = 0.0_dp
        if (.not. self%valid()) then
            matrix = ieee_value(0.0_dp,ieee_quiet_nan)
            return
        end if
        matrix(1:3,1:3) = self%r%as_matrix()
        matrix(1:3,4) = self%t
        matrix(4,4) = 1.0_dp
    end function rigid_as_matrix

    pure function cross3(a,b) result(c)
        real(dp), intent(in) :: a(3) !! first 3-vector
        real(dp), intent(in) :: b(3) !! second 3-vector
        real(dp) :: c(3)
        c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), &
            a(1)*b(2)-a(2)*b(1)]
    end function cross3

    ! Left Jacobian of SO(3) times a vector. Taylor polynomials preserve
    ! accuracy at zero rotation, without dividing by tiny angles.
    pure function so3_left_jacobian_times(w, v) result(t)
        real(dp), intent(in) :: w(3) !! rotation vector
        real(dp), intent(in) :: v(3) !! infinitesimal translation vector
        real(dp) :: t(3), theta, theta2, a, b
        theta = sqrt(sum(w*w))
        theta2 = theta*theta
        if (theta < small_angle) then
            a = 0.5_dp - theta2/24.0_dp + theta2**2/720.0_dp
            b = 1.0_dp/6.0_dp - theta2/120.0_dp + theta2**2/5040.0_dp
        else
            a = 2.0_dp*(sin(0.5_dp*theta)/theta)**2
            b = (theta-sin(theta))/(theta2*theta)
        end if
        t = v + a*cross3(w,v) + b*cross3(w,cross3(w,v))
    end function so3_left_jacobian_times

    pure function so3_left_jacobian_inverse_times(w,t) result(v)
        real(dp), intent(in) :: w(3) !! principal SO(3) rotation vector
        real(dp), intent(in) :: t(3) !! Cartesian translation
        real(dp) :: v(3), theta, theta2, c
        theta = sqrt(sum(w*w))
        theta2 = theta*theta
        if (theta < small_angle) then
            c = 1.0_dp/12.0_dp + theta2/720.0_dp + theta2**2/30240.0_dp
        else
            c = (1.0_dp - 0.5_dp*theta/tan(0.5_dp*theta))/theta2
        end if
        v = t - 0.5_dp*cross3(w,t) + c*cross3(w,cross3(w,t))
    end function so3_left_jacobian_inverse_times

    function rigid_from_exp_coords(exp_coords, status) result(tf)
        real(dp), intent(in) :: exp_coords(:) !! 6-vector (rx,ry,rz,vx,vy,vz)
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_type) :: tf
        type(rotation_type) :: r
        real(dp) :: t(3)
        integer :: code
        tf = invalid_rigid()
        call set_status(status, spatial_rigid_invalid_input)
        if (size(exp_coords) /= 6) return
        if (.not. all(ieee_is_finite(exp_coords))) return
        r = rotation_from_rotvec(exp_coords(1:3),status=code)
        if (code /= 0) return
        t = so3_left_jacobian_times(exp_coords(1:3),exp_coords(4:6))
        if (.not. all(ieee_is_finite(t))) return
        tf = rigid_from_components(t,r,status)
    end function rigid_from_exp_coords

    pure function rigid_as_exp_coords(self) result(coords)
        class(rigid_transform_type), intent(in) :: self !! rigid transform
        real(dp) :: coords(6)
        if (.not. self%valid()) then
            coords = ieee_value(0.0_dp,ieee_quiet_nan)
            return
        end if
        coords(1:3) = self%r%as_rotvec()
        coords(4:6) = so3_left_jacobian_inverse_times(coords(1:3),self%t)
    end function rigid_as_exp_coords

    pure function quat_product(a,b) result(q)
        real(dp), intent(in) :: a(4) !! left scalar-last quaternion
        real(dp), intent(in) :: b(4) !! right scalar-last quaternion
        real(dp) :: q(4)
        q(1:3) = a(4)*b(1:3) + b(4)*a(1:3) + cross3(a(1:3),b(1:3))
        q(4) = a(4)*b(4) - dot_product(a(1:3),b(1:3))
    end function quat_product

    pure function rigid_as_dual_quat(self, scalar_first) result(dq)
        class(rigid_transform_type), intent(in) :: self !! rigid transform
        logical, intent(in), optional :: scalar_first !! true for (w,x,y,z) within each quaternion
        real(dp) :: dq(8), q(4), d(4)
        logical :: first
        if (.not. self%valid()) then
            dq = ieee_value(0.0_dp,ieee_quiet_nan)
            return
        end if
        q = self%r%as_quat()
        d = 0.5_dp*quat_product([self%t,0.0_dp],q)
        first = .false.
        if (present(scalar_first)) first = scalar_first
        if (first) then
            dq = [q(4),q(1:3),d(4),d(1:3)]
        else
            dq = [q,d]
        end if
    end function rigid_as_dual_quat

    function rigid_from_dual_quat(dual_quat, scalar_first, status) result(tf)
        real(dp), intent(in) :: dual_quat(:) !! eight dual quaternion components
        logical, intent(in), optional :: scalar_first !! true for (w,x,y,z) within each quaternion
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_type) :: tf
        type(rotation_type) :: r
        real(dp) :: q(4), d(4), s, mag, tmp(4), t(3)
        logical :: first
        integer :: code
        tf = invalid_rigid()
        call set_status(status, spatial_rigid_invalid_input)
        if (size(dual_quat) /= 8) return
        if (.not. all(ieee_is_finite(dual_quat))) return
        first = .false.
        if (present(scalar_first)) first = scalar_first
        if (first) then
            q = [dual_quat(2:4),dual_quat(1)]
            d = [dual_quat(6:8),dual_quat(5)]
        else
            q = dual_quat(1:4)
            d = dual_quat(5:8)
        end if
        s = maxval(abs(q))
        if (s <= 0.0_dp) return
        q = q/s
        d = d/s
        mag = sqrt(sum(q*q))
        q = q/mag
        d = d/mag
        ! Enforce the unit-dual quaternion orthogonality constraint.
        d = d - dot_product(q,d)*q
        tmp = quat_product(d,[-q(1:3),q(4)])
        t = 2.0_dp*tmp(1:3)
        if (.not. all(ieee_is_finite(t))) return
        r = rotation_from_quat(q,status=code)
        if (code /= 0) return
        tf = rigid_from_components(t,r,status)
    end function rigid_from_dual_quat

    pure function rigid_apply_vector(self, vector, inverse) result(v)
        class(rigid_transform_type), intent(in) :: self !! rigid transform
        real(dp), intent(in) :: vector(:) !! Cartesian vector, length 3
        logical, intent(in), optional :: inverse !! use the inverse transform when true
        real(dp) :: v(3)
        logical :: backwards
        v = nan3()
        if (.not. self%valid() .or. size(vector) /= 3) return
        if (.not. all(ieee_is_finite(vector))) return
        backwards = .false.
        if (present(inverse)) backwards = inverse
        if (backwards) then
            v = self%r%apply(vector-self%t,inverse=.true.)
        else
            v = self%r%apply(vector)+self%t
        end if
    end function rigid_apply_vector

    pure function rigid_apply_points(self, points, inverse) result(v)
        class(rigid_transform_type), intent(in) :: self !! rigid transform
        real(dp), intent(in) :: points(:,:) !! point cloud, shape (n,3)
        logical, intent(in), optional :: inverse !! apply inverse transform when true
        real(dp), allocatable :: v(:,:)
        integer :: i
        allocate(v(size(points,1),3))
        v = ieee_value(0.0_dp,ieee_quiet_nan)
        if (size(points,2) /= 3 .or. .not. self%valid()) return
        do i=1,size(points,1)
            v(i,:) = self%apply(points(i,:),inverse)
        end do
    end function rigid_apply_points

    pure function rigid_inverse(self) result(tf)
        class(rigid_transform_type), intent(in) :: self !! rigid transform
        type(rigid_transform_type) :: tf
        tf = invalid_rigid()
        if (.not. self%valid()) return
        tf%r = self%r%inverse()
        tf%t = tf%r%apply(-self%t)
        tf%defined = all(ieee_is_finite(tf%t))
        if (.not. tf%defined) tf = invalid_rigid()
    end function rigid_inverse

    pure function rigid_compose(self, other) result(tf)
        class(rigid_transform_type), intent(in) :: self !! transform applied second
        type(rigid_transform_type), intent(in) :: other !! transform applied first
        type(rigid_transform_type) :: tf
        tf = invalid_rigid()
        if (.not. self%valid() .or. .not. other%valid()) return
        tf%r = self%r%compose(other%r)
        tf%t = self%r%apply(other%t)+self%t
        tf%defined = tf%r%valid() .and. all(ieee_is_finite(tf%t))
        if (.not. tf%defined) tf = invalid_rigid()
    end function rigid_compose

    function rigid_power(self, exponent, status) result(tf)
        class(rigid_transform_type), intent(in) :: self !! rigid transform
        real(dp), intent(in) :: exponent !! real exponent, may be fractional or negative
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_type) :: tf
        real(dp) :: x(6)
        tf = invalid_rigid()
        call set_status(status, spatial_rigid_invalid_input)
        if (.not. self%valid() .or. .not. ieee_is_finite(exponent)) return
        if (exponent == 0.0_dp) then
            tf = rigid_identity()
        else if (exponent == 1.0_dp) then
            tf = self
        else if (exponent == -1.0_dp) then
            tf = self%inverse()
        else
            x = self%as_exp_coords()*exponent
            if (.not. all(ieee_is_finite(x))) return
            tf = rigid_from_exp_coords(x)
        end if
        if (.not. tf%valid()) return
        call set_status(status, spatial_rigid_success)
    end function rigid_power

    pure logical function rigid_stack_valid(self) result(ok)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        integer :: i
        ok = .false.
        if (.not. self%defined .or. .not. allocated(self%values)) return
        if (size(self%values) == 0) return
        do i=1,size(self%values)
            if (.not. self%values(i)%valid()) return
        end do
        ok = .true.
    end function rigid_stack_valid

    pure integer function rigid_stack_size(self) result(n)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        n = 0
        if (allocated(self%values)) n = size(self%values)
    end function rigid_stack_size

    pure function rigid_stack_at(self, i) result(tf)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        integer, intent(in) :: i !! one-based index
        type(rigid_transform_type) :: tf
        tf = invalid_rigid()
        if (.not. self%valid()) return
        if (i < 1 .or. i > size(self%values)) return
        tf = self%values(i)
    end function rigid_stack_at

    function rigid_stack_from_components(translations, rotations, status) result(stack)
        real(dp), intent(in) :: translations(:,:) !! translations (n,3), n>=1
        type(rotation_type), intent(in) :: rotations(:) !! rotations, length n or 1
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_stack_type) :: stack
        integer :: i, n, nt, nr, it, ir
        call set_status(status, spatial_rigid_invalid_input)
        if (size(translations,2) /= 3) return
        nt = size(translations,1)
        nr = size(rotations)
        n = max(nt,nr)
        if (n < 1 .or. nt < 1 .or. nr < 1) return
        if ((nt /= n .and. nt /= 1) .or. (nr /= n .and. nr /= 1)) return
        allocate(stack%values(n))
        do i=1,n
            it = min(i,nt)
            ir = min(i,nr)
            stack%values(i) = rigid_from_components(translations(it,:),rotations(ir))
            if (.not. stack%values(i)%valid()) return
        end do
        stack%defined = .true.
        call set_status(status, spatial_rigid_success)
    end function rigid_stack_from_components

    function rigid_stack_from_matrix(matrices, status) result(stack)
        real(dp), intent(in) :: matrices(:,:,:) !! matrices with shape (n,4,4)
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_stack_type) :: stack
        integer :: i, n
        call set_status(status, spatial_rigid_invalid_input)
        n = size(matrices,1)
        if (n < 1 .or. size(matrices,2) /= 4 .or. size(matrices,3) /= 4) return
        allocate(stack%values(n))
        do i=1,n
            stack%values(i) = rigid_from_matrix(matrices(i,:,:))
            if (.not. stack%values(i)%valid()) return
        end do
        stack%defined = .true.
        call set_status(status, spatial_rigid_success)
    end function rigid_stack_from_matrix

    function rigid_stack_from_exp_coords(coords, status) result(stack)
        real(dp), intent(in) :: coords(:,:) !! exponential coords, shape (n,6)
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_stack_type) :: stack
        integer :: i, n
        call set_status(status, spatial_rigid_invalid_input)
        n = size(coords,1)
        if (n < 1 .or. size(coords,2) /= 6) return
        allocate(stack%values(n))
        do i=1,n
            stack%values(i) = rigid_from_exp_coords(coords(i,:))
            if (.not. stack%values(i)%valid()) return
        end do
        stack%defined = .true.
        call set_status(status, spatial_rigid_success)
    end function rigid_stack_from_exp_coords

    function rigid_stack_from_dual_quat(dqs, scalar_first, status) result(stack)
        real(dp), intent(in) :: dqs(:,:) !! dual quaternions, shape (n,8)
        logical, intent(in), optional :: scalar_first !! scalar-first order within both quaternions
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_stack_type) :: stack
        integer :: i, n
        call set_status(status, spatial_rigid_invalid_input)
        n = size(dqs,1)
        if (n < 1 .or. size(dqs,2) /= 8) return
        allocate(stack%values(n))
        do i=1,n
            stack%values(i) = rigid_from_dual_quat(dqs(i,:),scalar_first)
            if (.not. stack%values(i)%valid()) return
        end do
        stack%defined = .true.
        call set_status(status, spatial_rigid_success)
    end function rigid_stack_from_dual_quat

    pure function rigid_stack_as_matrix(self) result(matrices)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        real(dp), allocatable :: matrices(:,:,:)
        integer :: i
        allocate(matrices(self%size(),4,4))
        matrices = ieee_value(0.0_dp,ieee_quiet_nan)
        if (.not. self%valid()) return
        do i=1,self%size()
            matrices(i,:,:) = self%values(i)%as_matrix()
        end do
    end function rigid_stack_as_matrix

    pure function rigid_stack_as_exp_coords(self) result(coords)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        real(dp), allocatable :: coords(:,:)
        integer :: i
        allocate(coords(self%size(),6))
        coords = ieee_value(0.0_dp,ieee_quiet_nan)
        if (.not. self%valid()) return
        do i=1,self%size()
            coords(i,:) = self%values(i)%as_exp_coords()
        end do
    end function rigid_stack_as_exp_coords

    pure function rigid_stack_as_dual_quat(self, scalar_first) result(dqs)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        logical, intent(in), optional :: scalar_first !! quaternion component order
        real(dp), allocatable :: dqs(:,:)
        integer :: i
        allocate(dqs(self%size(),8))
        dqs = ieee_value(0.0_dp,ieee_quiet_nan)
        if (.not. self%valid()) return
        do i=1,self%size()
            dqs(i,:) = self%values(i)%as_dual_quat(scalar_first)
        end do
    end function rigid_stack_as_dual_quat

    pure function rigid_stack_translation(self) result(translations)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        real(dp), allocatable :: translations(:,:)
        integer :: i
        allocate(translations(self%size(),3))
        translations = ieee_value(0.0_dp,ieee_quiet_nan)
        if (.not. self%valid()) return
        do i=1,self%size()
            translations(i,:) = self%values(i)%translation()
        end do
    end function rigid_stack_translation

    pure function rigid_stack_apply_vector(self, vector, inverse) result(out)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        real(dp), intent(in) :: vector(:) !! common Cartesian vector, length 3
        logical, intent(in), optional :: inverse !! apply inverse transforms
        real(dp), allocatable :: out(:,:)
        integer :: i
        allocate(out(self%size(),3))
        out = ieee_value(0.0_dp,ieee_quiet_nan)
        if (.not. self%valid() .or. size(vector) /= 3) return
        do i=1,self%size()
            out(i,:) = self%values(i)%apply(vector,inverse)
        end do
    end function rigid_stack_apply_vector

    pure function rigid_stack_apply_points(self, points, inverse) result(out)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        real(dp), intent(in) :: points(:,:) !! Cartesian points, shape (m,3), m=n or 1, or n=1
        logical, intent(in), optional :: inverse !! apply inverse transforms
        real(dp), allocatable :: out(:,:)
        integer :: i, n, m, k, ip
        n = self%size()
        m = size(points,1)
        allocate(out(max(n,m),3))
        out = ieee_value(0.0_dp,ieee_quiet_nan)
        if (.not. self%valid() .or. size(points,2) /= 3 .or. m < 1) return
        if (n /= m .and. n /= 1 .and. m /= 1) return
        do i=1,max(n,m)
            k = min(i,n)
            ip = min(i,m)
            out(i,:) = self%values(k)%apply(points(ip,:),inverse)
        end do
    end function rigid_stack_apply_points

    pure function rigid_stack_inverse(self) result(stack)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        type(rigid_transform_stack_type) :: stack
        integer :: i
        if (.not. self%valid()) return
        allocate(stack%values(self%size()))
        do i=1,self%size()
            stack%values(i) = self%values(i)%inverse()
            if (.not. stack%values(i)%valid()) return
        end do
        stack%defined = .true.
    end function rigid_stack_inverse

    pure function rigid_stack_compose_stack(self, other) result(stack)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack applied second
        type(rigid_transform_stack_type), intent(in) :: other !! transform stack applied first
        type(rigid_transform_stack_type) :: stack
        integer :: i,n,m,l,is,io
        if (.not. self%valid() .or. .not. other%valid()) return
        n = self%size()
        m = other%size()
        if (n /= m .and. n /= 1 .and. m /= 1) return
        l = max(n,m)
        allocate(stack%values(l))
        do i=1,l
            is = min(i,n)
            io = min(i,m)
            stack%values(i) = self%values(is)%compose(other%values(io))
            if (.not. stack%values(i)%valid()) return
        end do
        stack%defined = .true.
    end function rigid_stack_compose_stack

    pure function rigid_stack_compose_single(self, other) result(stack)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack applied second
        type(rigid_transform_type), intent(in) :: other !! transform applied first
        type(rigid_transform_stack_type) :: stack
        integer :: i
        if (.not. self%valid() .or. .not. other%valid()) return
        allocate(stack%values(self%size()))
        do i=1,self%size()
            stack%values(i) = self%values(i)%compose(other)
            if (.not. stack%values(i)%valid()) return
        end do
        stack%defined = .true.
    end function rigid_stack_compose_single

    function rigid_stack_power(self, exponent, status) result(stack)
        class(rigid_transform_stack_type), intent(in) :: self !! transform stack
        real(dp), intent(in) :: exponent !! real transform exponent
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_stack_type) :: stack
        integer :: i
        call set_status(status, spatial_rigid_invalid_input)
        if (.not. self%valid() .or. .not. ieee_is_finite(exponent)) return
        allocate(stack%values(self%size()))
        do i=1,self%size()
            stack%values(i) = self%values(i)%power(exponent)
            if (.not. stack%values(i)%valid()) return
        end do
        stack%defined = .true.
        call set_status(status, spatial_rigid_success)
    end function rigid_stack_power

    function rigid_stack_mean(self, weights, status) result(tf)
        class(rigid_transform_stack_type), intent(in) :: self !! stack of transforms to average
        real(dp), intent(in), optional :: weights(:) !! nonnegative weights, length n
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rigid_transform_type) :: tf
        type(rotation_type), allocatable :: rotations(:)
        type(rotation_type) :: mean_rotation
        real(dp), allocatable :: w(:)
        real(dp) :: total, t(3), scale
        integer :: n, i, code
        tf = invalid_rigid()
        call set_status(status, spatial_rigid_invalid_input)
        if (.not. self%valid()) return
        n = self%size()
        allocate(w(n),rotations(n))
        w = 1.0_dp
        if (present(weights)) then
            if (size(weights) /= n) return
            if (.not. all(ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            w = weights
        end if
        scale = maxval(w)
        if (scale <= 0.0_dp) return
        w = w/scale
        total = sum(w)
        t = 0.0_dp
        do i=1,n
            rotations(i) = self%values(i)%rotation()
            t = t + (w(i)/total)*self%values(i)%t
        end do
        if (.not. all(ieee_is_finite(t))) return
        mean_rotation = rotation_mean(rotations, weights=w, status=code)
        if (code /= 0) return
        tf = rigid_from_components(t,mean_rotation,status)
    end function rigid_stack_mean

end module scifort_spatial_rigid_transform
