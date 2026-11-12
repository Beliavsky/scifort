! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Vectorized rotations / Rotation stacks. Quaternion representation and
! composition convention are inherited from scifort_spatial_transform.
module scifort_spatial_rotation_stack
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, &
        rotation_from_quat, rotation_from_rotvec, rotation_from_matrix, &
        rotation_from_mrp, rotation_from_euler, rotation_from_davenport, rotation_mean, &
        spatial_rotation_success, spatial_rotation_invalid_input
    implicit none
    private

    ! Explicit reproducible random stream. Park-Miller uses 64-bit arithmetic
    ! without signed overflow and is independent of the compiler's RANDOM_SEED.
    type, public :: rotation_rng_type
        private
        integer(int64) :: state = 0_int64
    contains
        procedure :: valid => rng_valid
        procedure :: uniform => rng_uniform
    end type rotation_rng_type

    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 16807_int64

    type, public :: rotation_stack_type
        private
        type(rotation_type), allocatable :: values(:)
    contains
        procedure :: valid => stack_valid
        procedure :: size => stack_size
        procedure :: at => stack_at
        procedure :: slice => stack_slice
        procedure :: as_quat => stack_as_quat
        procedure :: as_matrix => stack_as_matrix
        procedure :: as_rotvec => stack_as_rotvec
        procedure :: as_mrp => stack_as_mrp
        procedure :: as_euler => stack_as_euler
        procedure :: as_davenport => stack_as_davenport
        procedure, private :: stack_apply_vector
        procedure, private :: stack_apply_vectors
        generic :: apply => stack_apply_vector, stack_apply_vectors
        procedure :: inverse => stack_inverse
        procedure, private :: stack_compose_stack
        procedure, private :: stack_compose_rotation
        generic :: compose => stack_compose_stack, stack_compose_rotation
        procedure :: magnitude => stack_magnitude
        procedure :: mean => stack_mean
        procedure :: power => stack_power
        procedure :: approx_equal => stack_approx_equal
    end type rotation_stack_type

    public :: rotation_stack_identity, rotation_stack_from_rotations
    public :: rotation_stack_from_quat, rotation_stack_from_rotvec
    public :: rotation_stack_from_matrix, rotation_stack_from_mrp
    public :: rotation_stack_from_euler, rotation_stack_from_davenport, rotation_stack_concatenate
    public :: rotation_rng, rotation_random, rotation_stack_random

contains

    pure subroutine mark(status, value)
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        integer, intent(in) :: value !! output status code
        if (present(status)) status = value
    end subroutine mark

    function invalid_rotation() result(r)
        type(rotation_type) :: r
        r = rotation_from_quat([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
    end function invalid_rotation

    pure logical function stack_valid(self) result(ok)
        class(rotation_stack_type), intent(in) :: self !! stack to validate
        integer :: i
        ok = allocated(self%values)
        if (.not. ok) return
        do i = 1, size(self%values)
            if (.not. self%values(i)%valid()) then
                ok = .false.
                return
            end if
        end do
    end function stack_valid

    pure integer function stack_size(self) result(n)
        class(rotation_stack_type), intent(in) :: self !! stack to inspect
        n = 0
        if (allocated(self%values)) n = size(self%values)
    end function stack_size


    ! The classic Park-Miller generator has a 2**31-2 period. Two updates
    ! produce an approximately double-precision uniform variate. It is
    ! reproducible, not cryptographically secure, and not NumPy-bit-compatible.
    pure function rotation_rng(seed) result(rng)
        integer(int64), intent(in) :: seed !! arbitrary 64-bit reproducibility seed
        type(rotation_rng_type) :: rng
        rng%state = modulo(seed, rng_modulus - 1_int64) + 1_int64
    end function rotation_rng

    pure logical function rng_valid(self) result(ok)
        class(rotation_rng_type), intent(in) :: self !! random stream to inspect
        ok = self%state > 0_int64 .and. self%state < rng_modulus
    end function rng_valid

    function rng_uniform(self) result(u)
        class(rotation_rng_type), intent(inout) :: self !! stream advanced twice per sample
        real(dp) :: u
        integer(int64) :: first, second
        if (.not. self%valid()) then
            u = -1.0_dp
            return
        end if
        ! rng_multiplier*(rng_modulus-1) < 2**45: no signed overflow.
        self%state = modulo(rng_multiplier*self%state, rng_modulus)
        first = self%state - 1_int64
        self%state = modulo(rng_multiplier*self%state, rng_modulus)
        second = self%state - 1_int64
        u = (real(first,dp) + real(second,dp)/real(rng_modulus-1_int64,dp)) &
            / real(rng_modulus-1_int64,dp)
        u = min(u, 1.0_dp - epsilon(1.0_dp))
    end function rng_uniform

    ! The Shoemake transform maps three independent U[0,1) variates to a
    ! Haar-uniform rotation of SO(3), using scalar-last quaternions.
    function rotation_random(rng, status) result(r)
        type(rotation_rng_type), intent(inout), optional :: rng !! explicit random stream; otherwise RANDOM_NUMBER
        integer, intent(out), optional :: status !! zero on success, one on invalid RNG
        type(rotation_type) :: r
        real(dp) :: u(3), q(4), pi, a, b, x, y
        integer :: i, code
        r = invalid_rotation()
        call mark(status, spatial_rotation_invalid_input)
        if (present(rng)) then
            if (.not. rng%valid()) return
            do i = 1, 3
                u(i) = rng%uniform()
            end do
        else
            call random_number(u)
        end if
        x = sqrt(max(0.0_dp,1.0_dp-u(1)))
        y = sqrt(u(1))
        pi = acos(-1.0_dp)
        a = 2.0_dp*pi*u(2)
        b = 2.0_dp*pi*u(3)
        q = [x*sin(a), x*cos(a), y*sin(b), y*cos(b)]
        r = rotation_from_quat(q, status=code)
        if (code /= spatial_rotation_success) return
        call mark(status, spatial_rotation_success)
    end function rotation_random

    function rotation_stack_random(n, rng, status) result(stack)
        integer, intent(in) :: n !! number of independent Haar-uniform rotations; >= 0
        type(rotation_rng_type), intent(inout), optional :: rng !! reproducible stream, advanced 3 samples per rotation
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rotation_stack_type) :: stack
        integer :: i, code
        call mark(status, spatial_rotation_invalid_input)
        if (n < 0) return
        if (present(rng)) then
            if (.not. rng%valid()) return
        end if
        allocate(stack%values(n))
        do i = 1, n
            stack%values(i) = rotation_random(rng, code)
            if (code /= spatial_rotation_success) then
                deallocate(stack%values)
                return
            end if
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_stack_random

    function rotation_stack_identity(n, status) result(stack)
        integer, intent(in) :: n !! number of rotations, >= 0
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: stack
        integer :: i
        call mark(status, spatial_rotation_invalid_input)
        if (n < 0) return
        allocate(stack%values(n))
        do i = 1, n
            stack%values(i) = rotation_identity()
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_stack_identity

    function rotation_stack_from_rotations(rotations, status) result(stack)
        type(rotation_type), intent(in) :: rotations(:) !! input rotations, at least one
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: stack
        integer :: i
        call mark(status, spatial_rotation_invalid_input)
        if (size(rotations) == 0) return
        do i = 1, size(rotations)
            if (.not. rotations(i)%valid()) return
        end do
        allocate(stack%values(size(rotations)))
        stack%values = rotations
        call mark(status, spatial_rotation_success)
    end function rotation_stack_from_rotations

    function rotation_stack_from_quat(quats, scalar_first, status) result(stack)
        real(dp), intent(in) :: quats(:,:) !! shape (n,4), scalar-last by default
        logical, intent(in), optional :: scalar_first !! quaternion w first when true
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: stack
        integer :: i, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(quats,2) /= 4 .or. size(quats,1) < 1) return
        allocate(stack%values(size(quats,1)))
        do i = 1, size(quats,1)
            stack%values(i) = rotation_from_quat(quats(i,:), scalar_first, code)
            if (code /= spatial_rotation_success) then
                deallocate(stack%values)
                return
            end if
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_stack_from_quat

    function rotation_stack_from_rotvec(rotvecs, degrees, status) result(stack)
        real(dp), intent(in) :: rotvecs(:,:) !! shape (n,3), one rotation vector per row
        logical, intent(in), optional :: degrees !! true if vectors are in degrees
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: stack
        integer :: i, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(rotvecs,2) /= 3 .or. size(rotvecs,1) < 1) return
        allocate(stack%values(size(rotvecs,1)))
        do i = 1, size(rotvecs,1)
            stack%values(i) = rotation_from_rotvec(rotvecs(i,:), degrees, code)
            if (code /= spatial_rotation_success) then
                deallocate(stack%values)
                return
            end if
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_stack_from_rotvec

    function rotation_stack_from_matrix(matrices, status) result(stack)
        real(dp), intent(in) :: matrices(:,:,:) !! shape (n,3,3), matrix for each rotation
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: stack
        integer :: i, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(matrices,2) /= 3 .or. size(matrices,3) /= 3 .or. size(matrices,1) < 1) return
        allocate(stack%values(size(matrices,1)))
        do i = 1, size(matrices,1)
            stack%values(i) = rotation_from_matrix(matrices(i,:,:), code)
            if (code /= spatial_rotation_success) then
                deallocate(stack%values)
                return
            end if
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_stack_from_matrix

    function rotation_stack_from_mrp(mrps, status) result(stack)
        real(dp), intent(in) :: mrps(:,:) !! shape (n,3), MRP vectors
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: stack
        integer :: i, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(mrps,2) /= 3 .or. size(mrps,1) < 1) return
        allocate(stack%values(size(mrps,1)))
        do i = 1, size(mrps,1)
            stack%values(i) = rotation_from_mrp(mrps(i,:), code)
            if (code /= spatial_rotation_success) then
                deallocate(stack%values)
                return
            end if
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_stack_from_mrp

    function rotation_stack_from_euler(seq, angles, degrees, status) result(stack)
        character(len=*), intent(in) :: seq !! case-sensitive Euler axis sequence
        real(dp), intent(in) :: angles(:,:) !! shape (n,len(seq))
        logical, intent(in), optional :: degrees !! angles in degrees when true
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: stack
        integer :: i, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(angles,2) /= len(seq) .or. size(angles,1) < 1) return
        allocate(stack%values(size(angles,1)))
        do i = 1, size(angles,1)
            stack%values(i) = rotation_from_euler(seq, angles(i,:), degrees, code)
            if (code /= spatial_rotation_success) then
                deallocate(stack%values)
                return
            end if
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_stack_from_euler


    function rotation_stack_from_davenport(axes, order, angles, degrees, status) result(stack)
        real(dp), intent(in) :: axes(:,:) !! 1-3 input axes; shape (k,3)
        character(len=*), intent(in) :: order !! 'e'/'extrinsic' or 'i'/'intrinsic'
        real(dp), intent(in) :: angles(:,:) !! angles in shape (n,k)
        logical, intent(in), optional :: degrees !! true for degree arguments
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: stack
        integer :: i, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(angles,1) < 1 .or. size(angles,2) /= size(axes,1)) return
        allocate(stack%values(size(angles,1)))
        do i = 1, size(angles,1)
            stack%values(i) = rotation_from_davenport(axes,order,angles(i,:),degrees,code)
            if (code /= spatial_rotation_success) then
                deallocate(stack%values)
                return
            end if
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_stack_from_davenport

    function stack_at(self, index, status) result(r)
        class(rotation_stack_type), intent(in) :: self !! stack to inspect
        integer, intent(in) :: index !! one-based element index
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_type) :: r
        r = invalid_rotation()
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (index < 1 .or. index > self%size()) return
        r = self%values(index)
        call mark(status, spatial_rotation_success)
    end function stack_at

    function stack_slice(self, indices, status) result(out)
        class(rotation_stack_type), intent(in) :: self !! stack to subset
        integer, intent(in) :: indices(:) !! one-based indices, may repeat and reorder
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: out
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (any(indices < 1) .or. any(indices > self%size())) return
        allocate(out%values(size(indices)))
        out%values = self%values(indices)
        call mark(status, spatial_rotation_success)
    end function stack_slice

    function stack_as_quat(self, canonical, scalar_first) result(quats)
        class(rotation_stack_type), intent(in) :: self !! input stack
        logical, intent(in), optional :: canonical !! canonical signs when true
        logical, intent(in), optional :: scalar_first !! scalar first when true
        real(dp), allocatable :: quats(:,:) !! shape (n,4) on success, (0,4) if invalid
        integer :: i
        allocate(quats(0,4))
        if (.not. self%valid()) return
        deallocate(quats)
        allocate(quats(self%size(),4))
        do i = 1, self%size()
            quats(i,:) = self%values(i)%as_quat(canonical, scalar_first)
        end do
    end function stack_as_quat

    function stack_as_matrix(self) result(matrices)
        class(rotation_stack_type), intent(in) :: self !! input stack
        real(dp), allocatable :: matrices(:,:,:) !! shape (n,3,3) on success
        integer :: i
        allocate(matrices(0,3,3))
        if (.not. self%valid()) return
        deallocate(matrices)
        allocate(matrices(self%size(),3,3))
        do i = 1, self%size()
            matrices(i,:,:) = self%values(i)%as_matrix()
        end do
    end function stack_as_matrix

    function stack_as_rotvec(self, degrees) result(rotvecs)
        class(rotation_stack_type), intent(in) :: self !! input stack
        logical, intent(in), optional :: degrees !! vectors in degrees when true
        real(dp), allocatable :: rotvecs(:,:) !! shape (n,3) on success
        integer :: i
        allocate(rotvecs(0,3))
        if (.not. self%valid()) return
        deallocate(rotvecs)
        allocate(rotvecs(self%size(),3))
        do i = 1, self%size()
            rotvecs(i,:) = self%values(i)%as_rotvec(degrees)
        end do
    end function stack_as_rotvec

    function stack_as_mrp(self) result(mrps)
        class(rotation_stack_type), intent(in) :: self !! input stack
        real(dp), allocatable :: mrps(:,:) !! shape (n,3) on success
        integer :: i
        allocate(mrps(0,3))
        if (.not. self%valid()) return
        deallocate(mrps)
        allocate(mrps(self%size(),3))
        do i = 1, self%size()
            mrps(i,:) = self%values(i)%as_mrp()
        end do
    end function stack_as_mrp

    function stack_as_euler(self, seq, degrees, gimbal_lock, status) result(angles)
        class(rotation_stack_type), intent(in) :: self !! input stack
        character(len=*), intent(in) :: seq !! three axis letters
        logical, intent(in), optional :: degrees !! output degrees when true
        logical, allocatable, intent(out), optional :: gimbal_lock(:) !! lock flags per rotation
        integer, intent(out), optional :: status !! zero on success, one otherwise
        real(dp), allocatable :: angles(:,:) !! shape (n,3) on success
        logical :: locked
        integer :: i, code
        allocate(angles(0,3))
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        deallocate(angles)
        allocate(angles(self%size(),3))
        if (present(gimbal_lock)) allocate(gimbal_lock(self%size()))
        do i = 1, self%size()
            angles(i,:) = self%values(i)%as_euler(seq,degrees,locked,code)
            if (code /= spatial_rotation_success) then
                deallocate(angles)
                allocate(angles(0,3))
                if (present(gimbal_lock)) then
                    deallocate(gimbal_lock)
                    allocate(gimbal_lock(0))
                end if
                return
            end if
            if (present(gimbal_lock)) gimbal_lock(i) = locked
        end do
        call mark(status, spatial_rotation_success)
    end function stack_as_euler


    function stack_as_davenport(self, axes, order, degrees, gimbal_lock, status) result(angles)
        class(rotation_stack_type), intent(in) :: self !! rotation stack
        real(dp), intent(in) :: axes(:,:) !! three Davenport axes, shape (3,3)
        character(len=*), intent(in) :: order !! intrinsic or extrinsic sequence
        logical, intent(in), optional :: degrees !! true for degrees
        logical, allocatable, intent(out), optional :: gimbal_lock(:) !! one flag per orientation
        integer, intent(out), optional :: status !! zero on success, one otherwise
        real(dp), allocatable :: angles(:,:) !! shape (n,3), or (0,3) on error
        logical :: locked
        integer :: i, code
        allocate(angles(0,3))
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        deallocate(angles)
        allocate(angles(self%size(),3))
        if (present(gimbal_lock)) allocate(gimbal_lock(self%size()))
        do i = 1, self%size()
            angles(i,:) = self%values(i)%as_davenport(axes,order,degrees,locked,code)
            if (code /= spatial_rotation_success) then
                deallocate(angles)
                allocate(angles(0,3))
                if (present(gimbal_lock)) then
                    deallocate(gimbal_lock)
                    allocate(gimbal_lock(0))
                end if
                return
            end if
            if (present(gimbal_lock)) gimbal_lock(i) = locked
        end do
        call mark(status, spatial_rotation_success)
    end function stack_as_davenport

    ! Broadcast a common vector across every rotation.
    function stack_apply_vector(self, vector, inverse, status) result(out)
        class(rotation_stack_type), intent(in) :: self !! input rotation stack
        real(dp), intent(in) :: vector(:) !! single 3D vector
        logical, intent(in), optional :: inverse !! inverse rotations when true
        integer, intent(out), optional :: status !! zero on success, one otherwise
        real(dp), allocatable :: out(:,:) !! shape (n,3)
        integer :: i
        allocate(out(0,3))
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (size(vector) /= 3) return
        if (.not. all(ieee_is_finite(vector))) return
        deallocate(out)
        allocate(out(self%size(),3))
        do i = 1, self%size()
            out(i,:) = self%values(i)%apply(vector,inverse)
        end do
        call mark(status, spatial_rotation_success)
    end function stack_apply_vector

    ! NumPy/SciPy leading-axis broadcast: either both sizes match or one is 1.
    function stack_apply_vectors(self, vectors, inverse, status) result(out)
        class(rotation_stack_type), intent(in) :: self !! input rotation stack
        real(dp), intent(in) :: vectors(:,:) !! vectors in rows, shape (m,3)
        logical, intent(in), optional :: inverse !! inverse rotations when true
        integer, intent(out), optional :: status !! zero on success, one otherwise
        real(dp), allocatable :: out(:,:) !! shape (max(n,m),3)
        integer :: i, n, m, ir, iv
        real(dp) :: vector(3)
        type(rotation_type) :: r
        allocate(out(0,3))
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (size(vectors,2) /= 3) return
        if (.not. all(ieee_is_finite(vectors))) return
        n = self%size()
        m = size(vectors,1)
        if (n /= m .and. n /= 1 .and. m /= 1) return
        if (n == 0 .or. m == 0) return
        deallocate(out)
        allocate(out(max(n,m),3))
        do i = 1, size(out,1)
            ir = min(i,n)
            iv = min(i,m)
            vector = vectors(iv,:)
            r = self%values(ir)
            out(i,:) = r%apply(vector,inverse)
        end do
        call mark(status, spatial_rotation_success)
    end function stack_apply_vectors

    function stack_inverse(self, status) result(out)
        class(rotation_stack_type), intent(in) :: self !! input stack
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: out
        integer :: i
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        allocate(out%values(self%size()))
        do i = 1, self%size()
            out%values(i) = self%values(i)%inverse()
        end do
        call mark(status, spatial_rotation_success)
    end function stack_inverse

    function stack_compose_stack(self, other, status) result(out)
        class(rotation_stack_type), intent(in) :: self !! left rotation stack
        type(rotation_stack_type), intent(in) :: other !! right rotation stack
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: out
        integer :: i, n, m
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid() .or. .not. other%valid()) return
        n = self%size()
        m = other%size()
        if (n /= m .and. n /= 1 .and. m /= 1) return
        if (n == 0 .or. m == 0) return
        allocate(out%values(max(n,m)))
        do i = 1, size(out%values)
            out%values(i) = self%values(min(i,n))%compose(other%values(min(i,m)))
        end do
        call mark(status, spatial_rotation_success)
    end function stack_compose_stack

    function stack_compose_rotation(self, other, status) result(out)
        class(rotation_stack_type), intent(in) :: self !! left rotation stack
        type(rotation_type), intent(in) :: other !! one right rotation
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: out
        integer :: i
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid() .or. .not. other%valid()) return
        allocate(out%values(self%size()))
        do i = 1, self%size()
            out%values(i) = self%values(i)%compose(other)
        end do
        call mark(status, spatial_rotation_success)
    end function stack_compose_rotation

    function stack_magnitude(self) result(angles)
        class(rotation_stack_type), intent(in) :: self !! input rotation stack
        real(dp), allocatable :: angles(:) !! principal angles in radians
        integer :: i
        allocate(angles(0))
        if (.not. self%valid()) return
        deallocate(angles)
        allocate(angles(self%size()))
        do i = 1, self%size()
            angles(i) = self%values(i)%magnitude()
        end do
    end function stack_magnitude

    function stack_mean(self, weights, status) result(r)
        class(rotation_stack_type), intent(in) :: self !! stack to average
        real(dp), intent(in), optional :: weights(:) !! nonnegative weights
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_type) :: r
        r = invalid_rotation()
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (self%size() < 1) return
        r = rotation_mean(self%values, weights, status)
    end function stack_mean

    function stack_power(self, exponent, status) result(out)
        class(rotation_stack_type), intent(in) :: self !! stack to exponentiate
        real(dp), intent(in) :: exponent !! common exponent for all rotations
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: out
        integer :: i, code
        real(dp) :: vec(3)
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid() .or. .not. ieee_is_finite(exponent)) return
        allocate(out%values(self%size()))
        do i = 1, self%size()
            vec = self%values(i)%as_rotvec()
            if (.not. all(ieee_is_finite(vec*exponent))) then
                deallocate(out%values)
                return
            end if
            out%values(i) = rotation_from_rotvec(vec*exponent, status=code)
            if (code /= spatial_rotation_success) then
                deallocate(out%values)
                return
            end if
        end do
        call mark(status, spatial_rotation_success)
    end function stack_power

    function stack_approx_equal(self, other, atol, degrees, status) result(equal)
        class(rotation_stack_type), intent(in) :: self !! first rotation stack
        type(rotation_stack_type), intent(in) :: other !! second rotation stack
        real(dp), intent(in), optional :: atol !! maximum angular difference, default 1e-8 rad
        logical, intent(in), optional :: degrees !! atol given in degrees when true
        integer, intent(out), optional :: status !! zero on success, one otherwise
        logical, allocatable :: equal(:) !! one flag per broadcasted pair
        type(rotation_type) :: difference
        real(dp) :: threshold
        integer :: i, n, m
        allocate(equal(0))
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid() .or. .not. other%valid()) return
        n = self%size()
        m = other%size()
        if (n /= m .and. n /= 1 .and. m /= 1) return
        if (n == 0 .or. m == 0) return
        threshold = 1.0e-8_dp
        if (present(atol)) threshold = atol
        if (.not. ieee_is_finite(threshold) .or. threshold < 0.0_dp) return
        if (present(degrees)) then
            if (degrees) threshold = threshold*(acos(-1.0_dp)/180.0_dp)
        end if
        deallocate(equal)
        allocate(equal(max(n,m)))
        do i = 1, size(equal)
            difference = self%values(min(i,n))%compose(other%values(min(i,m))%inverse())
            equal(i) = difference%magnitude() <= threshold
        end do
        call mark(status, spatial_rotation_success)
    end function stack_approx_equal

    function rotation_stack_concatenate(stacks, status) result(out)
        type(rotation_stack_type), intent(in) :: stacks(:) !! stacks to concatenate
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_stack_type) :: out
        integer :: i, j, n
        call mark(status, spatial_rotation_invalid_input)
        n = 0
        do i = 1, size(stacks)
            if (.not. stacks(i)%valid()) return
            n = n + stacks(i)%size()
        end do
        allocate(out%values(n))
        j = 0
        do i = 1, size(stacks)
            n = stacks(i)%size()
            if (n == 0) cycle
            out%values(j+1:j+n) = stacks(i)%values
            j = j+n
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_stack_concatenate

end module scifort_spatial_rotation_stack
