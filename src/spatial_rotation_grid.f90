! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! A rank-two grid of 3D rotations. In SciPy >=1.17 rotations can have
! more than one batch dimension; this module implements the rank-two case.
! Its last representation dimension always holds quaternion/vector components.
module scifort_spatial_rotation_grid
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_quat, &
        rotation_from_rotvec, rotation_from_mrp, rotation_from_matrix, &
        rotation_from_euler, rotation_from_davenport, rotation_mean, spatial_rotation_success, &
        spatial_rotation_invalid_input
    use scifort_spatial_rotation_stack, only : rotation_stack_type, &
        rotation_stack_from_rotations, rotation_rng_type, rotation_random
    implicit none
    private

    type, public :: rotation_grid_type
        private
        type(rotation_type), allocatable :: values(:,:)
    contains
        procedure :: valid => grid_valid
        procedure :: shape => grid_shape
        procedure :: at => grid_at
        procedure :: as_quat => grid_as_quat
        procedure :: as_rotvec => grid_as_rotvec
        procedure :: as_matrix => grid_as_matrix
        procedure :: as_mrp => grid_as_mrp
        procedure :: as_euler => grid_as_euler
        procedure :: as_davenport => grid_as_davenport
        procedure, private :: grid_apply_vector
        procedure, private :: grid_apply_vectors
        generic :: apply => grid_apply_vector, grid_apply_vectors
        procedure :: inverse => grid_inverse
        procedure, private :: grid_compose_grid
        procedure, private :: grid_compose_rotation
        generic :: compose => grid_compose_grid, grid_compose_rotation
        procedure :: magnitude => grid_magnitude
        procedure :: mean => grid_mean
    end type rotation_grid_type

    public :: rotation_grid_from_quat, rotation_grid_from_rotvec
    public :: rotation_grid_from_matrix, rotation_grid_from_mrp, rotation_grid_from_euler
    public :: rotation_grid_from_rotations, rotation_grid_from_davenport, rotation_grid_random

contains

    pure subroutine mark(status, value)
        integer, intent(out), optional :: status !! zero on success, one on invalid inputs
        integer, intent(in) :: value !! integer status value
        if (present(status)) status = value
    end subroutine mark

    pure logical function grid_valid(self) result(ok)
        class(rotation_grid_type), intent(in) :: self !! rotation grid to validate
        integer :: i, j
        ok = allocated(self%values)
        if (.not. ok) return
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                if (.not. self%values(i,j)%valid()) then
                    ok = .false.
                    return
                end if
            end do
        end do
    end function grid_valid

    pure function grid_shape(self) result(dims)
        class(rotation_grid_type), intent(in) :: self !! rotation grid to inspect
        integer :: dims(2) !! batch dimension lengths, (0,0) when invalid
        dims = 0
        if (allocated(self%values)) dims = shape(self%values)
    end function grid_shape

    function rotation_grid_from_rotations(rotations, status) result(out)
        type(rotation_type), intent(in) :: rotations(:,:) !! scalar rotations in an (n,m) grid
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j
        call mark(status, spatial_rotation_invalid_input)
        if (any(shape(rotations) < 1)) return
        do j = 1, size(rotations,2)
            do i = 1, size(rotations,1)
                if (.not. rotations(i,j)%valid()) return
            end do
        end do
        out%values = rotations
        call mark(status, spatial_rotation_success)
    end function rotation_grid_from_rotations

    function rotation_grid_from_quat(quats, scalar_first, status) result(out)
        real(dp), intent(in) :: quats(:,:,:) !! quaternions, shape (n,m,4)
        logical, intent(in), optional :: scalar_first !! w-first when true
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(quats,3) /= 4 .or. size(quats,1) < 1 .or. size(quats,2) < 1) return
        allocate(out%values(size(quats,1),size(quats,2)))
        do j = 1, size(quats,2)
            do i = 1, size(quats,1)
                out%values(i,j) = rotation_from_quat(quats(i,j,:), scalar_first, code)
                if (code /= spatial_rotation_success) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_grid_from_quat

    function rotation_grid_from_rotvec(rotvecs, degrees, status) result(out)
        real(dp), intent(in) :: rotvecs(:,:,:) !! rotation vectors, shape (n,m,3)
        logical, intent(in), optional :: degrees !! angle unit flag; radians by default
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(rotvecs,3) /= 3 .or. size(rotvecs,1) < 1 .or. size(rotvecs,2) < 1) return
        allocate(out%values(size(rotvecs,1),size(rotvecs,2)))
        do j = 1, size(rotvecs,2)
            do i = 1, size(rotvecs,1)
                out%values(i,j) = rotation_from_rotvec(rotvecs(i,j,:), degrees, code)
                if (code /= spatial_rotation_success) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_grid_from_rotvec


    function rotation_grid_from_mrp(mrps, status) result(out)
        real(dp), intent(in) :: mrps(:,:,:) !! modified Rodrigues vectors, shape (n,m,3)
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(mrps,3) /= 3 .or. size(mrps,1) < 1 .or. size(mrps,2) < 1) return
        allocate(out%values(size(mrps,1),size(mrps,2)))
        do j = 1, size(mrps,2)
            do i = 1, size(mrps,1)
                out%values(i,j) = rotation_from_mrp(mrps(i,j,:), code)
                if (code /= spatial_rotation_success) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_grid_from_mrp

    function rotation_grid_from_matrix(matrices, status) result(out)
        real(dp), intent(in) :: matrices(:,:,:,:) !! proper orthogonal matrices, shape (n,m,3,3)
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j, code
        call mark(status, spatial_rotation_invalid_input)
        if (any(shape(matrices) < 1)) return
        if (size(matrices,3) /= 3 .or. size(matrices,4) /= 3) return
        allocate(out%values(size(matrices,1),size(matrices,2)))
        do j = 1, size(matrices,2)
            do i = 1, size(matrices,1)
                out%values(i,j) = rotation_from_matrix(matrices(i,j,:,:), code)
                if (code /= spatial_rotation_success) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_grid_from_matrix

    function rotation_grid_from_euler(seq, angles, degrees, status) result(out)
        character(len=*), intent(in) :: seq !! case-sensitive Euler sequence, length 1 to 3
        real(dp), intent(in) :: angles(:,:,:) !! Euler angles, shape (n,m,len(seq))
        logical, intent(in), optional :: degrees !! angles in degrees when true
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(angles,1) < 1 .or. size(angles,2) < 1) return
        if (size(angles,3) /= len(seq)) return
        allocate(out%values(size(angles,1),size(angles,2)))
        do j = 1, size(angles,2)
            do i = 1, size(angles,1)
                out%values(i,j) = rotation_from_euler(seq, angles(i,j,:), degrees, code)
                if (code /= spatial_rotation_success) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_grid_from_euler


    function rotation_grid_from_davenport(axes, order, angles, degrees, status) result(out)
        real(dp), intent(in) :: axes(:,:) !! Davenport axes (k,3), k=1..3
        character(len=*), intent(in) :: order !! intrinsic or extrinsic sequence
        real(dp), intent(in) :: angles(:,:,:) !! grid angles, shape (n,m,k)
        logical, intent(in), optional :: degrees !! true for degree arguments
        integer, intent(out), optional :: status !! zero on success, one otherwise
        type(rotation_grid_type) :: out
        integer :: i, j, code
        call mark(status, spatial_rotation_invalid_input)
        if (size(angles,1) < 1 .or. size(angles,2) < 1) return
        if (size(angles,3) /= size(axes,1)) return
        allocate(out%values(size(angles,1),size(angles,2)))
        do j = 1, size(angles,2)
            do i = 1, size(angles,1)
                out%values(i,j) = rotation_from_davenport(axes,order,angles(i,j,:),degrees,code)
                if (code /= spatial_rotation_success) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_grid_from_davenport

    function rotation_grid_random(dims, rng, status) result(out)
        integer, intent(in) :: dims(2) !! positive lengths for two batch axes
        type(rotation_rng_type), intent(inout), optional :: rng !! reproducible stream; otherwise intrinsic RNG
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j, code
        call mark(status, spatial_rotation_invalid_input)
        if (any(dims < 1)) return
        if (present(rng)) then
            if (.not. rng%valid()) return
        end if
        allocate(out%values(dims(1),dims(2)))
        ! Fortran column-major order; one stream shared by every element.
        do j = 1, dims(2)
            do i = 1, dims(1)
                out%values(i,j) = rotation_random(rng, code)
                if (code /= spatial_rotation_success) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function rotation_grid_random

    function grid_at(self, i, j, status) result(r)
        class(rotation_grid_type), intent(in) :: self !! grid to inspect
        integer, intent(in) :: i !! one-based first batch coordinate
        integer, intent(in) :: j !! one-based second batch coordinate
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_type) :: r
        r = rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (i < 1 .or. i > size(self%values,1)) return
        if (j < 1 .or. j > size(self%values,2)) return
        r = self%values(i,j)
        call mark(status, spatial_rotation_success)
    end function grid_at

    function grid_as_quat(self, canonical, scalar_first) result(q)
        class(rotation_grid_type), intent(in) :: self !! source rotation grid
        logical, intent(in), optional :: canonical !! canonical quaternion signs when true
        logical, intent(in), optional :: scalar_first !! return w first when true
        real(dp), allocatable :: q(:,:,:) !! shape (n,m,4); empty on invalid input
        integer :: i, j
        allocate(q(0,0,4))
        if (.not. self%valid()) return
        deallocate(q)
        allocate(q(size(self%values,1),size(self%values,2),4))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                q(i,j,:) = self%values(i,j)%as_quat(canonical,scalar_first)
            end do
        end do
    end function grid_as_quat

    function grid_as_rotvec(self, degrees) result(v)
        class(rotation_grid_type), intent(in) :: self !! source rotation grid
        logical, intent(in), optional :: degrees !! return degrees when true
        real(dp), allocatable :: v(:,:,:) !! shape (n,m,3); empty on invalid input
        integer :: i, j
        allocate(v(0,0,3))
        if (.not. self%valid()) return
        deallocate(v)
        allocate(v(size(self%values,1),size(self%values,2),3))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                v(i,j,:) = self%values(i,j)%as_rotvec(degrees)
            end do
        end do
    end function grid_as_rotvec

    function grid_as_matrix(self) result(m)
        class(rotation_grid_type), intent(in) :: self !! source rotation grid
        real(dp), allocatable :: m(:,:,:,:) !! shape (n,m,3,3); empty on invalid input
        integer :: i, j
        allocate(m(0,0,3,3))
        if (.not. self%valid()) return
        deallocate(m)
        allocate(m(size(self%values,1),size(self%values,2),3,3))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                m(i,j,:,:) = self%values(i,j)%as_matrix()
            end do
        end do
    end function grid_as_matrix


    function grid_as_mrp(self) result(v)
        class(rotation_grid_type), intent(in) :: self !! source rotation grid
        real(dp), allocatable :: v(:,:,:) !! modified Rodrigues vectors, shape (n,m,3)
        integer :: i, j
        allocate(v(0,0,3))
        if (.not. self%valid()) return
        deallocate(v)
        allocate(v(size(self%values,1),size(self%values,2),3))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                v(i,j,:) = self%values(i,j)%as_mrp()
            end do
        end do
    end function grid_as_mrp

    function grid_as_euler(self, seq, degrees, gimbal_lock, status) result(angles)
        class(rotation_grid_type), intent(in) :: self !! source rotation grid
        character(len=*), intent(in) :: seq !! case-sensitive three-axis Euler sequence
        logical, intent(in), optional :: degrees !! output degrees when true
        logical, allocatable, intent(out), optional :: gimbal_lock(:,:) !! shape (n,m), singularity flags
        integer, intent(out), optional :: status !! success=0; invalid=1
        real(dp), allocatable :: angles(:,:,:) !! shape (n,m,3) when valid
        integer :: i, j, code
        logical :: locked
        allocate(angles(0,0,3))
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        deallocate(angles)
        allocate(angles(size(self%values,1),size(self%values,2),3))
        if (present(gimbal_lock)) allocate(gimbal_lock(size(self%values,1),size(self%values,2)))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                angles(i,j,:) = self%values(i,j)%as_euler(seq, degrees, locked, code)
                if (code /= spatial_rotation_success) then
                    deallocate(angles)
                    allocate(angles(0,0,3))
                    if (present(gimbal_lock)) then
                        deallocate(gimbal_lock)
                        allocate(gimbal_lock(0,0))
                    end if
                    return
                end if
                if (present(gimbal_lock)) gimbal_lock(i,j) = locked
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function grid_as_euler


    function grid_as_davenport(self, axes, order, degrees, gimbal_lock, status) result(angles)
        class(rotation_grid_type), intent(in) :: self !! rotation grid
        real(dp), intent(in) :: axes(:,:) !! Davenport axes (3,3)
        character(len=*), intent(in) :: order !! intrinsic or extrinsic sequence
        logical, intent(in), optional :: degrees !! true for degrees
        logical, allocatable, intent(out), optional :: gimbal_lock(:,:) !! lock flags, shape (n,m)
        integer, intent(out), optional :: status !! zero on success, one otherwise
        real(dp), allocatable :: angles(:,:,:) !! shape (n,m,3), or (0,0,3) on error
        logical :: locked
        integer :: i, j, code
        allocate(angles(0,0,3))
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        deallocate(angles)
        allocate(angles(size(self%values,1),size(self%values,2),3))
        if (present(gimbal_lock)) allocate(gimbal_lock(size(self%values,1),size(self%values,2)))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                angles(i,j,:) = self%values(i,j)%as_davenport(axes,order,degrees,locked,code)
                if (code /= spatial_rotation_success) then
                    deallocate(angles)
                    allocate(angles(0,0,3))
                    if (present(gimbal_lock)) then
                        deallocate(gimbal_lock)
                        allocate(gimbal_lock(0,0))
                    end if
                    return
                end if
                if (present(gimbal_lock)) gimbal_lock(i,j) = locked
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function grid_as_davenport

    function grid_apply_vector(self, vector, inverse, status) result(v)
        class(rotation_grid_type), intent(in) :: self !! source rotation grid
        real(dp), intent(in) :: vector(:) !! one 3-component vector, broadcast everywhere
        logical, intent(in), optional :: inverse !! invert rotations before application
        integer, intent(out), optional :: status !! success=0; invalid=1
        real(dp), allocatable :: v(:,:,:) !! rotated vectors, shape (n,m,3)
        integer :: i, j
        allocate(v(0,0,3))
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (size(vector) /= 3) return
        if (.not. all(ieee_is_finite(vector))) return
        deallocate(v)
        allocate(v(size(self%values,1),size(self%values,2),3))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                v(i,j,:) = self%values(i,j)%apply(vector,inverse)
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function grid_apply_vector

    function grid_apply_vectors(self, vectors, inverse, status) result(v)
        class(rotation_grid_type), intent(in) :: self !! source rotation grid
        real(dp), intent(in) :: vectors(:,:,:) !! vectors, shape (k,l,3)
        logical, intent(in), optional :: inverse !! invert rotations before application
        integer, intent(out), optional :: status !! success=0; invalid=1
        real(dp), allocatable :: v(:,:,:) !! shape (max(n,k),max(m,l),3) on success
        integer :: i, j, n, m, k, l
        real(dp) :: point(3)
        type(rotation_type) :: r
        allocate(v(0,0,3))
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (size(vectors,3) /= 3) return
        if (.not. all(ieee_is_finite(vectors))) return
        n = size(self%values,1)
        m = size(self%values,2)
        k = size(vectors,1)
        l = size(vectors,2)
        if (k < 1 .or. l < 1) return
        if (n /= k .and. n /= 1 .and. k /= 1) return
        if (m /= l .and. m /= 1 .and. l /= 1) return
        deallocate(v)
        allocate(v(max(n,k),max(m,l),3))
        do j = 1, size(v,2)
            do i = 1, size(v,1)
                r = self%values(min(i,n),min(j,m))
                point = vectors(min(i,k),min(j,l),:)
                v(i,j,:) = r%apply(point,inverse)
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function grid_apply_vectors

    function grid_inverse(self, status) result(out)
        class(rotation_grid_type), intent(in) :: self !! source rotation grid
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        allocate(out%values(size(self%values,1),size(self%values,2)))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                out%values(i,j) = self%values(i,j)%inverse()
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function grid_inverse

    function grid_compose_grid(self, other, status) result(out)
        class(rotation_grid_type), intent(in) :: self !! left rotation grid
        type(rotation_grid_type), intent(in) :: other !! right rotation grid
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j, n, m, k, l
        type(rotation_type) :: first, second
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid() .or. .not. other%valid()) return
        n = size(self%values,1)
        m = size(self%values,2)
        k = size(other%values,1)
        l = size(other%values,2)
        if (n /= k .and. n /= 1 .and. k /= 1) return
        if (m /= l .and. m /= 1 .and. l /= 1) return
        allocate(out%values(max(n,k),max(m,l)))
        do j = 1, size(out%values,2)
            do i = 1, size(out%values,1)
                first = self%values(min(i,n),min(j,m))
                second = other%values(min(i,k),min(j,l))
                out%values(i,j) = first%compose(second)
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function grid_compose_grid

    function grid_compose_rotation(self, other, status) result(out)
        class(rotation_grid_type), intent(in) :: self !! left rotation grid
        type(rotation_type), intent(in) :: other !! one right rotation
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_grid_type) :: out
        integer :: i, j
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid() .or. .not. other%valid()) return
        allocate(out%values(size(self%values,1),size(self%values,2)))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                out%values(i,j) = self%values(i,j)%compose(other)
            end do
        end do
        call mark(status, spatial_rotation_success)
    end function grid_compose_rotation

    function grid_magnitude(self) result(angles)
        class(rotation_grid_type), intent(in) :: self !! rotation grid
        real(dp), allocatable :: angles(:,:) !! angular magnitudes in radians
        integer :: i, j
        allocate(angles(0,0))
        if (.not. self%valid()) return
        deallocate(angles)
        allocate(angles(size(self%values,1),size(self%values,2)))
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                angles(i,j) = self%values(i,j)%magnitude()
            end do
        end do
    end function grid_magnitude

    ! Weighted chordal L2 mean along one batch axis. The result is a
    ! rank-one rotation stack, matching scipy Rotation.mean(axis=...).
    function grid_mean(self, axis, weights, status) result(out)
        class(rotation_grid_type), intent(in) :: self !! source grid
        integer, intent(in) :: axis !! one-based reduction axis: 1 or 2
        real(dp), intent(in), optional :: weights(:) !! weights along reduced axis
        integer, intent(out), optional :: status !! success=0; invalid=1
        type(rotation_stack_type) :: out
        type(rotation_type), allocatable :: means(:)
        integer :: i, n, code
        call mark(status, spatial_rotation_invalid_input)
        if (.not. self%valid()) return
        if (axis /= 1 .and. axis /= 2) return
        if (present(weights)) then
            if (size(weights) /= size(self%values,axis)) return
            if (.not. all(ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            if (.not. any(weights > 0.0_dp)) return
        end if
        n = size(self%values,3-axis)
        allocate(means(n))
        do i = 1, n
            if (axis == 1) then
                means(i) = rotation_mean(self%values(:,i), weights, code)
            else
                means(i) = rotation_mean(self%values(i,:), weights, code)
            end if
            if (code /= spatial_rotation_success) return
        end do
        out = rotation_stack_from_rotations(means, status=code)
        if (code /= spatial_rotation_success) return
        call mark(status, spatial_rotation_success)
    end function grid_mean

end module scifort_spatial_rotation_grid
