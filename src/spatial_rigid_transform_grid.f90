! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rank-two batches of SE(3) rigid transforms (SciPy RigidTransform shape (n,m)).
! Active transforms, scalar-last dual quaternions, and NumPy-style broadcast
! of the two leading dimensions. All public indexing is one-based.
module scifort_spatial_rigid_transform_grid
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, rotation_mean
    use scifort_spatial_rotation_grid, only : rotation_grid_type
    use scifort_spatial_rigid_transform, only : rigid_transform_type, rigid_transform_stack_type, &
        rigid_identity, rigid_from_components, rigid_from_matrix, rigid_from_exp_coords, &
        rigid_from_dual_quat, rigid_stack_from_components, &
        spatial_rigid_success, spatial_rigid_invalid_input
    implicit none
    private

    type, public :: rigid_transform_grid_type
        private
        type(rigid_transform_type), allocatable :: values(:,:)
    contains
        procedure :: valid => grid_valid
        procedure :: shape => grid_shape
        procedure :: at => grid_at
        procedure :: row => grid_row
        procedure :: column => grid_column
        procedure :: translation => grid_translation
        procedure :: as_matrix => grid_as_matrix
        procedure :: as_exp_coords => grid_as_exp_coords
        procedure :: as_dual_quat => grid_as_dual_quat
        procedure, private :: grid_apply_vector
        procedure, private :: grid_apply_points
        procedure, private :: grid_apply_grid_points
        generic :: apply => grid_apply_vector, grid_apply_points, grid_apply_grid_points
        procedure :: inverse => grid_inverse
        procedure, private :: grid_compose_grid
        procedure, private :: grid_compose_stack
        procedure, private :: grid_compose_single
        generic :: compose => grid_compose_grid, grid_compose_stack, grid_compose_single
        procedure :: power => grid_power
        procedure :: mean => grid_mean_all
        procedure :: mean_axis => grid_mean_axis
    end type rigid_transform_grid_type

    public :: rigid_grid_identity, rigid_grid_from_transforms
    public :: rigid_grid_from_components, rigid_grid_from_translation, rigid_grid_from_rotation
    public :: rigid_grid_from_matrix, rigid_grid_from_exp_coords, rigid_grid_from_dual_quat

    interface rigid_grid_from_components
        module procedure grid_from_components_array
        module procedure grid_from_components_rotation_grid
    end interface rigid_grid_from_components

contains

    pure subroutine set_status(status, code)
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        integer, intent(in) :: code !! status code to return
        if (present(status)) status = code
    end subroutine set_status

    pure logical function grid_valid(self) result(ok)
        class(rigid_transform_grid_type), intent(in) :: self !! grid to inspect
        integer :: i, j
        ok = .false.
        if (.not. allocated(self%values)) return
        if (any(shape(self%values) < 1)) return
        do j = 1, size(self%values,2)
            do i = 1, size(self%values,1)
                if (.not. self%values(i,j)%valid()) return
            end do
        end do
        ok = .true.
    end function grid_valid

    pure function grid_shape(self) result(dims)
        class(rigid_transform_grid_type), intent(in) :: self !! grid to inspect
        integer :: dims(2) !! grid dimensions, or zeros if undefined
        dims = 0
        if (allocated(self%values)) dims = shape(self%values)
    end function grid_shape

    pure function grid_at(self, i, j) result(tf)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        integer, intent(in) :: i !! one-based row
        integer, intent(in) :: j !! one-based column
        type(rigid_transform_type) :: tf
        ! Return an invalid scalar for out-of-bounds or undefined grids.
        if (.not. self%valid()) then
            tf = invalid_scalar()
            return
        end if
        if (i < 1 .or. j < 1 .or. i > size(self%values,1) .or. j > size(self%values,2)) then
            tf = invalid_scalar()
            return
        end if
        tf = self%values(i,j)
    end function grid_at

    pure function invalid_scalar() result(tf)
        type(rigid_transform_type) :: tf
        ! Reuse the scalar invalid-result convention through an empty stack.
        type(rigid_transform_stack_type) :: empty
        tf = empty%at(1)
    end function invalid_scalar

    function rigid_grid_identity(dims, status) result(grid)
        integer, intent(in) :: dims(2) !! positive grid dimensions (n,m)
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: grid
        integer :: i, j
        call set_status(status, spatial_rigid_invalid_input)
        if (any(dims < 1)) return
        allocate(grid%values(dims(1),dims(2)))
        do j = 1, dims(2)
            do i = 1, dims(1)
                grid%values(i,j) = rigid_identity()
            end do
        end do
        call set_status(status, spatial_rigid_success)
    end function rigid_grid_identity

    function rigid_grid_from_transforms(transforms, status) result(grid)
        type(rigid_transform_type), intent(in) :: transforms(:,:) !! grid (n,m)
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: grid
        integer :: i, j
        call set_status(status, spatial_rigid_invalid_input)
        if (any(shape(transforms) < 1)) return
        do j = 1, size(transforms,2)
            do i = 1, size(transforms,1)
                if (.not. transforms(i,j)%valid()) return
            end do
        end do
        grid%values = transforms
        call set_status(status, spatial_rigid_success)
    end function rigid_grid_from_transforms

    function grid_from_components_array(translations, rotations, status) result(grid)
        real(dp), intent(in) :: translations(:,:,:) !! translations (n,m,3), broadcast axes permitted
        type(rotation_type), intent(in) :: rotations(:,:) !! rotation grid (n,m), broadcast axes permitted
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: grid
        integer :: n, m, i, j, k, l, code
        call set_status(status, spatial_rigid_invalid_input)
        if (size(translations,3) /= 3) return
        if (size(translations,1) < 1 .or. size(translations,2) < 1) return
        if (any(shape(rotations) < 1)) return
        n = max(size(translations,1), size(rotations,1))
        m = max(size(translations,2), size(rotations,2))
        if (.not. compatible(size(translations,1),n) .or. .not. compatible(size(translations,2),m)) return
        if (.not. compatible(size(rotations,1),n) .or. .not. compatible(size(rotations,2),m)) return
        allocate(grid%values(n,m))
        do j = 1, m
            l = min(j,size(translations,2))
            do i = 1, n
                k = min(i,size(translations,1))
                grid%values(i,j) = rigid_from_components(translations(k,l,:), &
                    rotations(min(i,size(rotations,1)),min(j,size(rotations,2))), status=code)
                if (code /= spatial_rigid_success) then
                    deallocate(grid%values)
                    return
                end if
            end do
        end do
        call set_status(status, spatial_rigid_success)
    end function grid_from_components_array

    function grid_from_components_rotation_grid(translations, rotations, status) result(grid)
        real(dp), intent(in) :: translations(:,:,:) !! translations with shape (n,m,3)
        type(rotation_grid_type), intent(in) :: rotations !! existing SciFort rotation grid
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: grid
        type(rotation_type), allocatable :: r(:,:)
        integer :: dims(2), i, j
        call set_status(status, spatial_rigid_invalid_input)
        if (.not. rotations%valid()) return
        dims = rotations%shape()
        allocate(r(dims(1),dims(2)))
        do j = 1, dims(2)
            do i = 1, dims(1)
                r(i,j) = rotations%at(i,j)
            end do
        end do
        grid = grid_from_components_array(translations,r,status)
    end function grid_from_components_rotation_grid

    function rigid_grid_from_translation(translations, status) result(grid)
        real(dp), intent(in) :: translations(:,:,:) !! translations (n,m,3)
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: grid
        type(rotation_type) :: identity_rotation(1,1)
        identity_rotation(1,1) = rotation_identity()
        grid = grid_from_components_array(translations,identity_rotation,status)
    end function rigid_grid_from_translation

    function rigid_grid_from_rotation(rotations, status) result(grid)
        type(rotation_grid_type), intent(in) :: rotations !! grid of rotations (n,m)
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: grid
        real(dp) :: translation(1,1,3)
        translation = 0.0_dp
        grid = grid_from_components_rotation_grid(translation,rotations,status)
    end function rigid_grid_from_rotation

    function rigid_grid_from_matrix(matrices, status) result(grid)
        real(dp), intent(in) :: matrices(:,:,:,:) !! homogeneous matrices (n,m,4,4)
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: grid
        integer :: i, j, code
        call set_status(status, spatial_rigid_invalid_input)
        if (size(matrices,3) /= 4 .or. size(matrices,4) /= 4) return
        if (size(matrices,1) < 1 .or. size(matrices,2) < 1) return
        allocate(grid%values(size(matrices,1),size(matrices,2)))
        do j = 1, size(matrices,2)
            do i = 1, size(matrices,1)
                grid%values(i,j) = rigid_from_matrix(matrices(i,j,:,:),status=code)
                if (code /= spatial_rigid_success) then
                    deallocate(grid%values)
                    return
                end if
            end do
        end do
        call set_status(status, spatial_rigid_success)
    end function rigid_grid_from_matrix

    function rigid_grid_from_exp_coords(coords, status) result(grid)
        real(dp), intent(in) :: coords(:,:,:) !! SE(3) exponentials (n,m,6)
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: grid
        integer :: i, j, code
        call set_status(status, spatial_rigid_invalid_input)
        if (size(coords,3) /= 6 .or. size(coords,1) < 1 .or. size(coords,2) < 1) return
        allocate(grid%values(size(coords,1),size(coords,2)))
        do j = 1, size(coords,2)
            do i = 1, size(coords,1)
                grid%values(i,j) = rigid_from_exp_coords(coords(i,j,:),status=code)
                if (code /= spatial_rigid_success) then
                    deallocate(grid%values)
                    return
                end if
            end do
        end do
        call set_status(status, spatial_rigid_success)
    end function rigid_grid_from_exp_coords

    function rigid_grid_from_dual_quat(duals, scalar_first, status) result(grid)
        real(dp), intent(in) :: duals(:,:,:) !! dual quaternions (n,m,8)
        logical, intent(in), optional :: scalar_first !! quaternion scalar-first representation
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: grid
        integer :: i, j, code
        call set_status(status, spatial_rigid_invalid_input)
        if (size(duals,3) /= 8 .or. size(duals,1) < 1 .or. size(duals,2) < 1) return
        allocate(grid%values(size(duals,1),size(duals,2)))
        do j = 1, size(duals,2)
            do i = 1, size(duals,1)
                grid%values(i,j) = rigid_from_dual_quat(duals(i,j,:),scalar_first,status=code)
                if (code /= spatial_rigid_success) then
                    deallocate(grid%values)
                    return
                end if
            end do
        end do
        call set_status(status, spatial_rigid_success)
    end function rigid_grid_from_dual_quat

    pure logical function compatible(s,n)
        integer, intent(in) :: s !! actual input dimension
        integer, intent(in) :: n !! broadcast result dimension
        compatible = (s == 1 .or. s == n)
    end function compatible

    pure function grid_translation(self) result(t)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        real(dp), allocatable :: t(:,:,:) !! translation array (n,m,3)
        integer :: i, j
        allocate(t(0,0,3))
        if (.not. self%valid()) return
        deallocate(t)
        allocate(t(size(self%values,1),size(self%values,2),3))
        do j=1,size(t,2)
            do i=1,size(t,1)
                t(i,j,:) = self%values(i,j)%translation()
            end do
        end do
    end function grid_translation

    pure function grid_as_matrix(self) result(m)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        real(dp), allocatable :: m(:,:,:,:) !! matrices (n,m,4,4)
        integer :: i, j
        allocate(m(0,0,4,4))
        if (.not. self%valid()) return
        deallocate(m)
        allocate(m(size(self%values,1),size(self%values,2),4,4))
        do j=1,size(m,2)
            do i=1,size(m,1)
                m(i,j,:,:) = self%values(i,j)%as_matrix()
            end do
        end do
    end function grid_as_matrix

    pure function grid_as_exp_coords(self) result(coords)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        real(dp), allocatable :: coords(:,:,:) !! SE(3) exponential coordinates (n,m,6)
        integer :: i, j
        allocate(coords(0,0,6))
        if (.not. self%valid()) return
        deallocate(coords)
        allocate(coords(size(self%values,1),size(self%values,2),6))
        do j=1,size(coords,2)
            do i=1,size(coords,1)
                coords(i,j,:) = self%values(i,j)%as_exp_coords()
            end do
        end do
    end function grid_as_exp_coords

    pure function grid_as_dual_quat(self, scalar_first) result(duals)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        logical, intent(in), optional :: scalar_first !! use scalar-first quaternions
        real(dp), allocatable :: duals(:,:,:) !! dual quaternions (n,m,8)
        integer :: i, j
        allocate(duals(0,0,8))
        if (.not. self%valid()) return
        deallocate(duals)
        allocate(duals(size(self%values,1),size(self%values,2),8))
        do j=1,size(duals,2)
            do i=1,size(duals,1)
                duals(i,j,:) = self%values(i,j)%as_dual_quat(scalar_first)
            end do
        end do
    end function grid_as_dual_quat

    function grid_row(self,i) result(stack)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        integer, intent(in) :: i !! one-based first dimension index
        type(rigid_transform_stack_type) :: stack
        real(dp), allocatable :: t(:,:)
        type(rotation_type), allocatable :: r(:)
        integer :: j, m
        if (.not. self%valid()) return
        if (i < 1 .or. i > size(self%values,1)) return
        m = size(self%values,2)
        allocate(t(m,3),r(m))
        do j=1,m
            t(j,:) = self%values(i,j)%translation()
            r(j) = self%values(i,j)%rotation()
        end do
        stack = rigid_stack_from_components(t,r)
    end function grid_row

    function grid_column(self,j) result(stack)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        integer, intent(in) :: j !! one-based second dimension index
        type(rigid_transform_stack_type) :: stack
        real(dp), allocatable :: t(:,:)
        type(rotation_type), allocatable :: r(:)
        integer :: i, n
        if (.not. self%valid()) return
        if (j < 1 .or. j > size(self%values,2)) return
        n = size(self%values,1)
        allocate(t(n,3),r(n))
        do i=1,n
            t(i,:) = self%values(i,j)%translation()
            r(i) = self%values(i,j)%rotation()
        end do
        stack = rigid_stack_from_components(t,r)
    end function grid_column

    pure function grid_apply_vector(self, point, inverse) result(out)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        real(dp), intent(in) :: point(:) !! common point, length 3
        logical, intent(in), optional :: inverse !! apply inverse transforms
        real(dp), allocatable :: out(:,:,:) !! transformed points (n,m,3)
        integer :: i, j
        allocate(out(0,0,3))
        if (.not. self%valid() .or. size(point) /= 3) return
        deallocate(out)
        allocate(out(size(self%values,1),size(self%values,2),3))
        do j=1,size(out,2)
            do i=1,size(out,1)
                out(i,j,:) = self%values(i,j)%apply(point,inverse)
            end do
        end do
    end function grid_apply_vector

    pure function grid_apply_points(self, points, inverse) result(out)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        real(dp), intent(in) :: points(:,:) !! points (k,3), broadcast with final grid dimension
        logical, intent(in), optional :: inverse !! apply inverse transforms
        real(dp), allocatable :: out(:,:,:) !! result (n,max(m,k),3)
        integer :: i, j, n, m, k, pj
        real(dp) :: xyz(3)
        allocate(out(0,0,3))
        if (.not. self%valid() .or. size(points,2) /= 3 .or. size(points,1) < 1) return
        n = size(self%values,1)
        m = size(self%values,2)
        k = size(points,1)
        if (.not. compatible(k,m) .and. .not. compatible(m,k)) return
        deallocate(out)
        allocate(out(n,max(m,k),3))
        do j=1,size(out,2)
            do i=1,n
                pj = min(j,k)
                xyz = points(pj,:)
                out(i,j,:) = self%values(i,min(j,m))%apply(xyz,inverse)
            end do
        end do
    end function grid_apply_points

    pure function grid_apply_grid_points(self, points, inverse) result(out)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        real(dp), intent(in) :: points(:,:,:) !! Cartesian points (a,b,3), broadcast axes
        logical, intent(in), optional :: inverse !! apply inverse transforms
        real(dp), allocatable :: out(:,:,:) !! transformed points (max(n,a), max(m,b),3)
        integer :: n, m, a, b, i, j, pi, pj
        real(dp) :: xyz(3)
        allocate(out(0,0,3))
        if (.not. self%valid()) return
        if (size(points,3) /= 3 .or. size(points,1) < 1 .or. size(points,2) < 1) return
        n = size(self%values,1)
        m = size(self%values,2)
        a = size(points,1)
        b = size(points,2)
        if (n /= a .and. n /= 1 .and. a /= 1) return
        if (m /= b .and. m /= 1 .and. b /= 1) return
        deallocate(out)
        allocate(out(max(n,a),max(m,b),3))
        do j=1,size(out,2)
            do i=1,size(out,1)
                pi = min(i,a)
                pj = min(j,b)
                xyz = points(pi,pj,:)
                out(i,j,:) = self%values(min(i,n),min(j,m))%apply(xyz,inverse)
            end do
        end do
    end function grid_apply_grid_points

    pure function grid_inverse(self) result(out)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        type(rigid_transform_grid_type) :: out
        integer :: i,j
        if (.not. self%valid()) return
        allocate(out%values(size(self%values,1),size(self%values,2)))
        do j=1,size(out%values,2)
            do i=1,size(out%values,1)
                out%values(i,j) = self%values(i,j)%inverse()
                if (.not. out%values(i,j)%valid()) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
    end function grid_inverse

    pure function grid_compose_grid(self, other) result(out)
        class(rigid_transform_grid_type), intent(in) :: self !! transform applied second
        type(rigid_transform_grid_type), intent(in) :: other !! transform applied first
        type(rigid_transform_grid_type) :: out
        integer :: a,b,c,d,i,j
        if (.not. self%valid() .or. .not. other%valid()) return
        a = size(self%values,1)
        b = size(self%values,2)
        c = size(other%values,1)
        d = size(other%values,2)
        if (a /= c .and. a /= 1 .and. c /= 1) return
        if (b /= d .and. b /= 1 .and. d /= 1) return
        allocate(out%values(max(a,c),max(b,d)))
        do j=1,size(out%values,2)
            do i=1,size(out%values,1)
                out%values(i,j) = self%values(min(i,a),min(j,b))%compose(other%values(min(i,c),min(j,d)))
                if (.not. out%values(i,j)%valid()) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
    end function grid_compose_grid

    pure function grid_compose_single(self, other) result(out)
        class(rigid_transform_grid_type), intent(in) :: self !! transform applied second
        type(rigid_transform_type), intent(in) :: other !! transform applied first
        type(rigid_transform_grid_type) :: out
        integer :: i,j
        if (.not. self%valid() .or. .not. other%valid()) return
        allocate(out%values(size(self%values,1),size(self%values,2)))
        do j=1,size(out%values,2)
            do i=1,size(out%values,1)
                out%values(i,j) = self%values(i,j)%compose(other)
                if (.not. out%values(i,j)%valid()) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
    end function grid_compose_single

    pure function grid_compose_stack(self, other) result(out)
        class(rigid_transform_grid_type), intent(in) :: self !! transform applied second
        type(rigid_transform_stack_type), intent(in) :: other !! rank-one transform stack, right aligned
        type(rigid_transform_grid_type) :: out
        integer :: i,j,n,m,k
        if (.not. self%valid() .or. .not. other%valid()) return
        n = size(self%values,1)
        m = size(self%values,2)
        k = other%size()
        if (m /= k .and. m /= 1 .and. k /= 1) return
        allocate(out%values(n,max(m,k)))
        do j=1,size(out%values,2)
            do i=1,n
                out%values(i,j) = self%values(i,min(j,m))%compose(other%at(min(j,k)))
                if (.not. out%values(i,j)%valid()) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
    end function grid_compose_stack

    function grid_power(self, exponent, status) result(out)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        real(dp), intent(in) :: exponent !! common exponent applied to every transform
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_grid_type) :: out
        integer :: i,j,code
        call set_status(status, spatial_rigid_invalid_input)
        if (.not. self%valid() .or. .not. ieee_is_finite(exponent)) return
        allocate(out%values(size(self%values,1),size(self%values,2)))
        do j=1,size(out%values,2)
            do i=1,size(out%values,1)
                out%values(i,j) = self%values(i,j)%power(exponent,status=code)
                if (code /= spatial_rigid_success) then
                    deallocate(out%values)
                    return
                end if
            end do
        end do
        call set_status(status, spatial_rigid_success)
    end function grid_power

    subroutine prepare_weights(self, weights, w, ok)
        class(rigid_transform_grid_type), intent(in) :: self !! input grid
        real(dp), intent(in), optional :: weights(:,:) !! finite nonnegative broadcast weights
        real(dp), allocatable, intent(out) :: w(:,:) !! scaled weights, shape grid
        logical, intent(out) :: ok !! validity flag
        integer :: n,m,i,j,a,b
        ok = .false.
        if (.not. self%valid()) return
        n = size(self%values,1)
        m = size(self%values,2)
        allocate(w(n,m))
        w = 1.0_dp
        if (present(weights)) then
            a = size(weights,1)
            b = size(weights,2)
            if (a < 1 .or. b < 1) return
            if (.not. compatible(a,n) .or. .not. compatible(b,m)) return
            if (.not. all(ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            do j=1,m
                do i=1,n
                    w(i,j) = weights(min(i,a),min(j,b))
                end do
            end do
        end if
        if (maxval(w) <= 0.0_dp) return
        w = w/maxval(w)
        ok = .true.
    end subroutine prepare_weights

    function grid_mean_all(self, weights, status) result(tf)
        class(rigid_transform_grid_type), intent(in) :: self !! grid of transforms to average
        real(dp), intent(in), optional :: weights(:,:) !! broadcast nonnegative weights
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_type) :: tf
        type(rotation_type), allocatable :: rotations(:)
        real(dp), allocatable :: w(:,:), wflat(:)
        real(dp) :: t(3), total
        integer :: i,j,k,n,m,code
        logical :: ok
        call set_status(status, spatial_rigid_invalid_input)
        tf = invalid_scalar()
        call prepare_weights(self,weights,w,ok)
        if (.not. ok) return
        n = size(w,1)
        m = size(w,2)
        allocate(rotations(n*m),wflat(n*m))
        total = sum(w)
        t = 0.0_dp
        k = 0
        do j=1,m
            do i=1,n
                k = k+1
                rotations(k) = self%values(i,j)%rotation()
                wflat(k) = w(i,j)
                t = t + (w(i,j)/total)*self%values(i,j)%translation()
            end do
        end do
        if (.not. all(ieee_is_finite(t))) return
        tf = rigid_from_components(t,rotation_mean(rotations,weights=wflat,status=code),status=status)
        if (code /= 0) then
            tf = invalid_scalar()
            call set_status(status, spatial_rigid_invalid_input)
        end if
    end function grid_mean_all

    function grid_mean_axis(self, axis, weights, status) result(stack)
        class(rigid_transform_grid_type), intent(in) :: self !! grid to reduce
        integer, intent(in) :: axis !! 1 averages rows, 2 averages columns
        real(dp), intent(in), optional :: weights(:,:) !! broadcast nonnegative weights
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(rigid_transform_stack_type) :: stack
        type(rotation_type), allocatable :: rotations(:), rots_out(:)
        real(dp), allocatable :: w(:,:), current(:), t(:,:)
        real(dp) :: total
        integer :: n,m,i,j,nt,ns,code
        logical :: ok
        call set_status(status, spatial_rigid_invalid_input)
        if (axis /= 1 .and. axis /= 2) return
        call prepare_weights(self,weights,w,ok)
        if (.not. ok) return
        n = size(w,1)
        m = size(w,2)
        nt = size(w,axis)
        ns = size(w,3-axis)
        allocate(rotations(nt),current(nt),rots_out(ns),t(ns,3))
        do j=1,ns
            do i=1,nt
                if (axis == 1) then
                    rotations(i) = self%values(i,j)%rotation()
                    current(i) = w(i,j)
                else
                    rotations(i) = self%values(j,i)%rotation()
                    current(i) = w(j,i)
                end if
            end do
            total = sum(current)
            if (total <= 0.0_dp) return
            t(j,:) = 0.0_dp
            do i=1,nt
                if (axis == 1) then
                    t(j,:) = t(j,:) + (current(i)/total)*self%values(i,j)%translation()
                else
                    t(j,:) = t(j,:) + (current(i)/total)*self%values(j,i)%translation()
                end if
            end do
            if (.not. all(ieee_is_finite(t(j,:)))) return
            rots_out(j) = rotation_mean(rotations,weights=current,status=code)
            if (code /= 0) return
        end do
        stack = rigid_stack_from_components(t,rots_out,status=status)
    end function grid_mean_axis

end module scifort_spatial_rigid_transform_grid
