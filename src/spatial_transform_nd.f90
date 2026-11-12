! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Runtime-rank rotation and SE(3) arrays with trailing component axes.
! A flat Fortran-order representation is used for outputs, accompanied by
! the batch shape, because Fortran cannot dynamically choose a result rank.
module scifort_spatial_transform_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, &
        rotation_from_quat, rotation_from_rotvec, rotation_from_matrix
    use scifort_spatial_rigid_transform, only : rigid_transform_type, rigid_identity, &
        rigid_from_matrix, rigid_from_exp_coords, rigid_from_dual_quat, &
        rigid_from_components
    implicit none
    private

    integer, parameter, public :: spatial_nd_success = 0
    integer, parameter, public :: spatial_nd_invalid_input = 1

    type, public :: rotation_nd_type
        private
        type(rotation_type), allocatable :: values(:)
        integer, allocatable :: dims(:)
    contains
        procedure :: valid => rotation_nd_valid
        procedure :: shape => rotation_nd_shape
        procedure :: size => rotation_nd_size
        procedure :: at => rotation_nd_at
        procedure :: as_quat => rotation_nd_as_quat
        procedure :: as_rotvec => rotation_nd_as_rotvec
        procedure :: as_matrix => rotation_nd_as_matrix
        procedure :: inverse => rotation_nd_inverse
        procedure :: compose => rotation_nd_compose
        procedure :: apply => rotation_nd_apply
    end type rotation_nd_type

    type, public :: rigid_nd_type
        private
        type(rigid_transform_type), allocatable :: values(:)
        integer, allocatable :: dims(:)
    contains
        procedure :: valid => rigid_nd_valid
        procedure :: shape => rigid_nd_shape
        procedure :: size => rigid_nd_size
        procedure :: at => rigid_nd_at
        procedure :: as_matrix => rigid_nd_as_matrix
        procedure :: as_exp_coords => rigid_nd_as_exp_coords
        procedure :: as_dual_quat => rigid_nd_as_dual_quat
        procedure :: translation => rigid_nd_translation
        procedure :: inverse => rigid_nd_inverse
        procedure :: compose => rigid_nd_compose
        procedure :: apply => rigid_nd_apply
        procedure :: power => rigid_nd_power
    end type rigid_nd_type

    public :: rotation_nd_from_quat, rotation_nd_from_rotvec, rotation_nd_from_matrix
    public :: rotation_nd_from_rotations, rotation_nd_identity
    public :: rigid_nd_from_matrix, rigid_nd_from_exp_coords, rigid_nd_from_dual_quat
    public :: rigid_nd_from_components, rigid_nd_from_transforms, rigid_nd_identity

contains

    pure subroutine set_status(status, value)
        integer, intent(out), optional :: status !! zero=success, nonzero=invalid input
        integer, intent(in) :: value !! requested status
        if (present(status)) status = value
    end subroutine set_status

    ! Last-axis flattening accepts any Fortran array rank 1 through 15.
    ! Input shape (*batch_dims, ncomp) becomes flat(product(batch_dims),ncomp).
    subroutine flatten_components(array, ncomp, flat, dims, ok)
        real(dp), intent(in) :: array(..) !! components in the final dimension
        integer, intent(in) :: ncomp !! number of trailing components
        real(dp), allocatable, intent(out) :: flat(:,:) !! Fortran-order component rows
        integer, allocatable, intent(out) :: dims(:) !! arbitrary batch dimensions
        logical, intent(out) :: ok !! validity flag
        integer, allocatable :: extents(:)
        ok = .false.
        select rank(array)
        rank(1)
            extents = shape(array)
            if (extents(1) /= ncomp) return
            dims = extents(:0)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(2)
            extents = shape(array)
            if (extents(2) /= ncomp) return
            dims = extents(:1)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(3)
            extents = shape(array)
            if (extents(3) /= ncomp) return
            dims = extents(:2)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(4)
            extents = shape(array)
            if (extents(4) /= ncomp) return
            dims = extents(:3)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(5)
            extents = shape(array)
            if (extents(5) /= ncomp) return
            dims = extents(:4)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(6)
            extents = shape(array)
            if (extents(6) /= ncomp) return
            dims = extents(:5)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(7)
            extents = shape(array)
            if (extents(7) /= ncomp) return
            dims = extents(:6)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(8)
            extents = shape(array)
            if (extents(8) /= ncomp) return
            dims = extents(:7)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(9)
            extents = shape(array)
            if (extents(9) /= ncomp) return
            dims = extents(:8)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(10)
            extents = shape(array)
            if (extents(10) /= ncomp) return
            dims = extents(:9)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(11)
            extents = shape(array)
            if (extents(11) /= ncomp) return
            dims = extents(:10)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(12)
            extents = shape(array)
            if (extents(12) /= ncomp) return
            dims = extents(:11)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(13)
            extents = shape(array)
            if (extents(13) /= ncomp) return
            dims = extents(:12)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(14)
            extents = shape(array)
            if (extents(14) /= ncomp) return
            dims = extents(:13)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank(15)
            extents = shape(array)
            if (extents(15) /= ncomp) return
            dims = extents(:14)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),ncomp])
        rank default
            return
        end select
        ok = .true.
    end subroutine flatten_components

    ! Last-two-axis flattening; matrices have batch axes followed by (n,n).
    subroutine flatten_matrices(array, n, flat, dims, ok)
        real(dp), intent(in) :: array(..) !! (*batch_dims,n,n) matrices
        integer, intent(in) :: n !! spatial matrix order (3 or 4)
        real(dp), allocatable, intent(out) :: flat(:,:) !! rows with n*n components
        integer, allocatable, intent(out) :: dims(:) !! batch dimensions
        logical, intent(out) :: ok !! validity flag
        integer, allocatable :: extents(:)
        ok = .false.
        select rank(array)
        rank(2)
            extents = shape(array)
            if (any(extents(1:2) /= n)) return
            dims = extents(:0)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(3)
            extents = shape(array)
            if (any(extents(2:3) /= n)) return
            dims = extents(:1)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(4)
            extents = shape(array)
            if (any(extents(3:4) /= n)) return
            dims = extents(:2)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(5)
            extents = shape(array)
            if (any(extents(4:5) /= n)) return
            dims = extents(:3)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(6)
            extents = shape(array)
            if (any(extents(5:6) /= n)) return
            dims = extents(:4)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(7)
            extents = shape(array)
            if (any(extents(6:7) /= n)) return
            dims = extents(:5)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(8)
            extents = shape(array)
            if (any(extents(7:8) /= n)) return
            dims = extents(:6)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(9)
            extents = shape(array)
            if (any(extents(8:9) /= n)) return
            dims = extents(:7)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(10)
            extents = shape(array)
            if (any(extents(9:10) /= n)) return
            dims = extents(:8)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(11)
            extents = shape(array)
            if (any(extents(10:11) /= n)) return
            dims = extents(:9)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(12)
            extents = shape(array)
            if (any(extents(11:12) /= n)) return
            dims = extents(:10)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(13)
            extents = shape(array)
            if (any(extents(12:13) /= n)) return
            dims = extents(:11)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(14)
            extents = shape(array)
            if (any(extents(13:14) /= n)) return
            dims = extents(:12)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank(15)
            extents = shape(array)
            if (any(extents(14:15) /= n)) return
            dims = extents(:13)
            if (any(dims < 1)) return
            flat = reshape(array,[product(dims),n*n])
        rank default
            return
        end select
        ok = .true.
    end subroutine flatten_matrices

    ! NumPy right-aligned broadcasting, returning a zero-length shape for
    ! a scalar batch. Invalid or empty input dimensions are rejected.
    function broadcast_shape(left, right, ok) result(dims)
        integer, intent(in) :: left(:) !! left batch shape
        integer, intent(in) :: right(:) !! right batch shape
        logical, intent(out) :: ok !! shape compatibility
        integer, allocatable :: dims(:) !! right-aligned broadcast dimensions
        integer :: i, j, a, b, nr
        ok = .false.
        if (any(left < 1) .or. any(right < 1)) return
        nr = max(size(left),size(right))
        allocate(dims(nr))
        do i = 1, nr
            a = 1
            b = 1
            j = i - (nr-size(left))
            if (j > 0) a = left(j)
            j = i - (nr-size(right))
            if (j > 0) b = right(j)
            if (a /= b .and. a /= 1 .and. b /= 1) then
                dims = 0
                return
            end if
            dims(i) = max(a,b)
        end do
        if (checked_count(dims)<1) then
            dims = 0
            return
        end if
        ok = .true.
    end function broadcast_shape

    pure integer function checked_count(dims) result(n)
        integer, intent(in) :: dims(:) !! positive batch lengths; empty shape means one scalar
        integer :: i
        n = 1
        do i=1,size(dims)
            if (dims(i)<1) then
                n = -1
                return
            end if
            if (n>huge(n)/dims(i)) then
                n = -1
                return
            end if
            n = n*dims(i)
        end do
    end function checked_count

    ! Convert one-based Fortran-order output index to an operand index.
    pure integer function source_index(output_index, output_dims, source_dims) result(index)
        integer, intent(in) :: output_index !! flat, one-based output index
        integer, intent(in) :: output_dims(:) !! broadcast result shape
        integer, intent(in) :: source_dims(:) !! original operand shape
        integer :: remainder, i, j, coord, stride
        remainder = output_index-1
        stride = 1
        index = 1
        do i=1,size(output_dims)
            coord = mod(remainder,output_dims(i))
            remainder = remainder/output_dims(i)
            j = i - (size(output_dims)-size(source_dims))
            if (j > 0) then
                if (source_dims(j) /= 1) index = index + coord*stride
                stride = stride*source_dims(j)
            end if
        end do
    end function source_index

    pure logical function valid_dims(dims, count) result(ok)
        integer, intent(in) :: dims(:) !! proposed batch dimensions
        integer, intent(in) :: count !! expected number of stored scalars
        ok = .false.
        if (size(dims) > 14 .or. count < 1) return
        ok = checked_count(dims)==count
    end function valid_dims

    function rotation_nd_from_rotations(rotations, batch_shape, status) result(out)
        type(rotation_type), intent(in) :: rotations(:) !! flat, Fortran-order rotations
        integer, intent(in) :: batch_shape(:) !! shape, rank 0..14
        integer, intent(out), optional :: status !! zero on success
        type(rotation_nd_type) :: out
        integer :: i
        call set_status(status,spatial_nd_invalid_input)
        if (.not. valid_dims(batch_shape,size(rotations))) return
        do i=1,size(rotations)
            if (.not. rotations(i)%valid()) return
        end do
        out%dims = batch_shape
        out%values = rotations
        call set_status(status,spatial_nd_success)
    end function rotation_nd_from_rotations

    function rotation_nd_identity(batch_shape, status) result(out)
        integer, intent(in) :: batch_shape(:) !! strictly positive axes, rank 0..14
        integer, intent(out), optional :: status !! zero on success
        type(rotation_nd_type) :: out
        integer :: i, n
        call set_status(status,spatial_nd_invalid_input)
        if (size(batch_shape)>14 .or. any(batch_shape<1)) return
        n = checked_count(batch_shape)
        if (n<1) return
        allocate(out%values(n))
        out%dims = batch_shape
        do i=1,n
            out%values(i) = rotation_identity()
        end do
        call set_status(status,spatial_nd_success)
    end function rotation_nd_identity

    function rotation_nd_from_quat(quats, scalar_first, status) result(out)
        real(dp), intent(in) :: quats(..) !! (*batch,4), rank 1..15
        logical, intent(in), optional :: scalar_first !! scalar component first when true
        integer, intent(out), optional :: status !! zero on success
        type(rotation_nd_type) :: out
        real(dp), allocatable :: flat(:,:)
        integer, allocatable :: dims(:)
        integer :: i, code
        logical :: ok
        call set_status(status,spatial_nd_invalid_input)
        call flatten_components(quats,4,flat,dims,ok)
        if (.not. ok) return
        allocate(out%values(size(flat,1)))
        do i=1,size(flat,1)
            out%values(i) = rotation_from_quat(flat(i,:),scalar_first,code)
            if (code/=0) then
                deallocate(out%values)
                return
            end if
        end do
        out%dims = dims
        call set_status(status,spatial_nd_success)
    end function rotation_nd_from_quat

    function rotation_nd_from_rotvec(rotvecs, degrees, status) result(out)
        real(dp), intent(in) :: rotvecs(..) !! (*batch,3), rank 1..15
        logical, intent(in), optional :: degrees !! input is in degrees
        integer, intent(out), optional :: status !! zero on success
        type(rotation_nd_type) :: out
        real(dp), allocatable :: flat(:,:)
        integer, allocatable :: dims(:)
        integer :: i, code
        logical :: ok
        call set_status(status,spatial_nd_invalid_input)
        call flatten_components(rotvecs,3,flat,dims,ok)
        if (.not. ok) return
        allocate(out%values(size(flat,1)))
        do i=1,size(flat,1)
            out%values(i) = rotation_from_rotvec(flat(i,:),degrees,code)
            if (code/=0) then
                deallocate(out%values)
                return
            end if
        end do
        out%dims = dims
        call set_status(status,spatial_nd_success)
    end function rotation_nd_from_rotvec

    function rotation_nd_from_matrix(matrices, status) result(out)
        real(dp), intent(in) :: matrices(..) !! (*batch,3,3), rank 2..15
        integer, intent(out), optional :: status !! zero on success
        type(rotation_nd_type) :: out
        real(dp), allocatable :: flat(:,:)
        integer, allocatable :: dims(:)
        integer :: i, code
        logical :: ok
        call set_status(status,spatial_nd_invalid_input)
        call flatten_matrices(matrices,3,flat,dims,ok)
        if (.not. ok) return
        allocate(out%values(size(flat,1)))
        do i=1,size(flat,1)
            out%values(i) = rotation_from_matrix(reshape(flat(i,:),[3,3]),code)
            if (code/=0) then
                deallocate(out%values)
                return
            end if
        end do
        out%dims = dims
        call set_status(status,spatial_nd_success)
    end function rotation_nd_from_matrix

    pure logical function rotation_nd_valid(self) result(ok)
        class(rotation_nd_type), intent(in) :: self !! rotation array
        integer :: i
        ok = .false.
        if (.not. allocated(self%dims) .or. .not. allocated(self%values)) return
        if (.not. valid_dims(self%dims,size(self%values))) return
        do i=1,size(self%values)
            if (.not. self%values(i)%valid()) return
        end do
        ok = .true.
    end function rotation_nd_valid

    pure function rotation_nd_shape(self) result(dims)
        class(rotation_nd_type), intent(in) :: self !! rotation array
        integer, allocatable :: dims(:) !! batch shape, empty for invalid arrays
        if (allocated(self%dims)) then
            dims = self%dims
        else
            allocate(dims(0))
        end if
    end function rotation_nd_shape

    pure integer function rotation_nd_size(self) result(n)
        class(rotation_nd_type), intent(in) :: self !! rotation array
        n = 0
        if (allocated(self%values)) n = size(self%values)
    end function rotation_nd_size

    pure integer function index_from_coords(coords, dims) result(i)
        integer, intent(in) :: coords(:) !! one-based multidimensional coordinates
        integer, intent(in) :: dims(:) !! batch dimensions
        integer :: j, stride
        i = 1
        stride = 1
        do j=1,size(dims)
            i = i + (coords(j)-1)*stride
            stride = stride*dims(j)
        end do
    end function index_from_coords

    function rotation_nd_at(self, indices) result(value)
        class(rotation_nd_type), intent(in) :: self !! rotation array
        integer, intent(in) :: indices(:) !! one-based batch indices
        type(rotation_type) :: value
        value = rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        if (.not. self%valid()) return
        if (size(indices)/=size(self%dims)) return
        if (any(indices<1) .or. any(indices>self%dims)) return
        value = self%values(index_from_coords(indices,self%dims))
    end function rotation_nd_at

    function rotation_nd_as_quat(self, scalar_first) result(quats)
        class(rotation_nd_type), intent(in) :: self !! rotation array
        logical, intent(in), optional :: scalar_first !! scalar-first convention
        real(dp), allocatable :: quats(:,:) !! (product(batch_shape),4)
        integer :: i
        allocate(quats(self%size(),4))
        if (.not. self%valid()) return
        do i=1,size(self%values)
            quats(i,:) = self%values(i)%as_quat(scalar_first=scalar_first)
        end do
    end function rotation_nd_as_quat

    function rotation_nd_as_rotvec(self, degrees) result(rotvecs)
        class(rotation_nd_type), intent(in) :: self !! rotation array
        logical, intent(in), optional :: degrees !! output angles in degrees
        real(dp), allocatable :: rotvecs(:,:) !! (product(batch_shape),3)
        integer :: i
        allocate(rotvecs(self%size(),3))
        if (.not. self%valid()) return
        do i=1,size(self%values)
            rotvecs(i,:) = self%values(i)%as_rotvec(degrees)
        end do
    end function rotation_nd_as_rotvec

    function rotation_nd_as_matrix(self) result(matrices)
        class(rotation_nd_type), intent(in) :: self !! rotation array
        real(dp), allocatable :: matrices(:,:,:) !! (product(batch_shape),3,3)
        integer :: i
        allocate(matrices(self%size(),3,3))
        if (.not. self%valid()) return
        do i=1,size(self%values)
            matrices(i,:,:) = self%values(i)%as_matrix()
        end do
    end function rotation_nd_as_matrix

    function rotation_nd_inverse(self) result(out)
        class(rotation_nd_type), intent(in) :: self !! rotation array
        type(rotation_nd_type) :: out
        integer :: i
        if (.not. self%valid()) return
        out%dims = self%dims
        allocate(out%values(size(self%values)))
        do i=1,size(self%values)
            out%values(i) = self%values(i)%inverse()
        end do
    end function rotation_nd_inverse

    function rotation_nd_compose(self, other, status) result(out)
        class(rotation_nd_type), intent(in) :: self !! outer rotations
        type(rotation_nd_type), intent(in) :: other !! inner rotations
        integer, intent(out), optional :: status !! zero on success
        type(rotation_nd_type) :: out
        integer, allocatable :: dims(:)
        integer :: i, j, k
        logical :: ok
        call set_status(status,spatial_nd_invalid_input)
        if (.not. self%valid() .or. .not. other%valid()) return
        dims = broadcast_shape(self%dims,other%dims,ok)
        if (.not. ok) return
        out%dims = dims
        allocate(out%values(checked_count(dims)))
        do i=1,size(out%values)
            j = source_index(i,dims,self%dims)
            k = source_index(i,dims,other%dims)
            out%values(i) = self%values(j)%compose(other%values(k))
            if (.not. out%values(i)%valid()) then
                deallocate(out%values,out%dims)
                return
            end if
        end do
        call set_status(status,spatial_nd_success)
    end function rotation_nd_compose

    function rotation_nd_apply(self, points, result_shape, status, inverse) result(mapped)
        class(rotation_nd_type), intent(in) :: self !! rotations (batch dimensions)
        real(dp), intent(in) :: points(..) !! points (*batch,3), rank 1..15
        integer, allocatable, intent(out), optional :: result_shape(:) !! broadcast shape
        integer, intent(out), optional :: status !! zero on success
        logical, intent(in), optional :: inverse !! apply inverse rotations
        real(dp), allocatable :: mapped(:,:) !! flattened output (n,3)
        real(dp), allocatable :: flat(:,:)
        integer, allocatable :: pdims(:), dims(:)
        integer :: i, j, k
        logical :: ok
        allocate(mapped(0,3))
        call set_status(status,spatial_nd_invalid_input)
        if (.not. self%valid()) return
        call flatten_components(points,3,flat,pdims,ok)
        if (.not. ok) return
        if (.not. all(ieee_is_finite(flat))) return
        dims = broadcast_shape(self%dims,pdims,ok)
        if (.not. ok) return
        deallocate(mapped)
        allocate(mapped(checked_count(dims),3))
        do i=1,size(mapped,1)
            j = source_index(i,dims,self%dims)
            k = source_index(i,dims,pdims)
            mapped(i,:) = self%values(j)%apply(flat(k,:),inverse)
        end do
        if (present(result_shape)) result_shape = dims
        call set_status(status,spatial_nd_success)
    end function rotation_nd_apply

    function rigid_nd_from_transforms(transforms, batch_shape, status) result(out)
        type(rigid_transform_type), intent(in) :: transforms(:) !! Fortran-order flat SE(3) transforms
        integer, intent(in) :: batch_shape(:) !! batch axes, rank 0..14
        integer, intent(out), optional :: status !! zero on success
        type(rigid_nd_type) :: out
        integer :: i
        call set_status(status,spatial_nd_invalid_input)
        if (.not. valid_dims(batch_shape,size(transforms))) return
        do i=1,size(transforms)
            if (.not. transforms(i)%valid()) return
        end do
        out%dims = batch_shape
        out%values = transforms
        call set_status(status,spatial_nd_success)
    end function rigid_nd_from_transforms

    function rigid_nd_identity(batch_shape, status) result(out)
        integer, intent(in) :: batch_shape(:) !! strictly positive dimensions, rank 0..14
        integer, intent(out), optional :: status !! zero on success
        type(rigid_nd_type) :: out
        integer :: n,i
        call set_status(status,spatial_nd_invalid_input)
        if (size(batch_shape)>14 .or. any(batch_shape<1)) return
        n = checked_count(batch_shape)
        if (n<1) return
        allocate(out%values(n))
        out%dims = batch_shape
        do i=1,n
            out%values(i) = rigid_identity()
        end do
        call set_status(status,spatial_nd_success)
    end function rigid_nd_identity

    function rigid_nd_from_matrix(matrices, status) result(out)
        real(dp), intent(in) :: matrices(..) !! (*batch,4,4), rank 2..15
        integer, intent(out), optional :: status !! zero on success
        type(rigid_nd_type) :: out
        real(dp), allocatable :: flat(:,:)
        integer, allocatable :: dims(:)
        integer :: i, code
        logical :: ok
        call set_status(status,spatial_nd_invalid_input)
        call flatten_matrices(matrices,4,flat,dims,ok)
        if (.not. ok) return
        allocate(out%values(size(flat,1)))
        do i=1,size(flat,1)
            out%values(i) = rigid_from_matrix(reshape(flat(i,:),[4,4]),code)
            if (code/=0) then
                deallocate(out%values)
                return
            end if
        end do
        out%dims = dims
        call set_status(status,spatial_nd_success)
    end function rigid_nd_from_matrix

    function rigid_nd_from_exp_coords(coords, status) result(out)
        real(dp), intent(in) :: coords(..) !! (*batch,6) rotation/translation twists
        integer, intent(out), optional :: status !! zero on success
        type(rigid_nd_type) :: out
        real(dp), allocatable :: flat(:,:)
        integer, allocatable :: dims(:)
        integer :: i, code
        logical :: ok
        call set_status(status,spatial_nd_invalid_input)
        call flatten_components(coords,6,flat,dims,ok)
        if (.not. ok) return
        allocate(out%values(size(flat,1)))
        do i=1,size(flat,1)
            out%values(i) = rigid_from_exp_coords(flat(i,:),code)
            if (code/=0) then
                deallocate(out%values)
                return
            end if
        end do
        out%dims = dims
        call set_status(status,spatial_nd_success)
    end function rigid_nd_from_exp_coords

    function rigid_nd_from_dual_quat(dual_quats, scalar_first, status) result(out)
        real(dp), intent(in) :: dual_quats(..) !! (*batch,8) dual quaternions
        logical, intent(in), optional :: scalar_first !! w-first input convention
        integer, intent(out), optional :: status !! zero on success
        type(rigid_nd_type) :: out
        real(dp), allocatable :: flat(:,:)
        integer, allocatable :: dims(:)
        integer :: i, code
        logical :: ok
        call set_status(status,spatial_nd_invalid_input)
        call flatten_components(dual_quats,8,flat,dims,ok)
        if (.not. ok) return
        allocate(out%values(size(flat,1)))
        do i=1,size(flat,1)
            out%values(i) = rigid_from_dual_quat(flat(i,:),scalar_first,code)
            if (code/=0) then
                deallocate(out%values)
                return
            end if
        end do
        out%dims = dims
        call set_status(status,spatial_nd_success)
    end function rigid_nd_from_dual_quat

    function rigid_nd_from_components(translations, rotations, status) result(out)
        real(dp), intent(in) :: translations(..) !! (*batch,3) translation vectors
        type(rotation_nd_type), intent(in) :: rotations !! rotations with compatible batch shape
        integer, intent(out), optional :: status !! zero on success
        type(rigid_nd_type) :: out
        real(dp), allocatable :: flat(:,:)
        integer, allocatable :: tdims(:), dims(:)
        integer :: i, j, k, code
        logical :: ok
        call set_status(status,spatial_nd_invalid_input)
        if (.not. rotations%valid()) return
        call flatten_components(translations,3,flat,tdims,ok)
        if (.not. ok) return
        if (.not. all(ieee_is_finite(flat))) return
        dims = broadcast_shape(tdims,rotations%dims,ok)
        if (.not. ok) return
        out%dims = dims
        allocate(out%values(checked_count(dims)))
        do i=1,size(out%values)
            j = source_index(i,dims,tdims)
            k = source_index(i,dims,rotations%dims)
            out%values(i) = rigid_from_components(flat(j,:),rotations%values(k),code)
            if (code/=0) then
                deallocate(out%values,out%dims)
                return
            end if
        end do
        call set_status(status,spatial_nd_success)
    end function rigid_nd_from_components

    pure logical function rigid_nd_valid(self) result(ok)
        class(rigid_nd_type), intent(in) :: self !! SE(3) array
        integer :: i
        ok = .false.
        if (.not. allocated(self%dims) .or. .not. allocated(self%values)) return
        if (.not. valid_dims(self%dims,size(self%values))) return
        do i=1,size(self%values)
            if (.not. self%values(i)%valid()) return
        end do
        ok = .true.
    end function rigid_nd_valid

    pure function rigid_nd_shape(self) result(dims)
        class(rigid_nd_type), intent(in) :: self !! SE(3) array
        integer, allocatable :: dims(:) !! batch shape, empty for invalid inputs
        if (allocated(self%dims)) then
            dims = self%dims
        else
            allocate(dims(0))
        end if
    end function rigid_nd_shape

    pure integer function rigid_nd_size(self) result(n)
        class(rigid_nd_type), intent(in) :: self !! SE(3) array
        n = 0
        if (allocated(self%values)) n = size(self%values)
    end function rigid_nd_size

    function rigid_nd_at(self, indices) result(value)
        class(rigid_nd_type), intent(in) :: self !! SE(3) array
        integer, intent(in) :: indices(:) !! one-based batch indices
        type(rigid_transform_type) :: value
        value = rigid_from_exp_coords([0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        if (.not. self%valid()) return
        if (size(indices)/=size(self%dims)) return
        if (any(indices<1) .or. any(indices>self%dims)) return
        value = self%values(index_from_coords(indices,self%dims))
    end function rigid_nd_at

    function rigid_nd_as_matrix(self) result(matrices)
        class(rigid_nd_type), intent(in) :: self !! SE(3) array
        real(dp), allocatable :: matrices(:,:,:) !! flattened (n,4,4)
        integer :: i
        allocate(matrices(self%size(),4,4))
        if (.not. self%valid()) return
        do i=1,size(self%values)
            matrices(i,:,:) = self%values(i)%as_matrix()
        end do
    end function rigid_nd_as_matrix

    function rigid_nd_as_exp_coords(self) result(coords)
        class(rigid_nd_type), intent(in) :: self !! SE(3) array
        real(dp), allocatable :: coords(:,:) !! flattened (n,6)
        integer :: i
        allocate(coords(self%size(),6))
        if (.not. self%valid()) return
        do i=1,size(self%values)
            coords(i,:) = self%values(i)%as_exp_coords()
        end do
    end function rigid_nd_as_exp_coords

    function rigid_nd_as_dual_quat(self, scalar_first) result(quats)
        class(rigid_nd_type), intent(in) :: self !! SE(3) array
        logical, intent(in), optional :: scalar_first !! w-first convention
        real(dp), allocatable :: quats(:,:) !! flattened (n,8)
        integer :: i
        allocate(quats(self%size(),8))
        if (.not. self%valid()) return
        do i=1,size(self%values)
            quats(i,:) = self%values(i)%as_dual_quat(scalar_first)
        end do
    end function rigid_nd_as_dual_quat

    function rigid_nd_translation(self) result(vecs)
        class(rigid_nd_type), intent(in) :: self !! SE(3) array
        real(dp), allocatable :: vecs(:,:) !! flattened (n,3)
        integer :: i
        allocate(vecs(self%size(),3))
        if (.not. self%valid()) return
        do i=1,size(self%values)
            vecs(i,:) = self%values(i)%translation()
        end do
    end function rigid_nd_translation

    function rigid_nd_inverse(self) result(out)
        class(rigid_nd_type), intent(in) :: self !! SE(3) array
        type(rigid_nd_type) :: out
        integer :: i
        if (.not. self%valid()) return
        out%dims = self%dims
        allocate(out%values(size(self%values)))
        do i=1,size(self%values)
            out%values(i) = self%values(i)%inverse()
        end do
    end function rigid_nd_inverse

    function rigid_nd_compose(self, other, status) result(out)
        class(rigid_nd_type), intent(in) :: self !! outer transformations
        type(rigid_nd_type), intent(in) :: other !! inner transformations
        integer, intent(out), optional :: status !! zero on success
        type(rigid_nd_type) :: out
        integer, allocatable :: dims(:)
        integer :: i, j, k
        logical :: ok
        call set_status(status,spatial_nd_invalid_input)
        if (.not. self%valid() .or. .not. other%valid()) return
        dims = broadcast_shape(self%dims,other%dims,ok)
        if (.not. ok) return
        out%dims = dims
        allocate(out%values(checked_count(dims)))
        do i=1,size(out%values)
            j = source_index(i,dims,self%dims)
            k = source_index(i,dims,other%dims)
            out%values(i) = self%values(j)%compose(other%values(k))
            if (.not. out%values(i)%valid()) then
                deallocate(out%values,out%dims)
                return
            end if
        end do
        call set_status(status,spatial_nd_success)
    end function rigid_nd_compose

    function rigid_nd_apply(self, points, result_shape, status, inverse) result(mapped)
        class(rigid_nd_type), intent(in) :: self !! SE(3) transformations
        real(dp), intent(in) :: points(..) !! points (*batch,3), rank 1..15
        integer, allocatable, intent(out), optional :: result_shape(:) !! broadcast shape
        integer, intent(out), optional :: status !! zero on success
        logical, intent(in), optional :: inverse !! apply inverse transformations
        real(dp), allocatable :: mapped(:,:) !! flattened (n,3)
        real(dp), allocatable :: flat(:,:)
        integer, allocatable :: pdims(:), dims(:)
        integer :: i, j, k
        logical :: ok
        allocate(mapped(0,3))
        call set_status(status,spatial_nd_invalid_input)
        if (.not. self%valid()) return
        call flatten_components(points,3,flat,pdims,ok)
        if (.not. ok) return
        if (.not. all(ieee_is_finite(flat))) return
        dims = broadcast_shape(self%dims,pdims,ok)
        if (.not. ok) return
        deallocate(mapped)
        allocate(mapped(checked_count(dims),3))
        do i=1,size(mapped,1)
            j = source_index(i,dims,self%dims)
            k = source_index(i,dims,pdims)
            mapped(i,:) = self%values(j)%apply(flat(k,:),inverse)
        end do
        if (present(result_shape)) result_shape = dims
        call set_status(status,spatial_nd_success)
    end function rigid_nd_apply

    function rigid_nd_power(self, exponent, status) result(out)
        class(rigid_nd_type), intent(in) :: self !! SE(3) transformations
        real(dp), intent(in) :: exponent !! common real exponent
        integer, intent(out), optional :: status !! zero on success
        type(rigid_nd_type) :: out
        integer :: i, code
        call set_status(status,spatial_nd_invalid_input)
        if (.not. self%valid() .or. .not. ieee_is_finite(exponent)) return
        out%dims = self%dims
        allocate(out%values(size(self%values)))
        do i=1,size(out%values)
            out%values(i) = self%values(i)%power(exponent,code)
            if (code/=0) then
                deallocate(out%values,out%dims)
                return
            end if
        end do
        call set_status(status,spatial_nd_success)
    end function rigid_nd_power

end module scifort_spatial_transform_nd
