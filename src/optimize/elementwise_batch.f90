! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Flattened NumPy-style trailing-axis broadcasting for scalar elementwise solvers.
! Used by all batch ranks, including the Fortran 2018 maximum rank of 15.
module scifort_optimize_elementwise_batch
    use scifort_kinds, only : dp
    implicit none
    private
    public :: elementwise_broadcast
    public :: elementwise_broadcast_optional
    public :: elementwise_flatten
    public :: elementwise_broadcast_arguments
    public :: elementwise_argument
    public :: elementwise_argument_from
    public :: elementwise_broadcast_argument_list
    public :: elementwise_typed_argument
    public :: elementwise_typed_argument_from
    public :: elementwise_broadcast_typed_argument_list
    public :: elementwise_pack_typed_arguments
    public :: elementwise_gather_typed_indices

    ! A real-valued objective argument with its own broadcast shape.  The
    ! components are private so callers construct valid snapshots using
    ! elementwise_argument_from(scalar_or_array).  One argument can be scalar
    ! while the next is a vector or a high-rank array.
    type :: elementwise_argument
        private
        real(dp), allocatable :: values(:)
        integer, allocatable :: dims(:)
    end type elementwise_argument

    ! Each scalar payload is separately owned. GNU Fortran 14.2 does not
    ! reliably reshape arrays of CLASS(*) (incorrect values or compiler ICE).
    ! Scalar polymorphic allocation followed by normal derived-type array
    ! assignment retains exact integer/logical/character/user-defined types.
    type :: elementwise_typed_value
        class(*), allocatable :: value !! one polymorphic value, with exact dynamic type
    end type elementwise_typed_value

    type :: elementwise_typed_argument
        type(elementwise_typed_value), allocatable :: values(:) !! one value per active lane
        integer, allocatable, private :: dims(:) !! original broadcast dimensions
    end type elementwise_typed_argument
contains

    function elementwise_box_value(value) result(item)
        class(*), intent(in) :: value !! typed scalar input
        type(elementwise_typed_value) :: item !! owned polymorphic copy
        allocate(item%value, source=value)
    end function elementwise_box_value

    function elementwise_typed_argument_from(arg) result(a)
        ! Flatten any rank in Fortran element order without using polymorphic
        ! RESHAPE (miscompiled by GNU Fortran 14.2 for CLASS(*) arrays).
        class(*), intent(in) :: arg(..) !! arbitrary type and rank
        type(elementwise_typed_argument) :: a !! flattened typed parameter
        integer :: pos
        integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

        a%dims = shape(arg)
        allocate(a%values(size(arg)))
        pos = 0
        select rank(arg)
        rank(0)
            a%values(1) = elementwise_box_value(arg)
        rank(1)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1))
            end do
        rank(2)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2))
            end do
            end do
        rank(3)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3))
            end do
            end do
            end do
        rank(4)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4))
            end do
            end do
            end do
            end do
        rank(5)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5))
            end do
            end do
            end do
            end do
            end do
        rank(6)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6))
            end do
            end do
            end do
            end do
            end do
            end do
        rank(7)
            do i7 = 1, size(arg, 7)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6, i7))
            end do
            end do
            end do
            end do
            end do
            end do
            end do
        rank(8)
            do i8 = 1, size(arg, 8)
            do i7 = 1, size(arg, 7)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6, i7, i8))
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
        rank(9)
            do i9 = 1, size(arg, 9)
            do i8 = 1, size(arg, 8)
            do i7 = 1, size(arg, 7)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6, i7, i8, i9))
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
        rank(10)
            do i10 = 1, size(arg, 10)
            do i9 = 1, size(arg, 9)
            do i8 = 1, size(arg, 8)
            do i7 = 1, size(arg, 7)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10))
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
        rank(11)
            do i11 = 1, size(arg, 11)
            do i10 = 1, size(arg, 10)
            do i9 = 1, size(arg, 9)
            do i8 = 1, size(arg, 8)
            do i7 = 1, size(arg, 7)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11))
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
        rank(12)
            do i12 = 1, size(arg, 12)
            do i11 = 1, size(arg, 11)
            do i10 = 1, size(arg, 10)
            do i9 = 1, size(arg, 9)
            do i8 = 1, size(arg, 8)
            do i7 = 1, size(arg, 7)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12))
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
        rank(13)
            do i13 = 1, size(arg, 13)
            do i12 = 1, size(arg, 12)
            do i11 = 1, size(arg, 11)
            do i10 = 1, size(arg, 10)
            do i9 = 1, size(arg, 9)
            do i8 = 1, size(arg, 8)
            do i7 = 1, size(arg, 7)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13))
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
        rank(14)
            do i14 = 1, size(arg, 14)
            do i13 = 1, size(arg, 13)
            do i12 = 1, size(arg, 12)
            do i11 = 1, size(arg, 11)
            do i10 = 1, size(arg, 10)
            do i9 = 1, size(arg, 9)
            do i8 = 1, size(arg, 8)
            do i7 = 1, size(arg, 7)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14))
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
        rank(15)
            do i15 = 1, size(arg, 15)
            do i14 = 1, size(arg, 14)
            do i13 = 1, size(arg, 13)
            do i12 = 1, size(arg, 12)
            do i11 = 1, size(arg, 11)
            do i10 = 1, size(arg, 10)
            do i9 = 1, size(arg, 9)
            do i8 = 1, size(arg, 8)
            do i7 = 1, size(arg, 7)
            do i6 = 1, size(arg, 6)
            do i5 = 1, size(arg, 5)
            do i4 = 1, size(arg, 4)
            do i3 = 1, size(arg, 3)
            do i2 = 1, size(arg, 2)
            do i1 = 1, size(arg, 1)
                pos = pos + 1
                a%values(pos) = elementwise_box_value(arg(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15))
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
            end do
        rank default
            error stop 'Unsupported argument rank'
        end select
    end function elementwise_typed_argument_from

    subroutine elementwise_broadcast_typed_argument_list(args, batch_shape, indices, valid)
        ! Map each source argument's independent broadcast dimensions onto the
        ! problem grid without converting any scalar to real(dp).
        type(elementwise_typed_argument), intent(in) :: args(:) !! typed parameters
        integer, intent(in) :: batch_shape(:) !! problem-grid dimensions
        integer, allocatable, intent(out) :: indices(:,:) !! (nargs, nproblems)
        logical, intent(inout) :: valid !! false for malformed input
        integer :: i, j

        if (.not. valid) return
        if (size(args) == 0) then
            valid = .false.
            return
        end if
        do i = 1, size(args)
            if (.not. allocated(args(i)%dims) .or. &
                    .not. allocated(args(i)%values)) then
                valid = .false.
                return
            end if
            if (size(args(i)%values) /= product(args(i)%dims)) then
                valid = .false.
                return
            end if
            do j = 1, size(args(i)%values)
                if (.not. allocated(args(i)%values(j)%value)) then
                    valid = .false.
                    return
                end if
            end do
        end do
        allocate(indices(size(args), product(batch_shape)))
        do i = 1, size(args)
            call elementwise_broadcast_indices(args(i)%dims, batch_shape, indices(i,:), valid)
            if (.not. valid) then
                deallocate(indices)
                return
            end if
        end do
    end subroutine elementwise_broadcast_typed_argument_list

    subroutine elementwise_broadcast_indices(arg_shape, batch_shape, indices, valid)
        integer, intent(in) :: arg_shape(:) !! source dimensions
        integer, intent(in) :: batch_shape(:) !! grid dimensions
        integer, intent(out) :: indices(:) !! one-based source element indices
        logical, intent(inout) :: valid !! false for incompatible trailing axes
        integer :: offset, i, d, rem, coord, ad, pos, stride

        if (.not. valid) return
        if (size(indices) /= product(batch_shape) .or. &
                size(arg_shape) > size(batch_shape)) then
            valid = .false.
            return
        end if
        offset = size(batch_shape) - size(arg_shape)
        do d = 1, size(arg_shape)
            if (arg_shape(d) /= 1 .and. arg_shape(d) /= batch_shape(d + offset)) then
                valid = .false.
                return
            end if
        end do
        do i = 1, size(indices)
            rem = i - 1
            pos = 1
            stride = 1
            do d = 1, size(batch_shape)
                coord = mod(rem, batch_shape(d))
                rem = rem / batch_shape(d)
                ad = d - offset
                if (ad >= 1) then
                    if (arg_shape(ad) /= 1) pos = pos + coord * stride
                    stride = stride * arg_shape(ad)
                end if
            end do
            indices(i) = pos
        end do
    end subroutine elementwise_broadcast_indices

    subroutine elementwise_gather_typed_indices(indices, lanes, selected)
        ! Gather active offsets into caller-owned storage, without allocating
        ! a temporary for the vector-subscripted rank-two array expression.
        ! The selected columns may repeat (initial endpoints or broadcasts).
        integer, intent(in) :: indices(:,:) !! (nargs, nproblems) source offsets
        integer, intent(in) :: lanes(:) !! active problem IDs, possibly repeated
        integer, intent(out) :: selected(:,:) !! (nargs, nactive) callback offsets
        integer :: k

        do k = 1, size(lanes)
            selected(:,k) = indices(:,lanes(k))
        end do
    end subroutine elementwise_gather_typed_indices

    subroutine elementwise_pack_typed_arguments(args, indices, packed)
        ! The solver selects a new active-lane index list at each evaluation.
        ! Copy each scalar through a helper: GNU Fortran 14.2 miscompiles
        ! direct nested polymorphic ALLOCATE(SOURCE=...) here.
        type(elementwise_typed_argument), intent(in) :: args(:) !! input snapshots
        integer, intent(in) :: indices(:,:) !! (nargs, active evaluations)
        type(elementwise_typed_argument), allocatable, intent(out) :: packed(:) !! callback values
        integer :: i, j

        allocate(packed(size(args)))
        do i = 1, size(args)
            packed(i)%dims = [size(indices, 2)]
            allocate(packed(i)%values(size(indices, 2)))
            do j = 1, size(indices, 2)
                packed(i)%values(j) = elementwise_box_value(args(i)%values(indices(i,j))%value)
            end do
        end do
    end subroutine elementwise_pack_typed_arguments

    function elementwise_argument_from(arg) result(a)
        real(dp), intent(in) :: arg(..) !! arbitrary-rank real objective argument
        type(elementwise_argument) :: a !! self-contained parameter and shape
        allocate(a%dims(rank(arg)))
        a%dims = shape(arg)
        call elementwise_flatten(arg, a%values)
    end function elementwise_argument_from

    subroutine elementwise_broadcast_argument_list(args, batch_shape, values, valid)
        ! Each argument has an independent trailing-axis broadcast shape.
        ! Each output row corresponds to one argument and each column to a
        ! lane.  The existing active-lane callback contract is unchanged.
        type(elementwise_argument), intent(in) :: args(:) !! independent arguments
        integer, intent(in) :: batch_shape(:) !! problem grid dimensions
        real(dp), allocatable, intent(out) :: values(:,:) !! (nargs, nproblems)
        logical, intent(inout) :: valid !! false for any invalid parameter
        integer :: i

        if (.not. valid) return
        if (size(args) == 0) then
            valid = .false.
            return
        end if
        ! Validate every argument before evaluating the objective.  The
        ! constructor always allocates both components, but default-created
        ! arguments must fail safely instead of accessing unallocated data.
        do i = 1, size(args)
            if (.not. allocated(args(i)%dims) .or. &
                    .not. allocated(args(i)%values)) then
                valid = .false.
                return
            end if
        end do
        allocate(values(size(args),product(batch_shape)))
        do i = 1, size(args)
            call elementwise_broadcast_flat(args(i)%dims, args(i)%values, &
                batch_shape, values(i,:), valid)
            if (.not. valid) then
                deallocate(values)
                return
            end if
        end do
    end subroutine elementwise_broadcast_argument_list

    subroutine elementwise_flatten(arg, values)
        ! Copy any rank (including scalar) in Fortran element order.
        ! Keep SELECT RANK in one place, never in the numerical controllers.
        real(dp), intent(in) :: arg(..) !! arbitrary-rank input
        real(dp), allocatable, intent(out) :: values(:) !! rank-one copy

        select rank(arg)
        rank(0)
            values = [arg]
        rank(1)
            values = reshape(arg, [size(arg)])
        rank(2)
            values = reshape(arg, [size(arg)])
        rank(3)
            values = reshape(arg, [size(arg)])
        rank(4)
            values = reshape(arg, [size(arg)])
        rank(5)
            values = reshape(arg, [size(arg)])
        rank(6)
            values = reshape(arg, [size(arg)])
        rank(7)
            values = reshape(arg, [size(arg)])
        rank(8)
            values = reshape(arg, [size(arg)])
        rank(9)
            values = reshape(arg, [size(arg)])
        rank(10)
            values = reshape(arg, [size(arg)])
        rank(11)
            values = reshape(arg, [size(arg)])
        rank(12)
            values = reshape(arg, [size(arg)])
        rank(13)
            values = reshape(arg, [size(arg)])
        rank(14)
            values = reshape(arg, [size(arg)])
        rank(15)
            values = reshape(arg, [size(arg)])
        rank default
            allocate(values(0))
        end select
    end subroutine elementwise_flatten

    subroutine elementwise_broadcast_optional(arg, batch_shape, values, valid)
        ! Leave an absent argument unallocated for forwarding to an optional dummy.
        ! A present argument is expanded to the flattened problem grid.
        real(dp), intent(in), optional :: arg(..) !! scalar or broadcastable array
        integer, intent(in) :: batch_shape(:) !! problem-grid dimensions
        real(dp), allocatable, intent(out) :: values(:) !! expanded values, or absent
        logical, intent(inout) :: valid !! cumulative broadcast validity

        if (.not. valid) return
        if (.not. present(arg)) return
        allocate(values(product(batch_shape)))
        values = 0.0_dp
        call elementwise_broadcast(arg, batch_shape, values, valid)
    end subroutine elementwise_broadcast_optional

    subroutine elementwise_broadcast(arg, batch_shape, values, valid)
        real(dp), intent(in), optional :: arg(..) !! scalar through rank-15 parameter
        integer, intent(in) :: batch_shape(:) !! result array shape
        real(dp), intent(inout) :: values(:) !! defaults, replaced when present
        logical, intent(inout) :: valid !! false for incompatible broadcasting
        real(dp), allocatable :: flat(:)
        integer, allocatable :: arg_shape(:)

        if (.not. valid) return
        if (size(values) /= product(batch_shape)) then
            valid = .false.
            return
        end if
        if (.not. present(arg)) return
        arg_shape = shape(arg)
        call elementwise_flatten(arg, flat)
        call elementwise_broadcast_flat(arg_shape, flat, batch_shape, values, valid)
    end subroutine elementwise_broadcast

    subroutine elementwise_broadcast_flat(arg_shape, flat, batch_shape, values, valid)
        ! Shared Fortran-order trailing-axis broadcast for already flattened
        ! arbitrary-rank data (also used for independently shaped args).
        integer, intent(in) :: arg_shape(:) !! source dimensions
        real(dp), intent(in) :: flat(:) !! source in Fortran array order
        integer, intent(in) :: batch_shape(:) !! destination grid dimensions
        real(dp), intent(inout) :: values(:) !! one value per problem
        logical, intent(inout) :: valid !! cumulative validity
        integer :: i, d, ad, offset, rem, coord, pos, stride

        if (.not. valid) return
        if (size(values) /= product(batch_shape)) then
            valid = .false.
            return
        end if
        if (size(flat) /= product(arg_shape)) then
            valid = .false.
            return
        end if
        if (size(arg_shape) == 0) then
            values = flat(1)
            return
        end if
        if (size(arg_shape) > size(batch_shape)) then
            valid = .false.
            return
        end if
        offset = size(batch_shape) - size(arg_shape)
        do d = 1, size(arg_shape)
            if (arg_shape(d) /= 1 .and. arg_shape(d) /= batch_shape(d + offset)) then
                valid = .false.
                return
            end if
        end do
        ! Singleton dimensions use coordinate 1.  The first Fortran array
        ! dimension varies fastest, also when a destination is noncontiguous.
        do i = 1, size(values)
            rem = i - 1
            pos = 1
            stride = 1
            do d = 1, size(batch_shape)
                coord = mod(rem, batch_shape(d))
                rem = rem / batch_shape(d)
                ad = d - offset
                if (ad >= 1) then
                    if (arg_shape(ad) /= 1) pos = pos + coord * stride
                    stride = stride * arg_shape(ad)
                end if
            end do
            values(i) = flat(pos)
        end do
    end subroutine elementwise_broadcast_flat

    subroutine elementwise_broadcast_arguments(args, batch_shape, values, valid)
        ! Broadcast an objective-parameter stack over the problem grid.
        ! The first axis counts parameters; remaining axes align to the
        ! trailing axes of the result grid (NumPy-style broadcasting).
        real(dp), intent(in) :: args(..) !! (nargs, optional grid axes)
        integer, intent(in) :: batch_shape(:) !! problem grid shape
        real(dp), allocatable, intent(out) :: values(:,:) !! (nargs, nproblems)
        logical, intent(inout) :: valid !! cumulative shape validity
        real(dp), allocatable :: flat(:)
        integer, allocatable :: dims(:), arg_grid(:)
        integer :: nargs, offset, d, j, rem, coord, pos, stride, ad

        if (.not. valid) return
        dims = shape(args)
        if (size(dims) == 0) then
            valid = .false.
            return
        end if
        nargs = dims(1)
        arg_grid = dims(2:)
        if (nargs < 1 .or. size(arg_grid) > size(batch_shape)) then
            valid = .false.
            return
        end if
        offset = size(batch_shape) - size(arg_grid)
        do d = 1, size(arg_grid)
            if (arg_grid(d) /= 1 .and. arg_grid(d) /= batch_shape(d + offset)) then
                valid = .false.
                return
            end if
        end do
        call elementwise_flatten(args, flat)
        allocate(values(nargs,product(batch_shape)))
        ! Fortran's first index varies fastest: the parameter axis remains
        ! contiguous inside each selected problem-grid cell.
        do j = 1, size(values,2)
            rem = j - 1
            pos = 1
            stride = 1
            do d = 1, size(batch_shape)
                coord = mod(rem, batch_shape(d))
                rem = rem / batch_shape(d)
                ad = d - offset
                if (ad >= 1) then
                    if (arg_grid(ad) /= 1) pos = pos + coord * stride
                    stride = stride * arg_grid(ad)
                end if
            end do
            values(:,j) = flat(nargs*(pos-1)+1:nargs*pos)
        end do
    end subroutine elementwise_broadcast_arguments

end module scifort_optimize_elementwise_batch
