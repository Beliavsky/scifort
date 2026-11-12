! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special_reductions_high_rank_weighted_support
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_special_reductions_low, only : low_logsumexp_signed => logsumexp_signed
    implicit none
    private

    public :: complex_log_from_phase
    public :: complex_nan
    public :: drop_dimension_shape
    public :: flatten_complex_weights
    public :: flatten_real_weights
    public :: weighted_complex_prepared
    public :: weighted_real_prepared

contains

    pure elemental function complex_log_from_phase(logabs, phase) result(value)
        real(dp), intent(in) :: logabs !! logarithm of a complex result magnitude
        complex(dp), intent(in) :: phase !! unit phase, or zero for an exact zero result
        complex(dp) :: value

        if (ieee_is_nan(logabs) .or. ieee_is_nan(real(phase, dp)) .or. &
                ieee_is_nan(aimag(phase))) then
            value = complex_nan()
        else if (abs(phase) == 0.0_dp) then
            value = cmplx(logabs, 0.0_dp, kind=dp)
        else
            value = cmplx(logabs, atan2(aimag(phase), real(phase, dp)), kind=dp)
        end if
    end function complex_log_from_phase

    pure elemental function complex_nan() result(value)
        complex(dp) :: value

        value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
    end function complex_nan

    pure function reduction_dims_mode_nd(rank_value, dims) result(mode)
        integer, intent(in) :: rank_value !! rank of the input array, from 1 through 15
        integer, intent(in) :: dims(:) !! unique Fortran dimensions selected for reduction
        integer :: mode

        integer :: i

        mode = 0
        do i = 1, size(dims)
            if (dims(i) < 1 .or. dims(i) > rank_value) then
                mode = -1
                return
            end if
            if (btest(mode, dims(i) - 1)) then
                mode = -1
                return
            end if
            mode = ibset(mode, dims(i) - 1)
        end do
    end function reduction_dims_mode_nd

    pure subroutine prepare_reduction_layout(xshape, dims, mode, outshape, strides, &
            reduction_size, valid)
        integer, intent(in) :: xshape(:) !! input extents in Fortran dimension order
        integer, intent(in) :: dims(:) !! unique dimensions selected for reduction
        integer, intent(out) :: mode !! bit mask identifying reduced dimensions
        integer, intent(out) :: outshape(size(xshape)) !! keepdims output extents
        integer, intent(out) :: strides(size(xshape)) !! column-major input strides
        integer, intent(out) :: reduction_size !! number of elements in each reduction group
        logical, intent(out) :: valid !! true when dimensions are unique and in range

        integer :: d

        mode = reduction_dims_mode_nd(size(xshape), dims)
        valid = mode >= 0
        if (.not. valid) then
            outshape = 0
            strides = 0
            reduction_size = 0
            return
        end if

        outshape = xshape
        strides(1) = 1
        do d = 2, size(xshape)
            strides(d) = strides(d - 1) * xshape(d - 1)
        end do

        reduction_size = 1
        do d = 1, size(xshape)
            if (btest(mode, d - 1)) then
                reduction_size = reduction_size * xshape(d)
                outshape(d) = 1
            end if
        end do
    end subroutine prepare_reduction_layout

    pure subroutine drop_dimension_shape(xshape, dim, outshape)
        integer, intent(in) :: xshape(:) !! input extents
        integer, intent(in) :: dim !! dimension omitted from the output shape
        integer, intent(out) :: outshape(size(xshape) - 1) !! extents with dim removed

        integer :: d
        integer :: j

        j = 0
        do d = 1, size(xshape)
            if (d == dim) cycle
            j = j + 1
            outshape(j) = xshape(d)
        end do
    end subroutine drop_dimension_shape

    pure subroutine reduction_group_indices(xshape, outshape, strides, mode, output_index, &
            indices, count)
        integer, intent(in) :: xshape(:) !! broadcast extents in Fortran dimension order
        integer, intent(in) :: outshape(:) !! keepdims output extents
        integer, intent(in) :: strides(:) !! column-major broadcast strides
        integer, intent(in) :: mode !! bit mask identifying reduced dimensions
        integer, intent(in) :: output_index !! one-based flattened keepdims group index
        integer, intent(out) :: indices(:) !! flattened broadcast indices for this group
        integer, intent(out) :: count !! number of populated entries in indices

        integer :: base
        integer :: coord
        integer :: d
        integer :: input_index
        integer :: q
        integer :: remainder
        integer :: reduction_size

        base = 0
        remainder = output_index - 1
        do d = 1, size(xshape)
            if (outshape(d) > 0) then
                coord = mod(remainder, outshape(d)) + 1
                remainder = remainder / outshape(d)
            else
                coord = 1
            end if
            if (.not. btest(mode, d - 1)) base = base + (coord - 1) * strides(d)
        end do

        reduction_size = 1
        do d = 1, size(xshape)
            if (btest(mode, d - 1)) reduction_size = reduction_size * xshape(d)
        end do

        count = reduction_size
        do q = 0, reduction_size - 1
            remainder = q
            input_index = base + 1
            do d = 1, size(xshape)
                if (.not. btest(mode, d - 1)) cycle
                coord = mod(remainder, xshape(d)) + 1
                remainder = remainder / xshape(d)
                input_index = input_index + (coord - 1) * strides(d)
            end do
            indices(q + 1) = input_index
        end do
    end subroutine reduction_group_indices

    pure subroutine prepare_broadcast_maps(xshape, bshape, joint_shape, xmap, bmap, valid)
        integer, intent(in) :: xshape(:) !! input extents
        integer, intent(in) :: bshape(:) !! weight extents before trailing alignment
        integer, intent(out) :: joint_shape(size(xshape)) !! NumPy-compatible joint extents
        integer, allocatable, intent(out) :: xmap(:) !! joint-index to input-index map
        integer, allocatable, intent(out) :: bmap(:) !! joint-index to weight-index map
        logical, intent(out) :: valid !! true when the shapes broadcast together

        integer, allocatable :: bstrides(:)
        integer :: bcoord
        integer :: bextent
        integer :: bj
        integer :: brank
        integer :: coord
        integer :: d
        integer :: j
        integer :: k
        integer :: remainder
        integer :: xcoord
        integer :: xindex
        integer :: xstrides(size(xshape))

        brank = size(bshape)
        valid = brank <= size(xshape)
        joint_shape = xshape
        if (.not. valid) then
            allocate(xmap(0), bmap(0))
            return
        end if

        do d = 1, size(xshape)
            bj = d - (size(xshape) - brank)
            if (bj > 0) then
                bextent = bshape(bj)
            else
                bextent = 1
            end if
            if (xshape(d) == bextent) then
                joint_shape(d) = xshape(d)
            else if (xshape(d) == 1) then
                joint_shape(d) = bextent
            else if (bextent == 1) then
                joint_shape(d) = xshape(d)
            else
                valid = .false.
                allocate(xmap(0), bmap(0))
                return
            end if
        end do

        xstrides(1) = 1
        do d = 2, size(xshape)
            xstrides(d) = xstrides(d - 1) * xshape(d - 1)
        end do
        allocate(bstrides(brank))
        if (brank > 0) then
            bstrides(1) = 1
            do j = 2, brank
                bstrides(j) = bstrides(j - 1) * bshape(j - 1)
            end do
        end if

        allocate(xmap(product(joint_shape)), bmap(product(joint_shape)))
        do k = 1, size(xmap)
            remainder = k - 1
            xindex = 1
            j = 1
            do d = 1, size(xshape)
                coord = mod(remainder, joint_shape(d)) + 1
                remainder = remainder / joint_shape(d)
                if (xshape(d) == 1) then
                    xcoord = 1
                else
                    xcoord = coord
                end if
                xindex = xindex + (xcoord - 1) * xstrides(d)

                bj = d - (size(xshape) - brank)
                if (bj > 0) then
                    if (bshape(bj) == 1) then
                        bcoord = 1
                    else
                        bcoord = coord
                    end if
                    j = j + (bcoord - 1) * bstrides(bj)
                end if
            end do
            xmap(k) = xindex
            bmap(k) = j
        end do
    end subroutine prepare_broadcast_maps

    pure subroutine weighted_real_prepared(x, xshape, dims, b, bshape, logabs, sign_value, &
            outshape, valid_dims)
        real(dp), intent(in) :: x(:) !! flattened real input values
        integer, intent(in) :: xshape(:) !! input extents
        integer, intent(in) :: dims(:) !! dimensions selected for reduction
        real(dp), intent(in) :: b(:) !! flattened real weights
        integer, intent(in) :: bshape(:) !! original weight extents
        real(dp), allocatable, intent(out) :: logabs(:) !! flattened keepdims log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:) !! flattened keepdims signs
        integer, intent(out) :: outshape(size(xshape)) !! keepdims result extents
        logical, intent(out) :: valid_dims !! true when dimensions are unique and in range

        real(dp), allocatable :: buffer(:)
        real(dp), allocatable :: weight_buffer(:)
        integer, allocatable :: bmap(:)
        integer :: count
        integer, allocatable :: group_indices(:)
        integer :: joint_shape(size(xshape))
        logical :: broadcast_ok
        integer :: mode
        integer :: o
        integer :: reduction_size
        integer :: strides(size(xshape))
        integer, allocatable :: xmap(:)

        if (reduction_dims_mode_nd(size(xshape), dims) < 0) then
            outshape = 0
            valid_dims = .false.
            allocate(logabs(0), sign_value(0))
            return
        end if

        call prepare_broadcast_maps(xshape, bshape, joint_shape, xmap, bmap, broadcast_ok)
        if (.not. broadcast_ok) then
            call prepare_reduction_layout(xshape, dims, mode, outshape, strides, &
                reduction_size, valid_dims)
            allocate(logabs(product(outshape)), sign_value(product(outshape)))
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(sign_value) > 0) sign_value = quiet_nan(0.0_dp)
            return
        end if

        call prepare_reduction_layout(joint_shape, dims, mode, outshape, strides, &
            reduction_size, valid_dims)
        allocate(logabs(product(outshape)), sign_value(product(outshape)))
        allocate(buffer(max(1, reduction_size)), weight_buffer(max(1, reduction_size)))
        allocate(group_indices(max(1, reduction_size)))
        do o = 1, size(logabs)
            call reduction_group_indices(joint_shape, outshape, strides, mode, o, &
                group_indices, count)
            if (count > 0) then
                buffer(1:count) = x(xmap(group_indices(1:count)))
                weight_buffer(1:count) = b(bmap(group_indices(1:count)))
            end if
            call low_logsumexp_signed(buffer(1:count), logabs(o), sign_value(o), &
                weight_buffer(1:count))
        end do
    end subroutine weighted_real_prepared

    pure subroutine weighted_complex_prepared(x, xshape, dims, b, bshape, logabs, phase, &
            outshape, valid_dims)
        complex(dp), intent(in) :: x(:) !! flattened complex input values
        integer, intent(in) :: xshape(:) !! input extents
        integer, intent(in) :: dims(:) !! dimensions selected for reduction
        complex(dp), intent(in) :: b(:) !! flattened complex weights
        integer, intent(in) :: bshape(:) !! original weight extents
        real(dp), allocatable, intent(out) :: logabs(:) !! flattened keepdims log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:) !! flattened keepdims unit phases
        integer, intent(out) :: outshape(size(xshape)) !! keepdims result extents
        logical, intent(out) :: valid_dims !! true when dimensions are unique and in range

        complex(dp), allocatable :: buffer(:)
        complex(dp), allocatable :: weight_buffer(:)
        integer, allocatable :: bmap(:)
        integer :: count
        integer, allocatable :: group_indices(:)
        integer :: joint_shape(size(xshape))
        logical :: broadcast_ok
        integer :: mode
        integer :: o
        integer :: reduction_size
        integer :: strides(size(xshape))
        integer, allocatable :: xmap(:)

        if (reduction_dims_mode_nd(size(xshape), dims) < 0) then
            outshape = 0
            valid_dims = .false.
            allocate(logabs(0), phase(0))
            return
        end if

        call prepare_broadcast_maps(xshape, bshape, joint_shape, xmap, bmap, broadcast_ok)
        if (.not. broadcast_ok) then
            call prepare_reduction_layout(xshape, dims, mode, outshape, strides, &
                reduction_size, valid_dims)
            allocate(logabs(product(outshape)), phase(product(outshape)))
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(phase) > 0) phase = complex_nan()
            return
        end if

        call prepare_reduction_layout(joint_shape, dims, mode, outshape, strides, &
            reduction_size, valid_dims)
        allocate(logabs(product(outshape)), phase(product(outshape)))
        allocate(buffer(max(1, reduction_size)), weight_buffer(max(1, reduction_size)))
        allocate(group_indices(max(1, reduction_size)))
        do o = 1, size(logabs)
            call reduction_group_indices(joint_shape, outshape, strides, mode, o, &
                group_indices, count)
            if (count > 0) then
                buffer(1:count) = x(xmap(group_indices(1:count)))
                weight_buffer(1:count) = b(bmap(group_indices(1:count)))
            end if
            call low_logsumexp_signed(buffer(1:count), logabs(o), phase(o), &
                weight_buffer(1:count))
        end do
    end subroutine weighted_complex_prepared

    pure subroutine flatten_real_weights(b, flat, bshape)
        real(dp), intent(in) :: b(..) !! scalar or array weights to flatten
        real(dp), allocatable, intent(out) :: flat(:) !! weights in Fortran element order
        integer, allocatable, intent(out) :: bshape(:) !! original weight extents

        select rank (b)
        rank (0)
            allocate(flat(1), bshape(0))
            flat(1) = b
        rank (1)
            allocate(bshape(1))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (2)
            allocate(bshape(2))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (3)
            allocate(bshape(3))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (4)
            allocate(bshape(4))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (5)
            allocate(bshape(5))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (6)
            allocate(bshape(6))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (7)
            allocate(bshape(7))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (8)
            allocate(bshape(8))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (9)
            allocate(bshape(9))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (10)
            allocate(bshape(10))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (11)
            allocate(bshape(11))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (12)
            allocate(bshape(12))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (13)
            allocate(bshape(13))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (14)
            allocate(bshape(14))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (15)
            allocate(bshape(15))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank default
            allocate(flat(0), bshape(16))
            bshape = -1
        end select
    end subroutine flatten_real_weights

    pure subroutine flatten_complex_weights(b, flat, bshape)
        complex(dp), intent(in) :: b(..) !! scalar or array weights to flatten
        complex(dp), allocatable, intent(out) :: flat(:) !! weights in Fortran element order
        integer, allocatable, intent(out) :: bshape(:) !! original weight extents

        select rank (b)
        rank (0)
            allocate(flat(1), bshape(0))
            flat(1) = b
        rank (1)
            allocate(bshape(1))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (2)
            allocate(bshape(2))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (3)
            allocate(bshape(3))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (4)
            allocate(bshape(4))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (5)
            allocate(bshape(5))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (6)
            allocate(bshape(6))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (7)
            allocate(bshape(7))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (8)
            allocate(bshape(8))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (9)
            allocate(bshape(9))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (10)
            allocate(bshape(10))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (11)
            allocate(bshape(11))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (12)
            allocate(bshape(12))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (13)
            allocate(bshape(13))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (14)
            allocate(bshape(14))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank (15)
            allocate(bshape(15))
            bshape = shape(b)
            allocate(flat(size(b)))
            if (size(b) > 0) flat = reshape(b, [size(b)])
        rank default
            allocate(flat(0), bshape(16))
            bshape = -1
        end select
    end subroutine flatten_complex_weights

end module scifort_special_reductions_high_rank_weighted_support
