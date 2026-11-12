! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special_reductions_high_rank
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_special_reductions_low, only : low_log_softmax => log_softmax, &
        low_logsumexp_signed => logsumexp_signed, low_softmax => softmax
    implicit none
    private

    public :: log_softmax
    public :: logsumexp
    public :: logsumexp_signed
    public :: logsumexp_keepdims
    public :: logsumexp_signed_keepdims
    public :: softmax

    interface logsumexp
        module procedure logsumexp_7d_all
        module procedure logsumexp_7d_dim
        module procedure logsumexp_complex_7d_all
        module procedure logsumexp_complex_7d_dim
        module procedure logsumexp_8d_all
        module procedure logsumexp_8d_dim
        module procedure logsumexp_complex_8d_all
        module procedure logsumexp_complex_8d_dim
        module procedure logsumexp_9d_all
        module procedure logsumexp_9d_dim
        module procedure logsumexp_complex_9d_all
        module procedure logsumexp_complex_9d_dim
        module procedure logsumexp_10d_all
        module procedure logsumexp_10d_dim
        module procedure logsumexp_complex_10d_all
        module procedure logsumexp_complex_10d_dim
        module procedure logsumexp_11d_all
        module procedure logsumexp_11d_dim
        module procedure logsumexp_complex_11d_all
        module procedure logsumexp_complex_11d_dim
        module procedure logsumexp_12d_all
        module procedure logsumexp_12d_dim
        module procedure logsumexp_complex_12d_all
        module procedure logsumexp_complex_12d_dim
        module procedure logsumexp_13d_all
        module procedure logsumexp_13d_dim
        module procedure logsumexp_complex_13d_all
        module procedure logsumexp_complex_13d_dim
        module procedure logsumexp_14d_all
        module procedure logsumexp_14d_dim
        module procedure logsumexp_complex_14d_all
        module procedure logsumexp_complex_14d_dim
        module procedure logsumexp_15d_all
        module procedure logsumexp_15d_dim
        module procedure logsumexp_complex_15d_all
        module procedure logsumexp_complex_15d_dim
    end interface logsumexp

    interface logsumexp_signed
        module procedure logsumexp_signed_7d_all
        module procedure logsumexp_signed_7d_dim
        module procedure logsumexp_signed_complex_7d_all
        module procedure logsumexp_signed_complex_7d_dim
        module procedure logsumexp_signed_8d_all
        module procedure logsumexp_signed_8d_dim
        module procedure logsumexp_signed_complex_8d_all
        module procedure logsumexp_signed_complex_8d_dim
        module procedure logsumexp_signed_9d_all
        module procedure logsumexp_signed_9d_dim
        module procedure logsumexp_signed_complex_9d_all
        module procedure logsumexp_signed_complex_9d_dim
        module procedure logsumexp_signed_10d_all
        module procedure logsumexp_signed_10d_dim
        module procedure logsumexp_signed_complex_10d_all
        module procedure logsumexp_signed_complex_10d_dim
        module procedure logsumexp_signed_11d_all
        module procedure logsumexp_signed_11d_dim
        module procedure logsumexp_signed_complex_11d_all
        module procedure logsumexp_signed_complex_11d_dim
        module procedure logsumexp_signed_12d_all
        module procedure logsumexp_signed_12d_dim
        module procedure logsumexp_signed_complex_12d_all
        module procedure logsumexp_signed_complex_12d_dim
        module procedure logsumexp_signed_13d_all
        module procedure logsumexp_signed_13d_dim
        module procedure logsumexp_signed_complex_13d_all
        module procedure logsumexp_signed_complex_13d_dim
        module procedure logsumexp_signed_14d_all
        module procedure logsumexp_signed_14d_dim
        module procedure logsumexp_signed_complex_14d_all
        module procedure logsumexp_signed_complex_14d_dim
        module procedure logsumexp_signed_15d_all
        module procedure logsumexp_signed_15d_dim
        module procedure logsumexp_signed_complex_15d_all
        module procedure logsumexp_signed_complex_15d_dim
    end interface logsumexp_signed

    interface logsumexp_keepdims
        module procedure logsumexp_keepdims_7d
        module procedure logsumexp_keepdims_complex_7d
        module procedure logsumexp_keepdims_8d
        module procedure logsumexp_keepdims_complex_8d
        module procedure logsumexp_keepdims_9d
        module procedure logsumexp_keepdims_complex_9d
        module procedure logsumexp_keepdims_10d
        module procedure logsumexp_keepdims_complex_10d
        module procedure logsumexp_keepdims_11d
        module procedure logsumexp_keepdims_complex_11d
        module procedure logsumexp_keepdims_12d
        module procedure logsumexp_keepdims_complex_12d
        module procedure logsumexp_keepdims_13d
        module procedure logsumexp_keepdims_complex_13d
        module procedure logsumexp_keepdims_14d
        module procedure logsumexp_keepdims_complex_14d
        module procedure logsumexp_keepdims_15d
        module procedure logsumexp_keepdims_complex_15d
    end interface logsumexp_keepdims

    interface logsumexp_signed_keepdims
        module procedure logsumexp_signed_keepdims_7d
        module procedure logsumexp_signed_keepdims_complex_7d
        module procedure logsumexp_signed_keepdims_8d
        module procedure logsumexp_signed_keepdims_complex_8d
        module procedure logsumexp_signed_keepdims_9d
        module procedure logsumexp_signed_keepdims_complex_9d
        module procedure logsumexp_signed_keepdims_10d
        module procedure logsumexp_signed_keepdims_complex_10d
        module procedure logsumexp_signed_keepdims_11d
        module procedure logsumexp_signed_keepdims_complex_11d
        module procedure logsumexp_signed_keepdims_12d
        module procedure logsumexp_signed_keepdims_complex_12d
        module procedure logsumexp_signed_keepdims_13d
        module procedure logsumexp_signed_keepdims_complex_13d
        module procedure logsumexp_signed_keepdims_14d
        module procedure logsumexp_signed_keepdims_complex_14d
        module procedure logsumexp_signed_keepdims_15d
        module procedure logsumexp_signed_keepdims_complex_15d
    end interface logsumexp_signed_keepdims

    interface softmax
        module procedure softmax_7d
        module procedure softmax_7d_dims
        module procedure softmax_complex_7d
        module procedure softmax_complex_7d_dims
        module procedure softmax_8d
        module procedure softmax_8d_dims
        module procedure softmax_complex_8d
        module procedure softmax_complex_8d_dims
        module procedure softmax_9d
        module procedure softmax_9d_dims
        module procedure softmax_complex_9d
        module procedure softmax_complex_9d_dims
        module procedure softmax_10d
        module procedure softmax_10d_dims
        module procedure softmax_complex_10d
        module procedure softmax_complex_10d_dims
        module procedure softmax_11d
        module procedure softmax_11d_dims
        module procedure softmax_complex_11d
        module procedure softmax_complex_11d_dims
        module procedure softmax_12d
        module procedure softmax_12d_dims
        module procedure softmax_complex_12d
        module procedure softmax_complex_12d_dims
        module procedure softmax_13d
        module procedure softmax_13d_dims
        module procedure softmax_complex_13d
        module procedure softmax_complex_13d_dims
        module procedure softmax_14d
        module procedure softmax_14d_dims
        module procedure softmax_complex_14d
        module procedure softmax_complex_14d_dims
        module procedure softmax_15d
        module procedure softmax_15d_dims
        module procedure softmax_complex_15d
        module procedure softmax_complex_15d_dims
    end interface softmax

    interface log_softmax
        module procedure log_softmax_7d
        module procedure log_softmax_7d_dims
        module procedure log_softmax_complex_7d
        module procedure log_softmax_complex_7d_dims
        module procedure log_softmax_8d
        module procedure log_softmax_8d_dims
        module procedure log_softmax_complex_8d
        module procedure log_softmax_complex_8d_dims
        module procedure log_softmax_9d
        module procedure log_softmax_9d_dims
        module procedure log_softmax_complex_9d
        module procedure log_softmax_complex_9d_dims
        module procedure log_softmax_10d
        module procedure log_softmax_10d_dims
        module procedure log_softmax_complex_10d
        module procedure log_softmax_complex_10d_dims
        module procedure log_softmax_11d
        module procedure log_softmax_11d_dims
        module procedure log_softmax_complex_11d
        module procedure log_softmax_complex_11d_dims
        module procedure log_softmax_12d
        module procedure log_softmax_12d_dims
        module procedure log_softmax_complex_12d
        module procedure log_softmax_complex_12d_dims
        module procedure log_softmax_13d
        module procedure log_softmax_13d_dims
        module procedure log_softmax_complex_13d
        module procedure log_softmax_complex_13d_dims
        module procedure log_softmax_14d
        module procedure log_softmax_14d_dims
        module procedure log_softmax_complex_14d
        module procedure log_softmax_complex_14d_dims
        module procedure log_softmax_15d
        module procedure log_softmax_15d_dims
        module procedure log_softmax_complex_15d
        module procedure log_softmax_complex_15d_dims
    end interface log_softmax

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

    ! Shared flattened implementation for ranks above six. The public API still
    ! uses fixed-rank wrappers because Fortran function result rank is static,
    ! while axis selection and group traversal live in these helpers.
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
        integer, intent(out) :: reduction_size !! number of input elements in each group
        logical, intent(out) :: valid !! true when all dimensions are unique and in range

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
        integer, intent(in) :: xshape(:) !! input extents in Fortran dimension order
        integer, intent(in) :: outshape(:) !! keepdims output extents
        integer, intent(in) :: strides(:) !! column-major input strides
        integer, intent(in) :: mode !! bit mask identifying reduced dimensions
        integer, intent(in) :: output_index !! one-based flattened keepdims group index
        integer, intent(out) :: indices(:) !! flattened input indices for this reduction group
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

    pure subroutine logsumexp_signed_keepdims_real_flat(x, xshape, dims, logabs, &
            sign_value, outshape, valid)
        real(dp), intent(in) :: x(:) !! flattened real input values
        integer, intent(in) :: xshape(:) !! original input extents
        integer, intent(in) :: dims(:) !! dimensions selected for reduction
        real(dp), allocatable, intent(out) :: logabs(:) !! flattened keepdims log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:) !! flattened keepdims signs
        integer, intent(out) :: outshape(size(xshape)) !! keepdims output extents
        logical, intent(out) :: valid !! true when the dimension list is valid

        real(dp), allocatable :: buffer(:)
        integer, allocatable :: indices(:)
        integer :: count
        integer :: mode
        integer :: o
        integer :: reduction_size
        integer :: strides(size(xshape))

        call prepare_reduction_layout(xshape, dims, mode, outshape, strides, reduction_size, valid)
        if (.not. valid) then
            allocate(logabs(0), sign_value(0))
            return
        end if

        allocate(logabs(product(outshape)), sign_value(product(outshape)))
        allocate(buffer(max(1, reduction_size)), indices(max(1, reduction_size)))
        do o = 1, size(logabs)
            call reduction_group_indices(xshape, outshape, strides, mode, o, indices, count)
            if (count > 0) buffer(1:count) = x(indices(1:count))
            call low_logsumexp_signed(buffer(1:count), logabs(o), sign_value(o))
        end do
    end subroutine logsumexp_signed_keepdims_real_flat

    pure subroutine logsumexp_signed_keepdims_complex_flat(x, xshape, dims, logabs, phase, &
            outshape, valid)
        complex(dp), intent(in) :: x(:) !! flattened complex input values
        integer, intent(in) :: xshape(:) !! original input extents
        integer, intent(in) :: dims(:) !! dimensions selected for reduction
        real(dp), allocatable, intent(out) :: logabs(:) !! flattened keepdims log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:) !! flattened keepdims unit phases
        integer, intent(out) :: outshape(size(xshape)) !! keepdims output extents
        logical, intent(out) :: valid !! true when the dimension list is valid

        complex(dp), allocatable :: buffer(:)
        integer, allocatable :: indices(:)
        integer :: count
        integer :: mode
        integer :: o
        integer :: reduction_size
        integer :: strides(size(xshape))

        call prepare_reduction_layout(xshape, dims, mode, outshape, strides, reduction_size, valid)
        if (.not. valid) then
            allocate(logabs(0), phase(0))
            return
        end if

        allocate(logabs(product(outshape)), phase(product(outshape)))
        allocate(buffer(max(1, reduction_size)), indices(max(1, reduction_size)))
        do o = 1, size(logabs)
            call reduction_group_indices(xshape, outshape, strides, mode, o, indices, count)
            if (count > 0) buffer(1:count) = x(indices(1:count))
            call low_logsumexp_signed(buffer(1:count), logabs(o), phase(o))
        end do
    end subroutine logsumexp_signed_keepdims_complex_flat

    pure subroutine normalize_real_flat_dims(x, xshape, dims, logarithmic, values, valid)
        real(dp), intent(in) :: x(:) !! flattened real input values
        integer, intent(in) :: xshape(:) !! original input extents
        integer, intent(in) :: dims(:) !! dimensions normalized jointly
        logical, intent(in) :: logarithmic !! return log-softmax when true
        real(dp), allocatable, intent(out) :: values(:) !! flattened normalized result
        logical, intent(out) :: valid !! true when the dimension list is valid

        real(dp), allocatable :: buffer(:)
        integer, allocatable :: indices(:)
        real(dp), allocatable :: reduced(:)
        integer :: count
        integer :: mode
        integer :: o
        integer :: outshape(size(xshape))
        integer :: reduction_size
        integer :: strides(size(xshape))

        call prepare_reduction_layout(xshape, dims, mode, outshape, strides, reduction_size, valid)
        allocate(values(size(x)))
        if (.not. valid) then
            values = quiet_nan(0.0_dp)
            return
        end if
        if (mode == 0) then
            if (logarithmic) then
                where (ieee_is_finite(x))
                    values = 0.0_dp
                elsewhere
                    values = quiet_nan(0.0_dp)
                end where
            else
                where (ieee_is_finite(x))
                    values = 1.0_dp
                elsewhere
                    values = quiet_nan(0.0_dp)
                end where
            end if
            return
        end if

        allocate(buffer(max(1, reduction_size)), indices(max(1, reduction_size)))
        do o = 1, product(outshape)
            call reduction_group_indices(xshape, outshape, strides, mode, o, indices, count)
            if (count > 0) buffer(1:count) = x(indices(1:count))
            if (logarithmic) then
                reduced = low_log_softmax(buffer(1:count))
            else
                reduced = low_softmax(buffer(1:count))
            end if
            if (count > 0) values(indices(1:count)) = reduced
        end do
    end subroutine normalize_real_flat_dims

    pure subroutine normalize_complex_flat_dims(x, xshape, dims, logarithmic, values, valid)
        complex(dp), intent(in) :: x(:) !! flattened complex input values
        integer, intent(in) :: xshape(:) !! original input extents
        integer, intent(in) :: dims(:) !! dimensions normalized jointly
        logical, intent(in) :: logarithmic !! return log-softmax when true
        complex(dp), allocatable, intent(out) :: values(:) !! flattened normalized result
        logical, intent(out) :: valid !! true when the dimension list is valid

        complex(dp), allocatable :: buffer(:)
        integer, allocatable :: indices(:)
        complex(dp), allocatable :: reduced(:)
        integer :: count
        integer :: mode
        integer :: o
        integer :: outshape(size(xshape))
        integer :: reduction_size
        integer :: strides(size(xshape))

        call prepare_reduction_layout(xshape, dims, mode, outshape, strides, reduction_size, valid)
        allocate(values(size(x)))
        if (.not. valid) then
            values = complex_nan()
            return
        end if
        if (mode == 0) then
            if (logarithmic) then
                where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                    values = cmplx(0.0_dp, 0.0_dp, kind=dp)
                elsewhere
                    values = complex_nan()
                end where
            else
                where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                    values = cmplx(1.0_dp, 0.0_dp, kind=dp)
                elsewhere
                    values = complex_nan()
                end where
            end if
            return
        end if

        allocate(buffer(max(1, reduction_size)), indices(max(1, reduction_size)))
        do o = 1, product(outshape)
            call reduction_group_indices(xshape, outshape, strides, mode, o, indices, count)
            if (count > 0) buffer(1:count) = x(indices(1:count))
            if (logarithmic) then
                reduced = low_log_softmax(buffer(1:count))
            else
                reduced = low_softmax(buffer(1:count))
            end if
            if (count > 0) values(indices(1:count)) = reduced
        end do
    end subroutine normalize_complex_flat_dims


    ! Rank-7 through rank-15 unweighted wrappers use the shared flattened helpers.
    pure subroutine logsumexp_signed_keepdims_7d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(7)

        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_7d

    pure function logsumexp_keepdims_7d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_7d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_7d

    pure function logsumexp_7d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call low_logsumexp_signed(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_7d_all

    pure function logsumexp_7d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        real(dp), allocatable :: values(:, :, :, :, :, :)

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        real(dp), allocatable :: signs_flattened(:, :, :, :, :, :)
        logical :: valid
        integer :: keep_shape(7)
        integer :: outshape(6)

        if (dim < 1 .or. dim > 7) then
        allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(values(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6)))
        if (size(values) == 0) return
        values = reshape(flat_logabs, shape(values))
        allocate(signs_flattened(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6)))
        signs_flattened = reshape(flat_sign, shape(signs_flattened))
        where (signs_flattened < 0.0_dp .or. ieee_is_nan(signs_flattened))
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_7d_dim

    pure subroutine logsumexp_signed_7d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_7d_all

    pure subroutine logsumexp_signed_7d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: keep_shape(7)
        integer :: outshape(6)

        if (dim < 1 .or. dim > 7) then
        allocate(logabs(0, 0, 0, 0, 0, 0), &
            & sign_value(0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_7d_dim

    pure function softmax_7d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7))

        if (present(dim)) then
            values = softmax_7d_dims(x, [dim])
        else
            values = softmax_7d_dims(x, [1, 2, 3, 4, 5, 6, 7])
        end if
    end function softmax_7d

    pure function softmax_7d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_7d_dims

    pure function log_softmax_7d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7))

        if (present(dim)) then
            values = log_softmax_7d_dims(x, [dim])
        else
            values = log_softmax_7d_dims(x, [1, 2, 3, 4, 5, 6, 7])
        end if
    end function log_softmax_7d

    pure function log_softmax_7d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_7d_dims

    pure subroutine logsumexp_signed_keepdims_complex_7d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(7)

        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_7d

    pure function logsumexp_keepdims_complex_7d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_7d(x, dims, logabs, phase)
        if (size(logabs) == 0) then
        allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_7d

    pure function logsumexp_complex_7d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values reduced over all dimensions
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_7d_all

    pure function logsumexp_complex_7d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        complex(dp), allocatable :: values(:, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :)

        if (dim < 1 .or. dim > 7) then
        allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_complex_7d_dim(x, dim, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6)))
        if (size(values) > 0) values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_7d_dim

    pure subroutine logsumexp_signed_complex_7d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the result magnitude
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_7d_all

    pure subroutine logsumexp_signed_complex_7d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: keep_shape(7)
        integer :: outshape(6)

        if (dim < 1 .or. dim > 7) then
        allocate(logabs(0, 0, 0, 0, 0, 0), &
            & phase(0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_phase, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_complex_7d_dim

    pure function softmax_complex_7d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7))

        if (present(dim)) then
            values = softmax_complex_7d_dims(x, [dim])
        else
            values = softmax_complex_7d_dims(x, [1, 2, 3, 4, 5, 6, 7])
        end if
    end function softmax_complex_7d

    pure function softmax_complex_7d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_complex_7d_dims

    pure function log_softmax_complex_7d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7))

        if (present(dim)) then
            values = log_softmax_complex_7d_dims(x, [dim])
        else
            values = log_softmax_complex_7d_dims(x, [1, 2, 3, 4, 5, 6, 7])
        end if
    end function log_softmax_complex_7d

    pure function log_softmax_complex_7d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! rank-7 complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_complex_7d_dims

    pure subroutine logsumexp_signed_keepdims_8d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(8)

        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_8d

    pure function logsumexp_keepdims_8d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_8d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_8d

    pure function logsumexp_8d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call low_logsumexp_signed(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_8d_all

    pure function logsumexp_8d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        real(dp), allocatable :: signs_flattened(:, :, :, :, :, :, :)
        logical :: valid
        integer :: keep_shape(8)
        integer :: outshape(7)

        if (dim < 1 .or. dim > 8) then
        allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(values(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7)))
        if (size(values) == 0) return
        values = reshape(flat_logabs, shape(values))
        allocate(signs_flattened(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), outshape(7)))
        signs_flattened = reshape(flat_sign, shape(signs_flattened))
        where (signs_flattened < 0.0_dp .or. ieee_is_nan(signs_flattened))
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_8d_dim

    pure subroutine logsumexp_signed_8d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_8d_all

    pure subroutine logsumexp_signed_8d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: keep_shape(8)
        integer :: outshape(7)

        if (dim < 1 .or. dim > 8) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0), &
            & sign_value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_8d_dim

    pure function softmax_8d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8))

        if (present(dim)) then
            values = softmax_8d_dims(x, [dim])
        else
            values = softmax_8d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8])
        end if
    end function softmax_8d

    pure function softmax_8d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_8d_dims

    pure function log_softmax_8d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8))

        if (present(dim)) then
            values = log_softmax_8d_dims(x, [dim])
        else
            values = log_softmax_8d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8])
        end if
    end function log_softmax_8d

    pure function log_softmax_8d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_8d_dims

    pure subroutine logsumexp_signed_keepdims_complex_8d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(8)

        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_8d

    pure function logsumexp_keepdims_complex_8d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_8d(x, dims, logabs, phase)
        if (size(logabs) == 0) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_8d

    pure function logsumexp_complex_8d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values reduced over all dimensions
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_8d_all

    pure function logsumexp_complex_8d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        complex(dp), allocatable :: values(:, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :)

        if (dim < 1 .or. dim > 8) then
        allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_complex_8d_dim(x, dim, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7)))
        if (size(values) > 0) values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_8d_dim

    pure subroutine logsumexp_signed_complex_8d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the result magnitude
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_8d_all

    pure subroutine logsumexp_signed_complex_8d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: keep_shape(8)
        integer :: outshape(7)

        if (dim < 1 .or. dim > 8) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0), &
            & phase(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_phase, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_complex_8d_dim

    pure function softmax_complex_8d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8))

        if (present(dim)) then
            values = softmax_complex_8d_dims(x, [dim])
        else
            values = softmax_complex_8d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8])
        end if
    end function softmax_complex_8d

    pure function softmax_complex_8d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_complex_8d_dims

    pure function log_softmax_complex_8d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8))

        if (present(dim)) then
            values = log_softmax_complex_8d_dims(x, [dim])
        else
            values = log_softmax_complex_8d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8])
        end if
    end function log_softmax_complex_8d

    pure function log_softmax_complex_8d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! rank-8 complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_complex_8d_dims

    pure subroutine logsumexp_signed_keepdims_9d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(9)

        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_9d

    pure function logsumexp_keepdims_9d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_9d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_9d

    pure function logsumexp_9d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call low_logsumexp_signed(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_9d_all

    pure function logsumexp_9d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        real(dp), allocatable :: signs_flattened(:, :, :, :, :, :, :, :)
        logical :: valid
        integer :: keep_shape(9)
        integer :: outshape(8)

        if (dim < 1 .or. dim > 9) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(values(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8)))
        if (size(values) == 0) return
        values = reshape(flat_logabs, shape(values))
        allocate(signs_flattened(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        signs_flattened = reshape(flat_sign, shape(signs_flattened))
        where (signs_flattened < 0.0_dp .or. ieee_is_nan(signs_flattened))
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_9d_dim

    pure subroutine logsumexp_signed_9d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_9d_all

    pure subroutine logsumexp_signed_9d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: keep_shape(9)
        integer :: outshape(8)

        if (dim < 1 .or. dim > 9) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0), &
            & sign_value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_9d_dim

    pure function softmax_9d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9))

        if (present(dim)) then
            values = softmax_9d_dims(x, [dim])
        else
            values = softmax_9d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9])
        end if
    end function softmax_9d

    pure function softmax_9d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_9d_dims

    pure function log_softmax_9d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9))

        if (present(dim)) then
            values = log_softmax_9d_dims(x, [dim])
        else
            values = log_softmax_9d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9])
        end if
    end function log_softmax_9d

    pure function log_softmax_9d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_9d_dims

    pure subroutine logsumexp_signed_keepdims_complex_9d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(9)

        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_9d

    pure function logsumexp_keepdims_complex_9d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_9d(x, dims, logabs, phase)
        if (size(logabs) == 0) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_9d

    pure function logsumexp_complex_9d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values reduced over all dimensions
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_9d_all

    pure function logsumexp_complex_9d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :)

        if (dim < 1 .or. dim > 9) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_complex_9d_dim(x, dim, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8)))
        if (size(values) > 0) values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_9d_dim

    pure subroutine logsumexp_signed_complex_9d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the result magnitude
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_9d_all

    pure subroutine logsumexp_signed_complex_9d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: keep_shape(9)
        integer :: outshape(8)

        if (dim < 1 .or. dim > 9) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0), &
            & phase(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_phase, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_complex_9d_dim

    pure function softmax_complex_9d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9))

        if (present(dim)) then
            values = softmax_complex_9d_dims(x, [dim])
        else
            values = softmax_complex_9d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9])
        end if
    end function softmax_complex_9d

    pure function softmax_complex_9d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_complex_9d_dims

    pure function log_softmax_complex_9d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9))

        if (present(dim)) then
            values = log_softmax_complex_9d_dims(x, [dim])
        else
            values = log_softmax_complex_9d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9])
        end if
    end function log_softmax_complex_9d

    pure function log_softmax_complex_9d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! rank-9 complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_complex_9d_dims

    pure subroutine logsumexp_signed_keepdims_10d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(10)

        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_10d

    pure function logsumexp_keepdims_10d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_10d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_10d

    pure function logsumexp_10d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call low_logsumexp_signed(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_10d_all

    pure function logsumexp_10d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        real(dp), allocatable :: signs_flattened(:, :, :, :, :, :, :, :, :)
        logical :: valid
        integer :: keep_shape(10)
        integer :: outshape(9)

        if (dim < 1 .or. dim > 10) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(values(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9)))
        if (size(values) == 0) return
        values = reshape(flat_logabs, shape(values))
        allocate(signs_flattened(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), outshape(9)))
        signs_flattened = reshape(flat_sign, shape(signs_flattened))
        where (signs_flattened < 0.0_dp .or. ieee_is_nan(signs_flattened))
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_10d_dim

    pure subroutine logsumexp_signed_10d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_10d_all

    pure subroutine logsumexp_signed_10d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: keep_shape(10)
        integer :: outshape(9)

        if (dim < 1 .or. dim > 10) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_10d_dim

    pure function softmax_10d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10))

        if (present(dim)) then
            values = softmax_10d_dims(x, [dim])
        else
            values = softmax_10d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        end if
    end function softmax_10d

    pure function softmax_10d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_10d_dims

    pure function log_softmax_10d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10))

        if (present(dim)) then
            values = log_softmax_10d_dims(x, [dim])
        else
            values = log_softmax_10d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        end if
    end function log_softmax_10d

    pure function log_softmax_10d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_10d_dims

    pure subroutine logsumexp_signed_keepdims_complex_10d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(10)

        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_10d

    pure function logsumexp_keepdims_complex_10d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_10d(x, dims, logabs, phase)
        if (size(logabs) == 0) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_10d

    pure function logsumexp_complex_10d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values reduced over all dimensions
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_10d_all

    pure function logsumexp_complex_10d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :)

        if (dim < 1 .or. dim > 10) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_complex_10d_dim(x, dim, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9)))
        if (size(values) > 0) values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_10d_dim

    pure subroutine logsumexp_signed_complex_10d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the result magnitude
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_10d_all

    pure subroutine logsumexp_signed_complex_10d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: keep_shape(10)
        integer :: outshape(9)

        if (dim < 1 .or. dim > 10) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & phase(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_phase, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_complex_10d_dim

    pure function softmax_complex_10d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10))

        if (present(dim)) then
            values = softmax_complex_10d_dims(x, [dim])
        else
            values = softmax_complex_10d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        end if
    end function softmax_complex_10d

    pure function softmax_complex_10d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_complex_10d_dims

    pure function log_softmax_complex_10d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10))

        if (present(dim)) then
            values = log_softmax_complex_10d_dims(x, [dim])
        else
            values = log_softmax_complex_10d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        end if
    end function log_softmax_complex_10d

    pure function log_softmax_complex_10d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! rank-10 complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_complex_10d_dims

    pure subroutine logsumexp_signed_keepdims_11d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(11)

        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_11d

    pure function logsumexp_keepdims_11d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_11d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_11d

    pure function logsumexp_11d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call low_logsumexp_signed(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_11d_all

    pure function logsumexp_11d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        real(dp), allocatable :: signs_flattened(:, :, :, :, :, :, :, :, :, :)
        logical :: valid
        integer :: keep_shape(11)
        integer :: outshape(10)

        if (dim < 1 .or. dim > 11) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(values(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10)))
        if (size(values) == 0) return
        values = reshape(flat_logabs, shape(values))
        allocate(signs_flattened(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), outshape(9), outshape(10)))
        signs_flattened = reshape(flat_sign, shape(signs_flattened))
        where (signs_flattened < 0.0_dp .or. ieee_is_nan(signs_flattened))
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_11d_dim

    pure subroutine logsumexp_signed_11d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_11d_all

    pure subroutine logsumexp_signed_11d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: keep_shape(11)
        integer :: outshape(10)

        if (dim < 1 .or. dim > 11) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_11d_dim

    pure function softmax_11d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11))

        if (present(dim)) then
            values = softmax_11d_dims(x, [dim])
        else
            values = softmax_11d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        end if
    end function softmax_11d

    pure function softmax_11d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_11d_dims

    pure function log_softmax_11d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11))

        if (present(dim)) then
            values = log_softmax_11d_dims(x, [dim])
        else
            values = log_softmax_11d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        end if
    end function log_softmax_11d

    pure function log_softmax_11d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_11d_dims

    pure subroutine logsumexp_signed_keepdims_complex_11d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(11)

        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_11d

    pure function logsumexp_keepdims_complex_11d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_11d(x, dims, logabs, phase)
        if (size(logabs) == 0) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10), size(logabs, 11)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_11d

    pure function logsumexp_complex_11d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values reduced over all dimensions
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_11d_all

    pure function logsumexp_complex_11d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :)

        if (dim < 1 .or. dim > 11) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_complex_11d_dim(x, dim, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10)))
        if (size(values) > 0) values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_11d_dim

    pure subroutine logsumexp_signed_complex_11d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the result magnitude
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_11d_all

    pure subroutine logsumexp_signed_complex_11d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: keep_shape(11)
        integer :: outshape(10)

        if (dim < 1 .or. dim > 11) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_phase, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_complex_11d_dim

    pure function softmax_complex_11d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11))

        if (present(dim)) then
            values = softmax_complex_11d_dims(x, [dim])
        else
            values = softmax_complex_11d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        end if
    end function softmax_complex_11d

    pure function softmax_complex_11d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_complex_11d_dims

    pure function log_softmax_complex_11d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11))

        if (present(dim)) then
            values = log_softmax_complex_11d_dims(x, [dim])
        else
            values = log_softmax_complex_11d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        end if
    end function log_softmax_complex_11d

    pure function log_softmax_complex_11d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! rank-11 complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_complex_11d_dims

    pure subroutine logsumexp_signed_keepdims_12d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(12)

        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_12d

    pure function logsumexp_keepdims_12d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_12d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_12d

    pure function logsumexp_12d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call low_logsumexp_signed(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_12d_all

    pure function logsumexp_12d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        real(dp), allocatable :: signs_flattened(:, :, :, :, :, :, :, :, :, :, :)
        logical :: valid
        integer :: keep_shape(12)
        integer :: outshape(11)

        if (dim < 1 .or. dim > 12) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(values(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)))
        if (size(values) == 0) return
        values = reshape(flat_logabs, shape(values))
        allocate(signs_flattened(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)))
        signs_flattened = reshape(flat_sign, shape(signs_flattened))
        where (signs_flattened < 0.0_dp .or. ieee_is_nan(signs_flattened))
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_12d_dim

    pure subroutine logsumexp_signed_12d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_12d_all

    pure subroutine logsumexp_signed_12d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: keep_shape(12)
        integer :: outshape(11)

        if (dim < 1 .or. dim > 12) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_12d_dim

    pure function softmax_12d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12))

        if (present(dim)) then
            values = softmax_12d_dims(x, [dim])
        else
            values = softmax_12d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        end if
    end function softmax_12d

    pure function softmax_12d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_12d_dims

    pure function log_softmax_12d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12))

        if (present(dim)) then
            values = log_softmax_12d_dims(x, [dim])
        else
            values = log_softmax_12d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        end if
    end function log_softmax_12d

    pure function log_softmax_12d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_12d_dims

    pure subroutine logsumexp_signed_keepdims_complex_12d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(12)

        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_12d

    pure function logsumexp_keepdims_complex_12d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_12d(x, dims, logabs, phase)
        if (size(logabs) == 0) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10), size(logabs, 11), size(logabs, 12)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_12d

    pure function logsumexp_complex_12d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values reduced over all dimensions
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_12d_all

    pure function logsumexp_complex_12d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :, :)

        if (dim < 1 .or. dim > 12) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_complex_12d_dim(x, dim, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10), size(logabs, 11)))
        if (size(values) > 0) values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_12d_dim

    pure subroutine logsumexp_signed_complex_12d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the result magnitude
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_12d_all

    pure subroutine logsumexp_signed_complex_12d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: keep_shape(12)
        integer :: outshape(11)

        if (dim < 1 .or. dim > 12) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_phase, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_complex_12d_dim

    pure function softmax_complex_12d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12))

        if (present(dim)) then
            values = softmax_complex_12d_dims(x, [dim])
        else
            values = softmax_complex_12d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        end if
    end function softmax_complex_12d

    pure function softmax_complex_12d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_complex_12d_dims

    pure function log_softmax_complex_12d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12))

        if (present(dim)) then
            values = log_softmax_complex_12d_dims(x, [dim])
        else
            values = log_softmax_complex_12d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        end if
    end function log_softmax_complex_12d

    pure function log_softmax_complex_12d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! rank-12 complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_complex_12d_dims

    pure subroutine logsumexp_signed_keepdims_13d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(13)

        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_13d

    pure function logsumexp_keepdims_13d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_13d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_13d

    pure function logsumexp_13d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call low_logsumexp_signed(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_13d_all

    pure function logsumexp_13d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        real(dp), allocatable :: signs_flattened(:, :, :, :, :, :, :, :, :, :, :, :)
        logical :: valid
        integer :: keep_shape(13)
        integer :: outshape(12)

        if (dim < 1 .or. dim > 13) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(values(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12)))
        if (size(values) == 0) return
        values = reshape(flat_logabs, shape(values))
        allocate(signs_flattened(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        signs_flattened = reshape(flat_sign, shape(signs_flattened))
        where (signs_flattened < 0.0_dp .or. ieee_is_nan(signs_flattened))
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_13d_dim

    pure subroutine logsumexp_signed_13d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_13d_all

    pure subroutine logsumexp_signed_13d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: keep_shape(13)
        integer :: outshape(12)

        if (dim < 1 .or. dim > 13) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_13d_dim

    pure function softmax_13d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13))

        if (present(dim)) then
            values = softmax_13d_dims(x, [dim])
        else
            values = softmax_13d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
        end if
    end function softmax_13d

    pure function softmax_13d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_13d_dims

    pure function log_softmax_13d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13))

        if (present(dim)) then
            values = log_softmax_13d_dims(x, [dim])
        else
            values = log_softmax_13d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
        end if
    end function log_softmax_13d

    pure function log_softmax_13d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_13d_dims

    pure subroutine logsumexp_signed_keepdims_complex_13d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(13)

        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_13d

    pure function logsumexp_keepdims_complex_13d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_13d(x, dims, logabs, phase)
        if (size(logabs) == 0) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10), size(logabs, 11), size(logabs, 12), size(logabs, 13)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_13d

    pure function logsumexp_complex_13d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values reduced over all dimensions
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_13d_all

    pure function logsumexp_complex_13d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :, :, :)

        if (dim < 1 .or. dim > 13) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_complex_13d_dim(x, dim, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10), size(logabs, 11), size(logabs, 12)))
        if (size(values) > 0) values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_13d_dim

    pure subroutine logsumexp_signed_complex_13d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the result magnitude
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_13d_all

    pure subroutine logsumexp_signed_complex_13d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: keep_shape(13)
        integer :: outshape(12)

        if (dim < 1 .or. dim > 13) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_phase, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_complex_13d_dim

    pure function softmax_complex_13d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13))

        if (present(dim)) then
            values = softmax_complex_13d_dims(x, [dim])
        else
            values = softmax_complex_13d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
        end if
    end function softmax_complex_13d

    pure function softmax_complex_13d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_complex_13d_dims

    pure function log_softmax_complex_13d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13))

        if (present(dim)) then
            values = log_softmax_complex_13d_dims(x, [dim])
        else
            values = log_softmax_complex_13d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
        end if
    end function log_softmax_complex_13d

    pure function log_softmax_complex_13d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-13 complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_complex_13d_dims

    pure subroutine logsumexp_signed_keepdims_14d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(14)

        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_14d

    pure function logsumexp_keepdims_14d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_14d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_14d

    pure function logsumexp_14d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call low_logsumexp_signed(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_14d_all

    pure function logsumexp_14d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        real(dp), allocatable :: signs_flattened(:, :, :, :, :, :, :, :, :, :, :, :, :)
        logical :: valid
        integer :: keep_shape(14)
        integer :: outshape(13)

        if (dim < 1 .or. dim > 14) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(values(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13)))
        if (size(values) == 0) return
        values = reshape(flat_logabs, shape(values))
        allocate(signs_flattened(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        signs_flattened = reshape(flat_sign, shape(signs_flattened))
        where (signs_flattened < 0.0_dp .or. ieee_is_nan(signs_flattened))
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_14d_dim

    pure subroutine logsumexp_signed_14d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_14d_all

    pure subroutine logsumexp_signed_14d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: keep_shape(14)
        integer :: outshape(13)

        if (dim < 1 .or. dim > 14) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_14d_dim

    pure function softmax_14d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14))

        if (present(dim)) then
            values = softmax_14d_dims(x, [dim])
        else
            values = softmax_14d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        end if
    end function softmax_14d

    pure function softmax_14d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_14d_dims

    pure function log_softmax_14d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14))

        if (present(dim)) then
            values = log_softmax_14d_dims(x, [dim])
        else
            values = log_softmax_14d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        end if
    end function log_softmax_14d

    pure function log_softmax_14d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_14d_dims

    pure subroutine logsumexp_signed_keepdims_complex_14d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(14)

        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_14d

    pure function logsumexp_keepdims_complex_14d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_14d(x, dims, logabs, phase)
        if (size(logabs) == 0) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10), size(logabs, 11), size(logabs, 12), size(logabs, 13), size(logabs, 14)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_14d

    pure function logsumexp_complex_14d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values reduced over all dimensions
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_14d_all

    pure function logsumexp_complex_14d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :)

        if (dim < 1 .or. dim > 14) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_complex_14d_dim(x, dim, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10), size(logabs, 11), size(logabs, 12), size(logabs, 13)))
        if (size(values) > 0) values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_14d_dim

    pure subroutine logsumexp_signed_complex_14d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the result magnitude
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_14d_all

    pure subroutine logsumexp_signed_complex_14d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: keep_shape(14)
        integer :: outshape(13)

        if (dim < 1 .or. dim > 14) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_phase, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_complex_14d_dim

    pure function softmax_complex_14d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14))

        if (present(dim)) then
            values = softmax_complex_14d_dims(x, [dim])
        else
            values = softmax_complex_14d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        end if
    end function softmax_complex_14d

    pure function softmax_complex_14d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_complex_14d_dims

    pure function log_softmax_complex_14d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14))

        if (present(dim)) then
            values = log_softmax_complex_14d_dims(x, [dim])
        else
            values = log_softmax_complex_14d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        end if
    end function log_softmax_complex_14d

    pure function log_softmax_complex_14d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-14 complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_complex_14d_dims

    pure subroutine logsumexp_signed_keepdims_15d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(15)

        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14), outshape(15)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14), outshape(15)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_15d

    pure function logsumexp_keepdims_15d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_15d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_15d

    pure function logsumexp_15d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call low_logsumexp_signed(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_15d_all

    pure function logsumexp_15d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        real(dp), allocatable :: signs_flattened(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
        logical :: valid
        integer :: keep_shape(15)
        integer :: outshape(14)

        if (dim < 1 .or. dim > 15) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(values(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14)))
        if (size(values) == 0) return
        values = reshape(flat_logabs, shape(values))
        allocate(signs_flattened(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        signs_flattened = reshape(flat_sign, shape(signs_flattened))
        where (signs_flattened < 0.0_dp .or. ieee_is_nan(signs_flattened))
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_15d_dim

    pure subroutine logsumexp_signed_15d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_15d_all

    pure subroutine logsumexp_signed_15d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: keep_shape(15)
        integer :: outshape(14)

        if (dim < 1 .or. dim > 15) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_real_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_sign, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14)), &
            & sign_value(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_15d_dim

    pure function softmax_15d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14), size(x, 15))

        if (present(dim)) then
            values = softmax_15d_dims(x, [dim])
        else
            values = softmax_15d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
        end if
    end function softmax_15d

    pure function softmax_15d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14), size(x, 15))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_15d_dims

    pure function log_softmax_15d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14), size(x, 15))

        if (present(dim)) then
            values = log_softmax_15d_dims(x, [dim])
        else
            values = log_softmax_15d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
        end if
    end function log_softmax_15d

    pure function log_softmax_15d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14), size(x, 15))

        real(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_real_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_15d_dims

    pure subroutine logsumexp_signed_keepdims_complex_15d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(15)

        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), dims, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14), outshape(15)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14), outshape(15)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_15d

    pure function logsumexp_keepdims_complex_15d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_15d(x, dims, logabs, phase)
        if (size(logabs) == 0) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10), size(logabs, 11), size(logabs, 12), size(logabs, 13), size(logabs, 14), &
            & size(logabs, 15)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_15d

    pure function logsumexp_complex_15d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, &
            & :, :, :) !! rank-15 complex values reduced over all dimensions
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_15d_all

    pure function logsumexp_complex_15d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        complex(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        if (dim < 1 .or. dim > 15) then
        allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_complex_15d_dim(x, dim, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), size(logabs, 7), size(logabs, 8), size(logabs, 9), &
            & size(logabs, 10), size(logabs, 11), size(logabs, 12), size(logabs, 13), size(logabs, 14)))
        if (size(values) > 0) values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_15d_dim

    pure subroutine logsumexp_signed_complex_15d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, &
            & :, :, :) !! rank-15 complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the result magnitude
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call low_logsumexp_signed(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_15d_all

    pure subroutine logsumexp_signed_complex_15d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: keep_shape(15)
        integer :: outshape(14)

        if (dim < 1 .or. dim > 15) then
        allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), &
            & phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call drop_dimension_shape(shape(x), dim, outshape)
        call logsumexp_signed_keepdims_complex_flat(reshape(x, [size(x)]), shape(x), [dim], &
            flat_logabs, flat_phase, keep_shape, valid)
        allocate(logabs(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14)), &
            & phase(outshape(1), outshape(2), outshape(3), outshape(4), outshape(5), outshape(6), &
            & outshape(7), outshape(8), outshape(9), outshape(10), outshape(11), outshape(12), &
            & outshape(13), outshape(14)))
        if (size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_complex_15d_dim

    pure function softmax_complex_15d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 complex values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14), size(x, 15))

        if (present(dim)) then
            values = softmax_complex_15d_dims(x, [dim])
        else
            values = softmax_complex_15d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
        end if
    end function softmax_complex_15d

    pure function softmax_complex_15d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14), size(x, 15))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .false., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function softmax_complex_15d_dims

    pure function log_softmax_complex_15d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 complex values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14), size(x, 15))

        if (present(dim)) then
            values = log_softmax_complex_15d_dims(x, [dim])
        else
            values = log_softmax_complex_15d_dims(x, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
        end if
    end function log_softmax_complex_15d

    pure function log_softmax_complex_15d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! rank-15 complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6), &
            & size(x, 7), size(x, 8), size(x, 9), size(x, 10), size(x, 11), size(x, 12), size(x, 13), &
            & size(x, 14), size(x, 15))

        complex(dp), allocatable :: flat_values(:)
        logical :: valid

        call normalize_complex_flat_dims(reshape(x, [size(x)]), shape(x), dims, .true., &
            flat_values, valid)
        if (size(values) > 0) values = reshape(flat_values, shape(values))
    end function log_softmax_complex_15d_dims

end module scifort_special_reductions_high_rank
