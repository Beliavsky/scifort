! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_weighted_engine
    use scifort_descriptive_lowrank, only : gmean_low => gmean, hmean_low => hmean, &
        pmean_low => pmean
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: mean_kind_gmean = 1
    integer, parameter, public :: mean_kind_hmean = 2
    integer, parameter, public :: mean_kind_pmean = 3

    public :: reduce_assumed_rank_weighted_mean

contains

    pure subroutine reduce_assumed_rank_weighted_mean(x, weights, axis, mean_kind, &
            reduced, out_shape, broadcast_shape, ok, p, nan_policy)
        real(dp), intent(in) :: x(..) !! observations with rank 4 to 15
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast observation weights
        integer, intent(in) :: axis !! Fortran axis to reduce
        integer, intent(in) :: mean_kind !! internal selector for geometric, harmonic, or power mean
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened reduction result
        integer, intent(out) :: out_shape(:) !! result shape with reduced axis removed
        integer, intent(out) :: broadcast_shape(:) !! joint broadcast shape before reduction
        logical, intent(out) :: ok !! true when rank and axis are valid
        real(dp), intent(in), optional :: p !! finite power for the power mean
        character(len=*), intent(in), optional :: nan_policy !! optional propagate/omit/raise policy

        integer :: coordinate
        integer, allocatable :: coordinates(:)
        integer :: d
        integer, allocatable :: dims_w(:)
        integer, allocatable :: dims_x(:)
        integer :: i
        integer :: linear
        integer :: n_out
        integer :: n_reduce
        integer :: offset
        integer :: out_dim
        integer :: rank_w
        integer :: rank_x
        integer :: source_coordinate
        integer :: source_index_w
        integer :: source_index_x
        integer, allocatable :: strides_w(:)
        integer, allocatable :: strides_x(:)
        logical :: valid_extent
        real(dp), allocatable :: flat_w(:)
        real(dp), allocatable :: flat_x(:)
        real(dp), allocatable :: slice_w(:)
        real(dp), allocatable :: slice_x(:)

        call flatten_assumed_real(x, flat_x, dims_x)
        call flatten_assumed_real(weights, flat_w, dims_w)
        rank_x = size(dims_x)
        rank_w = size(dims_w)

        ok = .false.
        if (rank_x < 4 .or. rank_x > 15) then
            allocate(reduced(0))
            out_shape = 0
            broadcast_shape = 0
            return
        end if
        if (size(out_shape) /= rank_x - 1 .or. size(broadcast_shape) /= rank_x) then
            allocate(reduced(0))
            out_shape = 0
            broadcast_shape = 0
            return
        end if
        if (axis < 1 .or. axis > rank_x) then
            allocate(reduced(0))
            out_shape = 0
            broadcast_shape = 0
            return
        end if
        if (mean_kind == mean_kind_pmean .and. .not. present(p)) then
            allocate(reduced(0))
            out_shape = 0
            broadcast_shape = 0
            return
        end if

        broadcast_shape = dims_x
        if (rank_w > rank_x) then
            call set_invalid_broadcast_result(dims_x, axis, reduced, out_shape, broadcast_shape)
            ok = .true.
            return
        end if

        offset = rank_x - rank_w
        do d = 1, rank_x
            if (d <= offset) cycle
            call broadcast_extent_pair(dims_x(d), dims_w(d - offset), &
                broadcast_shape(d), valid_extent)
            if (.not. valid_extent) then
                call set_invalid_broadcast_result(dims_x, axis, reduced, out_shape, &
                    broadcast_shape)
                ok = .true.
                return
            end if
        end do

        out_dim = 0
        do d = 1, rank_x
            if (d == axis) cycle
            out_dim = out_dim + 1
            out_shape(out_dim) = broadcast_shape(d)
        end do

        n_out = product(out_shape)
        n_reduce = broadcast_shape(axis)
        allocate(reduced(n_out))
        allocate(slice_x(n_reduce))
        allocate(slice_w(n_reduce))
        allocate(coordinates(rank_x))
        allocate(strides_x(rank_x))
        call make_strides(dims_x, strides_x)
        allocate(strides_w(rank_w))
        if (rank_w > 0) call make_strides(dims_w, strides_w)

        do i = 0, n_out - 1
            coordinates = 0
            linear = i
            out_dim = 0
            do d = 1, rank_x
                if (d == axis) cycle
                out_dim = out_dim + 1
                coordinate = mod(linear, out_shape(out_dim))
                linear = linear / out_shape(out_dim)
                coordinates(d) = coordinate
            end do

            do coordinate = 0, n_reduce - 1
                coordinates(axis) = coordinate
                source_index_x = 1
                do d = 1, rank_x
                    source_coordinate = coordinates(d)
                    if (dims_x(d) == 1) source_coordinate = 0
                    source_index_x = source_index_x + source_coordinate * strides_x(d)
                end do
                slice_x(coordinate + 1) = flat_x(source_index_x)

                source_index_w = 1
                if (rank_w > 0) then
                    do d = 1, rank_w
                        source_coordinate = coordinates(offset + d)
                        if (dims_w(d) == 1) source_coordinate = 0
                        source_index_w = source_index_w + source_coordinate * strides_w(d)
                    end do
                end if
                slice_w(coordinate + 1) = flat_w(source_index_w)
            end do

            select case (mean_kind)
            case (mean_kind_gmean)
                if (present(nan_policy)) then
                    reduced(i + 1) = gmean_low(slice_x, weights=slice_w, &
                        nan_policy=nan_policy)
                else
                    reduced(i + 1) = gmean_low(slice_x, weights=slice_w)
                end if
            case (mean_kind_hmean)
                if (present(nan_policy)) then
                    reduced(i + 1) = hmean_low(slice_x, weights=slice_w, &
                        nan_policy=nan_policy)
                else
                    reduced(i + 1) = hmean_low(slice_x, weights=slice_w)
                end if
            case (mean_kind_pmean)
                if (present(nan_policy)) then
                    reduced(i + 1) = pmean_low(slice_x, p, weights=slice_w, &
                        nan_policy=nan_policy)
                else
                    reduced(i + 1) = pmean_low(slice_x, p, weights=slice_w)
                end if
            case default
                reduced = quiet_nan(0.0_dp)
            end select
        end do
        ok = .true.
    end subroutine reduce_assumed_rank_weighted_mean

    pure subroutine flatten_assumed_real(a, flat, dims)
        real(dp), intent(in) :: a(..) !! scalar or array to flatten in Fortran element order
        real(dp), allocatable, intent(out) :: flat(:) !! flattened values
        integer, allocatable, intent(out) :: dims(:) !! source extents

        select rank (source => a)
        rank (0)
            allocate(dims(0))
            allocate(flat(1))
            flat(1) = source
        rank (1)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (2)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (3)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (4)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (5)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (6)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (7)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (8)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (9)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (10)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (11)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (12)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (13)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (14)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (15)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        end select
    end subroutine flatten_assumed_real

    pure subroutine broadcast_extent_pair(left, right, extent, valid)
        integer, intent(in) :: left !! first extent
        integer, intent(in) :: right !! second extent
        integer, intent(out) :: extent !! broadcast extent
        logical, intent(out) :: valid !! true when the extents are broadcast compatible

        if (left == right) then
            extent = left
            valid = .true.
        else if (left == 1) then
            extent = right
            valid = .true.
        else if (right == 1) then
            extent = left
            valid = .true.
        else
            extent = 0
            valid = .false.
        end if
    end subroutine broadcast_extent_pair

    pure subroutine make_strides(dims, strides)
        integer, intent(in) :: dims(:) !! source extents
        integer, intent(out) :: strides(:) !! Fortran column-major strides

        integer :: d

        if (size(dims) == 0) return
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
        end do
    end subroutine make_strides

    pure subroutine set_invalid_broadcast_result(source_shape, axis, reduced, &
            out_shape, broadcast_shape)
        integer, intent(in) :: source_shape(:) !! unbroadcast observation shape
        integer, intent(in) :: axis !! valid Fortran reduction axis
        real(dp), allocatable, intent(out) :: reduced(:) !! NaN-filled fallback result
        integer, intent(out) :: out_shape(:) !! fallback reduced shape
        integer, intent(out) :: broadcast_shape(:) !! fallback retained source shape

        integer :: d
        integer :: out_dim

        broadcast_shape = source_shape
        out_dim = 0
        do d = 1, size(source_shape)
            if (d == axis) cycle
            out_dim = out_dim + 1
            out_shape(out_dim) = source_shape(d)
        end do
        allocate(reduced(product(out_shape)))
        reduced = quiet_nan(0.0_dp)
    end subroutine set_invalid_broadcast_result

end module scifort_descriptive_weighted_engine
