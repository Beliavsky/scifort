! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Indexing and slicing for one-based N-dimensional COO sparse arrays.
! Axis-by-axis integer arrays use orthogonal (Cartesian) indexing; pointwise
! advanced indexing is available separately through coo_nd_gather.
! Unlike dense array slicing, no dense shape-product storage is allocated.
module scifort_sparse_coo_nd_index
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, valid_coo_nd
    implicit none
    private

    integer, parameter :: selector_all = 0
    integer, parameter :: selector_scalar = 1
    integer, parameter :: selector_slice = 2
    integer, parameter :: selector_array = 3
    integer(int64), parameter :: default_nnz_limit = 1000000_int64

    public :: coo_nd_selector, coo_nd_all, coo_nd_at, coo_nd_slice, coo_nd_indices
    public :: coo_nd_select, coo_nd_gather

    type :: coo_nd_selector
        private
        integer :: mode = selector_all
        integer :: first = 1
        integer :: last = 1
        integer :: stride = 1
        logical :: has_first = .false.
        logical :: has_last = .false.
        integer, allocatable :: indices(:)
    end type coo_nd_selector

    type :: selected_axis
        integer :: extent = 0
        integer :: first = 1
        integer :: last = 0
        integer :: stride = 1
    end type selected_axis

    type :: coordinate_matches
        integer, allocatable :: positions(:)
    end type coordinate_matches

contains

    pure function coo_nd_all() result(selection)
        type(coo_nd_selector) :: selection
        ! Default-initialized to the entire axis.
    end function coo_nd_all

    pure function coo_nd_at(index) result(selection)
        integer, intent(in) :: index !! One-based integer coordinate.
        type(coo_nd_selector) :: selection
        selection%mode = selector_scalar
        selection%first = index
    end function coo_nd_at

    pure function coo_nd_slice(first, last, stride) result(selection)
        integer, intent(in), optional :: first !! Inclusive one-based start; defaults to first/last element.
        integer, intent(in), optional :: last !! Inclusive one-based end; defaults to last/first element.
        integer, intent(in), optional :: stride !! Nonzero step; default 1, negative steps reverse the axis.
        type(coo_nd_selector) :: selection

        selection%mode = selector_slice
        if (present(first)) then
            selection%has_first = .true.
            selection%first = first
        end if
        if (present(last)) then
            selection%has_last = .true.
            selection%last = last
        end if
        if (present(stride)) selection%stride = stride
    end function coo_nd_slice

    pure function coo_nd_indices(indices) result(selection)
        integer, intent(in) :: indices(:) !! One-based selected coordinates, repeated values allowed.
        type(coo_nd_selector) :: selection

        selection%mode = selector_array
        selection%indices = indices
    end function coo_nd_indices

    pure subroutine prepare_axis(selection, extent, plan, ok)
        type(coo_nd_selector), intent(in) :: selection !! Index selector for a source axis.
        integer, intent(in) :: extent !! Size of the source axis.
        type(selected_axis), intent(out) :: plan !! Derived length and slice bounds.
        logical, intent(out) :: ok !! True when the selector is valid.
        integer(int64) :: first64, last64, stride64, n

        ok = .false.
        if (extent < 0) return
        select case (selection%mode)
        case (selector_all)
            plan%extent = extent
        case (selector_scalar)
            if (selection%first < 1 .or. selection%first > extent) return
            plan%extent = 1
            plan%first = selection%first
        case (selector_array)
            if (.not. allocated(selection%indices)) return
            if (any(selection%indices < 1) .or. any(selection%indices > extent)) return
            plan%extent = size(selection%indices)
        case (selector_slice)
            stride64 = int(selection%stride,int64)
            if (stride64 == 0_int64) return
            first64 = 1_int64
            last64 = int(extent,int64)
            if (stride64 < 0_int64) then
                first64 = int(extent,int64)
                last64 = 1_int64
            end if
            if (selection%has_first) first64 = int(selection%first,int64)
            if (selection%has_last) last64 = int(selection%last,int64)
            if (extent == 0 .and. .not. selection%has_first .and. .not. selection%has_last) then
                plan%extent = 0
                ok = .true.
                return
            end if
            if (first64 < 1_int64 .or. first64 > int(extent,int64)) return
            if (last64 < 1_int64 .or. last64 > int(extent,int64)) return
            plan%first = int(first64)
            plan%last = int(last64)
            plan%stride = selection%stride
            n = 0_int64
            if ((stride64 > 0 .and. first64 <= last64) .or. &
                (stride64 < 0 .and. first64 >= last64)) then
                n = 1_int64 + abs(last64-first64)/abs(stride64)
            end if
            if (n > int(huge(0),int64)) return
            plan%extent = int(n)
        case default
            return
        end select
        ok = .true.
    end subroutine prepare_axis

    subroutine locate_positions(selection, plan, coordinate, matches)
        type(coo_nd_selector), intent(in) :: selection !! Prepared axis selection.
        type(selected_axis), intent(in) :: plan !! Prepared axis information.
        integer, intent(in) :: coordinate !! One-based source coordinate.
        type(coordinate_matches), intent(inout) :: matches !! Matching destination coordinates.
        integer(int64) :: delta
        integer :: n, j, found

        select case (selection%mode)
        case (selector_all)
            allocate(matches%positions(1))
            matches%positions(1) = coordinate
        case (selector_scalar)
            n = 0
            if (coordinate == plan%first) n = 1
            allocate(matches%positions(n))
            if (n == 1) matches%positions(1) = 1
        case (selector_slice)
            n = 0
            found = 0
            if (plan%extent > 0) then
                delta = int(coordinate,int64) - int(plan%first,int64)
                if (mod(delta,int(plan%stride,int64)) == 0_int64) then
                    delta = delta/int(plan%stride,int64)
                    if (delta >= 0_int64 .and. delta < int(plan%extent,int64)) then
                        n = 1
                        found = int(delta)+1
                    end if
                end if
            end if
            allocate(matches%positions(n))
            if (n == 1) matches%positions(1) = found
        case (selector_array)
            n = count(selection%indices == coordinate)
            allocate(matches%positions(n))
            found = 0
            do j = 1, size(selection%indices)
                if (selection%indices(j) /= coordinate) cycle
                found = found+1
                matches%positions(found) = j
            end do
        end select
    end subroutine locate_positions

    subroutine entry_matches(selectors, plans, coordinate, matches, count_out, limit, ok)
        type(coo_nd_selector), intent(in) :: selectors(:) !! Axis selectors, one per dimension.
        type(selected_axis), intent(in) :: plans(:) !! Prepared selectors.
        integer, intent(in) :: coordinate(:) !! Original entry coordinates.
        type(coordinate_matches), intent(inout) :: matches(:) !! Per-axis destination matches.
        integer(int64), intent(out) :: count_out !! Replication count, bounded by limit.
        integer(int64), intent(in) :: limit !! Maximum generated entries.
        logical, intent(out) :: ok !! False when the count exceeds limit.
        integer :: ax, n

        count_out = 1_int64
        ok = .true.
        do ax = 1, size(selectors)
            if (allocated(matches(ax)%positions)) deallocate(matches(ax)%positions)
            call locate_positions(selectors(ax),plans(ax),coordinate(ax),matches(ax))
            n = size(matches(ax)%positions)
            if (n == 0) then
                count_out = 0_int64
                exit
            end if
            if (count_out > limit/int(n,int64)) then
                ok = .false.
                return
            end if
            count_out = count_out*int(n,int64)
        end do
    end subroutine entry_matches

    function coo_nd_select(a, selectors, status, keepdims, max_nnz) result(b)
        type(coo_nd_array), intent(in) :: a !! Source N-dimensional COO array.
        type(coo_nd_selector), intent(in) :: selectors(:) !! One selector per source axis (orthogonal indexing).
        integer, intent(out), optional :: status !! Zero on success, invalid input otherwise.
        logical, intent(in), optional :: keepdims !! Keep integer-indexed axes as size one (default false).
        integer(int64), intent(in), optional :: max_nnz !! Maximum output stored entries, default one million.
        type(coo_nd_array) :: b
        type(selected_axis), allocatable :: plans(:)
        type(coordinate_matches), allocatable :: matches(:)
        integer, allocatable :: axis_map(:), digits(:), result_dims(:)
        integer(int64) :: cap, copies, total
        integer :: ndim, nr, ax, k, i, j, out_index, output_axis
        logical :: valid, preserve_axes

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        ndim = size(a%shape)
        if (size(selectors) /= ndim) return
        cap = default_nnz_limit
        if (present(max_nnz)) cap = max_nnz
        if (cap < 0_int64) return
        cap = min(cap,int(huge(0),int64))
        preserve_axes = .false.
        if (present(keepdims)) preserve_axes = keepdims

        allocate(plans(ndim),matches(ndim),axis_map(ndim),digits(ndim))
        nr = 0
        do ax = 1, ndim
            call prepare_axis(selectors(ax),a%shape(ax),plans(ax),valid)
            if (.not. valid) return
            axis_map(ax) = 0
            if (.not. preserve_axes .and. selectors(ax)%mode == selector_scalar) cycle
            nr = nr+1
            axis_map(ax) = nr
        end do
        ! The base COO representation deliberately has no rank-zero type.
        ! Return a rank-one singleton when every axis is indexed by an integer.
        allocate(result_dims(max(nr,1)))
        result_dims = 1
        do ax = 1, ndim
            if (axis_map(ax) > 0) result_dims(axis_map(ax)) = plans(ax)%extent
        end do
        total = 0_int64
        if (all(result_dims > 0)) then
            do k = 1, size(a%data)
                call entry_matches(selectors,plans,a%coords(:,k),matches,copies,cap,valid)
                if (.not. valid) return
                if (copies > cap-total) return
                total = total+copies
            end do
        end if
        b%shape = result_dims
        allocate(b%coords(size(result_dims),int(total)),b%data(int(total)))
        out_index = 0
        do k = 1, size(a%data)
            if (out_index == int(total)) exit
            call entry_matches(selectors,plans,a%coords(:,k),matches,copies,cap,valid)
            if (.not. valid) error stop 'internal selection count failure'
            if (copies == 0_int64) cycle
            digits = 1
            do i = 1, int(copies)
                out_index = out_index+1
                b%coords(:,out_index) = 1
                do ax = 1, ndim
                    output_axis = axis_map(ax)
                    if (output_axis /= 0) then
                        b%coords(output_axis,out_index) = matches(ax)%positions(digits(ax))
                    end if
                end do
                b%data(out_index) = a%data(k)
                ! Rightmost selected axis varies fastest, as in C-order indexing.
                do j = ndim,1,-1
                    digits(j) = digits(j)+1
                    if (digits(j) <= size(matches(j)%positions)) exit
                    digits(j) = 1
                end do
            end do
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_select

    function coo_nd_gather(a, points, status, max_nnz) result(b)
        type(coo_nd_array), intent(in) :: a !! Source COO array, unsorted entries and duplicates allowed.
        integer, intent(in) :: points(:,:) !! One-based (ndim,npoints) point coordinates.
        integer, intent(out), optional :: status !! Zero on success.
        integer(int64), intent(in), optional :: max_nnz !! Maximum stored entries generated, default one million.
        type(coo_nd_array) :: b
        integer(int64) :: cap, total
        integer :: npoints, k, p, index_out

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (size(points,1) /= size(a%shape)) return
        npoints = size(points,2)
        do k = 1, size(a%shape)
            if (any(points(k,:) < 1) .or. any(points(k,:) > a%shape(k))) return
        end do
        cap = default_nnz_limit
        if (present(max_nnz)) cap = max_nnz
        if (cap < 0_int64) return
        cap = min(cap,int(huge(0),int64))
        total = 0_int64
        do p = 1, npoints
            do k = 1, size(a%data)
                if (.not. all(a%coords(:,k) == points(:,p))) cycle
                if (total == cap) return
                total = total+1_int64
            end do
        end do
        b%shape = [npoints]
        allocate(b%coords(1,int(total)),b%data(int(total)))
        index_out = 0
        do p = 1, npoints
            do k = 1, size(a%data)
                if (.not. all(a%coords(:,k) == points(:,p))) cycle
                index_out = index_out+1
                b%coords(1,index_out) = p
                b%data(index_out) = a%data(k)
            end do
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_gather

end module scifort_sparse_coo_nd_index
