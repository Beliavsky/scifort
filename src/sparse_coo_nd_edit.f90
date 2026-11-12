! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Boolean-mask indexing and transactional indexed updates for arbitrary-rank
! real COO arrays. Masks are flattened in C (last axis fastest) order. COO
! coordinates remain one-based; there is no dense intermediate for scatter.
module scifort_sparse_coo_nd_edit
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, valid_coo_nd, coo_nd_sum_duplicates
    use scifort_sparse_coo_nd_index, only : coo_nd_selector, coo_nd_all, &
        coo_nd_indices, coo_nd_select
    implicit none
    private

    integer(int64), parameter :: default_nnz_limit = 1000000_int64
    public :: coo_nd_mask_select, coo_nd_mask_axis, coo_nd_scatter, coo_nd_mask_assign

contains

    pure logical function valid_flat_mask(shape, mask_size) result(ok)
        integer, intent(in) :: shape(:) !! Source dimensions (each nonnegative).
        integer(int64), intent(in) :: mask_size !! Number of flat mask elements.
        integer(int64) :: product
        integer :: ax

        ok = .false.
        if (size(shape) < 1 .or. any(shape < 0)) return
        if (any(shape == 0)) then
            ok = (mask_size == 0_int64)
            return
        end if
        product = 1_int64
        do ax = 1, size(shape)
            if (product > mask_size/int(shape(ax),int64)) return
            product = product*int(shape(ax),int64)
        end do
        ok = (product == mask_size)
    end function valid_flat_mask

    pure integer(int64) function c_order_offset(shape, coordinates) result(index)
        integer, intent(in) :: shape(:) !! Dimensions of the COO array.
        integer, intent(in) :: coordinates(:) !! Valid one-based coordinates.
        integer :: ax

        index = int(coordinates(1),int64)
        do ax = 2, size(shape)
            index = (index-1_int64)*int(shape(ax),int64) + int(coordinates(ax),int64)
        end do
    end function c_order_offset

    function coo_nd_mask_select(a, mask, status, max_nnz) result(b)
        type(coo_nd_array), intent(in) :: a !! Source N-D COO array, including duplicates.
        logical, intent(in) :: mask(:) !! Full boolean mask, flattened in C order.
        integer, intent(out), optional :: status !! Zero on success; invalid input otherwise.
        integer(int64), intent(in), optional :: max_nnz !! Maximum stored output entries (default 1e6).
        type(coo_nd_array) :: b
        integer, allocatable :: position(:)
        integer(int64) :: limit, index, entries
        integer :: i, nselected, nout, k

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (.not. valid_flat_mask(a%shape,size(mask,kind=int64))) return
        limit = default_nnz_limit
        if (present(max_nnz)) limit = max_nnz
        if (limit < 0_int64) return
        limit = min(limit,int(huge(0),int64))

        allocate(position(size(mask)))
        nselected = 0
        do i = 1, size(mask)
            if (mask(i)) nselected = nselected+1
            position(i) = nselected
        end do
        entries = 0_int64
        do k = 1, size(a%data)
            index = c_order_offset(a%shape,a%coords(:,k))
            if (.not. mask(int(index))) cycle
            if (entries == limit) return
            entries = entries+1_int64
        end do
        b%shape = [nselected]
        allocate(b%coords(1,int(entries)),b%data(int(entries)))
        nout = 0
        do k = 1, size(a%data)
            index = c_order_offset(a%shape,a%coords(:,k))
            if (.not. mask(int(index))) cycle
            nout = nout+1
            b%coords(1,nout) = position(int(index))
            b%data(nout) = a%data(k)
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_mask_select

    function coo_nd_mask_axis(a, axis, mask, status, max_nnz) result(b)
        type(coo_nd_array), intent(in) :: a !! Source N-D COO array.
        integer, intent(in) :: axis !! One-based axis to compress by the mask.
        logical, intent(in) :: mask(:) !! Mask of exactly shape(axis) elements.
        integer, intent(out), optional :: status !! Zero on success.
        integer(int64), intent(in), optional :: max_nnz !! Maximum stored output entries.
        type(coo_nd_array) :: b
        type(coo_nd_selector), allocatable :: selectors(:)
        integer, allocatable :: selected(:)
        integer :: ax

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (axis < 1 .or. axis > size(a%shape)) return
        if (size(mask) /= a%shape(axis)) return
        allocate(selectors(size(a%shape)))
        do ax = 1, size(selectors)
            selectors(ax) = coo_nd_all()
        end do
        selected = pack([(ax,ax=1,size(mask))],mask)
        selectors(axis) = coo_nd_indices(selected)
        b = coo_nd_select(a,selectors,status,max_nnz=max_nnz)
    end function coo_nd_mask_axis

    pure logical function coordinates_less(coords, left, right) result(less)
        integer, intent(in) :: coords(:,:) !! Coordinates of update events.
        integer, intent(in) :: left !! First event index.
        integer, intent(in) :: right !! Second event index.
        integer :: ax

        less = .false.
        do ax = 1, size(coords,1)
            if (coords(ax,left) < coords(ax,right)) then
                less = .true.
                return
            end if
            if (coords(ax,left) > coords(ax,right)) return
        end do
    end function coordinates_less

    subroutine stable_coordinate_order(coords, order, scratch)
        integer, intent(in) :: coords(:,:) !! Coordinate keys, (ndim,nevents).
        integer, intent(out) :: order(:) !! Sorted permutation, stable for equal keys.
        integer, intent(out) :: scratch(:) !! Reusable sorting workspace.
        integer :: n, width, left, mid, right, i, j, k

        n = size(order)
        do i = 1, n
            order(i) = i
        end do
        width = 1
        do while (width < n)
            left = 1
            do while (left <= n)
                mid = left + min(width,n-left+1) - 1
                right = mid + min(width,n-mid)
                i = left
                j = mid+1
                do k = left, right
                    if (i > mid) then
                        scratch(k) = order(j)
                        j = j+1
                    else if (j > right) then
                        scratch(k) = order(i)
                        i = i+1
                    else if (coordinates_less(coords,order(j),order(i))) then
                        scratch(k) = order(j)
                        j = j+1
                    else
                        scratch(k) = order(i)
                        i = i+1
                    end if
                end do
                if (right == n) exit
                left = right+1
            end do
            order = scratch
            if (width >= (n+1)/2) exit
            width = width*2
        end do
    end subroutine stable_coordinate_order

    subroutine coo_nd_scatter(a, points, values, status, accumulate, keep_zeros, max_nnz)
        type(coo_nd_array), intent(inout) :: a !! Updated only on success; unchanged on failure.
        integer, intent(in) :: points(:,:) !! One-based (ndim,nupdates) update coordinates.
        real(dp), intent(in) :: values(:) !! One value per coordinate; finite values only.
        integer, intent(out), optional :: status !! Zero on success; invalid input otherwise.
        logical, intent(in), optional :: accumulate !! Add to existing values; default overwrite.
        logical, intent(in), optional :: keep_zeros !! Retain touched zeros; default false.
        integer(int64), intent(in), optional :: max_nnz !! Max resulting stored entries, default 1e6.
        type(coo_nd_array) :: original
        integer, allocatable :: event_coords(:,:), output_coords(:,:), order(:), scratch(:)
        real(dp), allocatable :: output_values(:)
        integer(int64) :: n_events64, limit
        integer :: ndim, n_original, n_updates, n_events, j, i, k, finish, nout, item
        real(dp) :: value
        logical :: is_add, retain_zero, touched

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        ndim = size(a%shape)
        n_updates = size(values)
        if (size(points,1) /= ndim .or. size(points,2) /= n_updates) return
        if (.not. all(ieee_is_finite(values))) return
        do i = 1, ndim
            if (any(points(i,:) < 1) .or. any(points(i,:) > a%shape(i))) return
        end do
        limit = default_nnz_limit
        if (present(max_nnz)) limit = max_nnz
        if (limit < 0_int64) return
        limit = min(limit,int(huge(0),int64))
        if (n_updates == 0) then
            if (present(status)) status = sparse_success
            return
        end if
        is_add = .false.
        if (present(accumulate)) is_add = accumulate
        retain_zero = .false.
        if (present(keep_zeros)) retain_zero = keep_zeros

        ! Coalesce the source in a temporary so invalid arithmetic never
        ! leaves a partly updated input. All subsequent work is transactional.
        original = a
        call coo_nd_sum_duplicates(original,j)
        if (j /= sparse_success) return
        n_original = size(original%data)
        n_events64 = int(n_original,int64) + int(n_updates,int64)
        if (n_events64 > int(huge(0),int64)) return
        n_events = int(n_events64)
        allocate(event_coords(ndim,n_events),order(n_events),scratch(n_events))
        event_coords(:,1:n_original) = original%coords
        event_coords(:,n_original+1:) = points
        call stable_coordinate_order(event_coords,order,scratch)
        allocate(output_coords(ndim,n_events),output_values(n_events))
        nout = 0
        k = 1
        do while (k <= n_events)
            finish = k
            value = 0.0_dp
            touched = .false.
            do
                item = order(finish)
                if (item <= n_original) then
                    value = original%data(item)
                else
                    touched = .true.
                    if (is_add) then
                        value = value + values(item-n_original)
                    else
                        value = values(item-n_original)
                    end if
                    if (.not. ieee_is_finite(value)) return
                end if
                if (finish == n_events) exit
                if (.not. all(event_coords(:,order(finish+1)) == event_coords(:,order(k)))) exit
                finish = finish+1
            end do
            if (.not. touched .or. retain_zero .or. value /= 0.0_dp) then
                if (int(nout,int64) == limit) return
                nout = nout+1
                output_coords(:,nout) = event_coords(:,order(k))
                output_values(nout) = value
            end if
            k = finish+1
        end do
        a%coords = output_coords(:,1:nout)
        a%data = output_values(1:nout)
        if (present(status)) status = sparse_success
    end subroutine coo_nd_scatter

    subroutine coo_nd_mask_assign(a, mask, value, status, accumulate, keep_zeros, max_nnz)
        type(coo_nd_array), intent(inout) :: a !! COO updated only on success.
        logical, intent(in) :: mask(:) !! Full boolean mask flattened in C order.
        real(dp), intent(in) :: value !! Scalar assigned or added to selected elements.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: accumulate !! Add scalar rather than overwrite.
        logical, intent(in), optional :: keep_zeros !! Retain touched explicit zeros.
        integer(int64), intent(in), optional :: max_nnz !! Bound on updates and resulting entries.
        integer, allocatable :: points(:,:)
        real(dp), allocatable :: values(:)
        integer(int64) :: limit, position
        integer :: nselected, index, ax, i, k

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (.not. valid_flat_mask(a%shape,size(mask,kind=int64))) return
        if (.not. ieee_is_finite(value)) return
        limit = default_nnz_limit
        if (present(max_nnz)) limit = max_nnz
        if (limit < 0_int64) return
        limit = min(limit,int(huge(0),int64))
        if (count(mask,kind=int64) > limit) return
        nselected = count(mask)
        if (nselected == 0) then
            if (present(status)) status = sparse_success
            return
        end if
        allocate(points(size(a%shape),nselected),values(nselected))
        values = value
        k = 0
        do index = 1, size(mask)
            if (.not. mask(index)) cycle
            k = k+1
            position = int(index,int64)-1_int64
            do ax = size(a%shape), 1, -1
                i = a%shape(ax)
                points(ax,k) = int(mod(position,int(i,int64)))+1
                position = position/int(i,int64)
            end do
        end do
        call coo_nd_scatter(a,points,values,status,accumulate,keep_zeros,max_nnz)
    end subroutine coo_nd_mask_assign

end module scifort_sparse_coo_nd_edit
