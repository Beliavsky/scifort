! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! N-dimensional coordinate sparse arrays, modeled after scipy.sparse.coo_array.
! Indices are one-based; coords(axis, entry) stores the coordinates of
! each entry. Duplicates and explicitly stored zeros are preserved until
! explicitly coalesced or eliminated. The existing two-dimensional COO,
! CSR and CSC interfaces remain unchanged.
module scifort_sparse_coo_nd
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo, only : coo_matrix, valid_coo, coo_from_triplets
    implicit none
    private

    public :: coo_nd_array, valid_coo_nd, coo_nd_from_triplets
    public :: coo_nd_from_coo, coo_nd_to_coo, coo_nd_get
    public :: coo_nd_sum_duplicates, coo_nd_eliminate_zeros, coo_nd_is_canonical
    public :: coo_nd_transpose, coo_nd_reshape, coo_nd_to_dense_flat
    public :: coo_nd_from_dense_flat, coo_nd_from_dense, coo_nd_to_dense

    ! The allocated components are deliberately public, just as the row,
    ! column and data arrays of scifort_sparse_coo%coo_matrix are public.
    type :: coo_nd_array
        integer, allocatable :: shape(:)  !! Length ndim; nonnegative extents.
        integer, allocatable :: coords(:,:) !! Shape (ndim,nnz); one-based coordinates.
        real(dp), allocatable :: data(:) !! Shape (nnz); duplicates/zeros allowed.
    end type coo_nd_array

contains

    pure logical function valid_coo_nd(a) result(ok)
        type(coo_nd_array), intent(in) :: a !! Candidate N-dimensional COO array.
        integer :: ax

        ok = .false.
        if (.not. allocated(a%shape)) return
        if (.not. allocated(a%coords)) return
        if (.not. allocated(a%data)) return
        if (size(a%shape) < 1) return
        if (any(a%shape < 0)) return
        if (size(a%data) == huge(0)) return
        if (size(a%coords,1) /= size(a%shape)) return
        if (size(a%coords,2) /= size(a%data)) return
        do ax = 1, size(a%shape)
            if (any(a%coords(ax,:) < 1)) return
            if (any(a%coords(ax,:) > a%shape(ax))) return
        end do
        if (.not. all(ieee_is_finite(a%data))) return
        ok = .true.
    end function valid_coo_nd

    function coo_nd_from_triplets(dims, coordinates, values, status) result(a)
        integer, intent(in) :: dims(:) !! Dimensions; zero extents allowed, rank >= 1.
        integer, intent(in) :: coordinates(:,:) !! One-based (ndim,nnz) coordinates.
        real(dp), intent(in) :: values(:) !! Finite values, including explicit zeros.
        integer, intent(out), optional :: status !! Zero on success, invalid input otherwise.
        type(coo_nd_array) :: a

        if (present(status)) status = sparse_invalid_input
        if (size(dims) < 1) return
        if (any(dims < 0)) return
        if (size(values) == huge(0)) return
        if (size(coordinates,1) /= size(dims)) return
        if (size(coordinates,2) /= size(values)) return
        if (.not. all(ieee_is_finite(values))) return
        block
            integer :: ax
            do ax = 1, size(dims)
                if (any(coordinates(ax,:) < 1)) return
                if (any(coordinates(ax,:) > dims(ax))) return
            end do
        end block
        a%shape = dims
        a%coords = coordinates
        a%data = values
        if (present(status)) status = sparse_success
    end function coo_nd_from_triplets

    function coo_nd_from_coo(a, status) result(b)
        type(coo_matrix), intent(in) :: a !! Existing two-dimensional COO input.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_nd_array) :: b
        integer, allocatable :: coordinates(:,:)
        integer :: n

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo(a)) return
        n = size(a%data)
        allocate(coordinates(2,n))
        coordinates(1,:) = a%row
        coordinates(2,:) = a%col
        b = coo_nd_from_triplets([a%n_rows,a%n_cols],coordinates,a%data,status)
    end function coo_nd_from_coo

    function coo_nd_to_coo(a, status) result(b)
        type(coo_nd_array), intent(in) :: a !! N-dimensional COO, must have rank two.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_matrix) :: b

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (size(a%shape) /= 2) return
        b = coo_from_triplets(a%shape(1),a%shape(2),a%coords(1,:), &
            a%coords(2,:),a%data,status)
    end function coo_nd_to_coo

    real(dp) function coo_nd_get(a, indices, status) result(value)
        type(coo_nd_array), intent(in) :: a !! N-dimensional COO input.
        integer, intent(in) :: indices(:) !! One-based coordinate of queried element.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: k

        value = 0.0_dp
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (size(indices) /= size(a%shape)) return
        if (any(indices < 1) .or. any(indices > a%shape)) return
        do k = 1, size(a%data)
            if (.not. all(a%coords(:,k) == indices)) cycle
            value = value + a%data(k)
            if (.not. ieee_is_finite(value)) then
                value = 0.0_dp
                return
            end if
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_get

    pure logical function coords_less(coords, i, j) result(less)
        integer, intent(in) :: coords(:,:) !! Coordinate storage (ndim,nnz).
        integer, intent(in) :: i !! First entry index.
        integer, intent(in) :: j !! Second entry index.
        integer :: ax

        less = .false.
        do ax = 1, size(coords,1)
            if (coords(ax,i) < coords(ax,j)) then
                less = .true.
                return
            end if
            if (coords(ax,i) > coords(ax,j)) return
        end do
    end function coords_less

    pure logical function coo_nd_is_canonical(a) result(ok)
        type(coo_nd_array), intent(in) :: a !! COO array to check.
        integer :: k

        ok = .false.
        if (.not. valid_coo_nd(a)) return
        do k = 2, size(a%data)
            if (.not. coords_less(a%coords,k-1,k)) return
        end do
        ok = .true.
    end function coo_nd_is_canonical

    subroutine coo_nd_sum_duplicates(a, status)
        type(coo_nd_array), intent(inout) :: a !! Array to canonicalize in C-order.
        integer, intent(out), optional :: status !! Zero on success; unchanged on failure.
        integer, allocatable :: index(:), scratch(:), canon_coords(:,:)
        real(dp), allocatable :: canon_data(:)
        integer :: n, ndim, width, left, mid, right, i, j, k, count_unique

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        n = size(a%data)
        ndim = size(a%shape)
        allocate(index(n),scratch(n),canon_coords(ndim,n),canon_data(n))
        do i = 1, n
            index(i) = i
        end do
        ! Stable, bottom-up mergesort, avoiding worst-case quadratic sorting.
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
                        scratch(k) = index(j)
                        j = j+1
                    else if (j > right) then
                        scratch(k) = index(i)
                        i = i+1
                    else if (coords_less(a%coords,index(j),index(i))) then
                        scratch(k) = index(j)
                        j = j+1
                    else
                        scratch(k) = index(i)
                        i = i+1
                    end if
                end do
                if (right == n) exit
                left = right+1
            end do
            index = scratch
            if (width >= (n+1)/2) exit
            width = width*2
        end do
        count_unique = 0
        do k = 1, n
            i = index(k)
            if (count_unique > 0) then
                if (all(canon_coords(:,count_unique) == a%coords(:,i))) then
                    canon_data(count_unique) = canon_data(count_unique) + a%data(i)
                    if (.not. ieee_is_finite(canon_data(count_unique))) return
                    cycle
                end if
            end if
            count_unique = count_unique+1
            canon_coords(:,count_unique) = a%coords(:,i)
            canon_data(count_unique) = a%data(i)
        end do
        a%coords = canon_coords(:,1:count_unique)
        a%data = canon_data(1:count_unique)
        if (present(status)) status = sparse_success
    end subroutine coo_nd_sum_duplicates

    subroutine coo_nd_eliminate_zeros(a, status)
        type(coo_nd_array), intent(inout) :: a !! COO array; duplicate positions retained.
        integer, intent(out), optional :: status !! Zero on success.
        logical, allocatable :: keep(:)
        integer, allocatable :: remaining_coords(:,:)
        integer :: k, j, remaining

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        keep = a%data /= 0.0_dp
        remaining = count(keep)
        allocate(remaining_coords(size(a%shape),remaining))
        j = 0
        do k = 1, size(a%data)
            if (.not. keep(k)) cycle
            j = j+1
            remaining_coords(:,j) = a%coords(:,k)
        end do
        a%coords = remaining_coords
        a%data = pack(a%data,keep)
        if (present(status)) status = sparse_success
    end subroutine coo_nd_eliminate_zeros

    function coo_nd_transpose(a, axes, status) result(b)
        type(coo_nd_array), intent(in) :: a !! COO input, of any positive rank.
        integer, intent(in), optional :: axes(:) !! One-based output axis order; default reversal.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_nd_array) :: b
        integer, allocatable :: permutation(:)
        integer :: ax, rank_a

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        rank_a = size(a%shape)
        allocate(permutation(rank_a))
        if (present(axes)) then
            if (size(axes) /= rank_a) return
            permutation = axes
        else
            permutation = [(rank_a+1-ax,ax=1,rank_a)]
        end if
        if (any(permutation < 1) .or. any(permutation > rank_a)) return
        do ax = 1, rank_a
            if (count(permutation == ax) /= 1) return
        end do
        b%shape = a%shape(permutation)
        b%coords = a%coords(permutation,:)
        b%data = a%data
        if (present(status)) status = sparse_success
    end function coo_nd_transpose

    ! Compute total number of cells without overflowing int64. It is valid
    ! to have enormous *sparse* arrays; this bound applies only to flat indices.
    pure subroutine extent_product(dims, total, ok)
        integer, intent(in) :: dims(:) !! Candidate nonnegative shape.
        integer(int64), intent(out) :: total !! Total elements, if representable.
        logical, intent(out) :: ok !! Success flag.
        integer :: ax

        total = 0_int64
        ok = .false.
        if (size(dims) < 1) return
        if (any(dims < 0)) return
        if (any(dims == 0)) then
            ok = .true.
            return
        end if
        total = 1_int64
        do ax = 1, size(dims)
            if (total > huge(total)/int(dims(ax),int64)) return
            total = total*int(dims(ax),int64)
        end do
        ok = .true.
    end subroutine extent_product

    pure logical function valid_order(order) result(ok)
        character(len=1), intent(in) :: order !! Storage order 'C' or 'F'.
        ok = (order == 'C' .or. order == 'c' .or. &
              order == 'F' .or. order == 'f')
    end function valid_order

    pure integer(int64) function flatten_coords(coordinates, dims, order) result(offset)
        integer, intent(in) :: coordinates(:) !! One-based element coordinate.
        integer, intent(in) :: dims(:) !! Extents whose product fits int64.
        character(len=1), intent(in) :: order !! 'C' or 'F', already validated.
        integer(int64) :: stride
        integer :: ax

        offset = 0_int64
        if (order == 'C' .or. order == 'c') then
            do ax = 1, size(dims)
                offset = offset*int(dims(ax),int64) + int(coordinates(ax)-1,int64)
            end do
        else
            stride = 1_int64
            do ax = 1, size(dims)
                offset = offset + stride*int(coordinates(ax)-1,int64)
                stride = stride*int(dims(ax),int64)
            end do
        end if
    end function flatten_coords

    pure subroutine unflatten_coords(offset, dims, order, coordinates)
        integer(int64), intent(in) :: offset !! Zero-based flattened index.
        integer, intent(in) :: dims(:) !! Nonzero extents whose product fits int64.
        character(len=1), intent(in) :: order !! 'C' or 'F', already validated.
        integer, intent(out) :: coordinates(:) !! Result, one-based coordinates.
        integer(int64) :: work
        integer :: ax

        work = offset
        if (order == 'C' .or. order == 'c') then
            do ax = size(dims),1,-1
                coordinates(ax) = int(mod(work,int(dims(ax),int64)))+1
                work = work/int(dims(ax),int64)
            end do
        else
            do ax = 1,size(dims)
                coordinates(ax) = int(mod(work,int(dims(ax),int64)))+1
                work = work/int(dims(ax),int64)
            end do
        end if
    end subroutine unflatten_coords

    function coo_nd_reshape(a, new_shape, order, status) result(b)
        type(coo_nd_array), intent(in) :: a !! COO input, duplicates and zeros retained.
        integer, intent(in) :: new_shape(:) !! Nonnegative new shape, equal element count.
        character(len=1), intent(in), optional :: order !! 'C' (default) or 'F'.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_nd_array) :: b
        integer(int64) :: nold, nnew, offset
        logical :: old_ok, new_ok
        integer :: k
        character(len=1) :: layout

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        layout = 'C'
        if (present(order)) layout = order
        if (.not. valid_order(layout)) return
        call extent_product(a%shape,nold,old_ok)
        call extent_product(new_shape,nnew,new_ok)
        if (.not. old_ok .or. .not. new_ok) return
        if (nold /= nnew) return
        b%shape = new_shape
        allocate(b%coords(size(new_shape),size(a%data)))
        b%data = a%data
        do k = 1, size(a%data)
            offset = flatten_coords(a%coords(:,k),a%shape,layout)
            call unflatten_coords(offset,new_shape,layout,b%coords(:,k))
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_reshape

    function coo_nd_to_dense_flat(a, status, order, max_elements) result(flat)
        type(coo_nd_array), intent(in) :: a !! COO input with possible duplicates.
        integer, intent(out), optional :: status !! Zero on success.
        character(len=1), intent(in), optional :: order !! 'C' default, or 'F'.
        integer(int64), intent(in), optional :: max_elements !! Dense allocation cap; default 10 million.
        real(dp), allocatable :: flat(:)
        integer(int64) :: total, location, limit
        logical :: ok
        character(len=1) :: layout
        integer :: k

        allocate(flat(0))
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        layout = 'C'
        if (present(order)) layout = order
        if (.not. valid_order(layout)) return
        limit = 10000000_int64
        if (present(max_elements)) limit = max_elements
        if (limit < 0_int64) return
        call extent_product(a%shape,total,ok)
        if (.not. ok) return
        if (total > min(int(huge(0),int64),limit)) return
        deallocate(flat)
        allocate(flat(int(total)))
        flat = 0.0_dp
        do k = 1, size(a%data)
            location = flatten_coords(a%coords(:,k),a%shape,layout) + 1_int64
            flat(int(location)) = flat(int(location)) + a%data(k)
            if (.not. ieee_is_finite(flat(int(location)))) then
                deallocate(flat)
                allocate(flat(0))
                return
            end if
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_to_dense_flat

    function coo_nd_from_dense_flat(dims, flat, status, order) result(a)
        integer, intent(in) :: dims(:) !! N-dimensional shape, total matches flat size.
        real(dp), intent(in) :: flat(:) !! Dense values in chosen flattening order.
        integer, intent(out), optional :: status !! Zero on success.
        character(len=1), intent(in), optional :: order !! 'C' default, or 'F'.
        type(coo_nd_array) :: a
        integer(int64) :: total
        logical :: ok
        integer, allocatable :: coordinates(:,:)
        real(dp), allocatable :: values(:)
        character(len=1) :: layout
        integer :: k, j, nnz, ierr

        if (present(status)) status = sparse_invalid_input
        layout = 'C'
        if (present(order)) layout = order
        if (.not. valid_order(layout)) return
        call extent_product(dims,total,ok)
        if (.not. ok) return
        if (total /= int(size(flat),int64)) return
        if (size(flat) == huge(0)) return
        if (.not. all(ieee_is_finite(flat))) return
        nnz = count(flat /= 0.0_dp)
        allocate(coordinates(size(dims),nnz),values(nnz))
        j = 0
        do k = 1, size(flat)
            if (flat(k) == 0.0_dp) cycle
            j = j+1
            call unflatten_coords(int(k-1,int64),dims,layout,coordinates(:,j))
            values(j) = flat(k)
        end do
        a = coo_nd_from_triplets(dims,coordinates,values,ierr)
        if (ierr /= sparse_success) return
        if (layout == 'F' .or. layout == 'f') then
            call coo_nd_sum_duplicates(a,ierr)
            if (ierr /= sparse_success) return
        end if
        if (present(status)) status = sparse_success
    end function coo_nd_from_dense_flat

    function coo_nd_from_dense(dense, status) result(a)
        real(dp), intent(in) :: dense(..) !! Dense assumed-rank real input (ranks 1 through 7).
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_nd_array) :: a
        real(dp), allocatable :: flat(:)
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        select rank(dense)
        rank(1)
            flat = reshape(dense,[size(dense)])
        rank(2)
            flat = reshape(dense,[size(dense)])
        rank(3)
            flat = reshape(dense,[size(dense)])
        rank(4)
            flat = reshape(dense,[size(dense)])
        rank(5)
            flat = reshape(dense,[size(dense)])
        rank(6)
            flat = reshape(dense,[size(dense)])
        rank(7)
            flat = reshape(dense,[size(dense)])
        rank default
            return
        end select
        a = coo_nd_from_dense_flat(shape(dense),flat,ierr,'F')
        if (present(status)) status = ierr
    end function coo_nd_from_dense

    subroutine coo_nd_to_dense(a, dense, status)
        type(coo_nd_array), intent(in) :: a !! COO array, ranks 1 through 7.
        real(dp), intent(out) :: dense(..) !! PREALLOCATED dense output, matching shape.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: flat(:)
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (rank(dense) /= size(a%shape)) return
        if (any(shape(dense) /= a%shape)) return
        flat = coo_nd_to_dense_flat(a,ierr,'F')
        if (ierr /= sparse_success) return
        select rank(dense)
        rank(1)
            dense = reshape(flat,shape(dense))
        rank(2)
            dense = reshape(flat,shape(dense))
        rank(3)
            dense = reshape(flat,shape(dense))
        rank(4)
            dense = reshape(flat,shape(dense))
        rank(5)
            dense = reshape(flat,shape(dense))
        rank(6)
            dense = reshape(flat,shape(dense))
        rank(7)
            dense = reshape(flat,shape(dense))
        rank default
            return
        end select
        if (present(status)) status = sparse_success
    end subroutine coo_nd_to_dense

end module scifort_sparse_coo_nd
