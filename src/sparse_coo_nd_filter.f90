! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Sparse Boolean predicates and coordinate-based masks for real N-D COO.
! A Boolean mask is represented by a real COO array whose nonzero entries
! mean .true.; unstored coordinates are always .false. No dense masks or
! dense tensor-sized lookup tables are constructed.
module scifort_sparse_coo_nd_filter
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, valid_coo_nd, &
        coo_nd_from_triplets, coo_nd_sum_duplicates
    implicit none
    private

    integer(int64), parameter :: default_nnz_limit = 1000000_int64
    public :: coo_nd_compare_scalar, coo_nd_count_nonzero
    public :: coo_nd_filter_sparse, coo_nd_filter_points, coo_nd_filter_value

contains

    pure logical function compare_value(x, op, threshold) result(match)
        real(dp), intent(in) :: x !! Aggregated coordinate value.
        character(*), intent(in) :: op !! One of >, >=, <, <=, ==, /=, abs>.
        real(dp), intent(in) :: threshold !! Finite real comparison value.

        select case (trim(op))
        case ('>')
            match = x > threshold
        case ('>=')
            match = x >= threshold
        case ('<')
            match = x < threshold
        case ('<=')
            match = x <= threshold
        case ('==')
            match = x == threshold
        case ('/=')
            match = x /= threshold
        case ('abs>')
            match = abs(x) > threshold
        case default
            match = .false.
        end select
    end function compare_value

    pure logical function supported_comparison(op) result(ok)
        character(*), intent(in) :: op !! Comparison operator.

        select case (trim(op))
        case ('>', '>=', '<', '<=', '==', '/=', 'abs>')
            ok = .true.
        case default
            ok = .false.
        end select
    end function supported_comparison

    function coo_nd_compare_scalar(a, op, threshold, status, max_nnz) result(mask)
        type(coo_nd_array), intent(in) :: a !! Input with possible duplicates and stored zeros.
        character(*), intent(in) :: op !! Comparison: >, >=, <, <=, ==, /=, abs>.
        real(dp), intent(in) :: threshold !! Scalar threshold; must be finite.
        integer, intent(out), optional :: status !! Zero if successful.
        integer(int64), intent(in), optional :: max_nnz !! Max true coordinates; default 1e6.
        type(coo_nd_array) :: mask
        type(coo_nd_array) :: canonical
        integer(int64) :: limit
        integer :: k, n, j, info

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (.not. ieee_is_finite(threshold)) return
        if (.not. supported_comparison(op)) return
        ! A COO mask has an implicit false background; reject predicates
        ! true at implicit zero rather than enumerating the dense complement.
        if (compare_value(0.0_dp, op, threshold)) return
        limit = default_nnz_limit
        if (present(max_nnz)) limit = max_nnz
        if (limit < 0_int64) return
        limit = min(limit, int(huge(0),int64))
        canonical = a
        call coo_nd_sum_duplicates(canonical,info)
        if (info /= sparse_success) return
        n = 0
        do k = 1, size(canonical%data)
            if (.not. compare_value(canonical%data(k), op, threshold)) cycle
            if (int(n,int64) == limit) return
            n = n+1
        end do
        mask%shape = a%shape
        allocate(mask%coords(size(a%shape),n),mask%data(n))
        j = 0
        do k = 1, size(canonical%data)
            if (.not. compare_value(canonical%data(k),op,threshold)) cycle
            j = j+1
            mask%coords(:,j) = canonical%coords(:,k)
        end do
        mask%data = 1.0_dp
        if (present(status)) status = sparse_success
    end function coo_nd_compare_scalar

    function coo_nd_count_nonzero(a, status) result(n)
        type(coo_nd_array), intent(in) :: a !! Input; duplicate coordinates are summed.
        integer, intent(out), optional :: status !! Zero on success.
        integer(int64) :: n !! Count of nonzero logical entries; -1 on invalid input.
        type(coo_nd_array) :: canonical
        integer :: info

        n = -1_int64
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        canonical = a
        call coo_nd_sum_duplicates(canonical,info)
        if (info /= sparse_success) return
        n = count(canonical%data /= 0.0_dp,kind=int64)
        if (present(status)) status = sparse_success
    end function coo_nd_count_nonzero

    pure integer function compare_point(point, coords, j) result(relation)
        integer, intent(in) :: point(:) !! One-based query coordinate tuple.
        integer, intent(in) :: coords(:,:) !! Sorted COO coordinate storage.
        integer, intent(in) :: j !! Column of coordinates to compare.
        integer :: ax

        relation = 0
        do ax = 1, size(point)
            if (point(ax) < coords(ax,j)) then
                relation = -1
                return
            else if (point(ax) > coords(ax,j)) then
                relation = 1
                return
            end if
        end do
    end function compare_point

    pure logical function mask_lookup(point, mask) result(found)
        integer, intent(in) :: point(:) !! One-based tuple.
        type(coo_nd_array), intent(in) :: mask !! Canonical sorted mask.
        integer :: low, high, mid, relation

        found = .false.
        low = 1
        high = size(mask%data)
        do while (low <= high)
            mid = low + (high-low)/2
            relation = compare_point(point, mask%coords, mid)
            if (relation < 0) then
                high = mid-1
            else if (relation > 0) then
                low = mid+1
            else
                found = mask%data(mid) /= 0.0_dp
                return
            end if
        end do
    end function mask_lookup

    function coo_nd_filter_sparse(a, mask, status, complement, max_nnz) result(b)
        type(coo_nd_array), intent(in) :: a !! Entries to retain/discard, including duplicates/zeros.
        type(coo_nd_array), intent(in) :: mask !! Real COO mask: summed nonzero=true, zero=false.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: complement !! Select coordinates absent/false in mask.
        integer(int64), intent(in), optional :: max_nnz !! Max stored output entries; default 1e6.
        type(coo_nd_array) :: b
        type(coo_nd_array) :: canonical_mask
        integer(int64) :: limit
        integer :: info, k, j, n
        logical :: invert, take

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a) .or. .not. valid_coo_nd(mask)) return
        if (size(a%shape) /= size(mask%shape)) return
        if (any(a%shape /= mask%shape)) return
        limit = default_nnz_limit
        if (present(max_nnz)) limit = max_nnz
        if (limit < 0_int64) return
        limit = min(limit, int(huge(0),int64))
        invert = .false.
        if (present(complement)) invert = complement
        canonical_mask = mask
        call coo_nd_sum_duplicates(canonical_mask,info)
        if (info /= sparse_success) return
        n = 0
        do k = 1, size(a%data)
            take = mask_lookup(a%coords(:,k),canonical_mask) .neqv. invert
            if (.not. take) cycle
            if (int(n,int64) == limit) return
            n = n+1
        end do
        b%shape = a%shape
        allocate(b%coords(size(a%shape),n),b%data(n))
        j = 0
        do k = 1, size(a%data)
            take = mask_lookup(a%coords(:,k),canonical_mask) .neqv. invert
            if (.not. take) cycle
            j = j+1
            b%coords(:,j) = a%coords(:,k)
            b%data(j) = a%data(k)
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_filter_sparse

    function coo_nd_filter_points(a, points, status, complement, max_nnz) result(b)
        type(coo_nd_array), intent(in) :: a !! Sparse source data.
        integer, intent(in) :: points(:,:) !! Coordinates (ndim,npoints), duplicates allowed.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: complement !! Invert coordinate-set membership.
        integer(int64), intent(in), optional :: max_nnz !! Max stored output entries.
        type(coo_nd_array) :: b, mask
        integer :: info

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (size(points,1) /= size(a%shape)) return
        mask = coo_nd_from_triplets(a%shape,points, &
            spread(1.0_dp,1,size(points,2)),info)
        if (info /= sparse_success) return
        b = coo_nd_filter_sparse(a,mask,status,complement,max_nnz)
    end function coo_nd_filter_points

    function coo_nd_filter_value(a, op, threshold, status, max_nnz) result(b)
        type(coo_nd_array), intent(in) :: a !! Sparse source data.
        character(*), intent(in) :: op !! Same comparison operators as compare_scalar.
        real(dp), intent(in) :: threshold !! Finite scalar threshold.
        integer, intent(out), optional :: status !! Zero on success.
        integer(int64), intent(in), optional :: max_nnz !! Max stored output entries.
        type(coo_nd_array) :: b, mask
        integer :: info

        if (present(status)) status = sparse_invalid_input
        mask = coo_nd_compare_scalar(a,op,threshold,info, &
            max_nnz=int(huge(0),int64))
        if (info /= sparse_success) return
        b = coo_nd_filter_sparse(a,mask,status,max_nnz=max_nnz)
    end function coo_nd_filter_value

end module scifort_sparse_coo_nd_filter
