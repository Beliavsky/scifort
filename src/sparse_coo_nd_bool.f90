! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Complement-aware sparse Boolean masks on arbitrary-rank COO grids.
! A mask stores a uniform implicit Boolean background, and the *exceptions*
! to that background as a lexicographically sorted coordinate set. This
! represents the complement of a sparse set without enumerating its entries.
module scifort_sparse_coo_nd_bool
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, valid_coo_nd, &
        coo_nd_from_triplets, coo_nd_sum_duplicates
    implicit none
    private

    integer(int64), parameter :: default_exception_limit = 1000000_int64
    integer, parameter :: operation_and = 1, operation_or = 2
    integer, parameter :: operation_xor = 3, operation_difference = 4

    type, public :: coo_nd_bool_mask
        integer, allocatable :: shape(:) !! Nonnegative extents; ndim >= 1.
        integer, allocatable :: coords(:,:) !! Sorted, unique (ndim,nexceptions), one-based.
        logical :: default_value = .false. !! Value at every coordinate not in coords.
    end type coo_nd_bool_mask

    public :: valid_coo_nd_bool_mask, coo_nd_bool_from_points
    public :: coo_nd_bool_from_coo, coo_nd_bool_to_coo
    public :: coo_nd_bool_compare_scalar, coo_nd_bool_not
    public :: coo_nd_bool_and, coo_nd_bool_or, coo_nd_bool_xor
    public :: coo_nd_bool_difference, coo_nd_bool_get
    public :: coo_nd_bool_count_true, coo_nd_bool_filter

contains

    pure integer function point_compare(point_a, point_b) result(relation)
        integer, intent(in) :: point_a(:) !! First coordinate tuple.
        integer, intent(in) :: point_b(:) !! Second coordinate tuple.
        integer :: ax

        relation = 0
        do ax = 1, size(point_a)
            if (point_a(ax) < point_b(ax)) then
                relation = -1
                return
            else if (point_a(ax) > point_b(ax)) then
                relation = 1
                return
            end if
        end do
    end function point_compare

    pure logical function valid_coo_nd_bool_mask(mask) result(ok)
        type(coo_nd_bool_mask), intent(in) :: mask !! Candidate sparse Boolean mask.
        integer :: k, ax

        ok = .false.
        if (.not. allocated(mask%shape)) return
        if (.not. allocated(mask%coords)) return
        if (size(mask%shape) == 0 .or. any(mask%shape < 0)) return
        if (size(mask%coords,1) /= size(mask%shape)) return
        do ax = 1, size(mask%shape)
            if (any(mask%coords(ax,:) < 1)) return
            if (any(mask%coords(ax,:) > mask%shape(ax))) return
        end do
        do k = 2, size(mask%coords,2)
            if (point_compare(mask%coords(:,k-1),mask%coords(:,k)) /= -1) return
        end do
        ok = .true.
    end function valid_coo_nd_bool_mask

    function coo_nd_bool_from_points(dims, points, status, default_value, max_exceptions) result(mask)
        integer, intent(in) :: dims(:) !! Nonnegative shape extents.
        integer, intent(in) :: points(:,:) !! One-based (ndim,npoints) exceptions; duplicates allowed.
        integer, intent(out), optional :: status !! sparse_success or sparse_invalid_input.
        logical, intent(in), optional :: default_value !! Uniform background; default false.
        integer(int64), intent(in), optional :: max_exceptions !! Maximum unique exceptions; default 1e6.
        type(coo_nd_bool_mask) :: mask
        type(coo_nd_array) :: work
        integer(int64) :: limit
        integer :: info

        if (present(status)) status = sparse_invalid_input
        limit = default_exception_limit
        if (present(max_exceptions)) limit = max_exceptions
        if (limit < 0_int64) return
        if (size(dims) < 1 .or. any(dims < 0)) return
        if (size(points,1) /= size(dims)) return
        if (size(points,2) == huge(0)) return
        work = coo_nd_from_triplets(dims,points, &
            spread(1.0_dp,1,size(points,2)),info)
        if (info /= sparse_success) return
        call coo_nd_sum_duplicates(work,info)
        if (info /= sparse_success) return
        if (int(size(work%data),int64) > limit) return
        mask%shape = dims
        mask%coords = work%coords
        mask%default_value = .false.
        if (present(default_value)) mask%default_value = default_value
        if (present(status)) status = sparse_success
    end function coo_nd_bool_from_points

    function coo_nd_bool_from_coo(a, status, complement, max_exceptions) result(mask)
        type(coo_nd_array), intent(in) :: a !! Real COO mask; duplicates are summed first.
        integer, intent(out), optional :: status !! sparse_success or sparse_invalid_input.
        logical, intent(in), optional :: complement !! True returns logical NOT of numerical nonzero mask.
        integer(int64), intent(in), optional :: max_exceptions !! Maximum stored exception coordinates.
        type(coo_nd_bool_mask) :: mask
        type(coo_nd_array) :: work
        integer(int64) :: limit
        integer :: n, k, j, info

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        limit = default_exception_limit
        if (present(max_exceptions)) limit = max_exceptions
        if (limit < 0_int64) return
        work = a
        call coo_nd_sum_duplicates(work,info)
        if (info /= sparse_success) return
        n = count(work%data /= 0.0_dp)
        if (int(n,int64) > limit) return
        mask%shape = a%shape
        allocate(mask%coords(size(a%shape),n))
        j = 0
        do k = 1, size(work%data)
            if (work%data(k) == 0.0_dp) cycle
            j = j+1
            mask%coords(:,j) = work%coords(:,k)
        end do
        mask%default_value = .false.
        if (present(complement)) mask%default_value = complement
        if (present(status)) status = sparse_success
    end function coo_nd_bool_from_coo

    function coo_nd_bool_to_coo(mask, status, max_nnz) result(a)
        type(coo_nd_bool_mask), intent(in) :: mask !! Boolean mask to materialize.
        integer, intent(out), optional :: status !! Reject implicit-true nonempty backgrounds.
        integer(int64), intent(in), optional :: max_nnz !! Maximum stored COO entries; default 1e6.
        type(coo_nd_array) :: a
        integer(int64) :: limit
        integer :: info

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(mask)) return
        limit = default_exception_limit
        if (present(max_nnz)) limit = max_nnz
        if (limit < 0_int64) return
        if (mask%default_value .and. all(mask%shape > 0)) return
        if (int(size(mask%coords,2),int64) > limit) return
        a = coo_nd_from_triplets(mask%shape,mask%coords, &
            spread(1.0_dp,1,size(mask%coords,2)),info)
        if (info /= sparse_success) return
        if (present(status)) status = sparse_success
    end function coo_nd_bool_to_coo

    function coo_nd_bool_not(mask, status) result(out)
        type(coo_nd_bool_mask), intent(in) :: mask !! Input Boolean mask.
        integer, intent(out), optional :: status !! sparse_success or sparse_invalid_input.
        type(coo_nd_bool_mask) :: out

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(mask)) return
        out = mask
        out%default_value = .not. mask%default_value
        if (present(status)) status = sparse_success
    end function coo_nd_bool_not

    pure logical function apply_binary(aval, bval, operation) result(value)
        logical, intent(in) :: aval !! Left Boolean value.
        logical, intent(in) :: bval !! Right Boolean value.
        integer, intent(in) :: operation !! Internal binary-operation identifier.

        select case (operation)
        case (operation_and)
            value = aval .and. bval
        case (operation_or)
            value = aval .or. bval
        case (operation_xor)
            value = aval .neqv. bval
        case (operation_difference)
            value = aval .and. .not. bval
        case default
            value = .false.
        end select
    end function apply_binary

    ! Merge the two exception sets. Outside their union the result is the
    ! Boolean operation of the two uniform defaults. A coordinate belongs
    ! to the result's exception set iff the operation differs there.
    function bool_binary(a,b,operation,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! Left sparse Boolean mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Right sparse Boolean mask.
        integer, intent(in) :: operation !! AND, OR, XOR or set difference.
        integer, intent(out), optional :: status !! sparse_success or sparse_invalid_input.
        integer(int64), intent(in), optional :: max_exceptions !! Maximum stored exceptions.
        type(coo_nd_bool_mask) :: out
        integer(int64) :: limit
        integer :: pass, i, j, na, nb, count_out, relation
        logical :: av, bv, default_out, current, use_a, use_b
        integer, allocatable :: coord(:)

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(a)) return
        if (.not. valid_coo_nd_bool_mask(b)) return
        if (size(a%shape) /= size(b%shape)) return
        if (any(a%shape /= b%shape)) return
        limit = default_exception_limit
        if (present(max_exceptions)) limit = max_exceptions
        if (limit < 0_int64) return
        na = size(a%coords,2)
        nb = size(b%coords,2)
        default_out = apply_binary(a%default_value,b%default_value,operation)
        allocate(coord(size(a%shape)))
        ! First pass counts; second fills exactly sized storage, with no
        ! allocation proportional to the full tensor extent.
        do pass = 1, 2
            i = 1
            j = 1
            count_out = 0
            do while (i <= na .or. j <= nb)
                use_a = .false.
                use_b = .false.
                if (i > na) then
                    use_b = .true.
                else if (j > nb) then
                    use_a = .true.
                else
                    relation = point_compare(a%coords(:,i),b%coords(:,j))
                    select case (relation)
                    case (-1)
                        use_a = .true.
                    case (0)
                        use_a = .true.
                        use_b = .true.
                    case (1)
                        use_b = .true.
                    end select
                end if
                av = a%default_value
                bv = b%default_value
                if (use_a) then
                    av = .not. av
                    coord = a%coords(:,i)
                    i = i+1
                end if
                if (use_b) then
                    bv = .not. bv
                    coord = b%coords(:,j)
                    j = j+1
                end if
                current = apply_binary(av,bv,operation)
                if (current .eqv. default_out) cycle
                if (pass == 1) then
                    if (int(count_out,int64) >= limit) return
                end if
                count_out = count_out+1
                if (pass == 2) out%coords(:,count_out) = coord
            end do
            if (pass == 1) then
                out%shape = a%shape
                out%default_value = default_out
                allocate(out%coords(size(a%shape),count_out))
            end if
        end do
        if (present(status)) status = sparse_success
    end function bool_binary

    function coo_nd_bool_and(a,b,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! Left mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Right mask.
        integer, intent(out), optional :: status !! Operation status.
        integer(int64), intent(in), optional :: max_exceptions !! Exception storage cap.
        type(coo_nd_bool_mask) :: out

        out = bool_binary(a,b,operation_and,status,max_exceptions)
    end function coo_nd_bool_and

    function coo_nd_bool_or(a,b,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! Left mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Right mask.
        integer, intent(out), optional :: status !! Operation status.
        integer(int64), intent(in), optional :: max_exceptions !! Exception storage cap.
        type(coo_nd_bool_mask) :: out

        out = bool_binary(a,b,operation_or,status,max_exceptions)
    end function coo_nd_bool_or

    function coo_nd_bool_xor(a,b,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! Left mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Right mask.
        integer, intent(out), optional :: status !! Operation status.
        integer(int64), intent(in), optional :: max_exceptions !! Exception storage cap.
        type(coo_nd_bool_mask) :: out

        out = bool_binary(a,b,operation_xor,status,max_exceptions)
    end function coo_nd_bool_xor

    function coo_nd_bool_difference(a,b,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! Left mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Right mask.
        integer, intent(out), optional :: status !! Operation status.
        integer(int64), intent(in), optional :: max_exceptions !! Exception storage cap.
        type(coo_nd_bool_mask) :: out

        out = bool_binary(a,b,operation_difference,status,max_exceptions)
    end function coo_nd_bool_difference

    pure logical function mask_lookup_validated(mask, point) result(value)
        type(coo_nd_bool_mask), intent(in) :: mask !! Already validated sparse Boolean mask.
        integer, intent(in) :: point(:) !! In-bounds coordinate tuple.
        integer :: lo, hi, mid, relation

        value = mask%default_value
        lo = 1
        hi = size(mask%coords,2)
        do while (lo <= hi)
            mid = lo+(hi-lo)/2
            relation = point_compare(point,mask%coords(:,mid))
            if (relation < 0) then
                hi = mid-1
            else if (relation > 0) then
                lo = mid+1
            else
                value = .not. value
                return
            end if
        end do
    end function mask_lookup_validated

    logical function coo_nd_bool_get(mask, point, status) result(value)
        type(coo_nd_bool_mask), intent(in) :: mask !! Boolean mask.
        integer, intent(in) :: point(:) !! One-based coordinate tuple.
        integer, intent(out), optional :: status !! Operation status.

        value = .false.
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(mask)) return
        if (size(point) /= size(mask%shape)) return
        if (any(point < 1) .or. any(point > mask%shape)) return
        value = mask_lookup_validated(mask,point)
        if (present(status)) status = sparse_success
    end function coo_nd_bool_get

    function coo_nd_bool_count_true(mask, status) result(ntrue)
        type(coo_nd_bool_mask), intent(in) :: mask !! Boolean mask.
        integer, intent(out), optional :: status !! Reject total count exceeding int64.
        integer(int64) :: ntrue !! True count or -1 on invalid input / 64-bit overflow.
        integer(int64) :: total, extent
        integer :: ax

        ntrue = -1_int64
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(mask)) return
        if (.not. mask%default_value) then
            ntrue = int(size(mask%coords,2),int64)
            if (present(status)) status = sparse_success
            return
        end if
        total = 0_int64
        if (all(mask%shape > 0)) then
            total = 1_int64
            do ax = 1, size(mask%shape)
                extent = int(mask%shape(ax),int64)
                if (total > huge(total)/extent) return
                total = total*extent
            end do
        end if
        ntrue = total-int(size(mask%coords,2),int64)
        if (present(status)) status = sparse_success
    end function coo_nd_bool_count_true

    pure logical function predicate(value,op,threshold) result(match)
        real(dp), intent(in) :: value !! Summed COO coordinate value.
        character(*), intent(in) :: op !! Comparison operator.
        real(dp), intent(in) :: threshold !! Comparison threshold.

        select case (trim(op))
        case ('>')
            match = value > threshold
        case ('>=')
            match = value >= threshold
        case ('<')
            match = value < threshold
        case ('<=')
            match = value <= threshold
        case ('==')
            match = value == threshold
        case ('/=')
            match = value /= threshold
        case ('abs>')
            match = abs(value) > threshold
        case default
            match = .false.
        end select
    end function predicate

    function coo_nd_bool_compare_scalar(a,op,threshold,status,max_exceptions) result(mask)
        type(coo_nd_array), intent(in) :: a !! Real COO array, possibly with duplicates/zeros.
        character(*), intent(in) :: op !! One of >, >=, <, <=, ==, /=, abs>.
        real(dp), intent(in) :: threshold !! Finite scalar threshold.
        integer, intent(out), optional :: status !! Operation status.
        integer(int64), intent(in), optional :: max_exceptions !! Exception storage cap.
        type(coo_nd_bool_mask) :: mask
        type(coo_nd_array) :: canonical
        integer(int64) :: limit
        integer :: k, j, n, info
        logical :: background

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (.not. ieee_is_finite(threshold)) return
        select case (trim(op))
        case ('>', '>=', '<', '<=', '==', '/=', 'abs>')
        case default
            return
        end select
        limit = default_exception_limit
        if (present(max_exceptions)) limit = max_exceptions
        if (limit < 0_int64) return
        canonical = a
        call coo_nd_sum_duplicates(canonical,info)
        if (info /= sparse_success) return
        background = predicate(0.0_dp,op,threshold)
        n = 0
        do k = 1, size(canonical%data)
            if (predicate(canonical%data(k),op,threshold) .eqv. background) cycle
            if (int(n,int64) >= limit) return
            n = n+1
        end do
        mask%shape = a%shape
        mask%default_value = background
        allocate(mask%coords(size(a%shape),n))
        j = 0
        do k = 1, size(canonical%data)
            if (predicate(canonical%data(k),op,threshold) .eqv. background) cycle
            j = j+1
            mask%coords(:,j) = canonical%coords(:,k)
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_bool_compare_scalar

    function coo_nd_bool_filter(a,mask,status,max_nnz) result(b)
        type(coo_nd_array), intent(in) :: a !! Real COO source; duplicates and zeros preserved.
        type(coo_nd_bool_mask), intent(in) :: mask !! Boolean mask with matching shape.
        integer, intent(out), optional :: status !! Operation status.
        integer(int64), intent(in), optional :: max_nnz !! Maximum output stored entries; default 1e6.
        type(coo_nd_array) :: b
        integer(int64) :: limit
        integer :: k, j, n
        logical :: selected

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        if (.not. valid_coo_nd_bool_mask(mask)) return
        if (size(a%shape) /= size(mask%shape)) return
        if (any(a%shape /= mask%shape)) return
        limit = default_exception_limit
        if (present(max_nnz)) limit = max_nnz
        if (limit < 0_int64) return
        n = 0
        do k = 1, size(a%data)
            selected = mask_lookup_validated(mask,a%coords(:,k))
            if (.not. selected) cycle
            if (int(n,int64) >= limit) return
            n = n+1
        end do
        b%shape = a%shape
        allocate(b%coords(size(a%shape),n),b%data(n))
        j = 0
        do k = 1, size(a%data)
            selected = mask_lookup_validated(mask,a%coords(:,k))
            if (.not. selected) cycle
            j = j+1
            b%coords(:,j) = a%coords(:,k)
            b%data(j) = a%data(k)
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_bool_filter

end module scifort_sparse_coo_nd_bool
