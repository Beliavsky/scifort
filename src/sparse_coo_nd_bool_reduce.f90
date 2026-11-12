! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Broadcasting and exact reductions for complement-aware N-D COO Boolean masks.
! No dense tensor is constructed. A sparse count array records a uniform
! int64 baseline and lexicographically sorted deviations from that baseline.
module scifort_sparse_coo_nd_bool_reduce
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, coo_nd_sum_duplicates
    use scifort_sparse_coo_nd_bool, only : coo_nd_bool_mask, valid_coo_nd_bool_mask, &
        coo_nd_bool_from_points, coo_nd_bool_and, coo_nd_bool_or, &
        coo_nd_bool_xor, coo_nd_bool_difference
    implicit none
    private

    integer(int64), parameter :: default_exception_limit = 1000000_int64

    type, public :: coo_nd_bool_count_array
        integer, allocatable :: shape(:) !! Output shape, always positive rank.
        integer, allocatable :: coords(:,:) !! Unique, sorted one-based coordinates.
        integer(int64), allocatable :: values(:) !! Exact counts differing from default_count.
        integer(int64) :: default_count = 0_int64 !! Value at all other coordinates.
    end type coo_nd_bool_count_array

    public :: coo_nd_bool_broadcast_to
    public :: coo_nd_bool_and_broadcast, coo_nd_bool_or_broadcast
    public :: coo_nd_bool_xor_broadcast, coo_nd_bool_difference_broadcast
    public :: coo_nd_bool_count_axis, coo_nd_bool_any_axis, coo_nd_bool_all_axis
    public :: coo_nd_bool_any, coo_nd_bool_all
    public :: valid_coo_nd_bool_count_array, coo_nd_bool_count_get

contains

    pure integer function point_compare(a,b) result(relation)
        integer, intent(in) :: a(:) !! First one-based coordinate.
        integer, intent(in) :: b(:) !! Second one-based coordinate.
        integer :: ax

        relation = 0
        do ax = 1, size(a)
            if (a(ax) < b(ax)) then
                relation = -1
                return
            else if (a(ax) > b(ax)) then
                relation = 1
                return
            end if
        end do
    end function point_compare

    pure logical function valid_coo_nd_bool_count_array(a) result(ok)
        type(coo_nd_bool_count_array), intent(in) :: a !! Candidate count array.
        integer :: ax, k

        ok = .false.
        if (.not. allocated(a%shape)) return
        if (.not. allocated(a%coords)) return
        if (.not. allocated(a%values)) return
        if (size(a%shape) < 1 .or. any(a%shape < 0)) return
        if (size(a%coords,1) /= size(a%shape)) return
        if (size(a%coords,2) /= size(a%values)) return
        if (a%default_count < 0_int64 .or. any(a%values < 0_int64)) return
        if (any(a%values == a%default_count)) return
        do ax = 1, size(a%shape)
            if (any(a%coords(ax,:) < 1)) return
            if (any(a%coords(ax,:) > a%shape(ax))) return
        end do
        do k = 2, size(a%values)
            if (point_compare(a%coords(:,k-1),a%coords(:,k)) /= -1) return
        end do
        ok = .true.
    end function valid_coo_nd_bool_count_array

    integer(int64) function coo_nd_bool_count_get(a, point, status) result(value)
        type(coo_nd_bool_count_array), intent(in) :: a !! Sparse int64 count array.
        integer, intent(in) :: point(:) !! One-based coordinate in the output shape.
        integer, intent(out), optional :: status !! Status code.
        integer :: lo, hi, mid, relation

        value = -1_int64
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_count_array(a)) return
        if (size(point) /= size(a%shape)) return
        if (any(point < 1) .or. any(point > a%shape)) return
        value = a%default_count
        lo = 1
        hi = size(a%values)
        do while (lo <= hi)
            mid = lo+(hi-lo)/2
            relation = point_compare(point,a%coords(:,mid))
            if (relation < 0) then
                hi = mid-1
            else if (relation > 0) then
                lo = mid+1
            else
                value = a%values(mid)
                exit
            end if
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_bool_count_get

    ! Broadcast dimensions are right-aligned as in NumPy. Extent one can
    ! broadcast to extent zero (the result is empty). No operation is
    ! permitted to silently truncate an axis or remove a leading axis.
    function coo_nd_bool_broadcast_to(a, target_shape, status, max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! Source sparse Boolean mask.
        integer, intent(in) :: target_shape(:) !! Requested output dimensions.
        integer, intent(out), optional :: status !! Status code.
        integer(int64), intent(in), optional :: max_exceptions !! Exception expansion cap.
        type(coo_nd_bool_mask) :: out
        integer(int64) :: limit, nout, repeat_count, remaining
        integer :: na, nt, shift, ax, sa, k, j, pos, dim
        integer, allocatable :: expanded(:,:)
        integer :: info

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(a)) return
        na = size(a%shape)
        nt = size(target_shape)
        if (nt < na .or. any(target_shape < 0)) return
        limit = default_exception_limit
        if (present(max_exceptions)) limit = max_exceptions
        if (limit < 0_int64) return
        shift = nt-na
        do ax = shift+1, nt
            sa = ax-shift
            if (a%shape(sa) /= 1 .and. a%shape(sa) /= target_shape(ax)) return
        end do
        if (nt == na) then
            if (all(target_shape == a%shape)) then
                if (int(size(a%coords,2),int64) > limit) return
                out = a
                if (present(status)) status = sparse_success
                return
            end if
        end if

        nout = int(size(a%coords,2),int64)
        ! A zero output extent produces an empty result irrespective of the
        ! nominal expansion factor along other axes.
        if (any(target_shape == 0)) nout = 0_int64
        if (nout > 0_int64) then
            do ax = 1, nt
                sa = ax-shift
                if (sa >= 1) then
                    if (a%shape(sa) /= 1) cycle
                end if
                dim = target_shape(ax)
                if (nout > limit/int(dim,int64)) return
                nout = nout*int(dim,int64)
            end do
        end if
        if (nout > limit .or. nout > int(huge(0),int64)) return
        allocate(expanded(nt,int(nout)))
        if (nout > 0_int64) then
            repeat_count = nout/int(size(a%coords,2),int64)
            pos = 0
            do k = 1, size(a%coords,2)
                do j = 0, int(repeat_count)-1
                    remaining = int(j,int64)
                    pos = pos+1
                    do ax = nt, 1, -1
                        sa = ax-shift
                        if (sa >= 1) then
                            if (a%shape(sa) /= 1) then
                                expanded(ax,pos) = a%coords(sa,k)
                                cycle
                            end if
                        end if
                        dim = target_shape(ax)
                        expanded(ax,pos) = int(mod(remaining,int(dim,int64)))+1
                        remaining = remaining/int(dim,int64)
                    end do
                end do
            end do
        end if
        out = coo_nd_bool_from_points(target_shape,expanded,info, &
            default_value=a%default_value,max_exceptions=limit)
        if (info /= sparse_success) return
        if (present(status)) status = sparse_success
    end function coo_nd_bool_broadcast_to

    ! Resolve common target shape. The right-aligned axes must either agree,
    ! or one of the extents must be one; 1 and 0 broadcast to 0.
    function common_broadcast_shape(a,b,ok) result(shape)
        integer, intent(in) :: a(:) !! First shape.
        integer, intent(in) :: b(:) !! Second shape.
        logical, intent(out) :: ok !! Whether shapes are broadcast-compatible.
        integer, allocatable :: shape(:)
        integer :: rank, ax, ia, ib, ea, eb

        ok = .false.
        rank = max(size(a),size(b))
        allocate(shape(rank))
        do ax = 1, rank
            ia = ax-(rank-size(a))
            ib = ax-(rank-size(b))
            ea = 1
            eb = 1
            if (ia >= 1) ea = a(ia)
            if (ib >= 1) eb = b(ib)
            if (ea == 1) then
                shape(ax) = eb
            else if (eb == 1 .or. ea == eb) then
                shape(ax) = ea
            else
                return
            end if
        end do
        ok = .true.
    end function common_broadcast_shape

    function broadcast_binary(a,b,operation,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! First mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Second mask.
        integer, intent(in) :: operation !! 1=AND, 2=OR, 3=XOR, 4=difference.
        integer, intent(out), optional :: status !! Status code.
        integer(int64), intent(in), optional :: max_exceptions !! Maximum intermediate/final exceptions.
        type(coo_nd_bool_mask) :: out
        type(coo_nd_bool_mask) :: aa, bb
        integer, allocatable :: dims(:)
        integer(int64) :: limit
        integer :: info
        logical :: ok

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(a)) return
        if (.not. valid_coo_nd_bool_mask(b)) return
        limit = default_exception_limit
        if (present(max_exceptions)) limit = max_exceptions
        if (limit < 0_int64) return
        dims = common_broadcast_shape(a%shape,b%shape,ok)
        if (.not. ok) return
        aa = coo_nd_bool_broadcast_to(a,dims,info,limit)
        if (info /= sparse_success) return
        bb = coo_nd_bool_broadcast_to(b,dims,info,limit)
        if (info /= sparse_success) return
        select case (operation)
        case (1)
            out = coo_nd_bool_and(aa,bb,info,limit)
        case (2)
            out = coo_nd_bool_or(aa,bb,info,limit)
        case (3)
            out = coo_nd_bool_xor(aa,bb,info,limit)
        case (4)
            out = coo_nd_bool_difference(aa,bb,info,limit)
        case default
            return
        end select
        if (info /= sparse_success) return
        if (present(status)) status = sparse_success
    end function broadcast_binary

    function coo_nd_bool_and_broadcast(a,b,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! First mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Second mask.
        integer, intent(out), optional :: status !! Status code.
        integer(int64), intent(in), optional :: max_exceptions !! Expansion cap.
        type(coo_nd_bool_mask) :: out

        out = broadcast_binary(a,b,1,status,max_exceptions)
    end function coo_nd_bool_and_broadcast

    function coo_nd_bool_or_broadcast(a,b,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! First mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Second mask.
        integer, intent(out), optional :: status !! Status code.
        integer(int64), intent(in), optional :: max_exceptions !! Expansion cap.
        type(coo_nd_bool_mask) :: out

        out = broadcast_binary(a,b,2,status,max_exceptions)
    end function coo_nd_bool_or_broadcast

    function coo_nd_bool_xor_broadcast(a,b,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! First mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Second mask.
        integer, intent(out), optional :: status !! Status code.
        integer(int64), intent(in), optional :: max_exceptions !! Expansion cap.
        type(coo_nd_bool_mask) :: out

        out = broadcast_binary(a,b,3,status,max_exceptions)
    end function coo_nd_bool_xor_broadcast

    function coo_nd_bool_difference_broadcast(a,b,status,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: a !! First mask.
        type(coo_nd_bool_mask), intent(in) :: b !! Second mask.
        integer, intent(out), optional :: status !! Status code.
        integer(int64), intent(in), optional :: max_exceptions !! Expansion cap.
        type(coo_nd_bool_mask) :: out

        out = broadcast_binary(a,b,4,status,max_exceptions)
    end function coo_nd_bool_difference_broadcast

    ! A rank-one reduction without keepdims is stored in shape [1], since
    ! N-D COO requires positive rank. For rank >=2, default is rank - 1.
    function reduced_shape(shape,axis,keepdims) result(dims)
        integer, intent(in) :: shape(:) !! Original shape.
        integer, intent(in) :: axis !! Fortran 1-based reduced axis.
        logical, intent(in) :: keepdims !! Retain singleton axis.
        integer, allocatable :: dims(:)
        integer :: rank

        rank = size(shape)
        if (keepdims) then
            dims = shape
            dims(axis) = 1
        else if (rank == 1) then
            dims = [1]
        else
            dims = [shape(:axis-1),shape(axis+1:)]
        end if
    end function reduced_shape

    function coo_nd_bool_count_axis(mask,axis,status,keepdims,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: mask !! Source Boolean mask.
        integer, intent(in) :: axis !! Reduced axis (1-based).
        integer, intent(out), optional :: status !! Status code.
        logical, intent(in), optional :: keepdims !! Whether to retain a singleton reduced axis.
        integer(int64), intent(in), optional :: max_exceptions !! Maximum reduced exception coordinates.
        type(coo_nd_bool_count_array) :: out
        type(coo_nd_array) :: reduced
        integer, allocatable :: dims(:), points(:,:)
        integer :: rank, k, info
        integer(int64) :: limit, n, axis_size
        logical :: keep

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(mask)) return
        rank = size(mask%shape)
        if (axis < 1 .or. axis > rank) return
        limit = default_exception_limit
        if (present(max_exceptions)) limit = max_exceptions
        if (limit < 0_int64) return
        keep = .false.
        if (present(keepdims)) keep = keepdims
        dims = reduced_shape(mask%shape,axis,keep)
        n = int(size(mask%coords,2),int64)
        if (n > int(huge(0),int64)) return
        allocate(points(size(dims),int(n)))
        if (n > 0_int64) then
            if (keep) then
                points = mask%coords
                points(axis,:) = 1
            else if (rank == 1) then
                points = 1
            else
                if (axis > 1) points(:axis-1,:) = mask%coords(:axis-1,:)
                if (axis < rank) points(axis:,:) = mask%coords(axis+1:,:)
            end if
        end if
        reduced = coo_nd_from_triplets(dims,points,spread(1.0_dp,1,int(n)),info)
        if (info /= sparse_success) return
        call coo_nd_sum_duplicates(reduced,info)
        if (info /= sparse_success) return
        if (int(size(reduced%data),int64) > limit) return
        out%shape = dims
        out%coords = reduced%coords
        out%default_count = 0_int64
        axis_size = int(mask%shape(axis),int64)
        if (mask%default_value) out%default_count = axis_size
        allocate(out%values(size(reduced%data)))
        do k = 1, size(reduced%data)
            n = int(reduced%data(k),int64) ! Exact: the reduced axis extent <= huge(default integer).
            if (mask%default_value) then
                out%values(k) = axis_size-n
            else
                out%values(k) = n
            end if
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_bool_count_axis

    function reduce_bool_axis(mask,axis,want_all,status,keepdims,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: mask !! Source Boolean mask.
        integer, intent(in) :: axis !! Reduced axis (1-based).
        logical, intent(in) :: want_all !! True for ALL, false for ANY.
        integer, intent(out), optional :: status !! Status code.
        logical, intent(in), optional :: keepdims !! Keep singleton reduced dimension.
        integer(int64), intent(in), optional :: max_exceptions !! Exception cap.
        type(coo_nd_bool_mask) :: out
        type(coo_nd_bool_count_array) :: counts
        integer, allocatable :: points(:,:)
        integer(int64) :: axis_size, limit
        integer :: info, i, j, n
        logical :: baseline, truth

        if (present(status)) status = sparse_invalid_input
        counts = coo_nd_bool_count_axis(mask,axis,info,keepdims,max_exceptions)
        if (info /= sparse_success) return
        limit = default_exception_limit
        if (present(max_exceptions)) limit = max_exceptions
        axis_size = int(mask%shape(axis),int64)
        if (want_all) then
            baseline = counts%default_count == axis_size
        else
            baseline = counts%default_count > 0_int64
        end if
        n = 0
        do i = 1, size(counts%values)
            if (want_all) then
                truth = counts%values(i) == axis_size
            else
                truth = counts%values(i) > 0_int64
            end if
            if (truth .eqv. baseline) cycle
            n = n+1
        end do
        if (int(n,int64) > limit) return
        allocate(points(size(counts%shape),n))
        j = 0
        do i = 1, size(counts%values)
            if (want_all) then
                truth = counts%values(i) == axis_size
            else
                truth = counts%values(i) > 0_int64
            end if
            if (truth .eqv. baseline) cycle
            j = j+1
            points(:,j) = counts%coords(:,i)
        end do
        out%shape = counts%shape
        out%coords = points
        out%default_value = baseline
        if (present(status)) status = sparse_success
    end function reduce_bool_axis

    function coo_nd_bool_any_axis(mask,axis,status,keepdims,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: mask !! Source mask.
        integer, intent(in) :: axis !! Reduced one-based axis.
        integer, intent(out), optional :: status !! Status code.
        logical, intent(in), optional :: keepdims !! Retain a singleton reduced dimension.
        integer(int64), intent(in), optional :: max_exceptions !! Intermediate/output exception cap.
        type(coo_nd_bool_mask) :: out

        out = reduce_bool_axis(mask,axis,.false.,status,keepdims,max_exceptions)
    end function coo_nd_bool_any_axis

    function coo_nd_bool_all_axis(mask,axis,status,keepdims,max_exceptions) result(out)
        type(coo_nd_bool_mask), intent(in) :: mask !! Source mask.
        integer, intent(in) :: axis !! Reduced one-based axis.
        integer, intent(out), optional :: status !! Status code.
        logical, intent(in), optional :: keepdims !! Retain a singleton reduced dimension.
        integer(int64), intent(in), optional :: max_exceptions !! Intermediate/output exception cap.
        type(coo_nd_bool_mask) :: out

        out = reduce_bool_axis(mask,axis,.true.,status,keepdims,max_exceptions)
    end function coo_nd_bool_all_axis

    ! Determine whether every coordinate is an exception without computing
    ! the full shape product; it may be greater than huge(int64).
    pure logical function every_point_is_exception(mask) result(complete)
        type(coo_nd_bool_mask), intent(in) :: mask !! Already validated mask.
        integer(int64) :: total, n
        integer :: axis

        complete = .false.
        if (any(mask%shape == 0)) return
        n = int(size(mask%coords,2),int64)
        total = 1_int64
        do axis = 1, size(mask%shape)
            if (total > n/int(mask%shape(axis),int64)) return
            total = total*int(mask%shape(axis),int64)
        end do
        complete = total == n
    end function every_point_is_exception

    logical function coo_nd_bool_any(mask,status) result(value)
        type(coo_nd_bool_mask), intent(in) :: mask !! Source mask.
        integer, intent(out), optional :: status !! Status code.

        value = .false.
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(mask)) return
        if (.not. any(mask%shape == 0)) then
            if (mask%default_value) then
                value = .not. every_point_is_exception(mask)
            else
                value = size(mask%coords,2) > 0
            end if
        end if
        if (present(status)) status = sparse_success
    end function coo_nd_bool_any

    logical function coo_nd_bool_all(mask,status) result(value)
        type(coo_nd_bool_mask), intent(in) :: mask !! Source mask.
        integer, intent(out), optional :: status !! Status code.

        value = .false.
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_bool_mask(mask)) return
        if (any(mask%shape == 0)) then
            value = .true.
        else if (mask%default_value) then
            value = size(mask%coords,2) == 0
        else
            value = every_point_is_exception(mask)
        end if
        if (present(status)) status = sparse_success
    end function coo_nd_bool_all

end module scifort_sparse_coo_nd_bool_reduce
