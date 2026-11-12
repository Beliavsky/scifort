! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Sparse reductions, NumPy-style broadcasting, and tensor contraction for
! the positive-rank N-dimensional COO representation. One-based coordinates;
! results are canonical COO unless otherwise noted. No dense intermediates.
module scifort_sparse_coo_nd_ops
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, valid_coo_nd, &
        coo_nd_sum_duplicates, coo_nd_eliminate_zeros
    implicit none
    private

    integer(int64), parameter :: default_nnz_limit = 1000000_int64

    public :: coo_nd_sum_all, coo_nd_sum_axes, coo_nd_sum_axis
    public :: coo_nd_broadcast_to, coo_nd_add, coo_nd_multiply
    public :: coo_nd_tensordot

contains

    real(dp) function coo_nd_sum_all(a, status) result(value)
        type(coo_nd_array), intent(in) :: a !! COO array whose entries are summed.
        integer, intent(out), optional :: status !! Zero on success; invalid input on overflow.
        integer :: i

        value = 0.0_dp
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        do i = 1, size(a%data)
            value = value + a%data(i)
            if (.not. ieee_is_finite(value)) then
                value = 0.0_dp
                return
            end if
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_sum_all

    pure subroutine axis_mask(rank_a, axes, mask, ok)
        integer, intent(in) :: rank_a !! Number of axes in source.
        integer, intent(in) :: axes(:) !! One-based axes; must not repeat.
        logical, intent(out) :: mask(:) !! True for a selected axis.
        logical, intent(out) :: ok !! False for an invalid axis list.
        integer :: j

        mask = .false.
        ok = .false.
        if (size(mask) /= rank_a) return
        do j = 1, size(axes)
            if (axes(j) < 1 .or. axes(j) > rank_a) return
            if (mask(axes(j))) return
            mask(axes(j)) = .true.
        end do
        ok = .true.
    end subroutine axis_mask

    function coo_nd_sum_axes(a, axes, keepdims, status) result(b)
        type(coo_nd_array), intent(in) :: a !! COO array to reduce.
        integer, intent(in) :: axes(:) !! Selected one-based axes; [] is identity.
        logical, intent(in), optional :: keepdims !! Keep reduced axes as length one.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_nd_array) :: b
        type(coo_nd_array) :: temp
        logical, allocatable :: selected(:)
        integer, allocatable :: remaining(:)
        logical :: valid, keep
        integer :: n, i, ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        n = size(a%shape)
        allocate(selected(n))
        call axis_mask(n, axes, selected, valid)
        if (.not. valid) return
        keep = .false.
        if (present(keepdims)) keep = keepdims
        if (all(selected) .and. .not. keep) return ! Scalar: use coo_nd_sum_all.
        if (keep) then
            temp%shape = a%shape
            where (selected) temp%shape = 1
            temp%coords = a%coords
            do i = 1, n
                if (selected(i)) temp%coords(i,:) = 1
            end do
        else
            remaining = pack([(i,i=1,n)],.not. selected)
            temp%shape = a%shape(remaining)
            temp%coords = a%coords(remaining,:)
        end if
        temp%data = a%data
        call coo_nd_sum_duplicates(temp,ierr)
        if (ierr /= sparse_success) return
        b = temp
        if (present(status)) status = sparse_success
    end function coo_nd_sum_axes

    function coo_nd_sum_axis(a, axis, keepdims, status) result(b)
        type(coo_nd_array), intent(in) :: a !! COO array to reduce.
        integer, intent(in) :: axis !! One-based reduced axis.
        logical, intent(in), optional :: keepdims !! Keep the axis with length one.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_nd_array) :: b

        b = coo_nd_sum_axes(a,[axis],keepdims,status)
    end function coo_nd_sum_axis

    pure subroutine broadcasting_shape(shape_a, shape_b, result_shape, ok)
        integer, intent(in) :: shape_a(:) !! First input shape.
        integer, intent(in) :: shape_b(:) !! Second input shape.
        integer, allocatable, intent(out) :: result_shape(:) !! Common broadcast shape.
        logical, intent(out) :: ok !! Compatible shapes.
        integer :: rank_a, rank_b, n, i, j, na, nb

        ok = .false.
        rank_a = size(shape_a)
        rank_b = size(shape_b)
        n = max(rank_a,rank_b)
        allocate(result_shape(n))
        do i = 1, n
            j = i - (n-rank_a)
            na = 1
            if (j > 0) na = shape_a(j)
            j = i - (n-rank_b)
            nb = 1
            if (j > 0) nb = shape_b(j)
            if (na == nb) then
                result_shape(i) = na
            else if (na == 1) then
                result_shape(i) = nb
            else if (nb == 1) then
                result_shape(i) = na
            else
                return
            end if
        end do
        ok = .true.
    end subroutine broadcasting_shape

    function coo_nd_broadcast_to(a, new_shape, status, max_nnz) result(b)
        type(coo_nd_array), intent(in) :: a !! Sparse source.
        integer, intent(in) :: new_shape(:) !! Explicit target shape; trailing axes align.
        integer, intent(out), optional :: status !! Zero on success, invalid on overflow.
        integer(int64), intent(in), optional :: max_nnz !! Max generated entries (default 1e6).
        type(coo_nd_array) :: b
        integer, allocatable :: factor(:), coordinates(:)
        integer(int64) :: cap, copies, count_out, t, left
        integer :: ra, rb, i, j, n, k, ix

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a)) return
        ra = size(a%shape)
        rb = size(new_shape)
        if (rb < ra) return
        if (any(new_shape < 0)) return
        cap = default_nnz_limit
        if (present(max_nnz)) cap = max_nnz
        if (cap < 0_int64) return
        cap = min(cap,int(huge(0),int64))
        allocate(factor(rb),coordinates(rb))
        do i = 1, rb
            j = i-(rb-ra)
            if (j < 1) then
                factor(i) = new_shape(i)
            else if (a%shape(j) == new_shape(i)) then
                factor(i) = 1
            else if (a%shape(j) == 1) then
                factor(i) = new_shape(i)
            else
                return
            end if
        end do
        ! An empty source has no entries to replicate, even when the shape's
        ! theoretical broadcast factor exceeds int64.
        n = size(a%data)
        copies = 1_int64
        if (n > 0) then
            if (any(factor == 0)) then
                copies = 0_int64
            else
                do i = 1, rb
                    if (copies > cap/int(factor(i),int64)) return
                    copies = copies*int(factor(i),int64)
                end do
            end if
        else
            copies = 0_int64
        end if
        count_out = int(n,int64)*copies
        if (count_out > cap) return
        b%shape = new_shape
        allocate(b%coords(rb,int(count_out)),b%data(int(count_out)))
        ix = 0
        do k = 1, n
            do t = 0_int64,copies-1_int64
                left = t
                do i = rb,1,-1
                    j = i-(rb-ra)
                    if (factor(i) > 1) then
                        coordinates(i) = int(mod(left,int(factor(i),int64)))+1
                        left = left/int(factor(i),int64)
                    else if (j > 0) then
                        coordinates(i) = a%coords(j,k)
                    else
                        coordinates(i) = 1
                    end if
                end do
                ix = ix+1
                b%coords(:,ix) = coordinates
                b%data(ix) = a%data(k)
            end do
        end do
        if (present(status)) status = sparse_success
    end function coo_nd_broadcast_to

    function coo_nd_add(a, b, status, max_nnz) result(c)
        type(coo_nd_array), intent(in) :: a !! Left input, with broadcasting.
        type(coo_nd_array), intent(in) :: b !! Right input, with broadcasting.
        integer, intent(out), optional :: status !! Zero on success.
        integer(int64), intent(in), optional :: max_nnz !! Entry cap before coalescing.
        type(coo_nd_array) :: c
        type(coo_nd_array) :: aa, bb, tmp
        integer, allocatable :: target(:)
        integer(int64) :: cap
        integer :: na, nb, ierr
        logical :: ok

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a) .or. .not. valid_coo_nd(b)) return
        call broadcasting_shape(a%shape,b%shape,target,ok)
        if (.not. ok) return
        cap = default_nnz_limit
        if (present(max_nnz)) cap = max_nnz
        if (cap < 0_int64) return
        aa = coo_nd_broadcast_to(a,target,ierr,cap)
        if (ierr /= sparse_success) return
        na = size(aa%data)
        bb = coo_nd_broadcast_to(b,target,ierr,cap-int(na,int64))
        if (ierr /= sparse_success) return
        nb = size(bb%data)
        tmp%shape = target
        allocate(tmp%coords(size(target),na+nb),tmp%data(na+nb))
        tmp%coords(:,1:na) = aa%coords
        tmp%coords(:,na+1:na+nb) = bb%coords
        tmp%data(1:na) = aa%data
        tmp%data(na+1:na+nb) = bb%data
        call coo_nd_sum_duplicates(tmp,ierr)
        if (ierr /= sparse_success) return
        c = tmp
        if (present(status)) status = sparse_success
    end function coo_nd_add

    pure integer function compare_columns(coords_a, col_a, coords_b, col_b, axes_a, axes_b) result(cmp)
        integer, intent(in) :: coords_a(:,:) !! Coordinates for first entry.
        integer, intent(in) :: col_a !! First entry index.
        integer, intent(in) :: coords_b(:,:) !! Coordinates for second entry.
        integer, intent(in) :: col_b !! Second entry index.
        integer, intent(in) :: axes_a(:) !! First axes to compare, in order.
        integer, intent(in) :: axes_b(:) !! Second axes to compare, in order.
        integer :: j, va, vb

        cmp = 0
        do j = 1, size(axes_a)
            va = coords_a(axes_a(j),col_a)
            vb = coords_b(axes_b(j),col_b)
            if (va < vb) then
                cmp = -1
                return
            end if
            if (va > vb) then
                cmp = 1
                return
            end if
        end do
    end function compare_columns

    function coo_nd_multiply(a, b, status, max_nnz) result(c)
        type(coo_nd_array), intent(in) :: a !! Left array, with broadcasting.
        type(coo_nd_array), intent(in) :: b !! Right array, with broadcasting.
        integer, intent(out), optional :: status !! Zero on success.
        integer(int64), intent(in), optional :: max_nnz !! Cap on expanded entries.
        type(coo_nd_array) :: c
        type(coo_nd_array) :: aa, bb, tmp
        integer, allocatable :: target(:), axes(:)
        integer(int64) :: cap
        integer :: ierr, ia, ib, k, cmp, n
        real(dp) :: prod
        logical :: ok

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a) .or. .not. valid_coo_nd(b)) return
        call broadcasting_shape(a%shape,b%shape,target,ok)
        if (.not. ok) return
        cap = default_nnz_limit
        if (present(max_nnz)) cap = max_nnz
        if (cap < 0_int64) return
        aa = coo_nd_broadcast_to(a,target,ierr,cap)
        if (ierr /= sparse_success) return
        bb = coo_nd_broadcast_to(b,target,ierr,cap)
        if (ierr /= sparse_success) return
        call coo_nd_sum_duplicates(aa,ierr)
        if (ierr /= sparse_success) return
        call coo_nd_sum_duplicates(bb,ierr)
        if (ierr /= sparse_success) return
        call coo_nd_eliminate_zeros(aa,ierr)
        call coo_nd_eliminate_zeros(bb,ierr)
        n = size(target)
        axes = [(k,k=1,n)]
        tmp%shape = target
        allocate(tmp%coords(n,min(size(aa%data),size(bb%data))))
        allocate(tmp%data(min(size(aa%data),size(bb%data))))
        ia = 1
        ib = 1
        k = 0
        do while (ia <= size(aa%data) .and. ib <= size(bb%data))
            cmp = compare_columns(aa%coords,ia,bb%coords,ib,axes,axes)
            if (cmp < 0) then
                ia = ia+1
            else if (cmp > 0) then
                ib = ib+1
            else
                prod = aa%data(ia)*bb%data(ib)
                if (.not. ieee_is_finite(prod)) return
                if (prod /= 0.0_dp) then
                    k = k+1
                    tmp%coords(:,k) = aa%coords(:,ia)
                    tmp%data(k) = prod
                end if
                ia = ia+1
                ib = ib+1
            end if
        end do
        c%shape = target
        c%coords = tmp%coords(:,1:k)
        c%data = tmp%data(1:k)
        if (present(status)) status = sparse_success
    end function coo_nd_multiply

    pure subroutine argsort_by_axes(coords, axes, perm)
        integer, intent(in) :: coords(:,:) !! Coordinate matrix.
        integer, intent(in) :: axes(:) !! Selected axes, in grouping order.
        integer, intent(out) :: perm(:) !! Index permutation in key order.
        integer, allocatable :: scratch(:)
        integer :: n, width, left, mid, right, i, j, k

        n = size(perm)
        perm = [(i,i=1,n)]
        allocate(scratch(n))
        width = 1
        do while (width < n)
            left = 1
            do while (left <= n)
                mid = left+min(width,n-left+1)-1
                right = mid+min(width,n-mid)
                i = left
                j = mid+1
                do k = left,right
                    if (i > mid) then
                        scratch(k) = perm(j)
                        j = j+1
                    else if (j > right) then
                        scratch(k) = perm(i)
                        i = i+1
                    else if (compare_columns(coords,perm(j),coords,perm(i),axes,axes) < 0) then
                        scratch(k) = perm(j)
                        j = j+1
                    else
                        scratch(k) = perm(i)
                        i = i+1
                    end if
                end do
                if (right == n) exit
                left = right+1
            end do
            perm = scratch
            if (width >= (n+1)/2) exit
            width = width*2
        end do
    end subroutine argsort_by_axes

    function coo_nd_tensordot(a, b, axes_a, axes_b, status, max_nnz) result(c)
        type(coo_nd_array), intent(in) :: a !! Left COO tensor.
        type(coo_nd_array), intent(in) :: b !! Right COO tensor.
        integer, intent(in) :: axes_a(:) !! One-based left axes to contract.
        integer, intent(in) :: axes_b(:) !! Matching right axes to contract, same order.
        integer, intent(out), optional :: status !! Zero on success; invalid on overflow.
        integer(int64), intent(in), optional :: max_nnz !! Candidate product cap (default 1e6).
        type(coo_nd_array) :: c
        type(coo_nd_array) :: aa, bb, tmp
        logical, allocatable :: mask_a(:), mask_b(:)
        integer, allocatable :: keep_a(:), keep_b(:), perm_a(:), perm_b(:), result_shape(:)
        integer(int64) :: cap, count_out, na_group, nb_group
        integer :: ierr, ra, rb, i, j, i2, j2, k, l, ip, jp, nout, cmp, ia, ib
        logical :: ok_a, ok_b
        real(dp) :: product_value

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd(a) .or. .not. valid_coo_nd(b)) return
        if (size(axes_a) /= size(axes_b)) return
        ra = size(a%shape)
        rb = size(b%shape)
        allocate(mask_a(ra),mask_b(rb))
        call axis_mask(ra,axes_a,mask_a,ok_a)
        call axis_mask(rb,axes_b,mask_b,ok_b)
        if (.not. ok_a .or. .not. ok_b) return
        do k = 1, size(axes_a)
            if (a%shape(axes_a(k)) /= b%shape(axes_b(k))) return
        end do
        keep_a = pack([(k,k=1,ra)],.not. mask_a)
        keep_b = pack([(k,k=1,rb)],.not. mask_b)
        if (size(keep_a)+size(keep_b) == 0) then
            ! Positive-rank storage: a complete contraction is a length-1 COO.
            result_shape = [1]
        else
            result_shape = [a%shape(keep_a),b%shape(keep_b)]
        end if
        cap = default_nnz_limit
        if (present(max_nnz)) cap = max_nnz
        if (cap < 0_int64) return
        cap = min(cap,int(huge(0),int64))
        aa = a
        bb = b
        call coo_nd_sum_duplicates(aa,ierr)
        if (ierr /= sparse_success) return
        call coo_nd_sum_duplicates(bb,ierr)
        if (ierr /= sparse_success) return
        call coo_nd_eliminate_zeros(aa,ierr)
        call coo_nd_eliminate_zeros(bb,ierr)
        allocate(perm_a(size(aa%data)),perm_b(size(bb%data)))
        call argsort_by_axes(aa%coords,axes_a,perm_a)
        call argsort_by_axes(bb%coords,axes_b,perm_b)
        ! Two passes through matching contracted keys: count and then produce
        ! the sparse outer products of each matching group. This avoids dense
        ! tensors and an O(nnz(A)*nnz(B)) all-pairs scan.
        count_out = 0_int64
        do l = 1, 2
            if (l == 2) then
                tmp%shape = result_shape
                allocate(tmp%coords(size(result_shape),int(count_out)))
                allocate(tmp%data(int(count_out)))
                nout = 0
            end if
            i = 1
            j = 1
            do while (i <= size(perm_a) .and. j <= size(perm_b))
                cmp = compare_columns(aa%coords,perm_a(i),bb%coords,perm_b(j),axes_a,axes_b)
                if (cmp < 0) then
                    i = i+1
                    cycle
                end if
                if (cmp > 0) then
                    j = j+1
                    cycle
                end if
                i2 = i+1
                do while (i2 <= size(perm_a))
                    if (compare_columns(aa%coords,perm_a(i),aa%coords,perm_a(i2),axes_a,axes_a) /= 0) exit
                    i2 = i2+1
                end do
                j2 = j+1
                do while (j2 <= size(perm_b))
                    if (compare_columns(bb%coords,perm_b(j),bb%coords,perm_b(j2),axes_b,axes_b) /= 0) exit
                    j2 = j2+1
                end do
                if (l == 1) then
                    na_group = int(i2-i,int64)
                    nb_group = int(j2-j,int64)
                    if (na_group > (cap-count_out)/nb_group) return
                    count_out = count_out + na_group*nb_group
                else
                    do ip = i,i2-1
                        ia = perm_a(ip)
                        do jp = j,j2-1
                            ib = perm_b(jp)
                            product_value = aa%data(ia)*bb%data(ib)
                            if (.not. ieee_is_finite(product_value)) return
                            nout = nout+1
                            if (size(keep_a)+size(keep_b) == 0) then
                                tmp%coords(1,nout) = 1
                            else
                                tmp%coords(:,nout) = [aa%coords(keep_a,ia),bb%coords(keep_b,ib)]
                            end if
                            tmp%data(nout) = product_value
                        end do
                    end do
                end if
                i = i2
                j = j2
            end do
        end do
        call coo_nd_sum_duplicates(tmp,ierr)
        if (ierr /= sparse_success) return
        call coo_nd_eliminate_zeros(tmp,ierr)
        c = tmp
        if (present(status)) status = sparse_success
    end function coo_nd_tensordot

end module scifort_sparse_coo_nd_ops
