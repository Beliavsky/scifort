! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Canonical rectangular compressed sparse row matrices (real binary64).
! Offsets, row numbers, and column numbers are one-based.  Explicit zeros
! are permitted; sorted column indices within each row are unique.
module scifort_sparse
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: sparse_success = 0
    integer, parameter, public :: sparse_invalid_input = 1
    public :: csr_matrix, valid_csr, csr_from_coo, csr_from_dense, csr_to_coo, csr_to_dense
    public :: csr_eye, csr_transpose, csr_matvec, csr_add, csr_matmul
    public :: csr_elementwise_multiply, csr_scale, csr_eliminate_zeros
    public :: csc_matrix, valid_csc, csr_to_csc, csc_to_csr, csc_from_coo
    public :: csc_from_dense, csc_to_dense, csc_matvec, csc_transpose
    public :: csr_diagonal, csr_diags, csr_setdiag, csr_triu, csr_tril
    public :: dia_matrix, valid_dia, dia_from_data, csr_to_dia, dia_to_csr
    public :: dia_from_dense, dia_to_dense, dia_matvec
    public :: bsr_matrix, valid_bsr, csr_to_bsr, bsr_to_csr
    public :: bsr_from_dense, bsr_to_dense, bsr_matvec
    public :: csr_hstack, csr_vstack, csr_block_diag, csr_kron, csr_kronsum
    public :: csr_get, csr_select, csr_slice
    public :: csr_bmat, csr_block, csr_set, csr_assign, csr_assign_coo
    public :: csr_delete_entries, csr_delete_rows, csr_delete_columns, csr_reshape

    interface csr_block
        module procedure csr_bmat
    end interface csr_block

    interface csr_assign
        module procedure csr_assign_dense, csr_assign_scalar, csr_assign_sparse
    end interface csr_assign

    type :: csr_matrix
        integer :: n_rows = 0
        integer :: n_cols = 0
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
        real(dp), allocatable :: data(:)
    end type csr_matrix

    ! Canonical CSC representation: colptr/indptr and row indices are one-based.
    ! Row indices within each column are strictly increasing; stored zeros valid.
    type :: csc_matrix
        integer :: n_rows = 0
        integer :: n_cols = 0
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
        real(dp), allocatable :: data(:)
    end type csc_matrix

    ! Canonical diagonal storage.  offsets are column-minus-row, sorted and
    ! distinct.  data(k,j) stores matrix entry (j-offsets(k),j), including
    ! zero-padded positions off the matrix.  All indices remain one-based.
    type :: dia_matrix
        integer :: n_rows = 0
        integer :: n_cols = 0
        integer, allocatable :: offsets(:)
        real(dp), allocatable :: data(:,:)
    end type dia_matrix

    ! Block CSR with rectangular, fixed-size blocks.  Entire matrix
    ! dimensions must be multiples of the block dimensions.  Block columns
    ! are sorted within each block row; data(:,:,k) is one dense block.
    type :: bsr_matrix
        integer :: n_rows = 0
        integer :: n_cols = 0
        integer :: block_rows = 1
        integer :: block_cols = 1
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
        real(dp), allocatable :: data(:,:,:)
    end type bsr_matrix

contains

    pure logical function valid_csr(a) result(ok)
        type(csr_matrix), intent(in) :: a !! Candidate one-based, canonical CSR matrix.
        integer :: i, k, nnz

        ok = .false.
        if (a%n_rows < 0 .or. a%n_cols < 0 .or. a%n_rows == huge(0)) return
        if (.not. allocated(a%indptr)) return
        if (.not. allocated(a%indices)) return
        if (.not. allocated(a%data)) return
        nnz = size(a%indices)
        if (nnz == huge(nnz)) return
        if (size(a%indptr) /= a%n_rows+1 .or. size(a%data) /= nnz) return
        if (a%indptr(1) /= 1 .or. a%indptr(a%n_rows+1) /= nnz+1) return
        if (any(a%indptr < 1) .or. any(a%indptr > nnz+1)) return
        if (any(a%indptr(2:) < a%indptr(:a%n_rows))) return
        if (any(a%indices < 1) .or. any(a%indices > a%n_cols)) return
        if (.not. all(ieee_is_finite(a%data))) return
        do i = 1, a%n_rows
            do k = a%indptr(i)+1, a%indptr(i+1)-1
                if (a%indices(k) <= a%indices(k-1)) return
            end do
        end do
        ok = .true.
    end function valid_csr

    function csr_from_coo(n_rows, n_cols, rows, columns, values, status) result(a)
        integer, intent(in) :: n_rows !! Number of matrix rows, >= 0.
        integer, intent(in) :: n_cols !! Number of matrix columns, >= 0.
        integer, intent(in) :: rows(:) !! One-based row indices of stored entries.
        integer, intent(in) :: columns(:) !! One-based column indices of stored entries.
        real(dp), intent(in) :: values(:) !! Finite values; duplicates are added together.
        integer, intent(out), optional :: status !! Zero if all coordinates and sums are valid.
        type(csr_matrix) :: a
        integer, allocatable :: ptr(:), cursor(:), cols(:), final_cols(:)
        real(dp), allocatable :: vals(:), final_vals(:)
        integer :: m, i, j, k, first, last, used
        real(dp) :: value

        if (present(status)) status = sparse_invalid_input
        m = size(rows)
        if (n_rows < 0 .or. n_rows == huge(n_rows) .or. n_cols < 0) return
        if (m == huge(m)) return
        if (size(columns) /= m .or. size(values) /= m) return
        if (any(rows < 1) .or. any(rows > n_rows)) return
        if (any(columns < 1) .or. any(columns > n_cols)) return
        if (.not. all(ieee_is_finite(values))) return

        allocate(ptr(n_rows+1), cursor(n_rows), cols(m), vals(m))
        ptr = 0
        ptr(1) = 1
        do k = 1, m
            ptr(rows(k)+1) = ptr(rows(k)+1) + 1
        end do
        do i = 1, n_rows
            ptr(i+1) = ptr(i) + ptr(i+1)
        end do
        cursor = ptr(:n_rows)
        do k = 1, m
            j = cursor(rows(k))
            cols(j) = columns(k)
            vals(j) = values(k)
            cursor(rows(k)) = j+1
        end do

        ! Stable ordering of duplicate values is not needed: their aggregate
        ! follows the column sort, and non-finite sums are rejected.
        allocate(final_cols(m), final_vals(m), a%indptr(n_rows+1))
        used = 0
        do i = 1, n_rows
            a%indptr(i) = used+1
            first = ptr(i)
            last = ptr(i+1)-1
            if (last < first) cycle
            call sort_pairs(cols(first:last), vals(first:last))
            k = first
            do while (k <= last)
                j = cols(k)
                value = vals(k)
                k = k+1
                do while (k <= last)
                    if (cols(k) /= j) exit
                    value = value + vals(k)
                    if (.not. ieee_is_finite(value)) then
                        a = csr_matrix()
                        return
                    end if
                    k = k+1
                end do
                used = used+1
                final_cols(used) = j
                final_vals(used) = value
            end do
        end do
        a%indptr(n_rows+1) = used+1
        a%n_rows = n_rows
        a%n_cols = n_cols
        a%indices = final_cols(:used)
        a%data = final_vals(:used)
        if (present(status)) status = sparse_success
    end function csr_from_coo

    function csr_from_dense(dense, status) result(a)
        real(dp), intent(in) :: dense(:,:) !! Finite matrix; implicit zeros not stored.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: a
        integer :: i, j, k, m, n, nnz

        if (present(status)) status = sparse_invalid_input
        m = size(dense,1)
        n = size(dense,2)
        if (m == huge(m)) return
        if (.not. all(ieee_is_finite(dense))) return
        if (count(dense /= 0.0_dp,kind=int64) >= int(huge(0),int64)) return
        nnz = count(dense /= 0.0_dp)
        allocate(a%indptr(m+1), a%indices(nnz), a%data(nnz))
        a%n_rows = m
        a%n_cols = n
        k = 0
        do i = 1, m
            a%indptr(i) = k+1
            do j = 1, n
                if (dense(i,j) == 0.0_dp) cycle
                k = k+1
                a%indices(k) = j
                a%data(k) = dense(i,j)
            end do
        end do
        a%indptr(m+1) = k+1
        if (present(status)) status = sparse_success
    end function csr_from_dense

    subroutine csr_to_coo(a, rows, columns, values, status)
        type(csr_matrix), intent(in) :: a !! Valid canonical CSR input.
        integer, allocatable, intent(out) :: rows(:) !! One-based stored row coordinates.
        integer, allocatable, intent(out) :: columns(:) !! One-based stored columns.
        real(dp), allocatable, intent(out) :: values(:) !! Stored values, including explicit zeros.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: i, k, m

        if (present(status)) status = sparse_invalid_input
        allocate(rows(0), columns(0), values(0))
        if (.not. valid_csr(a)) return
        m = size(a%data)
        deallocate(rows, columns, values)
        allocate(rows(m), columns(m), values(m))
        columns = a%indices
        values = a%data
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                rows(k) = i
            end do
        end do
        if (present(status)) status = sparse_success
    end subroutine csr_to_coo

    function csr_to_dense(a, status) result(dense)
        type(csr_matrix), intent(in) :: a !! Valid canonical CSR input.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: dense(:,:)
        integer :: i, k

        if (present(status)) status = sparse_invalid_input
        allocate(dense(0,0))
        if (.not. valid_csr(a)) return
        deallocate(dense)
        allocate(dense(a%n_rows,a%n_cols))
        dense = 0.0_dp
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                dense(i,a%indices(k)) = a%data(k)
            end do
        end do
        if (present(status)) status = sparse_success
    end function csr_to_dense

    function csr_eye(n_rows, n_cols, offset, status) result(a)
        integer, intent(in) :: n_rows !! Number of rows, >= 0.
        integer, intent(in), optional :: n_cols !! Number of columns, default n_rows.
        integer, intent(in), optional :: offset !! Diagonal above main (negative below), default 0.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: a
        integer :: i, j, n, k, nnz
        integer(int64) :: jj

        if (present(status)) status = sparse_invalid_input
        n = n_rows
        if (present(n_cols)) n = n_cols
        k = 0
        if (present(offset)) k = offset
        if (n_rows < 0 .or. n_rows == huge(n_rows) .or. n < 0) return
        nnz = 0
        do i = 1, n_rows
            jj = int(i,int64) + int(k,int64)
            if (jj >= 1_int64 .and. jj <= int(n,int64)) nnz = nnz+1
        end do
        a%n_rows = n_rows
        a%n_cols = n
        allocate(a%indptr(n_rows+1), a%indices(nnz), a%data(nnz))
        j = 0
        do i = 1, n_rows
            a%indptr(i) = j+1
            jj = int(i,int64) + int(k,int64)
            if (jj < 1_int64 .or. jj > int(n,int64)) cycle
            j = j+1
            a%indices(j) = int(jj)
            a%data(j) = 1.0_dp
        end do
        a%indptr(n_rows+1) = j+1
        if (present(status)) status = sparse_success
    end function csr_eye

    function csr_transpose(a, status) result(t)
        type(csr_matrix), intent(in) :: a !! Valid canonical CSR input.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: t
        integer, allocatable :: cursor(:)
        integer :: i, j, k, nnz, p

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        nnz = size(a%data)
        t%n_rows = a%n_cols
        t%n_cols = a%n_rows
        if (t%n_rows == huge(0)) then
            t = csr_matrix()
            return
        end if
        allocate(t%indptr(t%n_rows+1), t%indices(nnz), t%data(nnz), cursor(t%n_rows))
        t%indptr = 0
        t%indptr(1) = 1
        do k = 1, nnz
            t%indptr(a%indices(k)+1) = t%indptr(a%indices(k)+1)+1
        end do
        do i = 1, t%n_rows
            t%indptr(i+1) = t%indptr(i+1)+t%indptr(i)
        end do
        cursor = t%indptr(:t%n_rows)
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                j = a%indices(k)
                p = cursor(j)
                t%indices(p) = i
                t%data(p) = a%data(k)
                cursor(j) = p+1
            end do
        end do
        if (present(status)) status = sparse_success
    end function csr_transpose

    function csr_matvec(a, x, status) result(y)
        type(csr_matrix), intent(in) :: a !! Valid CSR left factor.
        real(dp), intent(in) :: x(:) !! Finite vector of length n_cols.
        integer, intent(out), optional :: status !! Zero on success, including finite products.
        real(dp), allocatable :: y(:)
        integer :: i, k
        real(dp) :: sum_value

        if (present(status)) status = sparse_invalid_input
        allocate(y(0))
        if (.not. valid_csr(a)) return
        if (size(x) /= a%n_cols) return
        if (.not. all(ieee_is_finite(x))) return
        deallocate(y)
        allocate(y(a%n_rows))
        y = 0.0_dp
        do i = 1, a%n_rows
            sum_value = 0.0_dp
            do k = a%indptr(i), a%indptr(i+1)-1
                sum_value = sum_value + a%data(k)*x(a%indices(k))
                if (.not. ieee_is_finite(sum_value)) then
                    deallocate(y)
                    allocate(y(0))
                    return
                end if
            end do
            y(i) = sum_value
        end do
        if (present(status)) status = sparse_success
    end function csr_matvec

    function csr_scale(a, factor, status) result(b)
        type(csr_matrix), intent(in) :: a !! Valid CSR input.
        real(dp), intent(in) :: factor !! Finite scalar multiplier.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: b

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (.not. ieee_is_finite(factor)) return
        b = a
        b%data = b%data*factor
        if (.not. all(ieee_is_finite(b%data))) then
            b = csr_matrix()
            return
        end if
        if (present(status)) status = sparse_success
    end function csr_scale

    subroutine csr_eliminate_zeros(a, status)
        type(csr_matrix), intent(inout) :: a !! Valid CSR, modified only if valid.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: new_cols(:), new_ptr(:)
        real(dp), allocatable :: new_vals(:)
        integer :: i, k, p, m

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        m = count(a%data /= 0.0_dp)
        allocate(new_ptr(a%n_rows+1), new_cols(m), new_vals(m))
        p = 0
        do i = 1, a%n_rows
            new_ptr(i) = p+1
            do k = a%indptr(i), a%indptr(i+1)-1
                if (a%data(k) == 0.0_dp) cycle
                p = p+1
                new_cols(p) = a%indices(k)
                new_vals(p) = a%data(k)
            end do
        end do
        new_ptr(a%n_rows+1) = p+1
        call move_alloc(new_ptr, a%indptr)
        call move_alloc(new_cols, a%indices)
        call move_alloc(new_vals, a%data)
        if (present(status)) status = sparse_success
    end subroutine csr_eliminate_zeros

    function csr_add(a, b, status, alpha, beta) result(c)
        type(csr_matrix), intent(in) :: a !! Valid CSR left addend.
        type(csr_matrix), intent(in) :: b !! Valid CSR right addend of same shape.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), intent(in), optional :: alpha !! Scale of a, default 1.
        real(dp), intent(in), optional :: beta !! Scale of b, default 1.
        type(csr_matrix) :: c
        integer, allocatable :: cols(:)
        real(dp), allocatable :: vals(:)
        integer :: i, ka, kb, j, used, na, nb
        real(dp) :: wa, wb, value

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a) .or. .not. valid_csr(b)) return
        if (a%n_rows /= b%n_rows .or. a%n_cols /= b%n_cols) return
        wa = 1.0_dp
        wb = 1.0_dp
        if (present(alpha)) wa = alpha
        if (present(beta)) wb = beta
        if (.not. ieee_is_finite(wa) .or. .not. ieee_is_finite(wb)) return
        na = size(a%data)
        nb = size(b%data)
        if (int(na,int64)+int(nb,int64) >= int(huge(0),int64)) return
        allocate(cols(na+nb), vals(na+nb), c%indptr(a%n_rows+1))
        used = 0
        do i = 1, a%n_rows
            c%indptr(i) = used+1
            ka = a%indptr(i)
            kb = b%indptr(i)
            do while (ka < a%indptr(i+1) .or. kb < b%indptr(i+1))
                if (ka == a%indptr(i+1)) then
                    j = b%indices(kb)
                    value = wb*b%data(kb)
                    kb = kb+1
                else if (kb == b%indptr(i+1)) then
                    j = a%indices(ka)
                    value = wa*a%data(ka)
                    ka = ka+1
                else if (a%indices(ka) < b%indices(kb)) then
                    j = a%indices(ka)
                    value = wa*a%data(ka)
                    ka = ka+1
                else if (a%indices(ka) > b%indices(kb)) then
                    j = b%indices(kb)
                    value = wb*b%data(kb)
                    kb = kb+1
                else
                    j = a%indices(ka)
                    value = wa*a%data(ka) + wb*b%data(kb)
                    ka = ka+1
                    kb = kb+1
                end if
                if (.not. ieee_is_finite(value)) then
                    c = csr_matrix()
                    return
                end if
                if (value == 0.0_dp) cycle
                used = used+1
                cols(used) = j
                vals(used) = value
            end do
        end do
        c%indptr(a%n_rows+1) = used+1
        c%n_rows = a%n_rows
        c%n_cols = a%n_cols
        c%indices = cols(:used)
        c%data = vals(:used)
        if (present(status)) status = sparse_success
    end function csr_add

    function csr_elementwise_multiply(a, b, status) result(c)
        type(csr_matrix), intent(in) :: a !! Valid CSR left factor.
        type(csr_matrix), intent(in) :: b !! Valid CSR right factor with same shape.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: c
        integer, allocatable :: cols(:)
        real(dp), allocatable :: vals(:)
        integer :: i, ka, kb, used
        real(dp) :: value

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a) .or. .not. valid_csr(b)) return
        if (a%n_rows /= b%n_rows .or. a%n_cols /= b%n_cols) return
        allocate(cols(min(size(a%data),size(b%data))), &
                 vals(min(size(a%data),size(b%data))), c%indptr(a%n_rows+1))
        used = 0
        do i = 1, a%n_rows
            c%indptr(i) = used+1
            ka = a%indptr(i)
            kb = b%indptr(i)
            do while (ka < a%indptr(i+1) .and. kb < b%indptr(i+1))
                if (a%indices(ka) < b%indices(kb)) then
                    ka = ka+1
                else if (a%indices(ka) > b%indices(kb)) then
                    kb = kb+1
                else
                    value = a%data(ka)*b%data(kb)
                    if (.not. ieee_is_finite(value)) then
                        c = csr_matrix()
                        return
                    end if
                    if (value /= 0.0_dp) then
                        used = used+1
                        cols(used) = a%indices(ka)
                        vals(used) = value
                    end if
                    ka = ka+1
                    kb = kb+1
                end if
            end do
        end do
        c%indptr(a%n_rows+1) = used+1
        c%n_rows = a%n_rows
        c%n_cols = a%n_cols
        c%indices = cols(:used)
        c%data = vals(:used)
        if (present(status)) status = sparse_success
    end function csr_elementwise_multiply

    function csr_matmul(a, b, status) result(c)
        type(csr_matrix), intent(in) :: a !! Valid CSR left matrix, m-by-k.
        type(csr_matrix), intent(in) :: b !! Valid CSR right matrix, k-by-n.
        integer, intent(out), optional :: status !! Zero on success; rejects nonfinite arithmetic.
        type(csr_matrix) :: c
        integer, allocatable :: marker(:), touched(:), counts(:)
        real(dp), allocatable :: workspace(:)
        integer :: i, pass, count_row, pos, k, j, t, used
        integer(int64) :: total

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a) .or. .not. valid_csr(b)) return
        if (a%n_cols /= b%n_rows) return
        if (b%n_cols == huge(0)) return
        allocate(marker(b%n_cols), touched(b%n_cols), workspace(b%n_cols))
        allocate(counts(a%n_rows))
        marker = 0
        total = 0_int64
        ! Two sparse-accumulator passes avoid storing O(nnz(a)*nnz(b)) products
        ! or allocating a dense m-by-n temporary.
        do pass = 1, 2
            if (pass == 2) then
                c%n_rows = a%n_rows
                c%n_cols = b%n_cols
                allocate(c%indptr(a%n_rows+1), c%indices(int(total)), c%data(int(total)))
                used = 0
            end if
            marker = 0
            do i = 1, a%n_rows
                if (pass == 2) c%indptr(i) = used+1
                count_row = 0
                do k = a%indptr(i), a%indptr(i+1)-1
                    do pos = b%indptr(a%indices(k)), b%indptr(a%indices(k)+1)-1
                        j = b%indices(pos)
                        if (marker(j) /= i) then
                            marker(j) = i
                            count_row = count_row+1
                            touched(count_row) = j
                            workspace(j) = a%data(k)*b%data(pos)
                        else
                            workspace(j) = workspace(j) + a%data(k)*b%data(pos)
                        end if
                        if (.not. ieee_is_finite(workspace(j))) then
                            c = csr_matrix()
                            return
                        end if
                    end do
                end do
                call sort_integers(touched(:count_row))
                t = 0
                do pos = 1, count_row
                    j = touched(pos)
                    if (workspace(j) == 0.0_dp) cycle
                    t = t+1
                    if (pass == 2) then
                        used = used+1
                        c%indices(used) = j
                        c%data(used) = workspace(j)
                    end if
                end do
                if (pass == 1) then
                    counts(i) = t
                    total = total+int(t,int64)
                    if (total >= int(huge(0),int64)) return
                else if (t /= counts(i)) then
                    c = csr_matrix()
                    return
                end if
            end do
        end do
        c%indptr(a%n_rows+1) = used+1
        if (present(status)) status = sparse_success
    end function csr_matmul


    ! Sparse assembly operations preserve canonical row order and stored zeros.
    ! Dimension and nnz arithmetic uses int64 before converting to default integers.
    ! General rectangular block assembly.  A default-constructed csr_matrix()
    ! denotes a missing/zero block; a valid 0-by-0 CSR is a present block.
    ! Each nonempty block row/column fixes its respective dimension.  For a
    ! completely missing block row/column, supply row_sizes/col_sizes.
    function csr_bmat(blocks, status, row_sizes, col_sizes) result(c)
        type(csr_matrix), intent(in) :: blocks(:,:) !! Block grid; csr_matrix() means an omitted zero block.
        integer, intent(out), optional :: status !! sparse_success or sparse_invalid_input.
        integer, intent(in), optional :: row_sizes(:) !! Explicit block-row dimensions, nonnegative.
        integer, intent(in), optional :: col_sizes(:) !! Explicit block-column dimensions, nonnegative.
        type(csr_matrix) :: c
        logical, allocatable :: occupied(:,:)
        integer, allocatable :: heights(:), widths(:)
        integer(int64) :: nrows, ncols, nnz
        integer :: nr, nc, br, bc, i, k, n, pos, local_row, column_offset, row_offset

        if (present(status)) status = sparse_invalid_input
        nr = size(blocks,1)
        nc = size(blocks,2)
        if (nr == 0 .and. nc == 0) then
            c = csr_from_coo(0,0,[integer::],[integer::],[real(dp)::])
            if (present(status)) status = sparse_success
            return
        end if
        if (nr == 0 .or. nc == 0) return
        allocate(occupied(nr,nc),heights(nr),widths(nc))
        occupied = .false.
        heights = -1
        widths = -1
        if (present(row_sizes)) then
            if (size(row_sizes) /= nr) return
            if (any(row_sizes < 0)) return
            heights = row_sizes
        end if
        if (present(col_sizes)) then
            if (size(col_sizes) /= nc) return
            if (any(col_sizes < 0)) return
            widths = col_sizes
        end if
        nnz = 0_int64
        do bc = 1, nc
            do br = 1, nr
                if (valid_csr(blocks(br,bc))) then
                    occupied(br,bc) = .true.
                    if (heights(br) >= 0) then
                        if (heights(br) /= blocks(br,bc)%n_rows) return
                    end if
                    if (widths(bc) >= 0) then
                        if (widths(bc) /= blocks(br,bc)%n_cols) return
                    end if
                    heights(br) = blocks(br,bc)%n_rows
                    widths(bc) = blocks(br,bc)%n_cols
                    nnz = nnz+int(size(blocks(br,bc)%data),int64)
                    if (nnz >= int(huge(0),int64)) return
                else
                    ! Distinguish an intentional missing block from malformed CSR.
                    if (blocks(br,bc)%n_rows /= 0 .or. blocks(br,bc)%n_cols /= 0) return
                    if (allocated(blocks(br,bc)%indptr)) return
                    if (allocated(blocks(br,bc)%indices)) return
                    if (allocated(blocks(br,bc)%data)) return
                end if
            end do
        end do
        if (any(heights < 0) .or. any(widths < 0)) return
        nrows = 0_int64
        ncols = 0_int64
        do br = 1, nr
            nrows = nrows+int(heights(br),int64)
            if (nrows >= int(huge(0),int64)) return
        end do
        do bc = 1, nc
            ncols = ncols+int(widths(bc),int64)
            if (ncols >= int(huge(0),int64)) return
        end do
        c%n_rows = int(nrows)
        c%n_cols = int(ncols)
        allocate(c%indptr(c%n_rows+1),c%indices(int(nnz)),c%data(int(nnz)))
        pos = 1
        row_offset = 0
        do br = 1, nr
            do local_row = 1, heights(br)
                c%indptr(row_offset+local_row) = pos
                column_offset = 0
                do bc = 1, nc
                    if (occupied(br,bc)) then
                        i = blocks(br,bc)%indptr(local_row)
                        k = blocks(br,bc)%indptr(local_row+1)-1
                        n = k-i+1
                        if (n > 0) then
                            c%indices(pos:pos+n-1) = blocks(br,bc)%indices(i:k)+column_offset
                            c%data(pos:pos+n-1) = blocks(br,bc)%data(i:k)
                            pos = pos+n
                        end if
                    end if
                    column_offset = column_offset+widths(bc)
                end do
            end do
            row_offset = row_offset+heights(br)
        end do
        c%indptr(c%n_rows+1) = pos
        if (present(status)) status = sparse_success
    end function csr_bmat

    ! Canonical sparse scatter assignment.  Entries are replaced, NOT added.
    ! Repeated (row,column) coordinates are resolved by last input occurrence.
    ! All validation precedes modification of a; errors leave a unchanged.
    subroutine csr_assign_coo(a, rows, columns, values, status)
        type(csr_matrix), intent(inout) :: a !! CSR matrix updated in-place on success only.
        integer, intent(in) :: rows(:) !! One-based row numbers for replacements.
        integer, intent(in) :: columns(:) !! One-based column numbers for replacements.
        real(dp), intent(in) :: values(:) !! Finite replacement values; explicit zeros stored.
        integer, intent(out), optional :: status !! Zero on successful assignment.
        type(csr_matrix) :: c
        integer(int64), allocatable :: keys(:)
        integer, allocatable :: order(:), ur(:), uc(:)
        real(dp), allocatable :: uv(:)
        integer :: m, k, t, u, nupdates, row, p, last, pass, used, old_col, new_col, pick
        integer(int64) :: max_nnz

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        m = size(rows)
        if (size(columns) /= m .or. size(values) /= m) return
        if (any(rows < 1) .or. any(rows > a%n_rows)) return
        if (any(columns < 1) .or. any(columns > a%n_cols)) return
        if (.not. all(ieee_is_finite(values))) return
        if (m == 0) then
            if (present(status)) status = sparse_success
            return
        end if
        max_nnz = int(size(a%data),int64)+int(m,int64)
        if (max_nnz >= int(huge(0),int64)) return
        allocate(keys(m),order(m),ur(m),uc(m),uv(m))
        do k = 1, m
            keys(k) = int(rows(k)-1,int64)*int(a%n_cols,int64)+int(columns(k),int64)
            order(k) = k
        end do
        call sort_sparse_assignments(keys,order)
        nupdates = 0
        do k = 1, m
            pick = order(k)
            if (k > 1) then
                if (keys(k) == keys(k-1)) then
                    ! This is the latest occurrence (order is increasing for ties).
                    uv(nupdates) = values(pick)
                    cycle
                end if
            end if
            nupdates = nupdates+1
            ur(nupdates) = rows(pick)
            uc(nupdates) = columns(pick)
            uv(nupdates) = values(pick)
        end do
        ! Two passes: count without overallocating, then fill canonical rows.
        do pass = 1, 2
            used = 0
            u = 1
            do row = 1, a%n_rows
                if (pass == 2) c%indptr(row) = used+1
                p = a%indptr(row)
                last = a%indptr(row+1)-1
                do
                    old_col = huge(0)
                    new_col = huge(0)
                    if (p <= last) old_col = a%indices(p)
                    if (u <= nupdates) then
                        if (ur(u) == row) new_col = uc(u)
                    end if
                    if (old_col == huge(0) .and. new_col == huge(0)) exit
                    used = used+1
                    if (pass == 2) then
                        if (new_col <= old_col) then
                            c%indices(used) = new_col
                            c%data(used) = uv(u)
                        else
                            c%indices(used) = old_col
                            c%data(used) = a%data(p)
                        end if
                    end if
                    if (new_col <= old_col) u = u+1
                    if (old_col <= new_col) p = p+1
                end do
            end do
            if (pass == 1) then
                c%n_rows = a%n_rows
                c%n_cols = a%n_cols
                allocate(c%indptr(a%n_rows+1),c%indices(used),c%data(used))
            end if
        end do
        c%indptr(a%n_rows+1) = used+1
        a = c
        if (present(status)) status = sparse_success
    end subroutine csr_assign_coo

    ! Sort assigned coordinates by (packed row-major key, original position).
    ! Tie ordering implements deterministic last-write-wins assignment.
    subroutine sort_sparse_assignments(keys,order)
        integer(int64), intent(inout) :: keys(:) !! Packed integer row-major keys.
        integer, intent(inout) :: order(:) !! Original occurrence indices.
        integer(int64) :: tmpkey
        integer :: n, i, last, tmporder

        n = size(keys)
        do i = n/2, 1, -1
            call sift_sparse_assignment(keys,order,i,n)
        end do
        do last = n, 2, -1
            tmpkey = keys(1)
            keys(1) = keys(last)
            keys(last) = tmpkey
            tmporder = order(1)
            order(1) = order(last)
            order(last) = tmporder
            call sift_sparse_assignment(keys,order,1,last-1)
        end do
    end subroutine sort_sparse_assignments

    subroutine sift_sparse_assignment(keys,order,first,last)
        integer(int64), intent(inout) :: keys(:) !! Heap keys.
        integer, intent(inout) :: order(:) !! Tie-breaking original positions.
        integer, intent(in) :: first !! Starting heap position.
        integer, intent(in) :: last !! Last valid heap position.
        integer :: root, child, tmporder
        integer(int64) :: tmpkey

        root = first
        do while (root <= last/2)
            child = 2*root
            if (child < last) then
                if (keys(child) < keys(child+1) .or. &
                    (keys(child) == keys(child+1) .and. order(child) < order(child+1))) child = child+1
            end if
            if (keys(root) > keys(child)) exit
            if (keys(root) == keys(child)) then
                if (order(root) >= order(child)) exit
            end if
            tmpkey = keys(root)
            keys(root) = keys(child)
            keys(child) = tmpkey
            tmporder = order(root)
            order(root) = order(child)
            order(child) = tmporder
            root = child
        end do
    end subroutine sift_sparse_assignment

    subroutine csr_set(a, row, column, value, status)
        type(csr_matrix), intent(inout) :: a !! CSR matrix modified only on success.
        integer, intent(in) :: row !! One-based row to change.
        integer, intent(in) :: column !! One-based column to change.
        real(dp), intent(in) :: value !! Finite replacement, zero kept explicitly.
        integer, intent(out), optional :: status !! Zero if updated.
        call csr_assign_coo(a,[row],[column],[value],status)
    end subroutine csr_set

    subroutine csr_assign_dense(a, rows, columns, values, status)
        type(csr_matrix), intent(inout) :: a !! CSR input/output, transactional on errors.
        integer, intent(in) :: rows(:) !! Selected one-based row indices, repetitions allowed.
        integer, intent(in) :: columns(:) !! Selected one-based column indices, repetitions allowed.
        real(dp), intent(in) :: values(:,:) !! Replacement shape (size(rows),size(columns)).
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: rr(:), cc(:)
        real(dp), allocatable :: vv(:)
        integer(int64) :: nm
        integer :: i, j, k, nr, nc

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (any(rows < 1) .or. any(rows > a%n_rows)) return
        if (any(columns < 1) .or. any(columns > a%n_cols)) return
        if (.not. all(ieee_is_finite(values))) return
        nr = size(rows)
        nc = size(columns)
        if (size(values,1) /= nr .or. size(values,2) /= nc) return
        nm = int(nr,int64)*int(nc,int64)
        if (nm >= int(huge(0),int64)) return
        allocate(rr(int(nm)),cc(int(nm)),vv(int(nm)))
        k = 0
        do j = 1, nc
            do i = 1, nr
                k = k+1
                rr(k) = rows(i)
                cc(k) = columns(j)
                vv(k) = values(i,j)
            end do
        end do
        call csr_assign_coo(a,rr,cc,vv,status)
    end subroutine csr_assign_dense

    subroutine csr_assign_scalar(a, rows, columns, value, status)
        type(csr_matrix), intent(inout) :: a !! CSR input/output, transactional on errors.
        integer, intent(in) :: rows(:) !! Selected one-based rows.
        integer, intent(in) :: columns(:) !! Selected one-based columns.
        real(dp), intent(in) :: value !! Scalar replacement, broadcast to selected entries.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: rr(:), cc(:)
        real(dp), allocatable :: vv(:)
        integer(int64) :: nm
        integer :: i, j, k, nr, nc

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (any(rows < 1) .or. any(rows > a%n_rows)) return
        if (any(columns < 1) .or. any(columns > a%n_cols)) return
        if (.not. ieee_is_finite(value)) return
        nr = size(rows)
        nc = size(columns)
        nm = int(nr,int64)*int(nc,int64)
        if (nm >= int(huge(0),int64)) return
        allocate(rr(int(nm)),cc(int(nm)),vv(int(nm)))
        k = 0
        do j = 1, nc
            do i = 1, nr
                k = k+1
                rr(k) = rows(i)
                cc(k) = columns(j)
                vv(k) = value
            end do
        end do
        call csr_assign_coo(a,rr,cc,vv,status)
    end subroutine csr_assign_scalar

    ! Overwrite an indexed submatrix from a sparse source without expanding
    ! the Cartesian selection.  Existing target entries in the selected
    ! rectangle are removed, including entries where the source is implicit
    ! zero.  Explicit source zero entries remain stored.  When indices repeat,
    ! the last source row and column occurrences determine the final value.
    ! Memory is O(nnz(a)+nnz(values)+size(rows)+size(columns)), never their
    ! Cartesian product or a%n_cols.  Validation precedes all changes.
    subroutine csr_assign_sparse(a, rows, columns, values, status)
        type(csr_matrix), intent(inout) :: a !! Target CSR, unchanged on failure.
        integer, intent(in) :: rows(:) !! One-based destination rows, repetitions allowed.
        integer, intent(in) :: columns(:) !! One-based destination columns, repetitions allowed.
        type(csr_matrix), intent(in) :: values !! CSR source of selected shape.
        integer, intent(out), optional :: status !! sparse_success or sparse_invalid_input.
        type(csr_matrix) :: updated
        integer, allocatable :: sel_rows(:), sel_cols(:), rr(:), cc(:)
        integer, allocatable :: keys(:), original(:)
        logical, allocatable :: last_row(:), last_col(:)
        real(dp), allocatable :: vv(:)
        integer(int64) :: capacity
        integer :: i, j, k, n, chosen, first, last, used, ierr
        logical :: select_row

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a) .or. .not. valid_csr(values)) return
        if (values%n_rows /= size(rows) .or. values%n_cols /= size(columns)) return
        if (any(rows < 1) .or. any(rows > a%n_rows)) return
        if (any(columns < 1) .or. any(columns > a%n_cols)) return
        if (size(rows) == 0 .or. size(columns) == 0) then
            if (present(status)) status = sparse_success
            return
        end if
        capacity = int(size(a%data),int64)+int(size(values%data),int64)
        if (capacity >= int(huge(0),int64)) return

        call sparse_sort_unique(rows,sel_rows)
        call sparse_sort_unique(columns,sel_cols)
        allocate(last_row(size(rows)), last_col(size(columns)))
        last_row = .false.
        last_col = .false.
        ! Mark the last occurrence of each selected destination row.
        keys = rows
        allocate(original(size(rows)))
        original = [(i,i=1,size(rows))]
        call sort_index_pairs(keys,original)
        i = 1
        do while (i <= size(keys))
            chosen = original(i)
            j = i+1
            do while (j <= size(keys))
                if (keys(j) /= keys(i)) exit
                chosen = max(chosen,original(j))
                j = j+1
            end do
            last_row(chosen) = .true.
            i = j
        end do
        deallocate(keys,original)
        keys = columns
        allocate(original(size(columns)))
        original = [(i,i=1,size(columns))]
        call sort_index_pairs(keys,original)
        i = 1
        do while (i <= size(keys))
            chosen = original(i)
            j = i+1
            do while (j <= size(keys))
                if (keys(j) /= keys(i)) exit
                chosen = max(chosen,original(j))
                j = j+1
            end do
            last_col(chosen) = .true.
            i = j
        end do

        allocate(rr(int(capacity)),cc(int(capacity)),vv(int(capacity)))
        used = 0
        do i = 1, a%n_rows
            first = sparse_lower_bound(sel_rows,i)
            select_row = first <= size(sel_rows)
            if (select_row) select_row = sel_rows(first) == i
            do k = a%indptr(i),a%indptr(i+1)-1
                if (select_row) then
                    j = sparse_lower_bound(sel_cols,a%indices(k))
                    if (j <= size(sel_cols)) then
                        if (sel_cols(j) == a%indices(k)) cycle
                    end if
                end if
                used = used+1
                rr(used) = i
                cc(used) = a%indices(k)
                vv(used) = a%data(k)
            end do
        end do
        do i = 1, values%n_rows
            if (.not. last_row(i)) cycle
            do k = values%indptr(i),values%indptr(i+1)-1
                n = values%indices(k)
                if (.not. last_col(n)) cycle
                used = used+1
                rr(used) = rows(i)
                cc(used) = columns(n)
                vv(used) = values%data(k)
            end do
        end do
        updated = csr_from_coo(a%n_rows,a%n_cols,rr(:used),cc(:used),vv(:used),ierr)
        if (ierr /= sparse_success) return
        a = updated
        if (present(status)) status = sparse_success
    end subroutine csr_assign_sparse

    ! Remove the specified stored coordinates without changing matrix shape.
    ! A missing coordinate is harmless; duplicate deletion requests are valid.
    subroutine csr_delete_entries(a, rows, columns, status)
        type(csr_matrix), intent(inout) :: a !! Modified only after successful validation.
        integer, intent(in) :: rows(:) !! One-based row coordinates, pairwise with columns.
        integer, intent(in) :: columns(:) !! One-based column coordinates.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: out
        integer(int64), allocatable :: keys(:)
        integer(int64) :: key
        integer, allocatable :: order(:)
        integer :: i, k, pos, pass, used, n, m

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        m = size(rows)
        if (size(columns) /= m) return
        if (any(rows < 1) .or. any(rows > a%n_rows)) return
        if (any(columns < 1) .or. any(columns > a%n_cols)) return
        if (m == 0) then
            if (present(status)) status = sparse_success
            return
        end if
        allocate(keys(m),order(m))
        do k = 1, m
            keys(k) = int(rows(k)-1,int64)*int(a%n_cols,int64)+int(columns(k),int64)
            order(k) = k
        end do
        call sort_sparse_assignments(keys,order)
        do pass = 1, 2
            used = 0
            pos = 1
            do i = 1, a%n_rows
                if (pass == 2) out%indptr(i) = used+1
                do k = a%indptr(i),a%indptr(i+1)-1
                    key = int(i-1,int64)*int(a%n_cols,int64)+int(a%indices(k),int64)
                    do while (pos <= m)
                        if (keys(pos) >= key) exit
                        pos = pos+1
                    end do
                    if (pos <= m) then
                        if (keys(pos) == key) cycle
                    end if
                    used = used+1
                    if (pass == 2) then
                        out%indices(used) = a%indices(k)
                        out%data(used) = a%data(k)
                    end if
                end do
            end do
            if (pass == 1) then
                out%n_rows = a%n_rows
                out%n_cols = a%n_cols
                allocate(out%indptr(a%n_rows+1),out%indices(used),out%data(used))
            end if
        end do
        out%indptr(out%n_rows+1) = used+1
        a = out
        if (present(status)) status = sparse_success
    end subroutine csr_delete_entries

    ! Remove complete rows; remaining rows retain their relative order.
    subroutine csr_delete_rows(a, rows, status)
        type(csr_matrix), intent(inout) :: a !! CSR matrix to resize.
        integer, intent(in) :: rows(:) !! One-based rows to remove, duplicates accepted.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: out
        integer, allocatable :: deleted(:)
        integer :: i, j, k, nnew, nnznew, row, pos, len

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (any(rows < 1) .or. any(rows > a%n_rows)) return
        call sparse_sort_unique(rows,deleted)
        if (size(deleted) == 0) then
            if (present(status)) status = sparse_success
            return
        end if
        nnew = a%n_rows-size(deleted)
        nnznew = size(a%data)
        do j = 1, size(deleted)
            i = deleted(j)
            nnznew = nnznew-(a%indptr(i+1)-a%indptr(i))
        end do
        out%n_rows = nnew
        out%n_cols = a%n_cols
        allocate(out%indptr(nnew+1),out%indices(nnznew),out%data(nnznew))
        row = 0
        pos = 1
        j = 1
        do i = 1, a%n_rows
            if (j <= size(deleted)) then
                if (i == deleted(j)) then
                    j = j+1
                    cycle
                end if
            end if
            row = row+1
            out%indptr(row) = pos
            len = a%indptr(i+1)-a%indptr(i)
            if (len > 0) then
                out%indices(pos:pos+len-1) = a%indices(a%indptr(i):a%indptr(i+1)-1)
                out%data(pos:pos+len-1) = a%data(a%indptr(i):a%indptr(i+1)-1)
                pos = pos+len
            end if
        end do
        out%indptr(nnew+1) = pos
        a = out
        if (present(status)) status = sparse_success
    end subroutine csr_delete_rows

    ! Remove complete columns and renumber the retained columns.  No dense
    ! column map is created, even for a matrix with very many logical columns.
    subroutine csr_delete_columns(a, columns, status)
        type(csr_matrix), intent(inout) :: a !! CSR matrix to resize.
        integer, intent(in) :: columns(:) !! One-based columns to remove, duplicates accepted.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: out
        integer, allocatable :: deleted(:)
        integer :: i, k, j, pass, used

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (any(columns < 1) .or. any(columns > a%n_cols)) return
        call sparse_sort_unique(columns,deleted)
        if (size(deleted) == 0) then
            if (present(status)) status = sparse_success
            return
        end if
        do pass = 1, 2
            used = 0
            do i = 1, a%n_rows
                if (pass == 2) out%indptr(i) = used+1
                do k = a%indptr(i),a%indptr(i+1)-1
                    j = sparse_lower_bound(deleted,a%indices(k))
                    if (j <= size(deleted)) then
                        if (deleted(j) == a%indices(k)) cycle
                    end if
                    used = used+1
                    if (pass == 2) then
                        out%indices(used) = a%indices(k)-(j-1)
                        out%data(used) = a%data(k)
                    end if
                end do
            end do
            if (pass == 1) then
                out%n_rows = a%n_rows
                out%n_cols = a%n_cols-size(deleted)
                allocate(out%indptr(a%n_rows+1),out%indices(used),out%data(used))
            end if
        end do
        out%indptr(a%n_rows+1) = used+1
        a = out
        if (present(status)) status = sparse_success
    end subroutine csr_delete_columns

    ! SciPy-compatible C-order reshape: equal row-major flattened elements
    ! before and after, with no dense temporary and stored zeros preserved.
    function csr_reshape(a, n_rows, n_cols, status) result(out)
        type(csr_matrix), intent(in) :: a !! Valid canonical CSR source.
        integer, intent(in) :: n_rows !! New number of rows, nonnegative.
        integer, intent(in) :: n_cols !! New number of columns, nonnegative.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: out
        integer(int64) :: flat
        integer :: i, k, row, col, used

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (n_rows < 0 .or. n_rows == huge(0) .or. n_cols < 0) return
        if (int(a%n_rows,int64)*int(a%n_cols,int64) /= &
            int(n_rows,int64)*int(n_cols,int64)) return
        out%n_rows = n_rows
        out%n_cols = n_cols
        allocate(out%indptr(n_rows+1),out%indices(size(a%data)),out%data(size(a%data)))
        out%indptr = 1
        if (size(a%data) > 0) then
            used = 0
            do i = 1, a%n_rows
                do k = a%indptr(i),a%indptr(i+1)-1
                    flat = int(i-1,int64)*int(a%n_cols,int64)+int(a%indices(k)-1,int64)
                    row = int(flat/int(n_cols,int64))+1
                    col = int(mod(flat,int(n_cols,int64)))+1
                    used = used+1
                    out%indices(used) = col
                    out%data(used) = a%data(k)
                    out%indptr(row+1) = out%indptr(row+1)+1
                end do
            end do
            do i = 1, n_rows
                out%indptr(i+1) = out%indptr(i+1)+out%indptr(i)-1
            end do
        end if
        if (present(status)) status = sparse_success
    end function csr_reshape

    subroutine sparse_sort_unique(values, unique_values)
        integer, intent(in) :: values(:) !! Possibly repeated integers.
        integer, allocatable, intent(out) :: unique_values(:) !! Sorted distinct integers.
        integer, allocatable :: work(:)
        integer :: i, n

        work = values
        call sort_integers(work)
        n = 0
        do i = 1, size(work)
            if (n > 0) then
                if (work(i) == work(n)) cycle
            end if
            n = n+1
            work(n) = work(i)
        end do
        unique_values = work(:n)
    end subroutine sparse_sort_unique

    pure integer function sparse_lower_bound(sorted, value) result(pos)
        integer, intent(in) :: sorted(:) !! Ascending sequence.
        integer, intent(in) :: value !! Element to find or insert.
        integer :: lo, hi, mid

        lo = 1
        hi = size(sorted)+1
        do while (lo < hi)
            mid = lo+(hi-lo)/2
            if (sorted(mid) < value) then
                lo = mid+1
            else
                hi = mid
            end if
        end do
        pos = lo
    end function sparse_lower_bound


    function csr_hstack(blocks, status) result(c)
        type(csr_matrix), intent(in) :: blocks(:) !! Nonempty array of CSR matrices sharing n_rows.
        integer, intent(out), optional :: status !! Zero on success; invalid shape/overflow otherwise.
        type(csr_matrix) :: c
        integer(int64) :: nc, nz
        integer :: t, i, k, first, last, pos, col_offset

        if (present(status)) status = sparse_invalid_input
        if (size(blocks) == 0) return
        if (.not. valid_csr(blocks(1))) return
        nc = 0_int64
        nz = 0_int64
        do t = 1, size(blocks)
            if (.not. valid_csr(blocks(t))) return
            if (blocks(t)%n_rows /= blocks(1)%n_rows) return
            nc = nc+int(blocks(t)%n_cols,int64)
            nz = nz+int(size(blocks(t)%data),int64)
            if (nc >= int(huge(0),int64) .or. nz >= int(huge(0),int64)) return
        end do
        c%n_rows = blocks(1)%n_rows
        c%n_cols = int(nc)
        allocate(c%indptr(c%n_rows+1), c%indices(int(nz)), c%data(int(nz)))
        pos = 1
        do i = 1, c%n_rows
            c%indptr(i) = pos
            col_offset = 0
            do t = 1, size(blocks)
                first = blocks(t)%indptr(i)
                last = blocks(t)%indptr(i+1)-1
                k = last-first+1
                if (k > 0) then
                    c%indices(pos:pos+k-1) = blocks(t)%indices(first:last)+col_offset
                    c%data(pos:pos+k-1) = blocks(t)%data(first:last)
                    pos = pos+k
                end if
                col_offset = col_offset+blocks(t)%n_cols
            end do
        end do
        c%indptr(c%n_rows+1) = pos
        if (present(status)) status = sparse_success
    end function csr_hstack

    function csr_vstack(blocks, status) result(c)
        type(csr_matrix), intent(in) :: blocks(:) !! Nonempty array of CSR matrices sharing n_cols.
        integer, intent(out), optional :: status !! Zero on success; invalid shape/overflow otherwise.
        type(csr_matrix) :: c
        integer(int64) :: nr, nz
        integer :: t, i, n, pos, row_offset

        if (present(status)) status = sparse_invalid_input
        if (size(blocks) == 0) return
        if (.not. valid_csr(blocks(1))) return
        nr = 0_int64
        nz = 0_int64
        do t = 1, size(blocks)
            if (.not. valid_csr(blocks(t))) return
            if (blocks(t)%n_cols /= blocks(1)%n_cols) return
            nr = nr+int(blocks(t)%n_rows,int64)
            nz = nz+int(size(blocks(t)%data),int64)
            if (nr >= int(huge(0),int64) .or. nz >= int(huge(0),int64)) return
        end do
        c%n_rows = int(nr)
        c%n_cols = blocks(1)%n_cols
        allocate(c%indptr(c%n_rows+1), c%indices(int(nz)), c%data(int(nz)))
        pos = 1
        row_offset = 0
        do t = 1, size(blocks)
            n = size(blocks(t)%data)
            do i = 1, blocks(t)%n_rows
                c%indptr(row_offset+i) = pos+blocks(t)%indptr(i)-1
            end do
            if (n > 0) then
                c%indices(pos:pos+n-1) = blocks(t)%indices
                c%data(pos:pos+n-1) = blocks(t)%data
                pos = pos+n
            end if
            row_offset = row_offset+blocks(t)%n_rows
        end do
        c%indptr(c%n_rows+1) = pos
        if (present(status)) status = sparse_success
    end function csr_vstack

    function csr_block_diag(blocks, status) result(c)
        type(csr_matrix), intent(in) :: blocks(:) !! CSR matrices placed along a block diagonal; may be empty.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: c
        integer(int64) :: nr, nc, nz
        integer :: t, i, n, pos, row_offset, col_offset

        if (present(status)) status = sparse_invalid_input
        nr = 0_int64
        nc = 0_int64
        nz = 0_int64
        do t = 1, size(blocks)
            if (.not. valid_csr(blocks(t))) return
            nr = nr+int(blocks(t)%n_rows,int64)
            nc = nc+int(blocks(t)%n_cols,int64)
            nz = nz+int(size(blocks(t)%data),int64)
            if (nr >= int(huge(0),int64) .or. nc >= int(huge(0),int64)) return
            if (nz >= int(huge(0),int64)) return
        end do
        c%n_rows = int(nr)
        c%n_cols = int(nc)
        allocate(c%indptr(c%n_rows+1), c%indices(int(nz)), c%data(int(nz)))
        pos = 1
        row_offset = 0
        col_offset = 0
        do t = 1, size(blocks)
            n = size(blocks(t)%data)
            do i = 1, blocks(t)%n_rows
                c%indptr(row_offset+i) = pos+blocks(t)%indptr(i)-1
            end do
            if (n > 0) then
                c%indices(pos:pos+n-1) = blocks(t)%indices+col_offset
                c%data(pos:pos+n-1) = blocks(t)%data
                pos = pos+n
            end if
            row_offset = row_offset+blocks(t)%n_rows
            col_offset = col_offset+blocks(t)%n_cols
        end do
        c%indptr(c%n_rows+1) = pos
        if (present(status)) status = sparse_success
    end function csr_block_diag

    function csr_kron(a, b, status) result(c)
        type(csr_matrix), intent(in) :: a !! First Kronecker factor, m-by-n.
        type(csr_matrix), intent(in) :: b !! Second Kronecker factor, p-by-q.
        integer, intent(out), optional :: status !! Zero on success, including finite products.
        type(csr_matrix) :: c
        integer(int64) :: nr, nc, nz
        integer :: i, j, ka, kb, pos, row
        real(dp) :: value

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a) .or. .not. valid_csr(b)) return
        nr = int(a%n_rows,int64)*int(b%n_rows,int64)
        nc = int(a%n_cols,int64)*int(b%n_cols,int64)
        nz = int(size(a%data),int64)*int(size(b%data),int64)
        if (nr >= int(huge(0),int64) .or. nc >= int(huge(0),int64)) return
        if (nz >= int(huge(0),int64)) return
        c%n_rows = int(nr)
        c%n_cols = int(nc)
        allocate(c%indptr(c%n_rows+1), c%indices(int(nz)), c%data(int(nz)))
        pos = 1
        row = 0
        do i = 1, a%n_rows
            do j = 1, b%n_rows
                row = row+1
                c%indptr(row) = pos
                do ka = a%indptr(i), a%indptr(i+1)-1
                    do kb = b%indptr(j), b%indptr(j+1)-1
                        value = a%data(ka)*b%data(kb)
                        if (.not. ieee_is_finite(value)) then
                            c = csr_matrix()
                            return
                        end if
                        c%indices(pos) = (a%indices(ka)-1)*b%n_cols+b%indices(kb)
                        c%data(pos) = value
                        pos = pos+1
                    end do
                end do
            end do
        end do
        c%indptr(c%n_rows+1) = pos
        if (present(status)) status = sparse_success
    end function csr_kron

    ! Kronecker sum B (x) I_m + I_n (x) A for square A(m,m), B(n,n).
    ! Build sorted CSR rows directly, merging the two component row streams.
    ! This avoids either Kronecker product temporary and any dense matrix.
    function csr_kronsum(a, b, status) result(c)
        type(csr_matrix), intent(in) :: a !! First square CSR operand, order m.
        type(csr_matrix), intent(in) :: b !! Second square CSR operand, order n.
        integer, intent(out), optional :: status !! sparse_success or invalid input.
        type(csr_matrix) :: c
        integer(int64) :: dimension, nnz_limit
        integer :: m, n, i, j, pa, pb, enda, endb, ia, ib, col, used, row
        real(dp) :: value

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a) .or. .not. valid_csr(b)) return
        if (a%n_rows /= a%n_cols .or. b%n_rows /= b%n_cols) return
        m = a%n_rows
        n = b%n_rows
        dimension = int(m,int64)*int(n,int64)
        nnz_limit = int(n,int64)*int(size(a%data),int64) + &
            int(m,int64)*int(size(b%data),int64)
        if (dimension >= int(huge(0),int64) .or. &
            nnz_limit >= int(huge(0),int64)) return
        c%n_rows = int(dimension)
        c%n_cols = int(dimension)
        allocate(c%indptr(c%n_rows+1),c%indices(int(nnz_limit)),c%data(int(nnz_limit)))
        used = 0
        row = 0
        do i = 1, n
            do j = 1, m
                row = row+1
                c%indptr(row) = used+1
                pa = a%indptr(j)
                enda = a%indptr(j+1)-1
                pb = b%indptr(i)
                endb = b%indptr(i+1)-1
                do while (pa <= enda .or. pb <= endb)
                    ia = huge(0)
                    ib = huge(0)
                    if (pa <= enda) ia = (i-1)*m+a%indices(pa)
                    if (pb <= endb) ib = (b%indices(pb)-1)*m+j
                    col = min(ia,ib)
                    if (ia == ib) then
                        value = a%data(pa)+b%data(pb)
                        pa = pa+1
                        pb = pb+1
                    else if (ia < ib) then
                        value = a%data(pa)
                        pa = pa+1
                    else
                        value = b%data(pb)
                        pb = pb+1
                    end if
                    if (.not. ieee_is_finite(value)) then
                        c = csr_matrix()
                        return
                    end if
                    used = used+1
                    c%indices(used) = col
                    c%data(used) = value
                end do
            end do
        end do
        c%indptr(c%n_rows+1) = used+1
        ! Equal diagonal coordinates from A and B are stored only once,
        ! even when the resulting sum is exactly zero.
        if (used < int(nnz_limit)) then
            c%indices = c%indices(:used)
            c%data = c%data(:used)
        end if
        if (present(status)) status = sparse_success
    end function csr_kronsum

    function csr_get(a, row, column, status) result(value)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(in) :: row !! One-based row, in [1,n_rows].
        integer, intent(in) :: column !! One-based column, in [1,n_cols].
        integer, intent(out), optional :: status !! Zero on success; invalid coordinates otherwise.
        real(dp) :: value
        integer :: low, high, middle

        value = 0.0_dp
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (row < 1 .or. row > a%n_rows) return
        if (column < 1 .or. column > a%n_cols) return
        low = a%indptr(row)
        high = a%indptr(row+1)-1
        do while (low <= high)
            middle = low+(high-low)/2
            if (a%indices(middle) < column) then
                low = middle+1
            else if (a%indices(middle) > column) then
                high = middle-1
            else
                value = a%data(middle)
                exit
            end if
        end do
        if (present(status)) status = sparse_success
    end function csr_get

    function csr_select(a, rows, columns, status) result(c)
        type(csr_matrix), intent(in) :: a !! Canonical CSR input matrix.
        integer, intent(in) :: rows(:) !! One-based row indices, reorderings/repetitions allowed.
        integer, intent(in) :: columns(:) !! One-based column indices, reorderings/repetitions allowed.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: c
        integer, allocatable :: sorted_columns(:), target_positions(:)
        integer(int64) :: total
        integer :: i, j, k, p, pos, n, first, column, last

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (size(rows) == huge(0) .or. size(columns) == huge(0)) return
        if (any(rows < 1) .or. any(rows > a%n_rows)) return
        if (any(columns < 1) .or. any(columns > a%n_cols)) return

        ! Sort the selected source columns with their output positions.  This
        ! needs O(number of selected columns) memory, not O(a%n_cols).
        sorted_columns = columns
        allocate(target_positions(size(columns)), c%indptr(size(rows)+1))
        target_positions = [(j,j=1,size(columns))]
        call sort_index_pairs(sorted_columns, target_positions)

        total = 0_int64
        do i = 1, size(rows)
            c%indptr(i) = int(total)+1
            do k = a%indptr(rows(i)), a%indptr(rows(i)+1)-1
                column = a%indices(k)
                call selected_column_range(sorted_columns, column, first, last)
                total = total+int(last-first+1,int64)
                if (total >= int(huge(0),int64)) then
                    c = csr_matrix()
                    return
                end if
            end do
        end do
        c%indptr(size(rows)+1) = int(total)+1
        c%n_rows = size(rows)
        c%n_cols = size(columns)
        allocate(c%indices(int(total)), c%data(int(total)))
        pos = 0
        do i = 1, size(rows)
            first = pos+1
            do k = a%indptr(rows(i)), a%indptr(rows(i)+1)-1
                column = a%indices(k)
                call selected_column_range(sorted_columns, column, j, last)
                do p = j, last
                    pos = pos+1
                    c%indices(pos) = target_positions(p)
                    c%data(pos) = a%data(k)
                end do
            end do
            n = pos-first+1
            if (n > 1) call sort_pairs(c%indices(first:pos), c%data(first:pos))
        end do
        if (present(status)) status = sparse_success
    end function csr_select

    subroutine selected_column_range(sorted, column, first, last)
        integer, intent(in) :: sorted(:) !! Ascending source-column keys, duplicates allowed.
        integer, intent(in) :: column !! Source-column index to locate.
        integer, intent(out) :: first !! First matching element (or insertion point).
        integer, intent(out) :: last !! Last matching element, possibly smaller than first.
        integer :: low, high, mid

        low = 1
        high = size(sorted)+1
        do while (low < high)
            mid = low+(high-low)/2
            if (sorted(mid) < column) then
                low = mid+1
            else
                high = mid
            end if
        end do
        first = low
        high = size(sorted)+1
        do while (low < high)
            mid = low+(high-low)/2
            if (sorted(mid) <= column) then
                low = mid+1
            else
                high = mid
            end if
        end do
        last = low-1
    end subroutine selected_column_range

    subroutine sort_index_pairs(keys, targets)
        integer, intent(inout) :: keys(:) !! Column keys sorted into ascending order.
        integer, intent(inout) :: targets(:) !! Original output-column positions.
        integer :: start, finish, key, target

        do start = size(keys)/2, 1, -1
            call sift_index_pair(keys, targets, start, size(keys))
        end do
        do finish = size(keys), 2, -1
            key = keys(1)
            keys(1) = keys(finish)
            keys(finish) = key
            target = targets(1)
            targets(1) = targets(finish)
            targets(finish) = target
            call sift_index_pair(keys, targets, 1, finish-1)
        end do
    end subroutine sort_index_pairs

    subroutine sift_index_pair(keys, targets, root, last)
        integer, intent(inout) :: keys(:) !! Heap keys.
        integer, intent(inout) :: targets(:) !! Corresponding output indices.
        integer, intent(in) :: root !! First heap position.
        integer, intent(in) :: last !! Last heap position.
        integer :: i, child, key, target

        i = root
        do while (i <= last/2)
            child = 2*i
            if (child < last) then
                if (keys(child) < keys(child+1)) child = child+1
            end if
            if (keys(i) >= keys(child)) exit
            key = keys(i)
            keys(i) = keys(child)
            keys(child) = key
            target = targets(i)
            targets(i) = targets(child)
            targets(child) = target
            i = child
        end do
    end subroutine sift_index_pair

    function csr_slice(a, row_first, row_last, col_first, col_last, status, row_step, col_step) result(c)
        type(csr_matrix), intent(in) :: a !! Canonical CSR input.
        integer, intent(in) :: row_first !! Inclusive Fortran row start, one-based.
        integer, intent(in) :: row_last !! Inclusive Fortran row end (may denote empty slice).
        integer, intent(in) :: col_first !! Inclusive Fortran column start, one-based.
        integer, intent(in) :: col_last !! Inclusive Fortran column end (may denote empty slice).
        integer, intent(out), optional :: status !! Zero on success.
        integer, intent(in), optional :: row_step !! Nonzero stride, default 1; negative reverses.
        integer, intent(in), optional :: col_step !! Nonzero stride, default 1; negative reverses.
        type(csr_matrix) :: c
        integer, allocatable :: rows(:), columns(:)
        integer :: rs, cs, info
        logical :: ok

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        rs = 1
        cs = 1
        if (present(row_step)) rs = row_step
        if (present(col_step)) cs = col_step
        call sparse_slice_indices(a%n_rows, row_first, row_last, rs, rows, ok)
        if (.not. ok) return
        call sparse_slice_indices(a%n_cols, col_first, col_last, cs, columns, ok)
        if (.not. ok) return
        c = csr_select(a, rows, columns, info)
        if (present(status)) status = info
    end function csr_slice

    subroutine sparse_slice_indices(n, first, last, stride, indices, ok)
        integer, intent(in) :: n !! Size of the selected axis.
        integer, intent(in) :: first !! Inclusive starting coordinate.
        integer, intent(in) :: last !! Inclusive stopping coordinate.
        integer, intent(in) :: stride !! Nonzero, signed traversal step.
        integer, allocatable, intent(out) :: indices(:) !! One-based selected positions.
        logical, intent(out) :: ok !! Whether slice coordinates are well formed.
        integer(int64) :: count64, k
        integer :: i

        ok = .false.
        allocate(indices(0))
        if (stride == 0) return
        if (stride > 0) then
            if (first < 1 .or. int(first,int64) > int(n,int64)+1_int64) return
            if (last < 0 .or. last > n) return
            count64 = 0_int64
            if (first <= last) count64 = 1_int64+(int(last,int64)-int(first,int64))/int(stride,int64)
        else
            if (first < 0 .or. first > n) return
            if (last < 1 .or. int(last,int64) > int(n,int64)+1_int64) return
            count64 = 0_int64
            if (first >= last) count64 = 1_int64+(int(first,int64)-int(last,int64))/(-int(stride,int64))
        end if
        if (count64 >= int(huge(0),int64)) return
        deallocate(indices)
        allocate(indices(int(count64)))
        do i = 1, size(indices)
            k = int(first,int64)+int(i-1,int64)*int(stride,int64)
            indices(i) = int(k)
        end do
        ok = .true.
    end subroutine sparse_slice_indices

    pure logical function valid_csc(a) result(ok)
        type(csc_matrix), intent(in) :: a !! Candidate canonical CSC matrix.
        integer :: j, k, nnz

        ok = .false.
        if (a%n_rows < 0 .or. a%n_cols < 0 .or. a%n_cols == huge(0)) return
        if (.not. allocated(a%indptr)) return
        if (.not. allocated(a%indices)) return
        if (.not. allocated(a%data)) return
        nnz = size(a%data)
        if (nnz == huge(0)) return
        if (size(a%indptr) /= a%n_cols+1 .or. size(a%indices) /= nnz) return
        if (a%indptr(1) /= 1 .or. a%indptr(a%n_cols+1) /= nnz+1) return
        if (any(a%indptr < 1) .or. any(a%indptr > nnz+1)) return
        if (any(a%indptr(2:) < a%indptr(:a%n_cols))) return
        if (any(a%indices < 1) .or. any(a%indices > a%n_rows)) return
        if (.not. all(ieee_is_finite(a%data))) return
        do j = 1, a%n_cols
            do k = a%indptr(j)+1, a%indptr(j+1)-1
                if (a%indices(k) <= a%indices(k-1)) return
            end do
        end do
        ok = .true.
    end function valid_csc

    function csr_to_csc(a, status) result(b)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix, m-by-n.
        integer, optional, intent(out) :: status !! Zero on success.
        type(csc_matrix) :: b
        integer, allocatable :: cursor(:)
        integer :: i, j, k, pos, nnz

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (a%n_cols == huge(0)) return
        nnz = size(a%data)
        b%n_rows = a%n_rows
        b%n_cols = a%n_cols
        allocate(b%indptr(a%n_cols+1), b%indices(nnz), b%data(nnz), cursor(a%n_cols))
        b%indptr = 0
        b%indptr(1) = 1
        do k = 1, nnz
            j = a%indices(k)
            b%indptr(j+1) = b%indptr(j+1)+1
        end do
        do j = 1, a%n_cols
            b%indptr(j+1) = b%indptr(j+1)+b%indptr(j)
        end do
        cursor = b%indptr(:a%n_cols)
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                j = a%indices(k)
                pos = cursor(j)
                b%indices(pos) = i
                b%data(pos) = a%data(k)
                cursor(j) = pos+1
            end do
        end do
        if (present(status)) status = sparse_success
    end function csr_to_csc

    function csc_to_csr(a, status) result(b)
        type(csc_matrix), intent(in) :: a !! Canonical CSC matrix, m-by-n.
        integer, optional, intent(out) :: status !! Zero on success.
        type(csr_matrix) :: b
        integer, allocatable :: cursor(:)
        integer :: i, j, k, pos, nnz

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csc(a)) return
        if (a%n_rows == huge(0)) return
        nnz = size(a%data)
        b%n_rows = a%n_rows
        b%n_cols = a%n_cols
        allocate(b%indptr(a%n_rows+1), b%indices(nnz), b%data(nnz), cursor(a%n_rows))
        b%indptr = 0
        b%indptr(1) = 1
        do k = 1, nnz
            i = a%indices(k)
            b%indptr(i+1) = b%indptr(i+1)+1
        end do
        do i = 1, a%n_rows
            b%indptr(i+1) = b%indptr(i+1)+b%indptr(i)
        end do
        cursor = b%indptr(:a%n_rows)
        do j = 1, a%n_cols
            do k = a%indptr(j), a%indptr(j+1)-1
                i = a%indices(k)
                pos = cursor(i)
                b%indices(pos) = j
                b%data(pos) = a%data(k)
                cursor(i) = pos+1
            end do
        end do
        if (present(status)) status = sparse_success
    end function csc_to_csr

    function csc_from_coo(n_rows, n_cols, rows, columns, values, status) result(a)
        integer, intent(in) :: n_rows !! Number of rows, >= 0.
        integer, intent(in) :: n_cols !! Number of columns, >= 0.
        integer, intent(in) :: rows(:) !! One-based row coordinates.
        integer, intent(in) :: columns(:) !! One-based column coordinates.
        real(dp), intent(in) :: values(:) !! Finite elements; duplicates summed.
        integer, optional, intent(out) :: status !! Zero on success.
        type(csc_matrix) :: a
        type(csr_matrix) :: b
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        b = csr_from_coo(n_rows, n_cols, rows, columns, values, ierr)
        if (ierr /= sparse_success) return
        a = csr_to_csc(b, ierr)
        if (present(status)) status = ierr
    end function csc_from_coo

    function csc_from_dense(dense, status) result(a)
        real(dp), intent(in) :: dense(:,:) !! Dense input, no implicit zero storage.
        integer, optional, intent(out) :: status !! Zero on success.
        type(csc_matrix) :: a
        type(csr_matrix) :: b
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        b = csr_from_dense(dense, ierr)
        if (ierr /= sparse_success) return
        a = csr_to_csc(b, ierr)
        if (present(status)) status = ierr
    end function csc_from_dense

    function csc_to_dense(a, status) result(dense)
        type(csc_matrix), intent(in) :: a !! Canonical CSC input.
        integer, optional, intent(out) :: status !! Zero on success.
        real(dp), allocatable :: dense(:,:)
        integer :: j, k

        if (present(status)) status = sparse_invalid_input
        allocate(dense(0,0))
        if (.not. valid_csc(a)) return
        deallocate(dense)
        allocate(dense(a%n_rows,a%n_cols))
        dense = 0.0_dp
        do j = 1, a%n_cols
            do k = a%indptr(j), a%indptr(j+1)-1
                dense(a%indices(k),j) = a%data(k)
            end do
        end do
        if (present(status)) status = sparse_success
    end function csc_to_dense

    function csc_matvec(a, x, status) result(y)
        type(csc_matrix), intent(in) :: a !! Canonical CSC matrix.
        real(dp), intent(in) :: x(:) !! Dense vector of length n_cols.
        integer, optional, intent(out) :: status !! Zero on success.
        real(dp), allocatable :: y(:)
        integer :: j, k

        if (present(status)) status = sparse_invalid_input
        allocate(y(0))
        if (.not. valid_csc(a)) return
        if (size(x) /= a%n_cols) return
        if (.not. all(ieee_is_finite(x))) return
        deallocate(y)
        allocate(y(a%n_rows))
        y = 0.0_dp
        do j = 1, a%n_cols
            do k = a%indptr(j), a%indptr(j+1)-1
                y(a%indices(k)) = y(a%indices(k)) + a%data(k)*x(j)
                if (.not. ieee_is_finite(y(a%indices(k)))) then
                    deallocate(y)
                    allocate(y(0))
                    return
                end if
            end do
        end do
        if (present(status)) status = sparse_success
    end function csc_matvec

    function csc_transpose(a, status) result(t)
        type(csc_matrix), intent(in) :: a !! Canonical CSC input.
        integer, optional, intent(out) :: status !! Zero on success.
        type(csc_matrix) :: t
        type(csr_matrix) :: b
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        b = csc_to_csr(a, ierr)
        if (ierr /= sparse_success) return
        ! CSR(A) has exactly the same physical layout as CSC(transpose(A)).
        t%n_rows = a%n_cols
        t%n_cols = a%n_rows
        call move_alloc(b%indptr, t%indptr)
        call move_alloc(b%indices, t%indices)
        call move_alloc(b%data, t%data)
        if (present(status)) status = sparse_success
    end function csc_transpose

    function csr_diagonal(a, offset, status) result(d)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, optional, intent(in) :: offset !! Column minus row, default zero.
        integer, optional, intent(out) :: status !! Zero on success.
        real(dp), allocatable :: d(:)
        integer(int64) :: i_first, i_last, j64, kk
        integer :: i, k, start_i

        if (present(status)) status = sparse_invalid_input
        allocate(d(0))
        if (.not. valid_csr(a)) return
        kk = 0_int64
        if (present(offset)) kk = int(offset,int64)
        i_first = max(1_int64,1_int64-kk)
        i_last = min(int(a%n_rows,int64),int(a%n_cols,int64)-kk)
        if (i_last < i_first) then
            if (present(status)) status = sparse_success
            return
        end if
        deallocate(d)
        allocate(d(int(i_last-i_first+1_int64)))
        d = 0.0_dp
        start_i = int(i_first)
        do i = start_i, int(i_last)
            j64 = int(i,int64)+kk
            do k = a%indptr(i), a%indptr(i+1)-1
                if (int(a%indices(k),int64) >= j64) then
                    if (int(a%indices(k),int64) == j64) d(i-start_i+1) = a%data(k)
                    exit
                end if
            end do
        end do
        if (present(status)) status = sparse_success
    end function csr_diagonal

    function csr_diags(n_rows, n_cols, offsets, diag_ptr, values, status) result(a)
        integer, intent(in) :: n_rows !! Number of rows, >= 0.
        integer, intent(in) :: n_cols !! Number of columns, >= 0.
        integer, intent(in) :: offsets(:) !! Unique column-minus-row offsets, any order.
        integer, intent(in) :: diag_ptr(:) !! One-based start of each packed diagonal.
        real(dp), intent(in) :: values(:) !! Diagonal entries, row-wise within each diagonal.
        integer, optional, intent(out) :: status !! Zero on success; invalid lengths rejected.
        type(csr_matrix) :: a
        integer, allocatable :: rows(:), cols(:)
        integer(int64) :: start_i, end_i, kk, len
        integer :: t, j, pos, n, ierr

        if (present(status)) status = sparse_invalid_input
        n = size(offsets)
        if (n_rows < 0 .or. n_rows == huge(0) .or. n_cols < 0) return
        if (n == huge(0) .or. size(values) == huge(0)) return
        if (size(diag_ptr) /= n+1) return
        if (diag_ptr(1) /= 1 .or. diag_ptr(n+1) /= size(values)+1) return
        if (any(diag_ptr < 1) .or. any(diag_ptr > size(values)+1)) return
        if (any(diag_ptr(2:) < diag_ptr(:n))) return
        if (.not. all(ieee_is_finite(values))) return
        do t = 1, n
            if (any(offsets(:t-1) == offsets(t))) return
            kk = int(offsets(t),int64)
            start_i = max(1_int64,1_int64-kk)
            end_i = min(int(n_rows,int64),int(n_cols,int64)-kk)
            len = max(0_int64,end_i-start_i+1_int64)
            if (int(diag_ptr(t+1)-diag_ptr(t),int64) /= len) return
        end do
        allocate(rows(size(values)), cols(size(values)))
        pos = 0
        do t = 1, n
            kk = int(offsets(t),int64)
            start_i = max(1_int64,1_int64-kk)
            do j = diag_ptr(t), diag_ptr(t+1)-1
                pos = pos+1
                rows(pos) = int(start_i+int(j-diag_ptr(t),int64))
                cols(pos) = int(int(rows(pos),int64)+kk)
            end do
        end do
        a = csr_from_coo(n_rows,n_cols,rows,cols,values,ierr)
        if (present(status)) status = ierr
    end function csr_diags

    subroutine csr_setdiag(a, values, offset, status)
        type(csr_matrix), intent(inout) :: a !! Matrix updated only on success.
        real(dp), intent(in) :: values(:) !! Entire target diagonal in increasing row order.
        integer, optional, intent(in) :: offset !! Column-minus-row, default zero.
        integer, optional, intent(out) :: status !! Zero on success.
        type(csr_matrix) :: b
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer(int64) :: start_i, end_i, kk, diag_len
        integer :: nnz, k, i, j, used, ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (.not. all(ieee_is_finite(values))) return
        kk = 0_int64
        if (present(offset)) kk = int(offset,int64)
        start_i = max(1_int64,1_int64-kk)
        end_i = min(int(a%n_rows,int64),int(a%n_cols,int64)-kk)
        diag_len = max(0_int64,end_i-start_i+1_int64)
        if (int(size(values),int64) /= diag_len) return
        nnz = size(a%data)
        if (int(nnz,int64)+diag_len >= int(huge(0),int64)) return
        allocate(rows(nnz+size(values)), cols(nnz+size(values)), vals(nnz+size(values)))
        used = 0
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                if (int(a%indices(k),int64) == int(i,int64)+kk) cycle
                used = used+1
                rows(used) = i
                cols(used) = a%indices(k)
                vals(used) = a%data(k)
            end do
        end do
        do j = 1, size(values)
            used = used+1
            rows(used) = int(start_i)+j-1
            cols(used) = int(int(rows(used),int64)+kk)
            vals(used) = values(j)
        end do
        b = csr_from_coo(a%n_rows,a%n_cols,rows(:used),cols(:used),vals(:used),ierr)
        if (ierr /= sparse_success) return
        a = b
        if (present(status)) status = sparse_success
    end subroutine csr_setdiag

    function csr_triu(a, offset, status) result(b)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, optional, intent(in) :: offset !! Keep column minus row >= offset.
        integer, optional, intent(out) :: status !! Zero on success.
        type(csr_matrix) :: b
        integer :: kk

        kk = 0
        if (present(offset)) kk = offset
        b = csr_triangular(a,kk,.true.,status)
    end function csr_triu

    function csr_tril(a, offset, status) result(b)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, optional, intent(in) :: offset !! Keep column minus row <= offset.
        integer, optional, intent(out) :: status !! Zero on success.
        type(csr_matrix) :: b
        integer :: kk

        kk = 0
        if (present(offset)) kk = offset
        b = csr_triangular(a,kk,.false.,status)
    end function csr_tril

    function csr_triangular(a, offset, upper, status) result(b)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(in) :: offset !! Diagonal offset.
        logical, intent(in) :: upper !! True for upper, false for lower triangle.
        integer, optional, intent(out) :: status !! Zero on success.
        type(csr_matrix) :: b
        integer :: i, k, nused
        integer, allocatable :: cols(:)
        real(dp), allocatable :: vals(:)
        integer(int64) :: delta

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        allocate(cols(size(a%data)), vals(size(a%data)), b%indptr(a%n_rows+1))
        nused = 0
        do i = 1, a%n_rows
            b%indptr(i) = nused+1
            do k = a%indptr(i), a%indptr(i+1)-1
                delta = int(a%indices(k),int64)-int(i,int64)
                if (upper .and. delta < int(offset,int64)) cycle
                if (.not. upper .and. delta > int(offset,int64)) cycle
                nused = nused+1
                cols(nused) = a%indices(k)
                vals(nused) = a%data(k)
            end do
        end do
        b%indptr(a%n_rows+1) = nused+1
        b%n_rows = a%n_rows
        b%n_cols = a%n_cols
        b%indices = cols(:nused)
        b%data = vals(:nused)
        if (present(status)) status = sparse_success
    end function csr_triangular

    ! DIA: offsets are column minus row; data(k,j) is aligned to column j.
    ! Padded values outside the matrix are ignored in numerical operations.
    pure logical function valid_dia(a) result(ok)
        type(dia_matrix), intent(in) :: a !! Candidate one-based DIA matrix.
        integer :: k

        ok = .false.
        if (a%n_rows < 0 .or. a%n_cols < 0) return
        if (.not. allocated(a%offsets) .or. .not. allocated(a%data)) return
        if (size(a%data,1) /= size(a%offsets)) return
        if (size(a%data,2) /= a%n_cols) return
        if (.not. all(ieee_is_finite(a%data))) return
        do k = 2, size(a%offsets)
            if (a%offsets(k) <= a%offsets(k-1)) return
        end do
        ok = .true.
    end function valid_dia

    function dia_from_data(n_rows, n_cols, offsets, data, status) result(a)
        integer, intent(in) :: n_rows !! Matrix row count, nonnegative.
        integer, intent(in) :: n_cols !! Matrix column count, nonnegative.
        integer, intent(in) :: offsets(:) !! Unique column-minus-row offsets, any order.
        real(dp), intent(in) :: data(:,:) !! Values aligned by column, shape (number of diagonals,n_cols).
        integer, intent(out), optional :: status !! 0 on successful construction.
        type(dia_matrix) :: a
        integer :: k, j, p, temp
        integer, allocatable :: order(:)

        if (present(status)) status = sparse_invalid_input
        if (n_rows < 0 .or. n_cols < 0) return
        if (size(data,1) /= size(offsets) .or. size(data,2) /= n_cols) return
        if (.not. all(ieee_is_finite(data))) return
        allocate(order(size(offsets)))
        do k = 1, size(order)
            order(k) = k
        end do
        ! Insertion sort of diagonal labels (not individual matrix entries).
        do k = 2, size(order)
            temp = order(k)
            j = k-1
            do while (j >= 1)
                if (offsets(order(j)) <= offsets(temp)) exit
                order(j+1) = order(j)
                j = j-1
            end do
            order(j+1) = temp
        end do
        do k = 2, size(order)
            if (offsets(order(k)) == offsets(order(k-1))) return
        end do
        a%n_rows = n_rows
        a%n_cols = n_cols
        allocate(a%offsets(size(offsets)),a%data(size(offsets),n_cols))
        do p = 1,size(order)
            a%offsets(p) = offsets(order(p))
            a%data(p,:) = data(order(p),:)
        end do
        if (present(status)) status = sparse_success
    end function dia_from_data

    function csr_to_dia(a, status) result(d)
        type(csr_matrix), intent(in) :: a !! Canonical one-based CSR input.
        integer, intent(out), optional :: status !! 0 on success; rejects excessive DIA storage.
        type(dia_matrix) :: d
        integer, allocatable :: offsets(:)
        real(dp), allocatable :: dummy(:)
        integer :: i, k, j, nd, n, delta, lo, hi, mid, position

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        n = size(a%data)
        allocate(offsets(n),dummy(n))
        k = 0
        do i = 1,a%n_rows
            do j = a%indptr(i),a%indptr(i+1)-1
                k = k+1
                offsets(k) = a%indices(j)-i
            end do
        end do
        dummy = 0.0_dp
        if (n > 1) call sort_pairs(offsets,dummy)
        nd = 0
        do k = 1,n
            if (nd > 0) then
                if (offsets(k) == offsets(nd)) cycle
            end if
            nd = nd+1
            offsets(nd) = offsets(k)
        end do
        if (int(nd,int64)*int(a%n_cols,int64) > int(huge(0),int64)) return
        d%n_rows = a%n_rows
        d%n_cols = a%n_cols
        d%offsets = offsets(:nd)
        allocate(d%data(nd,a%n_cols))
        d%data = 0.0_dp
        do i = 1,a%n_rows
            do k = a%indptr(i),a%indptr(i+1)-1
                j = a%indices(k)
                delta = j-i
                lo = 1
                hi = nd
                do while (lo <= hi)
                    mid = lo+(hi-lo)/2
                    if (d%offsets(mid) < delta) then
                        lo = mid+1
                    else if (d%offsets(mid) > delta) then
                        hi = mid-1
                    else
                        exit
                    end if
                end do
                position = mid
                d%data(position,j) = a%data(k)
            end do
        end do
        if (present(status)) status = sparse_success
    end function csr_to_dia

    function dia_to_csr(a, status) result(c)
        type(dia_matrix), intent(in) :: a !! Valid canonical DIA input.
        integer, intent(out), optional :: status !! 0 on success.
        type(csr_matrix) :: c
        integer :: i, j, k, used, nnz, n
        integer(int64) :: cnt, jj

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_dia(a)) return
        n = size(a%offsets)
        cnt = 0_int64
        do k = 1,n
            cnt = cnt + max(0_int64, min(int(a%n_rows,int64), &
                int(a%n_cols,int64)-int(a%offsets(k),int64)) &
                -max(1_int64, 1_int64-int(a%offsets(k),int64))+1_int64)
        end do
        if (cnt >= int(huge(0),int64)) return
        nnz = int(cnt)
        c%n_rows = a%n_rows
        c%n_cols = a%n_cols
        allocate(c%indptr(a%n_rows+1),c%indices(nnz),c%data(nnz))
        used = 0
        do i = 1,a%n_rows
            c%indptr(i) = used+1
            do k = 1,n
                jj = int(i,int64)+int(a%offsets(k),int64)
                if (jj < 1_int64 .or. jj > int(a%n_cols,int64)) cycle
                j = int(jj)
                used = used+1
                c%indices(used) = j
                c%data(used) = a%data(k,j)
            end do
        end do
        c%indptr(a%n_rows+1) = used+1
        if (present(status)) status = sparse_success
    end function dia_to_csr

    function dia_from_dense(dense, status) result(a)
        real(dp), intent(in) :: dense(:,:) !! Finite dense matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(dia_matrix) :: a
        type(csr_matrix) :: c
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        c = csr_from_dense(dense,ierr)
        if (ierr /= sparse_success) return
        a = csr_to_dia(c,ierr)
        if (present(status)) status = ierr
    end function dia_from_dense

    function dia_to_dense(a, status) result(out)
        type(dia_matrix), intent(in) :: a !! Valid canonical DIA input.
        integer, intent(out), optional :: status !! 0 on success.
        real(dp), allocatable :: out(:,:)
        integer :: k, j
        integer(int64) :: i

        if (present(status)) status = sparse_invalid_input
        allocate(out(0,0))
        if (.not. valid_dia(a)) return
        deallocate(out)
        allocate(out(a%n_rows,a%n_cols))
        out = 0.0_dp
        do k = 1,size(a%offsets)
            do j = 1,a%n_cols
                i = int(j,int64)-int(a%offsets(k),int64)
                if (i < 1_int64 .or. i > int(a%n_rows,int64)) cycle
                out(int(i),j) = a%data(k,j)
            end do
        end do
        if (present(status)) status = sparse_success
    end function dia_to_dense

    function dia_matvec(a, x, status) result(y)
        type(dia_matrix), intent(in) :: a !! Valid canonical DIA input.
        real(dp), intent(in) :: x(:) !! Vector of length n_cols.
        integer, intent(out), optional :: status !! 0 on success.
        real(dp), allocatable :: y(:)
        integer :: k, j
        integer(int64) :: i

        if (present(status)) status = sparse_invalid_input
        allocate(y(0))
        if (.not. valid_dia(a)) return
        if (size(x) /= a%n_cols .or. .not. all(ieee_is_finite(x))) return
        deallocate(y)
        allocate(y(a%n_rows))
        y = 0.0_dp
        do k = 1,size(a%offsets)
            do j = 1,a%n_cols
                i = int(j,int64)-int(a%offsets(k),int64)
                if (i < 1_int64 .or. i > int(a%n_rows,int64)) cycle
                y(int(i)) = y(int(i))+a%data(k,j)*x(j)
            end do
        end do
        if (.not. all(ieee_is_finite(y))) then
            deallocate(y)
            allocate(y(0))
            return
        end if
        if (present(status)) status = sparse_success
    end function dia_matvec

    pure logical function valid_bsr(a) result(ok)
        type(bsr_matrix), intent(in) :: a !! Candidate BSR matrix with fixed-size dense blocks.
        integer :: i, k, nnzb, nr, nc

        ok = .false.
        if (a%n_rows < 0 .or. a%n_cols < 0) return
        if (a%n_rows == huge(0) .or. a%n_cols == huge(0)) return
        if (a%block_rows < 1 .or. a%block_cols < 1) return
        if (mod(a%n_rows,a%block_rows) /= 0 .or. mod(a%n_cols,a%block_cols) /= 0) return
        if (.not. allocated(a%indptr) .or. .not. allocated(a%indices) .or. .not. allocated(a%data)) return
        nr = a%n_rows/a%block_rows
        nc = a%n_cols/a%block_cols
        nnzb = size(a%indices)
        if (nnzb == huge(nnzb)) return
        if (size(a%indptr) /= nr+1 .or. size(a%data,3) /= nnzb) return
        if (size(a%data,1) /= a%block_rows .or. size(a%data,2) /= a%block_cols) return
        if (a%indptr(1) /= 1 .or. a%indptr(nr+1) /= nnzb+1) return
        if (any(a%indptr < 1) .or. any(a%indptr > nnzb+1)) return
        if (any(a%indptr(2:) < a%indptr(:nr))) return
        if (any(a%indices < 1) .or. any(a%indices > nc)) return
        if (.not. all(ieee_is_finite(a%data))) return
        do i = 1,nr
            do k = a%indptr(i)+1,a%indptr(i+1)-1
                if (a%indices(k) <= a%indices(k-1)) return
            end do
        end do
        ok = .true.
    end function valid_bsr

    function csr_to_bsr(a, block_rows, block_cols, status) result(b)
        type(csr_matrix), intent(in) :: a !! Valid canonical CSR matrix.
        integer, intent(in) :: block_rows !! Number of rows per block, >0.
        integer, intent(in) :: block_cols !! Number of columns per block, >0.
        integer, intent(out), optional :: status !! 0 on success.
        type(bsr_matrix) :: b
        integer, allocatable :: marked(:), touched(:), map(:)
        real(dp), allocatable :: sort_values(:)
        integer :: nr, nc, i, r, k, j, used, block, nactive, col
        integer(int64) :: count_blocks

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        if (block_rows < 1 .or. block_cols < 1) return
        if (mod(a%n_rows,block_rows) /= 0 .or. mod(a%n_cols,block_cols) /= 0) return
        nr = a%n_rows/block_rows
        nc = a%n_cols/block_cols
        allocate(marked(nc),map(nc),touched(nc),sort_values(nc))
        marked = 0
        map = 0
        count_blocks = 0_int64
        do i = 1,nr
            nactive = 0
            do r = (i-1)*block_rows+1, i*block_rows
                do k = a%indptr(r),a%indptr(r+1)-1
                    j = (a%indices(k)-1)/block_cols+1
                    if (marked(j) == i) cycle
                    marked(j) = i
                    nactive = nactive+1
                end do
            end do
            count_blocks = count_blocks+int(nactive,int64)
            if (count_blocks >= int(huge(0),int64)) return
        end do
        if (count_blocks > int(huge(0),int64)/int(block_rows,int64)) return
        if (count_blocks*int(block_rows,int64) > int(huge(0),int64)/int(block_cols,int64)) return
        b%n_rows = a%n_rows
        b%n_cols = a%n_cols
        b%block_rows = block_rows
        b%block_cols = block_cols
        allocate(b%indptr(nr+1),b%indices(int(count_blocks)), &
            b%data(block_rows,block_cols,int(count_blocks)))
        b%data = 0.0_dp
        used = 0
        do i = 1,nr
            b%indptr(i) = used+1
            nactive = 0
            do r = (i-1)*block_rows+1,i*block_rows
                do k = a%indptr(r),a%indptr(r+1)-1
                    j = (a%indices(k)-1)/block_cols+1
                    if (map(j) /= 0) cycle
                    nactive = nactive+1
                    touched(nactive) = j
                    map(j) = -1
                end do
            end do
            if (nactive > 1) then
                sort_values(:nactive) = 0.0_dp
                call sort_pairs(touched(:nactive),sort_values(:nactive))
            end if
            do j = 1,nactive
                block = touched(j)
                b%indices(used+j) = block
                map(block) = used+j
            end do
            do r = (i-1)*block_rows+1,i*block_rows
                do k = a%indptr(r),a%indptr(r+1)-1
                    col = a%indices(k)
                    block = map((col-1)/block_cols+1)
                    b%data(mod(r-1,block_rows)+1,mod(col-1,block_cols)+1,block) = a%data(k)
                end do
            end do
            do j = 1,nactive
                map(touched(j)) = 0
            end do
            used = used+nactive
        end do
        b%indptr(nr+1) = used+1
        if (present(status)) status = sparse_success
    end function csr_to_bsr

    function bsr_to_csr(a, status) result(c)
        type(bsr_matrix), intent(in) :: a !! Canonical BSR input.
        integer, intent(out), optional :: status !! 0 on success.
        type(csr_matrix) :: c
        integer :: i, p, block, nnz, used, base, ir, ic
        integer(int64) :: total

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_bsr(a)) return
        total = int(size(a%indices),int64)*int(a%block_rows,int64)*int(a%block_cols,int64)
        if (total >= int(huge(0),int64)) return
        nnz = int(total)
        c%n_rows = a%n_rows
        c%n_cols = a%n_cols
        allocate(c%indptr(a%n_rows+1),c%indices(nnz),c%data(nnz))
        used = 0
        do i = 1,a%n_rows
            c%indptr(i) = used+1
            block = (i-1)/a%block_rows+1
            ir = mod(i-1,a%block_rows)+1
            do p = a%indptr(block),a%indptr(block+1)-1
                base = (a%indices(p)-1)*a%block_cols
                do ic = 1,a%block_cols
                    used = used+1
                    c%indices(used) = base+ic
                    c%data(used) = a%data(ir,ic,p)
                end do
            end do
        end do
        c%indptr(a%n_rows+1) = used+1
        if (present(status)) status = sparse_success
    end function bsr_to_csr

    function bsr_from_dense(dense, block_rows, block_cols, status) result(a)
        real(dp), intent(in) :: dense(:,:) !! Finite dense matrix.
        integer, intent(in) :: block_rows !! Rows in each block, >0.
        integer, intent(in) :: block_cols !! Columns in each block, >0.
        integer, intent(out), optional :: status !! 0 on success.
        type(bsr_matrix) :: a
        type(csr_matrix) :: c
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        c = csr_from_dense(dense,ierr)
        if (ierr /= sparse_success) return
        a = csr_to_bsr(c,block_rows,block_cols,ierr)
        if (present(status)) status = ierr
    end function bsr_from_dense

    function bsr_to_dense(a, status) result(out)
        type(bsr_matrix), intent(in) :: a !! Canonical BSR matrix.
        integer, intent(out), optional :: status !! 0 on success.
        real(dp), allocatable :: out(:,:)
        integer :: i, k, row_start, col_start

        if (present(status)) status = sparse_invalid_input
        allocate(out(0,0))
        if (.not. valid_bsr(a)) return
        deallocate(out)
        allocate(out(a%n_rows,a%n_cols))
        out = 0.0_dp
        do i = 1,a%n_rows/a%block_rows
            row_start = (i-1)*a%block_rows+1
            do k = a%indptr(i),a%indptr(i+1)-1
                col_start = (a%indices(k)-1)*a%block_cols+1
                out(row_start:row_start+a%block_rows-1, &
                    col_start:col_start+a%block_cols-1) = a%data(:,:,k)
            end do
        end do
        if (present(status)) status = sparse_success
    end function bsr_to_dense

    function bsr_matvec(a, x, status) result(y)
        type(bsr_matrix), intent(in) :: a !! Canonical BSR matrix.
        real(dp), intent(in) :: x(:) !! RHS vector, n_cols elements.
        integer, intent(out), optional :: status !! 0 on success.
        real(dp), allocatable :: y(:)
        integer :: i, k, row_start, col_start

        if (present(status)) status = sparse_invalid_input
        allocate(y(0))
        if (.not. valid_bsr(a)) return
        if (size(x) /= a%n_cols .or. .not. all(ieee_is_finite(x))) return
        deallocate(y)
        allocate(y(a%n_rows))
        y = 0.0_dp
        do i = 1,a%n_rows/a%block_rows
            row_start = (i-1)*a%block_rows+1
            do k = a%indptr(i),a%indptr(i+1)-1
                col_start = (a%indices(k)-1)*a%block_cols+1
                y(row_start:row_start+a%block_rows-1) = &
                    y(row_start:row_start+a%block_rows-1) + &
                    matmul(a%data(:,:,k),x(col_start:col_start+a%block_cols-1))
            end do
        end do
        if (.not. all(ieee_is_finite(y))) then
            deallocate(y)
            allocate(y(0))
            return
        end if
        if (present(status)) status = sparse_success
    end function bsr_matvec

    ! In-place heapsort of (column,value) pairs.  Unlike sorting a triplet
    ! list, only one CSR row's columns need sorting at any instant.
    subroutine sort_pairs(keys, values)
        integer, intent(inout) :: keys(:) !! Column keys, sorted increasing on return.
        real(dp), intent(inout) :: values(:) !! Values moved with their corresponding key.
        integer :: start, finish, key
        real(dp) :: value

        do start = size(keys)/2, 1, -1
            call sift_pair(keys, values, start, size(keys))
        end do
        do finish = size(keys), 2, -1
            key = keys(1)
            keys(1) = keys(finish)
            keys(finish) = key
            value = values(1)
            values(1) = values(finish)
            values(finish) = value
            call sift_pair(keys, values, 1, finish-1)
        end do
    end subroutine sort_pairs

    subroutine sift_pair(keys, values, root, last)
        integer, intent(inout) :: keys(:) !! Heap keys.
        real(dp), intent(inout) :: values(:) !! Aligned heap values.
        integer, intent(in) :: root !! Initial sift-down position.
        integer, intent(in) :: last !! Last position in heap.
        integer :: i, child, key
        real(dp) :: value

        i = root
        do while (i <= last/2)
            child = 2*i
            if (child < last) then
                if (keys(child) < keys(child+1)) child = child+1
            end if
            if (keys(i) >= keys(child)) exit
            key = keys(i)
            keys(i) = keys(child)
            keys(child) = key
            value = values(i)
            values(i) = values(child)
            values(child) = value
            i = child
        end do
    end subroutine sift_pair

    subroutine sort_integers(keys)
        integer, intent(inout) :: keys(:) !! Sorted increasingly in place.
        integer :: start, finish, key

        do start = size(keys)/2, 1, -1
            call sift_integer(keys, start, size(keys))
        end do
        do finish = size(keys), 2, -1
            key = keys(1)
            keys(1) = keys(finish)
            keys(finish) = key
            call sift_integer(keys, 1, finish-1)
        end do
    end subroutine sort_integers

    subroutine sift_integer(keys, root, last)
        integer, intent(inout) :: keys(:) !! Integer heap.
        integer, intent(in) :: root !! Initial sift-down position.
        integer, intent(in) :: last !! Last heap position.
        integer :: i, child, key

        i = root
        do while (i <= last/2)
            child = 2*i
            if (child < last) then
                if (keys(child) < keys(child+1)) child = child+1
            end if
            if (keys(i) >= keys(child)) exit
            key = keys(i)
            keys(i) = keys(child)
            keys(child) = key
            i = child
        end do
    end subroutine sift_integer

end module scifort_sparse
