! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Two-dimensional coordinate (COO) sparse matrices with one-based indices.
! Like scipy.sparse.coo_array, COO preserves unsorted coordinates, duplicate
! entries, and explicitly stored zeros until explicitly canonicalized.
module scifort_sparse_coo
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, valid_csr, valid_csc, &
        csr_from_coo, csr_to_coo, csc_from_coo, csc_to_csr, &
        sparse_success, sparse_invalid_input
    implicit none
    private

    public :: coo_matrix, valid_coo, coo_from_triplets, coo_from_dense
    public :: coo_from_csr, coo_from_csc, coo_to_csr, coo_to_csc
    public :: coo_to_dense, coo_matvec, coo_transpose
    public :: coo_sum_duplicates, coo_eliminate_zeros, coo_get
    public :: sparse_find, coo_find, csr_find, csc_find

    type :: coo_matrix
        integer :: n_rows = 0
        integer :: n_cols = 0
        integer, allocatable :: row(:) !! One-based row for each stored entry.
        integer, allocatable :: col(:) !! One-based column for each stored entry.
        real(dp), allocatable :: data(:) !! Values; repeats and zeros are allowed.
    end type coo_matrix

    interface sparse_find
        module procedure coo_find
        module procedure csr_find
        module procedure csc_find
    end interface sparse_find

contains

    pure logical function valid_coo(a) result(ok)
        type(coo_matrix), intent(in) :: a !! Candidate coordinate matrix.
        integer :: nnz

        ok = .false.
        if (a%n_rows < 0 .or. a%n_cols < 0) return
        if (a%n_rows == huge(0) .or. a%n_cols == huge(0)) return
        if (.not. allocated(a%row) .or. .not. allocated(a%col) .or. &
            .not. allocated(a%data)) return
        nnz = size(a%data)
        if (nnz == huge(nnz)) return
        if (size(a%row) /= nnz .or. size(a%col) /= nnz) return
        if (any(a%row < 1) .or. any(a%row > a%n_rows)) return
        if (any(a%col < 1) .or. any(a%col > a%n_cols)) return
        if (.not. all(ieee_is_finite(a%data))) return
        ok = .true.
    end function valid_coo

    function coo_from_triplets(n_rows, n_cols, rows, columns, values, status) result(a)
        integer, intent(in) :: n_rows !! Number of rows, nonnegative.
        integer, intent(in) :: n_cols !! Number of columns, nonnegative.
        integer, intent(in) :: rows(:) !! One-based coordinates, unsorted permitted.
        integer, intent(in) :: columns(:) !! One-based coordinates, duplicates permitted.
        real(dp), intent(in) :: values(:) !! Finite data; explicit zeros retained.
        integer, intent(out), optional :: status !! Zero on success; invalid input otherwise.
        type(coo_matrix) :: a

        if (present(status)) status = sparse_invalid_input
        if (n_rows < 0 .or. n_cols < 0) return
        if (n_rows == huge(0) .or. n_cols == huge(0)) return
        if (size(rows) /= size(columns) .or. size(rows) /= size(values)) return
        if (size(rows) == huge(0)) return
        if (any(rows < 1) .or. any(rows > n_rows)) return
        if (any(columns < 1) .or. any(columns > n_cols)) return
        if (.not. all(ieee_is_finite(values))) return
        a%n_rows = n_rows
        a%n_cols = n_cols
        a%row = rows
        a%col = columns
        a%data = values
        if (present(status)) status = sparse_success
    end function coo_from_triplets

    function coo_from_dense(dense, status) result(a)
        real(dp), intent(in) :: dense(:,:) !! Finite matrix; zeros are implicit.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_matrix) :: a
        integer :: i, j, k, nnz

        if (present(status)) status = sparse_invalid_input
        if (.not. all(ieee_is_finite(dense))) return
        if (size(dense,1) == huge(0) .or. size(dense,2) == huge(0)) return
        if (count(dense /= 0.0_dp, kind=int64) >= int(huge(0),int64)) return
        nnz = count(dense /= 0.0_dp)
        a%n_rows = size(dense,1)
        a%n_cols = size(dense,2)
        allocate(a%row(nnz), a%col(nnz), a%data(nnz))
        k = 0
        do i = 1, a%n_rows
            do j = 1, a%n_cols
                if (dense(i,j) == 0.0_dp) cycle
                k = k+1
                a%row(k) = i
                a%col(k) = j
                a%data(k) = dense(i,j)
            end do
        end do
        if (present(status)) status = sparse_success
    end function coo_from_dense

    function coo_from_csr(a, status) result(b)
        type(csr_matrix), intent(in) :: a !! Valid, canonical CSR input.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_matrix) :: b
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        call csr_to_coo(a,b%row,b%col,b%data,ierr)
        if (ierr /= sparse_success) return
        b%n_rows = a%n_rows
        b%n_cols = a%n_cols
        if (present(status)) status = sparse_success
    end function coo_from_csr

    function coo_from_csc(a, status) result(b)
        type(csc_matrix), intent(in) :: a !! Valid, canonical CSC input.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_matrix) :: b
        type(csr_matrix) :: work
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csc(a)) return
        work = csc_to_csr(a,ierr)
        if (ierr /= sparse_success) return
        b = coo_from_csr(work,ierr)
        if (present(status)) status = ierr
    end function coo_from_csc

    function coo_to_csr(a, status) result(b)
        type(coo_matrix), intent(in) :: a !! Valid COO; duplicates are summed.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: b
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo(a)) return
        b = csr_from_coo(a%n_rows,a%n_cols,a%row,a%col,a%data,ierr)
        if (present(status)) status = ierr
    end function coo_to_csr

    function coo_to_csc(a, status) result(b)
        type(coo_matrix), intent(in) :: a !! Valid COO; duplicates are summed.
        integer, intent(out), optional :: status !! Zero on success.
        type(csc_matrix) :: b
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo(a)) return
        b = csc_from_coo(a%n_rows,a%n_cols,a%row,a%col,a%data,ierr)
        if (present(status)) status = ierr
    end function coo_to_csc

    function coo_to_dense(a, status) result(dense)
        type(coo_matrix), intent(in) :: a !! Valid COO, including duplicate coordinates.
        integer, intent(out), optional :: status !! Zero on success; overflow rejected.
        real(dp), allocatable :: dense(:,:)
        integer :: k, i, j

        if (present(status)) status = sparse_invalid_input
        allocate(dense(0,0))
        if (.not. valid_coo(a)) return
        deallocate(dense)
        allocate(dense(a%n_rows,a%n_cols))
        dense = 0.0_dp
        do k = 1, size(a%data)
            i = a%row(k)
            j = a%col(k)
            dense(i,j) = dense(i,j) + a%data(k)
            if (.not. ieee_is_finite(dense(i,j))) then
                deallocate(dense)
                allocate(dense(0,0))
                return
            end if
        end do
        if (present(status)) status = sparse_success
    end function coo_to_dense

    function coo_matvec(a, x, status) result(y)
        type(coo_matrix), intent(in) :: a !! Valid COO input, duplicates supported.
        real(dp), intent(in) :: x(:) !! Finite vector of size n_cols.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: y(:)
        integer :: k, i

        if (present(status)) status = sparse_invalid_input
        allocate(y(0))
        if (.not. valid_coo(a)) return
        if (size(x) /= a%n_cols .or. .not. all(ieee_is_finite(x))) return
        deallocate(y)
        allocate(y(a%n_rows))
        y = 0.0_dp
        do k = 1, size(a%data)
            i = a%row(k)
            y(i) = y(i) + a%data(k)*x(a%col(k))
            if (.not. ieee_is_finite(y(i))) then
                deallocate(y)
                allocate(y(0))
                return
            end if
        end do
        if (present(status)) status = sparse_success
    end function coo_matvec

    function coo_transpose(a, status) result(b)
        type(coo_matrix), intent(in) :: a !! Valid COO input.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_matrix) :: b

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo(a)) return
        b%n_rows = a%n_cols
        b%n_cols = a%n_rows
        b%row = a%col
        b%col = a%row
        b%data = a%data
        if (present(status)) status = sparse_success
    end function coo_transpose

    subroutine coo_sum_duplicates(a, status)
        type(coo_matrix), intent(inout) :: a !! COO made canonical on success.
        integer, intent(out), optional :: status !! Zero on success; input unchanged on error.
        type(csr_matrix) :: b
        type(coo_matrix) :: canonical
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo(a)) return
        b = coo_to_csr(a,ierr)
        if (ierr /= sparse_success) return
        canonical = coo_from_csr(b,ierr)
        if (ierr /= sparse_success) return
        a = canonical
        if (present(status)) status = sparse_success
    end subroutine coo_sum_duplicates

    subroutine coo_eliminate_zeros(a, status)
        type(coo_matrix), intent(inout) :: a !! Remove stored zero entries only.
        integer, intent(out), optional :: status !! Zero on success.
        logical, allocatable :: keep(:)

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo(a)) return
        keep = a%data /= 0.0_dp
        a%row = pack(a%row,keep)
        a%col = pack(a%col,keep)
        a%data = pack(a%data,keep)
        if (present(status)) status = sparse_success
    end subroutine coo_eliminate_zeros

    real(dp) function coo_get(a, row, column, status) result(value)
        type(coo_matrix), intent(in) :: a !! Valid COO, including duplicate coordinates.
        integer, intent(in) :: row !! One-based row index.
        integer, intent(in) :: column !! One-based column index.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: k

        value = 0.0_dp
        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo(a)) return
        if (row < 1 .or. row > a%n_rows .or. column < 1 .or. column > a%n_cols) return
        do k = 1, size(a%data)
            if (a%row(k) /= row .or. a%col(k) /= column) cycle
            value = value + a%data(k)
            if (.not. ieee_is_finite(value)) then
                value = 0.0_dp
                return
            end if
        end do
        if (present(status)) status = sparse_success
    end function coo_get

    subroutine coo_find(a, rows, columns, values, status)
        type(coo_matrix), intent(in) :: a !! COO array, possibly with duplicates/zeros.
        integer, allocatable, intent(out) :: rows(:) !! Sorted one-based nonzero rows.
        integer, allocatable, intent(out) :: columns(:) !! Sorted one-based nonzero columns.
        real(dp), allocatable, intent(out) :: values(:) !! Summed nonzero values.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: b
        type(coo_matrix) :: canonical
        integer :: ierr
        logical, allocatable :: keep(:)

        if (present(status)) status = sparse_invalid_input
        allocate(rows(0),columns(0),values(0))
        b = coo_to_csr(a,ierr)
        if (ierr /= sparse_success) return
        canonical = coo_from_csr(b,ierr)
        if (ierr /= sparse_success) return
        keep = canonical%data /= 0.0_dp
        rows = pack(canonical%row,keep)
        columns = pack(canonical%col,keep)
        values = pack(canonical%data,keep)
        if (present(status)) status = sparse_success
    end subroutine coo_find

    subroutine csr_find(a, rows, columns, values, status)
        type(csr_matrix), intent(in) :: a !! CSR matrix.
        integer, allocatable, intent(out) :: rows(:) !! One-based nonzero rows.
        integer, allocatable, intent(out) :: columns(:) !! One-based nonzero columns.
        real(dp), allocatable, intent(out) :: values(:) !! Nonzero values.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_matrix) :: b
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        allocate(rows(0),columns(0),values(0))
        b = coo_from_csr(a,ierr)
        if (ierr /= sparse_success) return
        call coo_find(b,rows,columns,values,ierr)
        if (present(status)) status = ierr
    end subroutine csr_find

    subroutine csc_find(a, rows, columns, values, status)
        type(csc_matrix), intent(in) :: a !! CSC matrix.
        integer, allocatable, intent(out) :: rows(:) !! One-based nonzero rows.
        integer, allocatable, intent(out) :: columns(:) !! One-based nonzero columns.
        real(dp), allocatable, intent(out) :: values(:) !! Nonzero values.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_matrix) :: b
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        allocate(rows(0),columns(0),values(0))
        b = coo_from_csc(a,ierr)
        if (ierr /= sparse_success) return
        call coo_find(b,rows,columns,values,ierr)
        if (present(status)) status = ierr
    end subroutine csc_find

end module scifort_sparse_coo
