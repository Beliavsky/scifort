! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Complex(dp) canonical CSR and CSC matrices, one-based pointers and indices.
! Stored indices are strictly ascending within each row/column. Explicit
! zeros are retained, while duplicate COO coordinates are summed.
module scifort_sparse_complex
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input, csr_matrix, valid_csr
    use scifort_sparse_coo_nd_complex, only : coo_nd_complex_array, &
        valid_coo_nd_complex, coo_nd_complex_from_triplets, coo_nd_complex_sum_duplicates
    implicit none
    private

    public :: csr_complex_matrix, csc_complex_matrix
    public :: valid_csr_complex, valid_csc_complex
    public :: csr_complex_from_coo, csr_complex_from_coo_nd, csr_complex_to_coo_nd
    public :: csr_complex_from_dense, csr_complex_to_dense
    public :: csc_complex_from_coo, csc_complex_from_dense, csc_complex_to_dense
    public :: csr_complex_to_csc, csc_complex_to_csr
    public :: csr_complex_transpose, csr_complex_hermitian, csr_complex_conjugate
    public :: csc_complex_transpose, csc_complex_hermitian
    public :: csr_complex_matvec, csc_complex_matvec, csr_complex_matmat
    public :: csr_complex_from_real

    type :: csr_complex_matrix
        integer :: n_rows = 0
        integer :: n_cols = 0
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
        complex(dp), allocatable :: data(:)
    end type csr_complex_matrix

    type :: csc_complex_matrix
        integer :: n_rows = 0
        integer :: n_cols = 0
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
        complex(dp), allocatable :: data(:)
    end type csc_complex_matrix

contains

    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Complex number to check.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    pure logical function valid_csr_complex(a) result(ok)
        type(csr_complex_matrix), intent(in) :: a !! Candidate canonical CSR matrix.
        integer :: i, k, nnz

        ok = .false.
        if (a%n_rows < 0 .or. a%n_cols < 0 .or. a%n_rows == huge(0)) return
        if (.not. allocated(a%indptr) .or. .not. allocated(a%indices) .or. .not. allocated(a%data)) return
        nnz = size(a%data)
        if (nnz == huge(nnz)) return
        if (size(a%indptr) /= a%n_rows+1 .or. size(a%indices) /= nnz) return
        if (a%indptr(1) /= 1 .or. a%indptr(a%n_rows+1) /= nnz+1) return
        if (any(a%indptr < 1) .or. any(a%indptr > nnz+1)) return
        if (any(a%indptr(2:) < a%indptr(:a%n_rows))) return
        if (any(a%indices < 1) .or. any(a%indices > a%n_cols)) return
        if (.not. all(finite_complex(a%data))) return
        do i = 1, a%n_rows
            do k = a%indptr(i)+1, a%indptr(i+1)-1
                if (a%indices(k) <= a%indices(k-1)) return
            end do
        end do
        ok = .true.
    end function valid_csr_complex

    pure logical function valid_csc_complex(a) result(ok)
        type(csc_complex_matrix), intent(in) :: a !! Candidate canonical CSC matrix.
        integer :: i, k, nnz

        ok = .false.
        if (a%n_rows < 0 .or. a%n_cols < 0 .or. a%n_cols == huge(0)) return
        if (.not. allocated(a%indptr) .or. .not. allocated(a%indices) .or. .not. allocated(a%data)) return
        nnz = size(a%data)
        if (nnz == huge(nnz)) return
        if (size(a%indptr) /= a%n_cols+1 .or. size(a%indices) /= nnz) return
        if (a%indptr(1) /= 1 .or. a%indptr(a%n_cols+1) /= nnz+1) return
        if (any(a%indptr < 1) .or. any(a%indptr > nnz+1)) return
        if (any(a%indptr(2:) < a%indptr(:a%n_cols))) return
        if (any(a%indices < 1) .or. any(a%indices > a%n_rows)) return
        if (.not. all(finite_complex(a%data))) return
        do i = 1, a%n_cols
            do k = a%indptr(i)+1, a%indptr(i+1)-1
                if (a%indices(k) <= a%indices(k-1)) return
            end do
        end do
        ok = .true.
    end function valid_csc_complex

    function csr_complex_from_coo(n_rows, n_cols, rows, columns, values, status) result(a)
        integer, intent(in) :: n_rows !! Nonnegative number of rows.
        integer, intent(in) :: n_cols !! Nonnegative number of columns.
        integer, intent(in) :: rows(:) !! One-based COO row coordinates.
        integer, intent(in) :: columns(:) !! One-based COO column coordinates.
        complex(dp), intent(in) :: values(:) !! COO data; duplicate coordinates summed.
        integer, intent(out), optional :: status !! 0 on success, 1 on bad input or overflow.
        type(csr_complex_matrix) :: a
        type(coo_nd_complex_array) :: tmp
        integer, allocatable :: coordinates(:,:)
        integer :: n, i, p, info

        if (present(status)) status = sparse_invalid_input
        if (n_rows < 0 .or. n_cols < 0 .or. n_rows == huge(0)) return
        n = size(values)
        if (n == huge(0) .or. size(rows) /= n .or. size(columns) /= n) return
        allocate(coordinates(2,n))
        coordinates(1,:) = rows
        coordinates(2,:) = columns
        tmp = coo_nd_complex_from_triplets([n_rows,n_cols],coordinates,values,info)
        if (info /= sparse_success) return
        call coo_nd_complex_sum_duplicates(tmp,info)
        if (info /= sparse_success) return
        n = size(tmp%data)
        a%n_rows = n_rows
        a%n_cols = n_cols
        allocate(a%indptr(n_rows+1), a%indices(n), a%data(n))
        a%indices = tmp%coords(2,:)
        a%data = tmp%data
        p = 1
        do i = 1, n_rows
            a%indptr(i) = p
            do while (p <= n)
                if (tmp%coords(1,p) /= i) exit
                p = p+1
            end do
        end do
        a%indptr(n_rows+1) = n+1
        if (present(status)) status = sparse_success
    end function csr_complex_from_coo

    function csr_complex_from_coo_nd(a, status) result(b)
        type(coo_nd_complex_array), intent(in) :: a !! Rank-two complex COO input.
        integer, intent(out), optional :: status !! 0 on success.
        type(csr_complex_matrix) :: b

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo_nd_complex(a)) return
        if (size(a%shape) /= 2) return
        b = csr_complex_from_coo(a%shape(1),a%shape(2),a%coords(1,:),a%coords(2,:),a%data,status)
    end function csr_complex_from_coo_nd

    function csr_complex_to_coo_nd(a, status) result(b)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR input.
        integer, intent(out), optional :: status !! 0 on success.
        type(coo_nd_complex_array) :: b
        integer, allocatable :: coords(:,:)
        integer :: i, k

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr_complex(a)) return
        allocate(coords(2,size(a%data)))
        coords(2,:) = a%indices
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                coords(1,k) = i
            end do
        end do
        b = coo_nd_complex_from_triplets([a%n_rows,a%n_cols],coords,a%data,status)
    end function csr_complex_to_coo_nd

    function csr_complex_from_dense(dense, status) result(a)
        complex(dp), intent(in) :: dense(:,:) !! Finite complex dense matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csr_complex_matrix) :: a
        integer :: m, n, nnz, i, j, p

        if (present(status)) status = sparse_invalid_input
        m = size(dense,1)
        n = size(dense,2)
        if (m == huge(0)) return
        if (.not. all(finite_complex(dense))) return
        if (count(dense /= cmplx(0.0_dp,0.0_dp,dp),kind=int64) >= int(huge(0),int64)) return
        nnz = count(dense /= cmplx(0.0_dp,0.0_dp,dp))
        a%n_rows = m
        a%n_cols = n
        allocate(a%indptr(m+1),a%indices(nnz),a%data(nnz))
        p = 1
        do i = 1, m
            a%indptr(i) = p
            do j = 1, n
                if (dense(i,j) == cmplx(0.0_dp,0.0_dp,dp)) cycle
                a%indices(p) = j
                a%data(p) = dense(i,j)
                p = p+1
            end do
        end do
        a%indptr(m+1) = p
        if (present(status)) status = sparse_success
    end function csr_complex_from_dense

    function csr_complex_to_dense(a, status, max_elements) result(dense)
        type(csr_complex_matrix), intent(in) :: a !! Valid canonical CSR input.
        integer, intent(out), optional :: status !! 0 on success.
        integer(int64), intent(in), optional :: max_elements !! Optional dense allocation limit.
        complex(dp), allocatable :: dense(:,:)
        integer :: i, k

        if (present(status)) status = sparse_invalid_input
        allocate(dense(0,0))
        if (.not. valid_csr_complex(a)) return
        if (present(max_elements)) then
            if (max_elements < 0_int64) return
            if (int(a%n_rows,int64)*int(a%n_cols,int64) > max_elements) return
        end if
        deallocate(dense)
        allocate(dense(a%n_rows,a%n_cols))
        dense = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                dense(i,a%indices(k)) = a%data(k)
            end do
        end do
        if (present(status)) status = sparse_success
    end function csr_complex_to_dense

    function csr_complex_to_csc(a, status) result(b)
        type(csr_complex_matrix), intent(in) :: a !! Canonical complex CSR matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csc_complex_matrix) :: b
        integer, allocatable :: cursor(:)
        integer :: i, j, k, p, nnz

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr_complex(a)) return
        if (a%n_cols == huge(0)) return
        nnz = size(a%data)
        b%n_rows = a%n_rows
        b%n_cols = a%n_cols
        allocate(b%indptr(a%n_cols+1),b%indices(nnz),b%data(nnz),cursor(a%n_cols))
        b%indptr = 0
        do k = 1, nnz
            j = a%indices(k)
            b%indptr(j+1) = b%indptr(j+1)+1
        end do
        b%indptr(1) = 1
        do j = 1, a%n_cols
            b%indptr(j+1) = b%indptr(j+1)+b%indptr(j)
        end do
        cursor = b%indptr(:a%n_cols)
        do i = 1, a%n_rows
            do k = a%indptr(i),a%indptr(i+1)-1
                j = a%indices(k)
                p = cursor(j)
                b%indices(p) = i
                b%data(p) = a%data(k)
                cursor(j) = p+1
            end do
        end do
        if (present(status)) status = sparse_success
    end function csr_complex_to_csc

    function csc_complex_to_csr(a, status) result(b)
        type(csc_complex_matrix), intent(in) :: a !! Canonical complex CSC matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csr_complex_matrix) :: b
        integer, allocatable :: cursor(:)
        integer :: i, j, k, p, nnz

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csc_complex(a)) return
        if (a%n_rows == huge(0)) return
        nnz = size(a%data)
        b%n_rows = a%n_rows
        b%n_cols = a%n_cols
        allocate(b%indptr(a%n_rows+1),b%indices(nnz),b%data(nnz),cursor(a%n_rows))
        b%indptr = 0
        do k = 1, nnz
            i = a%indices(k)
            b%indptr(i+1) = b%indptr(i+1)+1
        end do
        b%indptr(1) = 1
        do i = 1, a%n_rows
            b%indptr(i+1) = b%indptr(i+1)+b%indptr(i)
        end do
        cursor = b%indptr(:a%n_rows)
        do j = 1, a%n_cols
            do k = a%indptr(j),a%indptr(j+1)-1
                i = a%indices(k)
                p = cursor(i)
                b%indices(p) = j
                b%data(p) = a%data(k)
                cursor(i) = p+1
            end do
        end do
        if (present(status)) status = sparse_success
    end function csc_complex_to_csr

    function csc_complex_from_coo(n_rows, n_cols, rows, columns, values, status) result(b)
        integer, intent(in) :: n_rows !! Number of rows.
        integer, intent(in) :: n_cols !! Number of columns.
        integer, intent(in) :: rows(:) !! One-based row coordinates.
        integer, intent(in) :: columns(:) !! One-based column coordinates.
        complex(dp), intent(in) :: values(:) !! Values; duplicates summed.
        integer, intent(out), optional :: status !! 0 on success.
        type(csc_complex_matrix) :: b
        type(csr_complex_matrix) :: a
        integer :: info

        if (present(status)) status = sparse_invalid_input
        a = csr_complex_from_coo(n_rows,n_cols,rows,columns,values,info)
        if (info /= sparse_success) return
        b = csr_complex_to_csc(a,status)
    end function csc_complex_from_coo

    function csc_complex_from_dense(dense, status) result(b)
        complex(dp), intent(in) :: dense(:,:) !! Finite dense matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csc_complex_matrix) :: b
        type(csr_complex_matrix) :: a
        integer :: info

        if (present(status)) status = sparse_invalid_input
        a = csr_complex_from_dense(dense,info)
        if (info /= sparse_success) return
        b = csr_complex_to_csc(a,status)
    end function csc_complex_from_dense

    function csc_complex_to_dense(a, status, max_elements) result(dense)
        type(csc_complex_matrix), intent(in) :: a !! Canonical CSC matrix.
        integer, intent(out), optional :: status !! 0 on success.
        integer(int64), intent(in), optional :: max_elements !! Optional allocation cap.
        complex(dp), allocatable :: dense(:,:)
        integer :: j, k

        if (present(status)) status = sparse_invalid_input
        allocate(dense(0,0))
        if (.not. valid_csc_complex(a)) return
        if (present(max_elements)) then
            if (max_elements < 0_int64) return
            if (int(a%n_rows,int64)*int(a%n_cols,int64) > max_elements) return
        end if
        deallocate(dense)
        allocate(dense(a%n_rows,a%n_cols))
        dense = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1, a%n_cols
            do k = a%indptr(j),a%indptr(j+1)-1
                dense(a%indices(k),j) = a%data(k)
            end do
        end do
        if (present(status)) status = sparse_success
    end function csc_complex_to_dense

    function csr_complex_transpose(a, status) result(t)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csr_complex_matrix) :: t
        type(csc_complex_matrix) :: b
        integer :: info

        if (present(status)) status = sparse_invalid_input
        b = csr_complex_to_csc(a,info)
        if (info /= sparse_success) return
        t%n_rows = a%n_cols
        t%n_cols = a%n_rows
        t%indptr = b%indptr
        t%indices = b%indices
        t%data = b%data
        if (present(status)) status = sparse_success
    end function csr_complex_transpose

    function csr_complex_hermitian(a, status) result(t)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csr_complex_matrix) :: t
        integer :: info

        t = csr_complex_transpose(a,info)
        if (present(status)) status = info
        if (info == sparse_success) t%data = conjg(t%data)
    end function csr_complex_hermitian

    function csr_complex_conjugate(a, status) result(t)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csr_complex_matrix) :: t

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr_complex(a)) return
        t = a
        t%data = conjg(t%data)
        if (present(status)) status = sparse_success
    end function csr_complex_conjugate

    function csc_complex_transpose(a, status) result(t)
        type(csc_complex_matrix), intent(in) :: a !! Canonical CSC matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csc_complex_matrix) :: t
        type(csr_complex_matrix) :: b
        integer :: info

        if (present(status)) status = sparse_invalid_input
        b = csc_complex_to_csr(a,info)
        if (info /= sparse_success) return
        t%n_rows = a%n_cols
        t%n_cols = a%n_rows
        t%indptr = b%indptr
        t%indices = b%indices
        t%data = b%data
        if (present(status)) status = sparse_success
    end function csc_complex_transpose

    function csc_complex_hermitian(a, status) result(t)
        type(csc_complex_matrix), intent(in) :: a !! Canonical CSC matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csc_complex_matrix) :: t
        integer :: info

        t = csc_complex_transpose(a,info)
        if (info == sparse_success) t%data = conjg(t%data)
        if (present(status)) status = info
    end function csc_complex_hermitian

    function csr_complex_matvec(a, x, status) result(y)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR coefficient matrix.
        complex(dp), intent(in) :: x(:) !! Finite vector of length n_cols.
        integer, intent(out), optional :: status !! 0 on success.
        complex(dp), allocatable :: y(:)
        complex(dp) :: acc
        integer :: i, k

        if (present(status)) status = sparse_invalid_input
        allocate(y(0))
        if (.not. valid_csr_complex(a)) return
        if (size(x) /= a%n_cols) return
        if (.not. all(finite_complex(x))) return
        deallocate(y)
        allocate(y(a%n_rows))
        y = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1, a%n_rows
            acc = cmplx(0.0_dp,0.0_dp,dp)
            do k = a%indptr(i),a%indptr(i+1)-1
                acc = acc + a%data(k)*x(a%indices(k))
                if (.not. finite_complex(acc)) then
                    deallocate(y)
                    allocate(y(0))
                    return
                end if
            end do
            y(i) = acc
        end do
        if (present(status)) status = sparse_success
    end function csr_complex_matvec

    function csc_complex_matvec(a, x, status) result(y)
        type(csc_complex_matrix), intent(in) :: a !! Canonical CSC matrix.
        complex(dp), intent(in) :: x(:) !! Vector of length n_cols.
        integer, intent(out), optional :: status !! 0 on success.
        complex(dp), allocatable :: y(:)
        integer :: j, k, i

        if (present(status)) status = sparse_invalid_input
        allocate(y(0))
        if (.not. valid_csc_complex(a)) return
        if (size(x) /= a%n_cols) return
        if (.not. all(finite_complex(x))) return
        deallocate(y)
        allocate(y(a%n_rows))
        y = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1, a%n_cols
            do k = a%indptr(j),a%indptr(j+1)-1
                i = a%indices(k)
                y(i) = y(i) + a%data(k)*x(j)
                if (.not. finite_complex(y(i))) then
                    deallocate(y)
                    allocate(y(0))
                    return
                end if
            end do
        end do
        if (present(status)) status = sparse_success
    end function csc_complex_matvec

    function csr_complex_matmat(a, x, status) result(y)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR matrix.
        complex(dp), intent(in) :: x(:,:) !! Dense right factor (n_cols,k).
        integer, intent(out), optional :: status !! 0 on success.
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: col(:)
        integer :: j, info

        if (present(status)) status = sparse_invalid_input
        allocate(y(0,0))
        if (.not. valid_csr_complex(a)) return
        if (size(x,1) /= a%n_cols) return
        if (.not. all(finite_complex(x))) return
        deallocate(y)
        allocate(y(a%n_rows,size(x,2)))
        do j = 1, size(x,2)
            col = csr_complex_matvec(a,x(:,j),info)
            if (info /= sparse_success) then
                deallocate(y)
                allocate(y(0,0))
                return
            end if
            y(:,j) = col
        end do
        if (present(status)) status = sparse_success
    end function csr_complex_matmat

    function csr_complex_from_real(a, status) result(b)
        type(csr_matrix), intent(in) :: a !! Canonical real CSR matrix.
        integer, intent(out), optional :: status !! 0 on success.
        type(csr_complex_matrix) :: b

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(a)) return
        b%n_rows = a%n_rows
        b%n_cols = a%n_cols
        b%indptr = a%indptr
        b%indices = a%indices
        b%data = cmplx(a%data,0.0_dp,dp)
        if (present(status)) status = sparse_success
    end function csr_complex_from_real

end module scifort_sparse_complex
