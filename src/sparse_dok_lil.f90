! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Incrementally editable, one-based sparse DOK and LIL storage.
! DOK uses open addressing with tombstones; LIL keeps sorted columns in
! dynamically sized row buffers.  Assignment of zero removes a stored entry,
! while conversion from CSR/COO retains explicitly stored zeros.
module scifort_sparse_dok_lil
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, valid_csr, &
        sparse_success, sparse_invalid_input, csr_from_dense
    use scifort_sparse_coo, only : coo_matrix, valid_coo, coo_from_triplets, &
        coo_from_csr, coo_from_csc, coo_to_csr, coo_to_csc, coo_transpose
    implicit none
    private

    public :: dok_matrix, lil_matrix, valid_dok, valid_lil
    public :: dok_empty, lil_empty, dok_nnz, lil_nnz
    public :: dok_set, lil_set, dok_get, lil_get, dok_add, lil_add
    public :: dok_delete, lil_delete, dok_clear, lil_clear
    public :: dok_from_dense, lil_from_dense, dok_to_dense, lil_to_dense
    public :: dok_from_coo, lil_from_coo, dok_to_coo, lil_to_coo
    public :: dok_from_csr, lil_from_csr, dok_to_csr, lil_to_csr
    public :: dok_from_csc, lil_from_csc, dok_to_csc, lil_to_csc
    public :: dok_from_lil, lil_from_dok, dok_matvec, lil_matvec
    public :: dok_transpose, lil_transpose, dok_eliminate_zeros, lil_eliminate_zeros

    type :: dok_matrix
        private
        integer, public :: n_rows = 0
        integer, public :: n_cols = 0
        integer :: used = 0
        integer :: tombstones = 0
        integer, allocatable :: row(:)
        integer, allocatable :: col(:)
        integer, allocatable :: state(:)  ! 0 empty, 1 occupied, 2 deleted
        real(dp), allocatable :: data(:)
    end type dok_matrix

    type, private :: lil_row
        integer :: used = 0
        integer, allocatable :: col(:)
        real(dp), allocatable :: data(:)
    end type lil_row

    type :: lil_matrix
        private
        integer, public :: n_rows = 0
        integer, public :: n_cols = 0
        type(lil_row), allocatable :: rows(:)
    end type lil_matrix

contains

    pure logical function dok_ready(a) result(ok)
        type(dok_matrix), intent(in) :: a !! Candidate initialized DOK matrix.
        integer :: cap

        ok = .false.
        if (a%n_rows < 0 .or. a%n_cols < 0) return
        if (a%n_rows == huge(0) .or. a%n_cols == huge(0)) return
        if (.not. allocated(a%state)) return
        if (.not. allocated(a%row)) return
        if (.not. allocated(a%col)) return
        if (.not. allocated(a%data)) return
        cap = size(a%state)
        if (cap < 16) return
        if (size(a%row) /= cap .or. size(a%col) /= cap) return
        if (size(a%data) /= cap) return
        if (a%used < 0 .or. a%tombstones < 0) return
        if (int(a%used,int64)+int(a%tombstones,int64) > int(cap,int64)) return
        ok = .true.
    end function dok_ready

    pure integer function dok_hash(i, j, cap) result(slot)
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        integer, intent(in) :: cap !! Positive hash table capacity.
        integer(int64) :: key

        ! The products fit int64 for the supported default-integer indices.
        key = int(i-1,int64)*1315423911_int64 + int(j-1,int64)*2654435761_int64
        slot = int(modulo(key, int(cap,int64))) + 1
    end function dok_hash

    pure subroutine dok_slot(a, i, j, slot, found)
        type(dok_matrix), intent(in) :: a !! Initialized DOK matrix.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        integer, intent(out) :: slot !! Matching slot or available insertion slot.
        logical, intent(out) :: found !! Whether an existing key matched.
        integer :: first_deleted, k, idx, cap

        cap = size(a%state)
        first_deleted = 0
        found = .false.
        slot = 0
        idx = dok_hash(i,j,cap)
        do k = 1, cap
            select case (a%state(idx))
            case (0)
                slot = idx
                if (first_deleted /= 0) slot = first_deleted
                return
            case (1)
                if (a%row(idx) == i .and. a%col(idx) == j) then
                    slot = idx
                    found = .true.
                    return
                end if
            case (2)
                if (first_deleted == 0) first_deleted = idx
            end select
            idx = idx + 1
            if (idx > cap) idx = 1
        end do
        slot = first_deleted
    end subroutine dok_slot

    pure logical function valid_dok(a) result(ok)
        type(dok_matrix), intent(in) :: a !! Candidate hashed sparse matrix.
        integer :: k, slot, n_used, n_deleted
        logical :: found

        ok = .false.
        if (.not. dok_ready(a)) return
        if (any(a%state < 0) .or. any(a%state > 2)) return
        n_used = count(a%state == 1)
        n_deleted = count(a%state == 2)
        if (n_used /= a%used .or. n_deleted /= a%tombstones) return
        if (n_used == size(a%state)) return
        do k = 1, size(a%state)
            if (a%state(k) /= 1) cycle
            if (a%row(k) < 1 .or. a%row(k) > a%n_rows) return
            if (a%col(k) < 1 .or. a%col(k) > a%n_cols) return
            if (.not. ieee_is_finite(a%data(k))) return
            call dok_slot(a, a%row(k), a%col(k), slot, found)
            if (.not. found .or. slot /= k) return
        end do
        ok = .true.
    end function valid_dok

    function dok_empty(n_rows, n_cols, status) result(a)
        integer, intent(in) :: n_rows !! Nonnegative number of rows.
        integer, intent(in) :: n_cols !! Nonnegative number of columns.
        integer, intent(out), optional :: status !! Zero on success, one otherwise.
        type(dok_matrix) :: a

        if (present(status)) status = sparse_invalid_input
        if (n_rows < 0 .or. n_cols < 0) return
        if (n_rows == huge(0) .or. n_cols == huge(0)) return
        a%n_rows = n_rows
        a%n_cols = n_cols
        allocate(a%row(16), a%col(16), a%state(16), a%data(16))
        a%state = 0
        a%row = 0
        a%col = 0
        a%data = 0.0_dp
        if (present(status)) status = sparse_success
    end function dok_empty

    subroutine dok_rehash(a, new_capacity)
        type(dok_matrix), intent(inout) :: a !! DOK matrix to resize in place.
        integer, intent(in) :: new_capacity !! New capacity exceeding live entries.
        type(dok_matrix) :: b
        integer :: k, slot
        logical :: found

        b%n_rows = a%n_rows
        b%n_cols = a%n_cols
        allocate(b%row(new_capacity), b%col(new_capacity), &
            b%state(new_capacity), b%data(new_capacity))
        b%row = 0
        b%col = 0
        b%state = 0
        b%data = 0.0_dp
        do k = 1, size(a%state)
            if (a%state(k) /= 1) cycle
            call dok_slot(b, a%row(k), a%col(k), slot, found)
            b%state(slot) = 1
            b%row(slot) = a%row(k)
            b%col(slot) = a%col(k)
            b%data(slot) = a%data(k)
            b%used = b%used + 1
        end do
        a = b
    end subroutine dok_rehash

    subroutine dok_put(a, i, j, value, keep_zero, status)
        type(dok_matrix), intent(inout) :: a !! Mutable, initialized DOK matrix.
        integer, intent(in) :: i !! One-based row, within matrix bounds.
        integer, intent(in) :: j !! One-based column, within matrix bounds.
        real(dp), intent(in) :: value !! Finite value to assign.
        logical, intent(in) :: keep_zero !! Retain explicit zero during conversion.
        integer, intent(out) :: status !! Success or invalid input.
        integer :: slot, cap, new_cap
        logical :: found

        status = sparse_invalid_input
        if (.not. dok_ready(a)) return
        if (i < 1 .or. i > a%n_rows .or. j < 1 .or. j > a%n_cols) return
        if (.not. ieee_is_finite(value)) return
        call dok_slot(a, i, j, slot, found)
        if (found) then
            if (value == 0.0_dp .and. .not. keep_zero) then
                a%state(slot) = 2
                a%used = a%used - 1
                a%tombstones = a%tombstones + 1
                if (a%used == 0) then
                    a%state = 0
                    a%tombstones = 0
                end if
            else
                a%data(slot) = value
            end if
            status = sparse_success
            return
        end if
        if (value == 0.0_dp .and. .not. keep_zero) then
            status = sparse_success
            return
        end if
        cap = size(a%state)
        if (10_int64*(int(a%used,int64)+int(a%tombstones,int64)+1_int64) >= &
            7_int64*int(cap,int64)) then
            new_cap = cap
            if (int(a%used,int64)+1_int64 >= int(cap,int64)/2_int64) then
                if (cap > huge(0)/2) return
                new_cap = 2*cap
            end if
            call dok_rehash(a, new_cap)
            call dok_slot(a, i, j, slot, found)
        end if
        if (slot == 0) return
        if (a%state(slot) == 2) a%tombstones = a%tombstones-1
        a%state(slot) = 1
        a%row(slot) = i
        a%col(slot) = j
        a%data(slot) = value
        a%used = a%used+1
        status = sparse_success
    end subroutine dok_put

    subroutine dok_set(a, i, j, value, status)
        type(dok_matrix), intent(inout) :: a !! Matrix modified in place.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        real(dp), intent(in) :: value !! Finite replacement; zero deletes entry.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: ierr

        call dok_put(a,i,j,value,.false.,ierr)
        if (present(status)) status = ierr
    end subroutine dok_set

    function dok_get(a, i, j, status) result(value)
        type(dok_matrix), intent(in) :: a !! Initialized DOK matrix.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp) :: value
        integer :: slot
        logical :: found

        value = 0.0_dp
        if (present(status)) status = sparse_invalid_input
        if (.not. dok_ready(a)) return
        if (i < 1 .or. i > a%n_rows .or. j < 1 .or. j > a%n_cols) return
        call dok_slot(a, i, j, slot, found)
        if (found) value = a%data(slot)
        if (present(status)) status = sparse_success
    end function dok_get

    subroutine dok_add(a, i, j, value, status)
        type(dok_matrix), intent(inout) :: a !! Matrix receiving additive update.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        real(dp), intent(in) :: value !! Finite increment.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: ierr
        real(dp) :: old, updated

        if (present(status)) status = sparse_invalid_input
        if (.not. ieee_is_finite(value)) return
        old = dok_get(a,i,j,ierr)
        if (ierr /= sparse_success) return
        updated = old + value
        if (.not. ieee_is_finite(updated)) return
        call dok_set(a,i,j,updated,ierr)
        if (present(status)) status = ierr
    end subroutine dok_add

    subroutine dok_delete(a, i, j, status)
        type(dok_matrix), intent(inout) :: a !! Matrix from which an entry is removed.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        integer, intent(out), optional :: status !! Zero on success (even if absent).
        integer :: slot
        logical :: found

        if (present(status)) status = sparse_invalid_input
        if (.not. dok_ready(a)) return
        if (i < 1 .or. i > a%n_rows .or. j < 1 .or. j > a%n_cols) return
        call dok_slot(a,i,j,slot,found)
        if (found) then
            a%state(slot) = 2
            a%used = a%used-1
            a%tombstones = a%tombstones+1
            if (a%used == 0) then
                a%state = 0
                a%tombstones = 0
            end if
        end if
        if (present(status)) status = sparse_success
    end subroutine dok_delete

    subroutine dok_clear(a, status)
        type(dok_matrix), intent(inout) :: a !! Matrix emptied without changing shape.
        integer, intent(out), optional :: status !! Zero on success.

        if (present(status)) status = sparse_invalid_input
        if (.not. dok_ready(a)) return
        a%state = 0
        a%used = 0
        a%tombstones = 0
        if (present(status)) status = sparse_success
    end subroutine dok_clear

    subroutine dok_eliminate_zeros(a, status)
        type(dok_matrix), intent(inout) :: a !! Remove stored zeros inherited from conversion.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: k

        if (present(status)) status = sparse_invalid_input
        if (.not. dok_ready(a)) return
        do k = 1, size(a%state)
            if (a%state(k) /= 1) cycle
            if (a%data(k) /= 0.0_dp) cycle
            a%state(k) = 2
            a%used = a%used-1
            a%tombstones = a%tombstones+1
        end do
        if (a%used == 0) then
            a%state = 0
            a%tombstones = 0
        end if
        if (present(status)) status = sparse_success
    end subroutine dok_eliminate_zeros

    integer function dok_nnz(a, status) result(nnz)
        type(dok_matrix), intent(in) :: a !! Candidate DOK matrix.
        integer, intent(out), optional :: status !! Zero on success.

        nnz = -1
        if (present(status)) status = sparse_invalid_input
        if (.not. dok_ready(a)) return
        nnz = a%used
        if (present(status)) status = sparse_success
    end function dok_nnz

    function dok_from_csr(csr, status) result(a)
        type(csr_matrix), intent(in) :: csr !! Canonical CSR input with stored zeros.
        integer, intent(out), optional :: status !! Zero on success.
        type(dok_matrix) :: a
        integer :: i, k, ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(csr)) return
        a = dok_empty(csr%n_rows,csr%n_cols,ierr)
        if (ierr /= sparse_success) return
        do i = 1, csr%n_rows
            do k = csr%indptr(i), csr%indptr(i+1)-1
                call dok_put(a,i,csr%indices(k),csr%data(k),.true.,ierr)
                if (ierr /= sparse_success) then
                    a = dok_matrix()
                    return
                end if
            end do
        end do
        if (present(status)) status = sparse_success
    end function dok_from_csr

    function dok_from_coo(coo, status) result(a)
        type(coo_matrix), intent(in) :: coo !! COO input; duplicates are summed.
        integer, intent(out), optional :: status !! Zero on success.
        type(dok_matrix) :: a
        type(csr_matrix) :: csr
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo(coo)) return
        csr = coo_to_csr(coo,ierr)
        if (ierr /= sparse_success) return
        a = dok_from_csr(csr,ierr)
        if (present(status)) status = ierr
    end function dok_from_coo

    function dok_from_csc(csc, status) result(a)
        type(csc_matrix), intent(in) :: csc !! Canonical CSC input.
        integer, intent(out), optional :: status !! Zero on success.
        type(dok_matrix) :: a
        type(coo_matrix) :: coo
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        coo = coo_from_csc(csc,ierr)
        if (ierr /= sparse_success) return
        a = dok_from_coo(coo,ierr)
        if (present(status)) status = ierr
    end function dok_from_csc

    function dok_from_dense(dense, status) result(a)
        real(dp), intent(in) :: dense(:,:) !! Finite dense matrix; implicit zeros skipped.
        integer, intent(out), optional :: status !! Zero on success.
        type(dok_matrix) :: a
        integer :: i, j, ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. all(ieee_is_finite(dense))) return
        a = dok_empty(size(dense,1),size(dense,2),ierr)
        if (ierr /= sparse_success) return
        do i = 1, a%n_rows
            do j = 1, a%n_cols
                if (dense(i,j) == 0.0_dp) cycle
                call dok_set(a,i,j,dense(i,j),ierr)
                if (ierr /= sparse_success) then
                    a = dok_matrix()
                    return
                end if
            end do
        end do
        if (present(status)) status = sparse_success
    end function dok_from_dense

    function dok_to_coo(a, status) result(coo)
        type(dok_matrix), intent(in) :: a !! Valid DOK input.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_matrix) :: coo
        integer, allocatable :: rr(:), cc(:)
        real(dp), allocatable :: vv(:)
        integer :: k, n, ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_dok(a)) return
        allocate(rr(a%used),cc(a%used),vv(a%used))
        n = 0
        do k = 1, size(a%state)
            if (a%state(k) /= 1) cycle
            n = n+1
            rr(n) = a%row(k)
            cc(n) = a%col(k)
            vv(n) = a%data(k)
        end do
        coo = coo_from_triplets(a%n_rows,a%n_cols,rr,cc,vv,ierr)
        if (present(status)) status = ierr
    end function dok_to_coo

    function dok_to_csr(a, status) result(csr)
        type(dok_matrix), intent(in) :: a !! Valid DOK input.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: csr
        type(coo_matrix) :: coo
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        coo = dok_to_coo(a,ierr)
        if (ierr /= sparse_success) return
        csr = coo_to_csr(coo,ierr)
        if (present(status)) status = ierr
    end function dok_to_csr

    function dok_to_csc(a, status) result(csc)
        type(dok_matrix), intent(in) :: a !! Valid DOK input.
        integer, intent(out), optional :: status !! Zero on success.
        type(csc_matrix) :: csc
        type(coo_matrix) :: coo
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        coo = dok_to_coo(a,ierr)
        if (ierr /= sparse_success) return
        csc = coo_to_csc(coo,ierr)
        if (present(status)) status = ierr
    end function dok_to_csc

    function dok_to_dense(a, status) result(dense)
        type(dok_matrix), intent(in) :: a !! Valid DOK input.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: dense(:,:)
        integer :: k

        if (present(status)) status = sparse_invalid_input
        allocate(dense(0,0))
        if (.not. valid_dok(a)) return
        deallocate(dense)
        allocate(dense(a%n_rows,a%n_cols))
        dense = 0.0_dp
        do k = 1, size(a%state)
            if (a%state(k) /= 1) cycle
            dense(a%row(k),a%col(k)) = a%data(k)
        end do
        if (present(status)) status = sparse_success
    end function dok_to_dense

    function dok_matvec(a, x, status) result(y)
        type(dok_matrix), intent(in) :: a !! Valid DOK input.
        real(dp), intent(in) :: x(:) !! Finite vector of length n_cols.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: y(:)
        integer :: k, i

        if (present(status)) status = sparse_invalid_input
        allocate(y(0))
        if (.not. valid_dok(a)) return
        if (size(x) /= a%n_cols) return
        if (.not. all(ieee_is_finite(x))) return
        deallocate(y)
        allocate(y(a%n_rows))
        y = 0.0_dp
        do k = 1, size(a%state)
            if (a%state(k) /= 1) cycle
            i = a%row(k)
            y(i) = y(i) + a%data(k)*x(a%col(k))
            if (.not. ieee_is_finite(y(i))) then
                deallocate(y)
                allocate(y(0))
                return
            end if
        end do
        if (present(status)) status = sparse_success
    end function dok_matvec

    function dok_transpose(a, status) result(b)
        type(dok_matrix), intent(in) :: a !! Valid DOK input.
        integer, intent(out), optional :: status !! Zero on success.
        type(dok_matrix) :: b
        integer :: k, ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_dok(a)) return
        b = dok_empty(a%n_cols,a%n_rows,ierr)
        if (ierr /= sparse_success) return
        do k = 1, size(a%state)
            if (a%state(k) /= 1) cycle
            call dok_put(b,a%col(k),a%row(k),a%data(k),.true.,ierr)
            if (ierr /= sparse_success) then
                b = dok_matrix()
                return
            end if
        end do
        if (present(status)) status = sparse_success
    end function dok_transpose

    pure logical function lil_ready(a) result(ok)
        type(lil_matrix), intent(in) :: a !! Candidate initialized LIL matrix.

        ok = .false.
        if (a%n_rows < 0 .or. a%n_cols < 0) return
        if (a%n_rows == huge(0) .or. a%n_cols == huge(0)) return
        if (.not. allocated(a%rows)) return
        if (size(a%rows) /= a%n_rows) return
        ok = .true.
    end function lil_ready

    pure logical function valid_lil(a) result(ok)
        type(lil_matrix), intent(in) :: a !! Candidate sorted-row LIL matrix.
        integer :: i, nnz

        ok = .false.
        if (.not. lil_ready(a)) return
        do i = 1, a%n_rows
            nnz = a%rows(i)%used
            if (nnz < 0) return
            if (allocated(a%rows(i)%col) .neqv. allocated(a%rows(i)%data)) return
            if (.not. allocated(a%rows(i)%col)) then
                if (nnz /= 0) return
                cycle
            end if
            if (size(a%rows(i)%col) /= size(a%rows(i)%data)) return
            if (nnz > size(a%rows(i)%col)) return
            if (any(a%rows(i)%col(:nnz) < 1)) return
            if (any(a%rows(i)%col(:nnz) > a%n_cols)) return
            if (nnz > 1) then
                if (any(a%rows(i)%col(2:nnz) <= a%rows(i)%col(:nnz-1))) return
            end if
            if (.not. all(ieee_is_finite(a%rows(i)%data(:nnz)))) return
        end do
        ok = .true.
    end function valid_lil

    function lil_empty(n_rows, n_cols, status) result(a)
        integer, intent(in) :: n_rows !! Nonnegative number of rows.
        integer, intent(in) :: n_cols !! Nonnegative number of columns.
        integer, intent(out), optional :: status !! Zero on success.
        type(lil_matrix) :: a

        if (present(status)) status = sparse_invalid_input
        if (n_rows < 0 .or. n_cols < 0) return
        if (n_rows == huge(0) .or. n_cols == huge(0)) return
        a%n_rows = n_rows
        a%n_cols = n_cols
        allocate(a%rows(n_rows))
        if (present(status)) status = sparse_success
    end function lil_empty

    pure subroutine lil_slot(row, j, pos, found)
        type(lil_row), intent(in) :: row !! Sorted row to search.
        integer, intent(in) :: j !! One-based target column.
        integer, intent(out) :: pos !! Found position or insertion position.
        logical, intent(out) :: found !! Whether target column exists.
        integer :: lo, hi, mid

        lo = 1
        hi = row%used + 1
        do while (lo < hi)
            mid = lo + (hi-lo)/2
            if (row%col(mid) < j) then
                lo = mid+1
            else
                hi = mid
            end if
        end do
        pos = lo
        found = .false.
        if (pos <= row%used) found = row%col(pos) == j
    end subroutine lil_slot

    subroutine lil_grow(row)
        type(lil_row), intent(inout) :: row !! Mutable LIL row to enlarge.
        integer, allocatable :: new_col(:)
        real(dp), allocatable :: new_data(:)
        integer :: cap, next_cap

        cap = 0
        if (allocated(row%col)) cap = size(row%col)
        if (cap == 0) then
            next_cap = 4
        else
            next_cap = min(huge(0)-1,2*cap)
        end if
        allocate(new_col(next_cap),new_data(next_cap))
        if (row%used > 0) then
            new_col(:row%used) = row%col(:row%used)
            new_data(:row%used) = row%data(:row%used)
        end if
        call move_alloc(new_col,row%col)
        call move_alloc(new_data,row%data)
    end subroutine lil_grow

    subroutine lil_put(a, i, j, value, keep_zero, status)
        type(lil_matrix), intent(inout) :: a !! LIL matrix updated in place.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        real(dp), intent(in) :: value !! Finite replacement value.
        logical, intent(in) :: keep_zero !! Keep explicitly stored zero from conversion.
        integer, intent(out) :: status !! Zero on success.
        integer :: pos, k
        logical :: found

        status = sparse_invalid_input
        if (.not. lil_ready(a)) return
        if (i < 1 .or. i > a%n_rows .or. j < 1 .or. j > a%n_cols) return
        if (.not. ieee_is_finite(value)) return
        call lil_slot(a%rows(i),j,pos,found)
        if (found) then
            if (value == 0.0_dp .and. .not. keep_zero) then
                do k = pos, a%rows(i)%used-1
                    a%rows(i)%col(k) = a%rows(i)%col(k+1)
                    a%rows(i)%data(k) = a%rows(i)%data(k+1)
                end do
                a%rows(i)%used = a%rows(i)%used-1
            else
                a%rows(i)%data(pos) = value
            end if
            status = sparse_success
            return
        end if
        if (value == 0.0_dp .and. .not. keep_zero) then
            status = sparse_success
            return
        end if
        if (a%rows(i)%used == huge(0)-1) return
        if (.not. allocated(a%rows(i)%col)) then
            call lil_grow(a%rows(i))
        else if (a%rows(i)%used == size(a%rows(i)%col)) then
            if (size(a%rows(i)%col) > (huge(0)-1)/2) return
            call lil_grow(a%rows(i))
        end if
        do k = a%rows(i)%used, pos, -1
            a%rows(i)%col(k+1) = a%rows(i)%col(k)
            a%rows(i)%data(k+1) = a%rows(i)%data(k)
        end do
        a%rows(i)%col(pos) = j
        a%rows(i)%data(pos) = value
        a%rows(i)%used = a%rows(i)%used + 1
        status = sparse_success
    end subroutine lil_put

    subroutine lil_set(a, i, j, value, status)
        type(lil_matrix), intent(inout) :: a !! Matrix modified in place.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        real(dp), intent(in) :: value !! Finite replacement; zero deletes the key.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: ierr

        call lil_put(a,i,j,value,.false.,ierr)
        if (present(status)) status = ierr
    end subroutine lil_set

    function lil_get(a, i, j, status) result(value)
        type(lil_matrix), intent(in) :: a !! Initialized LIL matrix.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp) :: value
        integer :: pos
        logical :: found

        value = 0.0_dp
        if (present(status)) status = sparse_invalid_input
        if (.not. lil_ready(a)) return
        if (i < 1 .or. i > a%n_rows .or. j < 1 .or. j > a%n_cols) return
        call lil_slot(a%rows(i),j,pos,found)
        if (found) value = a%rows(i)%data(pos)
        if (present(status)) status = sparse_success
    end function lil_get

    subroutine lil_add(a, i, j, value, status)
        type(lil_matrix), intent(inout) :: a !! Matrix receiving additive update.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        real(dp), intent(in) :: value !! Finite increment.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: ierr
        real(dp) :: old, updated

        if (present(status)) status = sparse_invalid_input
        if (.not. ieee_is_finite(value)) return
        old = lil_get(a,i,j,ierr)
        if (ierr /= sparse_success) return
        updated = old + value
        if (.not. ieee_is_finite(updated)) return
        call lil_set(a,i,j,updated,ierr)
        if (present(status)) status = ierr
    end subroutine lil_add

    subroutine lil_delete(a, i, j, status)
        type(lil_matrix), intent(inout) :: a !! Matrix from which an entry is deleted.
        integer, intent(in) :: i !! One-based row coordinate.
        integer, intent(in) :: j !! One-based column coordinate.
        integer, intent(out), optional :: status !! Zero on success (also if absent).
        integer :: pos, k
        logical :: found

        if (present(status)) status = sparse_invalid_input
        if (.not. lil_ready(a)) return
        if (i < 1 .or. i > a%n_rows .or. j < 1 .or. j > a%n_cols) return
        call lil_slot(a%rows(i),j,pos,found)
        if (found) then
            do k = pos, a%rows(i)%used-1
                a%rows(i)%col(k) = a%rows(i)%col(k+1)
                a%rows(i)%data(k) = a%rows(i)%data(k+1)
            end do
            a%rows(i)%used = a%rows(i)%used-1
        end if
        if (present(status)) status = sparse_success
    end subroutine lil_delete

    subroutine lil_clear(a, status)
        type(lil_matrix), intent(inout) :: a !! Empty all rows, preserving shape.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: i

        if (present(status)) status = sparse_invalid_input
        if (.not. lil_ready(a)) return
        do i = 1, a%n_rows
            a%rows(i)%used = 0
        end do
        if (present(status)) status = sparse_success
    end subroutine lil_clear

    subroutine lil_eliminate_zeros(a, status)
        type(lil_matrix), intent(inout) :: a !! Remove explicit zeros after conversion.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: i, k, n

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_lil(a)) return
        do i = 1, a%n_rows
            n = 0
            do k = 1, a%rows(i)%used
                if (a%rows(i)%data(k) == 0.0_dp) cycle
                n = n+1
                a%rows(i)%col(n) = a%rows(i)%col(k)
                a%rows(i)%data(n) = a%rows(i)%data(k)
            end do
            a%rows(i)%used = n
        end do
        if (present(status)) status = sparse_success
    end subroutine lil_eliminate_zeros

    integer function lil_nnz(a, status) result(nnz)
        type(lil_matrix), intent(in) :: a !! Candidate LIL matrix.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: i
        integer(int64) :: total

        nnz = -1
        if (present(status)) status = sparse_invalid_input
        if (.not. lil_ready(a)) return
        total = 0_int64
        do i = 1, a%n_rows
            total = total + int(a%rows(i)%used,int64)
        end do
        if (total >= int(huge(0),int64)) return
        nnz = int(total)
        if (present(status)) status = sparse_success
    end function lil_nnz

    function lil_from_csr(csr, status) result(a)
        type(csr_matrix), intent(in) :: csr !! Canonical CSR input, zeros retained.
        integer, intent(out), optional :: status !! Zero on success.
        type(lil_matrix) :: a
        integer :: i, k, first, last, n, ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_csr(csr)) return
        a = lil_empty(csr%n_rows,csr%n_cols,ierr)
        if (ierr /= sparse_success) return
        do i = 1, a%n_rows
            first = csr%indptr(i)
            last = csr%indptr(i+1)-1
            n = last-first+1
            if (n == 0) cycle
            allocate(a%rows(i)%col(n),a%rows(i)%data(n))
            a%rows(i)%used = n
            do k = 1, n
                a%rows(i)%col(k) = csr%indices(first+k-1)
                a%rows(i)%data(k) = csr%data(first+k-1)
            end do
        end do
        if (present(status)) status = sparse_success
    end function lil_from_csr

    function lil_from_coo(coo, status) result(a)
        type(coo_matrix), intent(in) :: coo !! COO with duplicates and explicit zeros.
        integer, intent(out), optional :: status !! Zero on success.
        type(lil_matrix) :: a
        type(csr_matrix) :: csr
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_coo(coo)) return
        csr = coo_to_csr(coo,ierr)
        if (ierr /= sparse_success) return
        a = lil_from_csr(csr,ierr)
        if (present(status)) status = ierr
    end function lil_from_coo

    function lil_from_csc(csc, status) result(a)
        type(csc_matrix), intent(in) :: csc !! Canonical CSC input.
        integer, intent(out), optional :: status !! Zero on success.
        type(lil_matrix) :: a
        type(coo_matrix) :: coo
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        coo = coo_from_csc(csc,ierr)
        if (ierr /= sparse_success) return
        a = lil_from_coo(coo,ierr)
        if (present(status)) status = ierr
    end function lil_from_csc

    function lil_from_dense(dense, status) result(a)
        real(dp), intent(in) :: dense(:,:) !! Finite dense matrix; zeros implicit.
        integer, intent(out), optional :: status !! Zero on success.
        type(lil_matrix) :: a
        type(csr_matrix) :: csr
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        csr = csr_from_dense(dense,ierr)
        if (ierr /= sparse_success) return
        a = lil_from_csr(csr,ierr)
        if (present(status)) status = ierr
    end function lil_from_dense

    function lil_to_csr(a, status) result(csr)
        type(lil_matrix), intent(in) :: a !! Valid sorted-row LIL input.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: csr
        integer :: nnz, ierr, i, first, last, k

        if (present(status)) status = sparse_invalid_input
        if (.not. valid_lil(a)) return
        nnz = lil_nnz(a,ierr)
        if (ierr /= sparse_success) return
        allocate(csr%indptr(a%n_rows+1),csr%indices(nnz),csr%data(nnz))
        csr%n_rows = a%n_rows
        csr%n_cols = a%n_cols
        k = 1
        do i = 1, a%n_rows
            csr%indptr(i) = k
            first = k
            last = k + a%rows(i)%used - 1
            if (last >= first) then
                csr%indices(first:last) = a%rows(i)%col(:a%rows(i)%used)
                csr%data(first:last) = a%rows(i)%data(:a%rows(i)%used)
            end if
            k = last+1
        end do
        csr%indptr(a%n_rows+1) = k
        if (present(status)) status = sparse_success
    end function lil_to_csr

    function lil_to_coo(a, status) result(coo)
        type(lil_matrix), intent(in) :: a !! Valid sorted-row LIL input.
        integer, intent(out), optional :: status !! Zero on success.
        type(coo_matrix) :: coo
        type(csr_matrix) :: csr
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        csr = lil_to_csr(a,ierr)
        if (ierr /= sparse_success) return
        coo = coo_from_csr(csr,ierr)
        if (present(status)) status = ierr
    end function lil_to_coo

    function lil_to_csc(a, status) result(csc)
        type(lil_matrix), intent(in) :: a !! Valid LIL input.
        integer, intent(out), optional :: status !! Zero on success.
        type(csc_matrix) :: csc
        type(coo_matrix) :: coo
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        coo = lil_to_coo(a,ierr)
        if (ierr /= sparse_success) return
        csc = coo_to_csc(coo,ierr)
        if (present(status)) status = ierr
    end function lil_to_csc

    function lil_to_dense(a, status) result(dense)
        type(lil_matrix), intent(in) :: a !! Valid LIL matrix.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: dense(:,:)
        integer :: i, k

        if (present(status)) status = sparse_invalid_input
        allocate(dense(0,0))
        if (.not. valid_lil(a)) return
        deallocate(dense)
        allocate(dense(a%n_rows,a%n_cols))
        dense = 0.0_dp
        do i = 1, a%n_rows
            do k = 1, a%rows(i)%used
                dense(i,a%rows(i)%col(k)) = a%rows(i)%data(k)
            end do
        end do
        if (present(status)) status = sparse_success
    end function lil_to_dense

    function lil_matvec(a, x, status) result(y)
        type(lil_matrix), intent(in) :: a !! Valid LIL matrix.
        real(dp), intent(in) :: x(:) !! Finite vector of length n_cols.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: y(:)
        integer :: i, k

        if (present(status)) status = sparse_invalid_input
        allocate(y(0))
        if (.not. valid_lil(a)) return
        if (size(x) /= a%n_cols .or. .not. all(ieee_is_finite(x))) return
        deallocate(y)
        allocate(y(a%n_rows))
        y = 0.0_dp
        do i = 1, a%n_rows
            do k = 1, a%rows(i)%used
                y(i) = y(i) + a%rows(i)%data(k)*x(a%rows(i)%col(k))
                if (.not. ieee_is_finite(y(i))) then
                    deallocate(y)
                    allocate(y(0))
                    return
                end if
            end do
        end do
        if (present(status)) status = sparse_success
    end function lil_matvec

    function lil_transpose(a, status) result(b)
        type(lil_matrix), intent(in) :: a !! Valid LIL matrix.
        integer, intent(out), optional :: status !! Zero on success.
        type(lil_matrix) :: b
        type(coo_matrix) :: coo, trans
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        coo = lil_to_coo(a,ierr)
        if (ierr /= sparse_success) return
        trans = coo_transpose(coo,ierr)
        if (ierr /= sparse_success) return
        b = lil_from_coo(trans,ierr)
        if (present(status)) status = ierr
    end function lil_transpose

    function dok_from_lil(a, status) result(b)
        type(lil_matrix), intent(in) :: a !! LIL matrix to convert to DOK.
        integer, intent(out), optional :: status !! Zero on success.
        type(dok_matrix) :: b
        type(csr_matrix) :: csr
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        csr = lil_to_csr(a,ierr)
        if (ierr /= sparse_success) return
        b = dok_from_csr(csr,ierr)
        if (present(status)) status = ierr
    end function dok_from_lil

    function lil_from_dok(a, status) result(b)
        type(dok_matrix), intent(in) :: a !! DOK matrix to convert to LIL.
        integer, intent(out), optional :: status !! Zero on success.
        type(lil_matrix) :: b
        type(csr_matrix) :: csr
        integer :: ierr

        if (present(status)) status = sparse_invalid_input
        csr = dok_to_csr(a,ierr)
        if (ierr /= sparse_success) return
        b = lil_from_csr(csr,ierr)
        if (present(status)) status = ierr
    end function lil_from_dok

end module scifort_sparse_dok_lil
