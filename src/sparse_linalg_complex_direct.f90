! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Complex pivoted sparse LU factorization and direct solves for CSR and CSC.
! One-based canonical CSR factors store unit-lower L below the diagonal and
! U on/above it.  Partial row pivoting is applied without column reordering.
! The method is an unblocked sparse Gaussian elimination, not a SuperLU clone.
module scifort_sparse_linalg_complex_direct
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, csc_complex_to_csr
    implicit none
    private

    integer, parameter, public :: complex_direct_success = 0
    integer, parameter, public :: complex_direct_invalid_input = -1
    integer, parameter, public :: complex_direct_breakdown = -2
    integer, parameter, public :: complex_direct_singular = -4

    type, public :: splu_complex_factor
        type(csr_complex_matrix) :: lu !! Combined, pivoted sparse L and U factors.
        integer, allocatable :: row_perm(:) !! P(i,row_perm(i)) = 1; P*A=L*U.
    end type splu_complex_factor

    public :: splu_complex, valid_splu_complex, splu_complex_solve, spsolve_complex

    interface splu_complex
        module procedure splu_complex_csr
        module procedure splu_complex_csc
    end interface splu_complex
    interface splu_complex_solve
        module procedure splu_complex_solve_vector
        module procedure splu_complex_solve_matrix
    end interface splu_complex_solve
    interface spsolve_complex
        module procedure spsolve_complex_csr_vector
        module procedure spsolve_complex_csr_matrix
        module procedure spsolve_complex_csc_vector
        module procedure spsolve_complex_csc_matrix
    end interface spsolve_complex

    type :: mutable_complex_row
        integer, allocatable :: columns(:)
        complex(dp), allocatable :: values(:)
        integer :: used = 0
    end type mutable_complex_row

contains

    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Candidate complex value.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    pure integer function find_column(row, column) result(pos)
        type(mutable_complex_row), intent(in) :: row !! Sorted row with used entries.
        integer, intent(in) :: column !! One-based column being sought.
        integer :: left, right, middle

        left = 1
        right = row%used
        do while (left <= right)
            middle = left + (right-left)/2
            if (row%columns(middle) == column) then
                pos = middle
                return
            else if (row%columns(middle) < column) then
                left = middle+1
            else
                right = middle-1
            end if
        end do
        pos = -left
    end function find_column

    pure complex(dp) function row_value(row, column) result(value)
        type(mutable_complex_row), intent(in) :: row !! Sorted sparse factor row.
        integer, intent(in) :: column !! One-based column.
        integer :: pos

        value = cmplx(0.0_dp,0.0_dp,dp)
        pos = find_column(row,column)
        if (pos > 0) value = row%values(pos)
    end function row_value

    subroutine set_row_value(row, column, value, info)
        type(mutable_complex_row), intent(inout) :: row !! Row modified in place.
        integer, intent(in) :: column !! One-based column to replace or insert.
        complex(dp), intent(in) :: value !! Replacement, including explicit zero.
        integer, intent(out) :: info !! Zero, or negative for overflow/allocation failure.
        integer :: p, new_capacity, current_capacity, stat
        integer, allocatable :: new_columns(:)
        complex(dp), allocatable :: new_values(:)

        info = complex_direct_success
        p = find_column(row,column)
        if (p > 0) then
            row%values(p) = value
            return
        end if
        if (row%used == huge(0)) then
            info = complex_direct_breakdown
            return
        end if
        current_capacity = size(row%columns)
        if (row%used == current_capacity) then
            if (current_capacity > (huge(0)-1)/2) then
                info = complex_direct_breakdown
                return
            end if
            new_capacity = max(4,2*current_capacity+1)
            allocate(new_columns(new_capacity),new_values(new_capacity),stat=stat)
            if (stat /= 0) then
                info = complex_direct_breakdown
                return
            end if
            new_columns(:row%used) = row%columns(:row%used)
            new_values(:row%used) = row%values(:row%used)
            call move_alloc(new_columns,row%columns)
            call move_alloc(new_values,row%values)
        end if
        p = -p
        if (p <= row%used) then
            row%columns(p+1:row%used+1) = row%columns(p:row%used)
            row%values(p+1:row%used+1) = row%values(p:row%used)
        end if
        row%columns(p) = column
        row%values(p) = value
        row%used = row%used+1
    end subroutine set_row_value

    subroutine splu_complex_csc(a, factor, info, pivot_tol)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC matrix.
        type(splu_complex_factor), intent(out) :: factor !! Pivoted LU on success.
        integer, intent(out) :: info !! 0 success; -1 input; -2 breakdown; -4 singular.
        real(dp), intent(in), optional :: pivot_tol !! Absolute minimum pivot magnitude, default 0.
        type(csr_complex_matrix) :: csr
        integer :: ierr

        info = complex_direct_invalid_input
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,ierr)
        if (ierr /= 0) return
        call splu_complex_csr(csr,factor,info,pivot_tol)
    end subroutine splu_complex_csc

    subroutine splu_complex_csr(a, factor, info, pivot_tol)
        type(csr_complex_matrix), intent(in) :: a !! Square canonical CSR matrix.
        type(splu_complex_factor), intent(out) :: factor !! Pivoted LU on success.
        integer, intent(out) :: info !! 0 success; -1 input; -2 breakdown; -4 singular.
        real(dp), intent(in), optional :: pivot_tol !! Absolute minimum pivot magnitude, default 0.
        type(mutable_complex_row), allocatable :: rows(:)
        type(mutable_complex_row) :: temporary
        integer :: n, i, j, k, p, pos, ierr, nnz, offset, count, permutation
        integer(int64) :: nstored
        real(dp) :: tol, best, magnitude
        complex(dp) :: pivot, multiplier, new_value, value

        info = complex_direct_invalid_input
        if (.not. valid_csr_complex(a)) return
        n = a%n_rows
        if (a%n_cols /= n) return
        tol = 0.0_dp
        if (present(pivot_tol)) tol = pivot_tol
        if (.not. ieee_is_finite(tol) .or. tol < 0.0_dp) return
        allocate(rows(n),factor%row_perm(n))
        factor%row_perm = [(i,i=1,n)]
        do i = 1, n
            count = a%indptr(i+1)-a%indptr(i)
            allocate(rows(i)%columns(count),rows(i)%values(count))
            rows(i)%used = count
            rows(i)%columns = a%indices(a%indptr(i):a%indptr(i+1)-1)
            rows(i)%values = a%data(a%indptr(i):a%indptr(i+1)-1)
        end do
        do k = 1,n
            p = k
            best = abs(row_value(rows(k),k))
            do i = k+1,n
                magnitude = abs(row_value(rows(i),k))
                if (magnitude > best) then
                    p = i
                    best = magnitude
                end if
            end do
            if (best <= tol) then
                factor = splu_complex_factor()
                info = complex_direct_singular
                return
            end if
            if (.not. ieee_is_finite(best)) then
                factor = splu_complex_factor()
                info = complex_direct_breakdown
                return
            end if
            if (p /= k) then
                temporary = rows(k)
                rows(k) = rows(p)
                rows(p) = temporary
                permutation = factor%row_perm(k)
                factor%row_perm(k) = factor%row_perm(p)
                factor%row_perm(p) = permutation
            end if
            pivot = row_value(rows(k),k)
            do i = k+1,n
                value = row_value(rows(i),k)
                if (value == cmplx(0.0_dp,0.0_dp,dp)) cycle
                multiplier = value/pivot
                if (.not. finite_complex(multiplier)) then
                    factor = splu_complex_factor()
                    info = complex_direct_breakdown
                    return
                end if
                call set_row_value(rows(i),k,multiplier,ierr)
                if (ierr /= 0) then
                    factor = splu_complex_factor()
                    info = ierr
                    return
                end if
                do pos = 1,rows(k)%used
                    j = rows(k)%columns(pos)
                    if (j <= k) cycle
                    new_value = row_value(rows(i),j) - multiplier*rows(k)%values(pos)
                    if (.not. finite_complex(new_value)) then
                        factor = splu_complex_factor()
                        info = complex_direct_breakdown
                        return
                    end if
                    call set_row_value(rows(i),j,new_value,ierr)
                    if (ierr /= 0) then
                        factor = splu_complex_factor()
                        info = ierr
                        return
                    end if
                end do
            end do
        end do
        nstored = 0_int64
        do i = 1,n
            nstored = nstored + int(rows(i)%used,int64)
            if (nstored >= int(huge(0),int64)) then
                factor = splu_complex_factor()
                info = complex_direct_breakdown
                return
            end if
        end do
        nnz = int(nstored)
        factor%lu%n_rows = n
        factor%lu%n_cols = n
        allocate(factor%lu%indptr(n+1),factor%lu%indices(nnz),factor%lu%data(nnz))
        offset = 1
        do i = 1,n
            factor%lu%indptr(i) = offset
            count = rows(i)%used
            factor%lu%indices(offset:offset+count-1) = rows(i)%columns(:count)
            factor%lu%data(offset:offset+count-1) = rows(i)%values(:count)
            offset = offset+count
        end do
        factor%lu%indptr(n+1) = offset
        info = complex_direct_success
    end subroutine splu_complex_csr

    pure logical function valid_splu_complex(factor) result(ok)
        type(splu_complex_factor), intent(in) :: factor !! Candidate LU factors.
        logical, allocatable :: seen(:)
        logical :: has_pivot
        integer :: n, i, k

        ok = .false.
        if (.not. valid_csr_complex(factor%lu)) return
        n = factor%lu%n_rows
        if (factor%lu%n_cols /= n) return
        if (.not. allocated(factor%row_perm)) return
        if (size(factor%row_perm) /= n) return
        if (any(factor%row_perm < 1) .or. any(factor%row_perm > n)) return
        allocate(seen(n))
        seen = .false.
        do i = 1,n
            if (seen(factor%row_perm(i))) return
            seen(factor%row_perm(i)) = .true.
            has_pivot = .false.
            do k = factor%lu%indptr(i),factor%lu%indptr(i+1)-1
                if (factor%lu%indices(k) /= i) cycle
                if (factor%lu%data(k) /= cmplx(0.0_dp,0.0_dp,dp)) has_pivot = .true.
            end do
            if (.not. has_pivot) return
        end do
        ok = .true.
    end function valid_splu_complex

    subroutine splu_complex_solve_vector(factor, b, x, info, trans)
        type(splu_complex_factor), intent(in) :: factor !! Valid pivoted LU factors.
        complex(dp), intent(in) :: b(:) !! Right-hand-side vector of size n.
        complex(dp), allocatable, intent(out) :: x(:) !! Solution, or zero-size on failure.
        integer, intent(out) :: info !! 0 success; -1 input; -2 breakdown.
        character(len=*), intent(in), optional :: trans !! 'N', 'T', or 'H'; default 'N'.
        complex(dp), allocatable :: block(:,:), rhs(:,:)

        allocate(rhs(size(b),1))
        rhs(:,1) = b
        call splu_complex_solve_matrix(factor,rhs,block,info,trans)
        if (info == complex_direct_success) then
            x = block(:,1)
        else
            allocate(x(0))
        end if
    end subroutine splu_complex_solve_vector

    subroutine splu_complex_solve_matrix(factor, b, x, info, trans)
        type(splu_complex_factor), intent(in) :: factor !! Valid pivoted LU factors.
        complex(dp), intent(in) :: b(:,:) !! Multiple RHS, shape (n,nrhs).
        complex(dp), allocatable, intent(out) :: x(:,:) !! Solution, or (0,0) on failure.
        integer, intent(out) :: info !! 0 success; -1 input; -2 breakdown.
        character(len=*), intent(in), optional :: trans !! 'N', 'T', or 'H'; default 'N'.
        complex(dp), allocatable :: work(:,:)
        complex(dp) :: pivot, coefficient
        character(len=1) :: operation
        integer :: n, i, j, k
        logical :: adjoint

        info = complex_direct_invalid_input
        allocate(x(0,0))
        if (.not. valid_splu_complex(factor)) return
        n = factor%lu%n_rows
        if (size(b,1) /= n) return
        if (.not. all(finite_complex(b))) return
        operation = 'N'
        if (present(trans)) then
            if (len_trim(trans) /= 1) return
            operation = trans(1:1)
        end if
        if (operation /= 'N' .and. operation /= 'T' .and. operation /= 'H') return
        adjoint = operation == 'H'
        deallocate(x)
        allocate(x(n,size(b,2)))
        if (operation == 'N') then
            do i = 1,n
                x(i,:) = b(factor%row_perm(i),:)
            end do
            ! L*y=P*b; stored multipliers are strictly below the diagonal.
            do i = 1,n
                do k = factor%lu%indptr(i),factor%lu%indptr(i+1)-1
                    j = factor%lu%indices(k)
                    if (j >= i) exit
                    x(i,:) = x(i,:) - factor%lu%data(k)*x(j,:)
                end do
            end do
            ! U*x=y.
            do i = n,1,-1
                pivot = cmplx(0.0_dp,0.0_dp,dp)
                do k = factor%lu%indptr(i),factor%lu%indptr(i+1)-1
                    j = factor%lu%indices(k)
                    if (j == i) pivot = factor%lu%data(k)
                    if (j > i) x(i,:) = x(i,:) - factor%lu%data(k)*x(j,:)
                end do
                x(i,:) = x(i,:)/pivot
            end do
        else
            ! A^T=U^T*L^T*P and A^H=U^H*L^H*P.
            ! Forward solve U^(T/H)*z=b via row-wise scatter updates.
            x = b
            do i = 1,n
                pivot = cmplx(0.0_dp,0.0_dp,dp)
                do k = factor%lu%indptr(i),factor%lu%indptr(i+1)-1
                    if (factor%lu%indices(k) == i) then
                        pivot = factor%lu%data(k)
                        exit
                    end if
                end do
                if (adjoint) pivot = conjg(pivot)
                x(i,:) = x(i,:)/pivot
                do k = factor%lu%indptr(i),factor%lu%indptr(i+1)-1
                    j = factor%lu%indices(k)
                    if (j <= i) cycle
                    coefficient = factor%lu%data(k)
                    if (adjoint) coefficient = conjg(coefficient)
                    x(j,:) = x(j,:) - coefficient*x(i,:)
                end do
            end do
            ! Back solve L^(T/H)*y=z, then undo P's row permutation.
            do i = n,1,-1
                do k = factor%lu%indptr(i),factor%lu%indptr(i+1)-1
                    j = factor%lu%indices(k)
                    if (j >= i) exit
                    coefficient = factor%lu%data(k)
                    if (adjoint) coefficient = conjg(coefficient)
                    x(j,:) = x(j,:) - coefficient*x(i,:)
                end do
            end do
            work = x
            do i = 1,n
                x(factor%row_perm(i),:) = work(i,:)
            end do
        end if
        if (.not. all(finite_complex(x))) then
            deallocate(x)
            allocate(x(0,0))
            info = complex_direct_breakdown
            return
        end if
        info = complex_direct_success
    end subroutine splu_complex_solve_matrix

    subroutine spsolve_complex_csr_vector(a, b, x, info, pivot_tol)
        type(csr_complex_matrix), intent(in) :: a !! Square complex CSR matrix.
        complex(dp), intent(in) :: b(:) !! Dense RHS vector.
        complex(dp), allocatable, intent(out) :: x(:) !! Direct solution on success.
        integer, intent(out) :: info !! 0 success; negative on invalid/singular input.
        real(dp), intent(in), optional :: pivot_tol !! Absolute minimum pivot magnitude.
        type(splu_complex_factor) :: factor

        call splu_complex(a,factor,info,pivot_tol)
        if (info /= complex_direct_success) then
            allocate(x(0))
            return
        end if
        call splu_complex_solve(factor,b,x,info)
    end subroutine spsolve_complex_csr_vector

    subroutine spsolve_complex_csr_matrix(a, b, x, info, pivot_tol)
        type(csr_complex_matrix), intent(in) :: a !! Square complex CSR matrix.
        complex(dp), intent(in) :: b(:,:) !! Dense RHS block.
        complex(dp), allocatable, intent(out) :: x(:,:) !! Direct solution on success.
        integer, intent(out) :: info !! 0 success; negative on invalid/singular input.
        real(dp), intent(in), optional :: pivot_tol !! Absolute minimum pivot magnitude.
        type(splu_complex_factor) :: factor

        call splu_complex(a,factor,info,pivot_tol)
        if (info /= complex_direct_success) then
            allocate(x(0,0))
            return
        end if
        call splu_complex_solve(factor,b,x,info)
    end subroutine spsolve_complex_csr_matrix

    subroutine spsolve_complex_csc_vector(a, b, x, info, pivot_tol)
        type(csc_complex_matrix), intent(in) :: a !! Square complex CSC matrix.
        complex(dp), intent(in) :: b(:) !! Dense RHS vector.
        complex(dp), allocatable, intent(out) :: x(:) !! Direct solution on success.
        integer, intent(out) :: info !! 0 success; negative on invalid/singular input.
        real(dp), intent(in), optional :: pivot_tol !! Absolute minimum pivot magnitude.
        type(splu_complex_factor) :: factor

        call splu_complex(a,factor,info,pivot_tol)
        if (info /= complex_direct_success) then
            allocate(x(0))
            return
        end if
        call splu_complex_solve(factor,b,x,info)
    end subroutine spsolve_complex_csc_vector

    subroutine spsolve_complex_csc_matrix(a, b, x, info, pivot_tol)
        type(csc_complex_matrix), intent(in) :: a !! Square complex CSC matrix.
        complex(dp), intent(in) :: b(:,:) !! Dense RHS block.
        complex(dp), allocatable, intent(out) :: x(:,:) !! Direct solution on success.
        integer, intent(out) :: info !! 0 success; negative on invalid/singular input.
        real(dp), intent(in), optional :: pivot_tol !! Absolute minimum pivot magnitude.
        type(splu_complex_factor) :: factor

        call splu_complex(a,factor,info,pivot_tol)
        if (info /= complex_direct_success) then
            allocate(x(0,0))
            return
        end if
        call splu_complex_solve(factor,b,x,info)
    end subroutine spsolve_complex_csc_matrix

end module scifort_sparse_linalg_complex_direct
