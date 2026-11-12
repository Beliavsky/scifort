! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Sparse LU factorization for solve_ivp constant CSC Jacobians.  The implementation
! performs Gaussian elimination with partial row pivoting while storing each
! factor row as sorted sparse column/value arrays.  Fill-in is inserted as needed.

module scifort_integrate_sparse_lu
    use scifort_kinds, only : dp
    use scifort_integrate_sparse_jacobian, only : ivp_sparse_matrix, &
        ivp_complex_sparse_matrix
    implicit none
    private

    type :: real_sparse_row
        integer, allocatable :: col(:)
        real(dp), allocatable :: value(:)
    end type real_sparse_row

    type :: complex_sparse_row
        integer, allocatable :: col(:)
        complex(dp), allocatable :: value(:)
    end type complex_sparse_row

    type, public :: ivp_sparse_lu_real
        integer :: n = 0
        integer, allocatable :: piv(:)
        type(real_sparse_row), allocatable :: row(:)
    contains
        procedure :: nnz => real_sparse_lu_nnz
        procedure :: solve => real_sparse_lu_solve
    end type ivp_sparse_lu_real

    type, public :: ivp_sparse_lu_complex
        integer :: n = 0
        integer, allocatable :: piv(:)
        type(complex_sparse_row), allocatable :: row(:)
    contains
        procedure :: nnz => complex_sparse_lu_nnz
        procedure :: solve => complex_sparse_lu_solve
    end type ivp_sparse_lu_complex

    interface sparse_lu_factor
        module procedure factor_real_csc_real
        module procedure factor_real_csc_complex
        module procedure factor_complex_csc_complex
    end interface sparse_lu_factor

    public :: sparse_lu_factor

contains

    subroutine factor_real_csc_real(jac, diagonal, jac_scale, lu, info)
        type(ivp_sparse_matrix), intent(in) :: jac !! square real CSC matrix to factor
        real(dp), intent(in) :: diagonal !! scalar diagonal contribution to the matrix
        real(dp), intent(in) :: jac_scale !! multiplier applied to every stored Jacobian entry
        type(ivp_sparse_lu_real), intent(out) :: lu !! sparse LU factors and row pivots
        integer, intent(out) :: info !! zero on success; pivot column on numerical singularity

        integer :: i
        integer :: j
        integer :: k

        call initialize_real_lu(lu, jac%nrows)
        do i = 1, jac%nrows
            call real_row_add(lu%row(i), i, diagonal)
        end do
        do j = 1, jac%ncols
            do k = jac%indptr(j), jac%indptr(j + 1) - 1
                i = jac%indices(k)
                call real_row_add(lu%row(i), j, jac_scale * jac%data(k))
            end do
        end do
        call factor_real_rows(lu, info)
    end subroutine factor_real_csc_real

    subroutine factor_real_csc_complex(jac, diagonal, jac_scale, lu, info)
        type(ivp_sparse_matrix), intent(in) :: jac !! square real CSC matrix to factor
        complex(dp), intent(in) :: diagonal !! complex diagonal contribution to the matrix
        complex(dp), intent(in) :: jac_scale !! multiplier applied to Jacobian entries
        type(ivp_sparse_lu_complex), intent(out) :: lu !! sparse complex LU factors and pivots
        integer, intent(out) :: info !! zero on success; pivot column on numerical singularity

        integer :: i
        integer :: j
        integer :: k

        call initialize_complex_lu(lu, jac%nrows)
        do i = 1, jac%nrows
            call complex_row_add(lu%row(i), i, diagonal)
        end do
        do j = 1, jac%ncols
            do k = jac%indptr(j), jac%indptr(j + 1) - 1
                i = jac%indices(k)
                call complex_row_add(lu%row(i), j, jac_scale * jac%data(k))
            end do
        end do
        call factor_complex_rows(lu, info)
    end subroutine factor_real_csc_complex

    subroutine factor_complex_csc_complex(jac, diagonal, jac_scale, lu, info)
        type(ivp_complex_sparse_matrix), intent(in) :: jac !! square complex CSC matrix to factor
        complex(dp), intent(in) :: diagonal !! complex diagonal contribution to the matrix
        complex(dp), intent(in) :: jac_scale !! multiplier applied to Jacobian entries
        type(ivp_sparse_lu_complex), intent(out) :: lu !! sparse complex LU factors and pivots
        integer, intent(out) :: info !! zero on success; pivot column on numerical singularity

        integer :: i
        integer :: j
        integer :: k

        call initialize_complex_lu(lu, jac%nrows)
        do i = 1, jac%nrows
            call complex_row_add(lu%row(i), i, diagonal)
        end do
        do j = 1, jac%ncols
            do k = jac%indptr(j), jac%indptr(j + 1) - 1
                i = jac%indices(k)
                call complex_row_add(lu%row(i), j, jac_scale * jac%data(k))
            end do
        end do
        call factor_complex_rows(lu, info)
    end subroutine factor_complex_csc_complex

    subroutine initialize_real_lu(lu, n)
        type(ivp_sparse_lu_real), intent(out) :: lu !! sparse factor object to initialize
        integer, intent(in) :: n !! square matrix order, nonnegative

        integer :: i

        lu%n = n
        allocate(lu%piv(n), lu%row(n))
        lu%piv = [(i, i = 1, n)]
        do i = 1, n
            allocate(lu%row(i)%col(0), lu%row(i)%value(0))
        end do
    end subroutine initialize_real_lu

    subroutine initialize_complex_lu(lu, n)
        type(ivp_sparse_lu_complex), intent(out) :: lu !! sparse complex factor to initialize
        integer, intent(in) :: n !! square matrix order, nonnegative

        integer :: i

        lu%n = n
        allocate(lu%piv(n), lu%row(n))
        lu%piv = [(i, i = 1, n)]
        do i = 1, n
            allocate(lu%row(i)%col(0), lu%row(i)%value(0))
        end do
    end subroutine initialize_complex_lu

    subroutine factor_real_rows(lu, info)
        type(ivp_sparse_lu_real), intent(inout) :: lu !! sparse matrix overwritten by LU factors
        integer, intent(out) :: info !! zero on success; pivot column on numerical singularity

        real(dp) :: aik
        real(dp) :: delta
        real(dp) :: multiplier
        real(dp) :: pivot_abs
        real(dp) :: pivot_value
        integer :: i
        integer :: j
        integer :: k
        integer :: p
        integer :: pivot

        info = 0
        do k = 1, lu%n
            pivot = k
            pivot_abs = 0.0_dp
            do i = k, lu%n
                aik = real_row_get(lu%row(i), k)
                if (abs(aik) > pivot_abs) then
                    pivot_abs = abs(aik)
                    pivot = i
                end if
            end do
            if (pivot_abs <= epsilon(1.0_dp) * max(1.0_dp, pivot_abs)) then
                info = k
                return
            end if
            lu%piv(k) = pivot
            if (pivot /= k) call swap_real_rows(lu%row(k), lu%row(pivot))
            pivot_value = real_row_get(lu%row(k), k)
            do i = k + 1, lu%n
                aik = real_row_get(lu%row(i), k)
                if (aik == 0.0_dp) cycle
                multiplier = aik / pivot_value
                call real_row_set(lu%row(i), k, multiplier)
                do p = 1, size(lu%row(k)%col)
                    j = lu%row(k)%col(p)
                    if (j <= k) cycle
                    delta = -multiplier * lu%row(k)%value(p)
                    call real_row_add(lu%row(i), j, delta)
                end do
            end do
        end do
    end subroutine factor_real_rows

    subroutine factor_complex_rows(lu, info)
        type(ivp_sparse_lu_complex), intent(inout) :: lu !! sparse matrix overwritten by LU factors
        integer, intent(out) :: info !! zero on success; pivot column on numerical singularity

        complex(dp) :: aik
        complex(dp) :: delta
        complex(dp) :: multiplier
        complex(dp) :: pivot_value
        real(dp) :: pivot_abs
        integer :: i
        integer :: j
        integer :: k
        integer :: p
        integer :: pivot

        info = 0
        do k = 1, lu%n
            pivot = k
            pivot_abs = 0.0_dp
            do i = k, lu%n
                aik = complex_row_get(lu%row(i), k)
                if (abs(aik) > pivot_abs) then
                    pivot_abs = abs(aik)
                    pivot = i
                end if
            end do
            if (pivot_abs <= epsilon(1.0_dp) * max(1.0_dp, pivot_abs)) then
                info = k
                return
            end if
            lu%piv(k) = pivot
            if (pivot /= k) call swap_complex_rows(lu%row(k), lu%row(pivot))
            pivot_value = complex_row_get(lu%row(k), k)
            do i = k + 1, lu%n
                aik = complex_row_get(lu%row(i), k)
                if (aik == cmplx(0.0_dp, 0.0_dp, kind=dp)) cycle
                multiplier = aik / pivot_value
                call complex_row_set(lu%row(i), k, multiplier)
                do p = 1, size(lu%row(k)%col)
                    j = lu%row(k)%col(p)
                    if (j <= k) cycle
                    delta = -multiplier * lu%row(k)%value(p)
                    call complex_row_add(lu%row(i), j, delta)
                end do
            end do
        end do
    end subroutine factor_complex_rows

    subroutine real_sparse_lu_solve(self, b, x, info)
        class(ivp_sparse_lu_real), intent(in) :: self !! factored sparse real matrix
        real(dp), intent(in) :: b(:) !! right-hand side, length n
        real(dp), intent(out) :: x(:) !! solution, length n
        integer, intent(out) :: info !! zero on success; diagonal index on singularity

        real(dp) :: diagonal
        real(dp) :: tmp
        integer :: i
        integer :: j
        integer :: k
        integer :: p

        info = 0
        if (size(b) /= self%n .or. size(x) /= self%n) then
            info = -1
            return
        end if
        x = b
        do k = 1, self%n
            if (self%piv(k) /= k) then
                tmp = x(k)
                x(k) = x(self%piv(k))
                x(self%piv(k)) = tmp
            end if
        end do
        do i = 2, self%n
            do p = 1, size(self%row(i)%col)
                j = self%row(i)%col(p)
                if (j >= i) exit
                x(i) = x(i) - self%row(i)%value(p) * x(j)
            end do
        end do
        do i = self%n, 1, -1
            diagonal = 0.0_dp
            do p = 1, size(self%row(i)%col)
                j = self%row(i)%col(p)
                if (j == i) then
                    diagonal = self%row(i)%value(p)
                else if (j > i) then
                    x(i) = x(i) - self%row(i)%value(p) * x(j)
                end if
            end do
            if (diagonal == 0.0_dp) then
                info = i
                return
            end if
            x(i) = x(i) / diagonal
        end do
    end subroutine real_sparse_lu_solve

    subroutine complex_sparse_lu_solve(self, b, x, info)
        class(ivp_sparse_lu_complex), intent(in) :: self !! factored sparse complex matrix
        complex(dp), intent(in) :: b(:) !! complex right-hand side, length n
        complex(dp), intent(out) :: x(:) !! complex solution, length n
        integer, intent(out) :: info !! zero on success; diagonal index on singularity

        complex(dp) :: diagonal
        complex(dp) :: tmp
        integer :: i
        integer :: j
        integer :: k
        integer :: p

        info = 0
        if (size(b) /= self%n .or. size(x) /= self%n) then
            info = -1
            return
        end if
        x = b
        do k = 1, self%n
            if (self%piv(k) /= k) then
                tmp = x(k)
                x(k) = x(self%piv(k))
                x(self%piv(k)) = tmp
            end if
        end do
        do i = 2, self%n
            do p = 1, size(self%row(i)%col)
                j = self%row(i)%col(p)
                if (j >= i) exit
                x(i) = x(i) - self%row(i)%value(p) * x(j)
            end do
        end do
        do i = self%n, 1, -1
            diagonal = cmplx(0.0_dp, 0.0_dp, kind=dp)
            do p = 1, size(self%row(i)%col)
                j = self%row(i)%col(p)
                if (j == i) then
                    diagonal = self%row(i)%value(p)
                else if (j > i) then
                    x(i) = x(i) - self%row(i)%value(p) * x(j)
                end if
            end do
            if (diagonal == cmplx(0.0_dp, 0.0_dp, kind=dp)) then
                info = i
                return
            end if
            x(i) = x(i) / diagonal
        end do
    end subroutine complex_sparse_lu_solve

    integer function real_sparse_lu_nnz(self) result(value)
        class(ivp_sparse_lu_real), intent(in) :: self !! sparse real LU factor object
        integer :: i

        value = 0
        do i = 1, self%n
            value = value + size(self%row(i)%col)
        end do
    end function real_sparse_lu_nnz

    integer function complex_sparse_lu_nnz(self) result(value)
        class(ivp_sparse_lu_complex), intent(in) :: self !! sparse complex LU factor object
        integer :: i

        value = 0
        do i = 1, self%n
            value = value + size(self%row(i)%col)
        end do
    end function complex_sparse_lu_nnz

    real(dp) function real_row_get(row, column) result(value)
        type(real_sparse_row), intent(in) :: row !! sorted sparse row to query
        integer, intent(in) :: column !! one-based column index to query
        integer :: p

        value = 0.0_dp
        do p = 1, size(row%col)
            if (row%col(p) == column) then
                value = row%value(p)
                return
            end if
            if (row%col(p) > column) return
        end do
    end function real_row_get

    complex(dp) function complex_row_get(row, column) result(value)
        type(complex_sparse_row), intent(in) :: row !! sorted sparse complex row to query
        integer, intent(in) :: column !! one-based column index to query
        integer :: p

        value = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do p = 1, size(row%col)
            if (row%col(p) == column) then
                value = row%value(p)
                return
            end if
            if (row%col(p) > column) return
        end do
    end function complex_row_get

    subroutine real_row_set(row, column, value)
        type(real_sparse_row), intent(inout) :: row !! sparse row to update
        integer, intent(in) :: column !! one-based column index to set
        real(dp), intent(in) :: value !! replacement value at the requested column

        integer :: p

        do p = 1, size(row%col)
            if (row%col(p) == column) then
                row%value(p) = value
                return
            end if
            if (row%col(p) > column) exit
        end do
        call real_row_insert(row, column, value)
    end subroutine real_row_set

    subroutine complex_row_set(row, column, value)
        type(complex_sparse_row), intent(inout) :: row !! sparse complex row to update
        integer, intent(in) :: column !! one-based column index to set
        complex(dp), intent(in) :: value !! replacement value at the requested column

        integer :: p

        do p = 1, size(row%col)
            if (row%col(p) == column) then
                row%value(p) = value
                return
            end if
            if (row%col(p) > column) exit
        end do
        call complex_row_insert(row, column, value)
    end subroutine complex_row_set

    subroutine real_row_add(row, column, delta)
        type(real_sparse_row), intent(inout) :: row !! sparse row receiving an increment
        integer, intent(in) :: column !! one-based column index to update
        real(dp), intent(in) :: delta !! value added to the requested entry

        integer :: p

        if (delta == 0.0_dp) return
        do p = 1, size(row%col)
            if (row%col(p) == column) then
                row%value(p) = row%value(p) + delta
                return
            end if
            if (row%col(p) > column) exit
        end do
        call real_row_insert(row, column, delta)
    end subroutine real_row_add

    subroutine complex_row_add(row, column, delta)
        type(complex_sparse_row), intent(inout) :: row !! sparse complex row receiving increment
        integer, intent(in) :: column !! one-based column index to update
        complex(dp), intent(in) :: delta !! complex value added to the requested entry

        integer :: p

        if (delta == cmplx(0.0_dp, 0.0_dp, kind=dp)) return
        do p = 1, size(row%col)
            if (row%col(p) == column) then
                row%value(p) = row%value(p) + delta
                return
            end if
            if (row%col(p) > column) exit
        end do
        call complex_row_insert(row, column, delta)
    end subroutine complex_row_add

    subroutine real_row_insert(row, column, value)
        type(real_sparse_row), intent(inout) :: row !! sorted sparse row receiving a new entry
        integer, intent(in) :: column !! one-based column index to insert
        real(dp), intent(in) :: value !! value of the new entry

        integer, allocatable :: new_col(:)
        real(dp), allocatable :: new_value(:)
        integer :: n
        integer :: pos

        n = size(row%col)
        pos = 1
        do while (pos <= n)
            if (row%col(pos) > column) exit
            pos = pos + 1
        end do
        allocate(new_col(n + 1), new_value(n + 1))
        if (pos > 1) then
            new_col(1:pos - 1) = row%col(1:pos - 1)
            new_value(1:pos - 1) = row%value(1:pos - 1)
        end if
        new_col(pos) = column
        new_value(pos) = value
        if (pos <= n) then
            new_col(pos + 1:n + 1) = row%col(pos:n)
            new_value(pos + 1:n + 1) = row%value(pos:n)
        end if
        call move_alloc(new_col, row%col)
        call move_alloc(new_value, row%value)
    end subroutine real_row_insert

    subroutine complex_row_insert(row, column, value)
        type(complex_sparse_row), intent(inout) :: row !! sorted sparse row receiving a new entry
        integer, intent(in) :: column !! one-based column index to insert
        complex(dp), intent(in) :: value !! complex value of the new entry

        integer, allocatable :: new_col(:)
        complex(dp), allocatable :: new_value(:)
        integer :: n
        integer :: pos

        n = size(row%col)
        pos = 1
        do while (pos <= n)
            if (row%col(pos) > column) exit
            pos = pos + 1
        end do
        allocate(new_col(n + 1), new_value(n + 1))
        if (pos > 1) then
            new_col(1:pos - 1) = row%col(1:pos - 1)
            new_value(1:pos - 1) = row%value(1:pos - 1)
        end if
        new_col(pos) = column
        new_value(pos) = value
        if (pos <= n) then
            new_col(pos + 1:n + 1) = row%col(pos:n)
            new_value(pos + 1:n + 1) = row%value(pos:n)
        end if
        call move_alloc(new_col, row%col)
        call move_alloc(new_value, row%value)
    end subroutine complex_row_insert

    subroutine swap_real_rows(a, b)
        type(real_sparse_row), intent(inout) :: a !! first sparse row exchanged with b
        type(real_sparse_row), intent(inout) :: b !! second sparse row exchanged with a

        integer, allocatable :: col_tmp(:)
        real(dp), allocatable :: value_tmp(:)

        call move_alloc(a%col, col_tmp)
        call move_alloc(a%value, value_tmp)
        call move_alloc(b%col, a%col)
        call move_alloc(b%value, a%value)
        call move_alloc(col_tmp, b%col)
        call move_alloc(value_tmp, b%value)
    end subroutine swap_real_rows

    subroutine swap_complex_rows(a, b)
        type(complex_sparse_row), intent(inout) :: a !! first sparse complex row exchanged with b
        type(complex_sparse_row), intent(inout) :: b !! second sparse complex row exchanged with a

        integer, allocatable :: col_tmp(:)
        complex(dp), allocatable :: value_tmp(:)

        call move_alloc(a%col, col_tmp)
        call move_alloc(a%value, value_tmp)
        call move_alloc(b%col, a%col)
        call move_alloc(b%value, a%value)
        call move_alloc(col_tmp, b%col)
        call move_alloc(value_tmp, b%value)
    end subroutine swap_complex_rows

end module scifort_integrate_sparse_lu
