! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Incomplete LU with zero fill for canonical, one-based, complex CSR matrices.
! The LU factors retain the input's structural pattern, with unit diagonal L.
module scifort_sparse_complex_ilu
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, valid_csr_complex
    implicit none
    private

    integer, parameter, public :: complex_ilu_success = 0
    integer, parameter, public :: complex_ilu_invalid_input = -1
    integer, parameter, public :: complex_ilu_breakdown = -2
    integer, parameter, public :: complex_ilu_bad_diagonal = -3

    type, public :: complex_ilu0_factor
        type(csr_complex_matrix) :: lu        !! Lower multipliers and upper triangular matrix.
        integer, allocatable :: diag_ptr(:) !! One-based diagonal slot in each row.
    end type complex_ilu0_factor

    public :: ilu0_complex, valid_ilu0_complex, ilu0_complex_solve
    public :: ilu0_complex_apply

contains

    pure elemental logical function finite_complex_ilu(z) result(ok)
        complex(dp), intent(in) :: z !! Number to validate.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex_ilu

    pure integer function row_position(a, row, column) result(position)
        type(csr_complex_matrix), intent(in) :: a !! Sorted canonical CSR matrix.
        integer, intent(in) :: row !! One-based row index.
        integer, intent(in) :: column !! One-based column index.
        integer :: low, high, middle

        position = 0
        low = a%indptr(row)
        high = a%indptr(row+1)-1
        do while (low <= high)
            middle = low+(high-low)/2
            if (a%indices(middle) == column) then
                position = middle
                return
            else if (a%indices(middle) < column) then
                low = middle+1
            else
                high = middle-1
            end if
        end do
    end function row_position

    pure logical function valid_ilu0_complex(factor) result(ok)
        type(complex_ilu0_factor), intent(in) :: factor !! Candidate zero-fill LU.
        integer :: i, p

        ok = .false.
        if (.not. valid_csr_complex(factor%lu)) return
        if (factor%lu%n_rows /= factor%lu%n_cols) return
        if (.not. allocated(factor%diag_ptr)) return
        if (size(factor%diag_ptr) /= factor%lu%n_rows) return
        do i = 1, factor%lu%n_rows
            p = factor%diag_ptr(i)
            if (p < factor%lu%indptr(i) .or. p >= factor%lu%indptr(i+1)) return
            if (factor%lu%indices(p) /= i) return
            if (abs(factor%lu%data(p)) <= tiny(1.0_dp)) return
        end do
        ok = .true.
    end function valid_ilu0_complex

    subroutine ilu0_complex(a, factor, info, diagonal_shift)
        type(csr_complex_matrix), intent(in) :: a !! Canonical square complex CSR input.
        type(complex_ilu0_factor), intent(out) :: factor !! Factor or invalid partial factor on failure.
        integer, intent(out) :: info !! 0 success; -1 invalid, -2 nonfinite, -3 bad pivot.
        real(dp), intent(in), optional :: diagonal_shift !! Nonnegative real shift added to the stored diagonal.
        integer :: n, i, p, j, q, target
        complex(dp) :: multiplier, updated
        real(dp) :: shift

        info = complex_ilu_invalid_input
        if (.not. valid_csr_complex(a)) return
        n = a%n_rows
        if (n /= a%n_cols) return
        shift = 0.0_dp
        if (present(diagonal_shift)) shift = diagonal_shift
        if (.not. ieee_is_finite(shift) .or. shift < 0.0_dp) return
        allocate(factor%diag_ptr(n))
        do i = 1, n
            factor%diag_ptr(i) = row_position(a,i,i)
            if (factor%diag_ptr(i) == 0) then
                info = complex_ilu_bad_diagonal
                return
            end if
        end do
        factor%lu = a
        do i = 1, n
            p = factor%diag_ptr(i)
            factor%lu%data(p) = factor%lu%data(p)+cmplx(shift,0.0_dp,dp)
            if (.not. finite_complex_ilu(factor%lu%data(p))) then
                info = complex_ilu_breakdown
                return
            end if
        end do
        do i = 1, n
            do p = factor%lu%indptr(i), factor%diag_ptr(i)-1
                j = factor%lu%indices(p)
                if (abs(factor%lu%data(factor%diag_ptr(j))) <= tiny(1.0_dp)) then
                    info = complex_ilu_bad_diagonal
                    return
                end if
                multiplier = factor%lu%data(p)/factor%lu%data(factor%diag_ptr(j))
                if (.not. finite_complex_ilu(multiplier)) then
                    info = complex_ilu_breakdown
                    return
                end if
                factor%lu%data(p) = multiplier
                do q = factor%diag_ptr(j)+1, factor%lu%indptr(j+1)-1
                    target = row_position(factor%lu,i,factor%lu%indices(q))
                    if (target == 0) cycle
                    updated = factor%lu%data(target)-multiplier*factor%lu%data(q)
                    if (.not. finite_complex_ilu(updated)) then
                        info = complex_ilu_breakdown
                        return
                    end if
                    factor%lu%data(target) = updated
                end do
            end do
            if (abs(factor%lu%data(factor%diag_ptr(i))) <= tiny(1.0_dp)) then
                info = complex_ilu_bad_diagonal
                return
            end if
        end do
        info = complex_ilu_success
    end subroutine ilu0_complex

    ! Fast triangular application for internal solvers. The factor was checked
    ! once by their prepare() routine. The adjoint path solves (LU)^H*x=b.
    pure subroutine ilu0_complex_apply(factor, b, x, adjoint)
        type(complex_ilu0_factor), intent(in) :: factor !! Prevalidated ILU(0) factor.
        complex(dp), intent(in) :: b(:) !! Right-hand side, size n.
        complex(dp), intent(out) :: x(:) !! Solution vector, size n.
        logical, intent(in), optional :: adjoint !! Solve with conjugate transpose if true.
        integer :: i, j, p, n
        logical :: transposed

        transposed = .false.
        if (present(adjoint)) transposed = adjoint
        n = factor%lu%n_rows
        x = b
        if (.not. transposed) then
            ! L*y=b, with unit diagonal L; then U*x=y.
            do i = 1, n
                do p = factor%lu%indptr(i), factor%diag_ptr(i)-1
                    j = factor%lu%indices(p)
                    x(i) = x(i)-factor%lu%data(p)*x(j)
                end do
            end do
            do i = n, 1, -1
                do p = factor%diag_ptr(i)+1, factor%lu%indptr(i+1)-1
                    j = factor%lu%indices(p)
                    x(i) = x(i)-factor%lu%data(p)*x(j)
                end do
                x(i) = x(i)/factor%lu%data(factor%diag_ptr(i))
            end do
        else
            ! U^H*y=b: scatter updates along columns of U^H.
            do i = 1, n
                x(i) = x(i)/conjg(factor%lu%data(factor%diag_ptr(i)))
                do p = factor%diag_ptr(i)+1, factor%lu%indptr(i+1)-1
                    j = factor%lu%indices(p)
                    x(j) = x(j)-conjg(factor%lu%data(p))*x(i)
                end do
            end do
            ! L^H*x=y: unit diagonal, descending scatter elimination.
            do i = n, 1, -1
                do p = factor%lu%indptr(i), factor%diag_ptr(i)-1
                    j = factor%lu%indices(p)
                    x(j) = x(j)-conjg(factor%lu%data(p))*x(i)
                end do
            end do
        end if
    end subroutine ilu0_complex_apply

    subroutine ilu0_complex_solve(factor, b, x, info, adjoint)
        type(complex_ilu0_factor), intent(in) :: factor !! Complete factorization.
        complex(dp), intent(in) :: b(:) !! Right-hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Solution, empty on invalid input.
        integer, intent(out) :: info !! Zero success, negative on invalid input/nonfinite result.
        logical, intent(in), optional :: adjoint !! Solve (LU)^H when present and true.

        info = complex_ilu_invalid_input
        allocate(x(0))
        if (.not. valid_ilu0_complex(factor)) return
        if (size(b) /= factor%lu%n_rows) return
        if (.not. all(finite_complex_ilu(b))) return
        deallocate(x)
        allocate(x(size(b)))
        call ilu0_complex_apply(factor,b,x,adjoint)
        if (.not. all(finite_complex_ilu(x))) then
            info = complex_ilu_breakdown
            return
        end if
        info = complex_ilu_success
    end subroutine ilu0_complex_solve

end module scifort_sparse_complex_ilu
