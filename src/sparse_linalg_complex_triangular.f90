! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dependency-free complex binary64 triangular solves for one-based canonical
! CSR/CSC matrices.  The transpose/adjoint modes are extensions of SciPy's
! sparse spsolve_triangular interface; the input is never modified.
module scifort_sparse_linalg_complex_triangular
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, csc_complex_to_csr
    implicit none
    private

    integer, parameter, public :: complex_triangular_success = 0
    integer, parameter, public :: complex_triangular_invalid_input = -1
    integer, parameter, public :: complex_triangular_breakdown = -2
    integer, parameter, public :: complex_triangular_bad_diagonal = -3

    public :: spsolve_triangular_complex

    interface spsolve_triangular_complex
        module procedure triangular_csr_vector, triangular_csr_matrix
        module procedure triangular_csc_vector, triangular_csc_matrix
    end interface spsolve_triangular_complex

contains

    pure elemental logical function finite_complex(value) result(ok)
        complex(dp), intent(in) :: value !! Complex candidate to validate.
        ok = ieee_is_finite(real(value,dp)) .and. ieee_is_finite(aimag(value))
    end function finite_complex

    subroutine triangular_csr_vector(a, b, x, info, lower, unit_diagonal, trans)
        type(csr_complex_matrix), intent(in) :: a !! Square canonical CSR coefficient matrix.
        complex(dp), intent(in) :: b(:) !! Finite RHS vector of length n.
        complex(dp), allocatable, intent(out) :: x(:) !! Length n on success; empty on failure.
        integer, intent(out) :: info !! 0 success, -1 invalid, -2 overflow, -3 bad diagonal.
        logical, intent(in), optional :: lower !! Lower triangle if true, default true.
        logical, intent(in), optional :: unit_diagonal !! Treat diagonal as one, default false.
        character(len=*), intent(in), optional :: trans !! N/T/H/C for normal/transpose/adjoint.
        complex(dp), allocatable :: rhs(:,:), solution(:,:)

        allocate(rhs(size(b),1))
        rhs(:,1) = b
        call triangular_csr_matrix(a,rhs,solution,info,lower,unit_diagonal,trans)
        if (info == complex_triangular_success) then
            x = solution(:,1)
        else
            allocate(x(0))
        end if
    end subroutine triangular_csr_vector

    subroutine triangular_csc_vector(a, b, x, info, lower, unit_diagonal, trans)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC coefficient matrix.
        complex(dp), intent(in) :: b(:) !! Finite RHS vector of length n.
        complex(dp), allocatable, intent(out) :: x(:) !! Length n on success; empty on failure.
        integer, intent(out) :: info !! 0 success, -1 invalid, -2 overflow, -3 bad diagonal.
        logical, intent(in), optional :: lower !! Lower triangle if true, default true.
        logical, intent(in), optional :: unit_diagonal !! Treat diagonal as one, default false.
        character(len=*), intent(in), optional :: trans !! N/T/H/C for normal/transpose/adjoint.
        complex(dp), allocatable :: rhs(:,:), solution(:,:)

        allocate(rhs(size(b),1))
        rhs(:,1) = b
        call triangular_csc_matrix(a,rhs,solution,info,lower,unit_diagonal,trans)
        if (info == complex_triangular_success) then
            x = solution(:,1)
        else
            allocate(x(0))
        end if
    end subroutine triangular_csc_vector

    subroutine triangular_csc_matrix(a, b, x, info, lower, unit_diagonal, trans)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC coefficient matrix.
        complex(dp), intent(in) :: b(:,:) !! Finite RHS block (n,nrhs).
        complex(dp), allocatable, intent(out) :: x(:,:) !! Shape (n,nrhs) on success.
        integer, intent(out) :: info !! 0 success, -1 invalid, -2 overflow, -3 bad diagonal.
        logical, intent(in), optional :: lower !! Lower triangle if true, default true.
        logical, intent(in), optional :: unit_diagonal !! Treat diagonal as one, default false.
        character(len=*), intent(in), optional :: trans !! N/T/H/C for normal/transpose/adjoint.
        type(csr_complex_matrix) :: csr
        integer :: status

        info = complex_triangular_invalid_input
        allocate(x(0,0))
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        deallocate(x)
        call triangular_csr_matrix(csr,b,x,info,lower,unit_diagonal,trans)
    end subroutine triangular_csc_matrix

    subroutine triangular_csr_matrix(a, b, x, info, lower, unit_diagonal, trans)
        type(csr_complex_matrix), intent(in) :: a !! Square canonical CSR coefficient matrix.
        complex(dp), intent(in) :: b(:,:) !! Finite RHS block (n,nrhs).
        complex(dp), allocatable, intent(out) :: x(:,:) !! Shape (n,nrhs) on success.
        integer, intent(out) :: info !! 0 success, -1 invalid, -2 overflow, -3 bad diagonal.
        logical, intent(in), optional :: lower !! Lower triangle if true, default true.
        logical, intent(in), optional :: unit_diagonal !! Treat diagonal as one, default false.
        character(len=*), intent(in), optional :: trans !! N/T/H/C for normal/transpose/adjoint.
        integer :: n, i, j, p, first, last, step
        logical :: is_lower, is_unit, is_adj, is_trans, has_diagonal
        character(len=1) :: operation
        complex(dp) :: diagonal, coefficient

        info = complex_triangular_invalid_input
        allocate(x(0,0))
        if (.not. valid_csr_complex(a)) return
        n = a%n_rows
        if (a%n_cols /= n .or. size(b,1) /= n) return
        if (.not. all(finite_complex(b))) return
        is_lower = .true.
        if (present(lower)) is_lower = lower
        is_unit = .false.
        if (present(unit_diagonal)) is_unit = unit_diagonal
        operation = 'N'
        if (present(trans)) then
            if (len_trim(trans) /= 1) return
            operation = trans(1:1)
        end if
        if (operation /= 'N' .and. operation /= 'T' .and. &
            operation /= 'H' .and. operation /= 'C') return
        is_trans = operation /= 'N'
        is_adj = operation == 'H' .or. operation == 'C'

        ! Reject non-triangular nonzero entries; stored zeros are harmless.
        ! Check every pivot before modifying the result (including nrhs=0).
        do i = 1,n
            diagonal = cmplx(1.0_dp,0.0_dp,dp)
            has_diagonal = is_unit
            do p = a%indptr(i), a%indptr(i+1)-1
                j = a%indices(p)
                if ((is_lower .and. j > i) .or. (.not. is_lower .and. j < i)) then
                    if (a%data(p) /= cmplx(0.0_dp,0.0_dp,dp)) return
                end if
                if (j == i .and. .not. is_unit) then
                    diagonal = a%data(p)
                    has_diagonal = .true.
                end if
            end do
            if (.not. has_diagonal .or. diagonal == cmplx(0.0_dp,0.0_dp,dp)) then
                info = complex_triangular_bad_diagonal
                return
            end if
        end do

        if (is_lower .neqv. is_trans) then
            first = 1
            last = n
            step = 1
        else
            first = n
            last = 1
            step = -1
        end if
        deallocate(x)
        allocate(x(n,size(b,2)))
        x = b
        do i = first,last,step
            diagonal = cmplx(1.0_dp,0.0_dp,dp)
            do p = a%indptr(i),a%indptr(i+1)-1
                j = a%indices(p)
                if (j == i) then
                    if (.not. is_unit) diagonal = a%data(p)
                else if (.not. is_trans) then
                    if ((is_lower .and. j < i) .or. (.not. is_lower .and. j > i)) then
                        x(i,:) = x(i,:) - a%data(p)*x(j,:)
                    end if
                end if
            end do
            if (is_adj) diagonal = conjg(diagonal)
            x(i,:) = x(i,:)/diagonal
            if (.not. all(finite_complex(x(i,:)))) exit
            if (is_trans) then
                do p = a%indptr(i),a%indptr(i+1)-1
                    j = a%indices(p)
                    if ((is_lower .and. j < i) .or. (.not. is_lower .and. j > i)) then
                        coefficient = a%data(p)
                        if (is_adj) coefficient = conjg(coefficient)
                        x(j,:) = x(j,:) - coefficient*x(i,:)
                        if (.not. all(finite_complex(x(j,:)))) exit
                    end if
                end do
            end if
        end do
        if (.not. all(finite_complex(x))) then
            deallocate(x)
            allocate(x(0,0))
            info = complex_triangular_breakdown
            return
        end if
        info = complex_triangular_success
    end subroutine triangular_csr_matrix

end module scifort_sparse_linalg_complex_triangular
