! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Determinant and inverse of dense square real or complex matrices.
! Singular matrices have determinant zero; inv reports a singular status.
module scifort_linalg_det_inv
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_lu, only : lu_indices
    use scifort_linalg_dense, only : lu_factor, lu_solve
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: det, inv
    interface det
        module procedure det_real
        module procedure det_complex
    end interface
    interface inv
        module procedure inv_real
        module procedure inv_complex
    end interface
contains
    subroutine det_real(a, value, status, check_finite)
        real(dp), intent(in) :: a(:, :) !! square coefficient matrix, possibly singular or empty
        real(dp), intent(out) :: value !! determinant, exactly zero for singular matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        real(dp), allocatable :: l(:, :), u(:, :)
        integer, allocatable :: perm(:)
        integer :: n, i, j, parity

        value = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (size(a,2) /= n) return
        if (n == 0) then
            value = 1.0_dp
            status = linalg_status_success
            return
        end if
        allocate(l(n,n),u(n,n),perm(n))
        call lu_indices(a,perm,l,u,status,check_finite)
        if (status /= linalg_status_success) return
        ! The parity of the one-based inverse row permutation determines sign.
        parity = 0
        do i = 1, n-1
            do j = i+1, n
                if (perm(j) < perm(i)) parity = 1-parity
            end do
        end do
        value = 1.0_dp
        if (parity == 1) value = -value
        do i = 1, n
            if (u(i,i) == 0.0_dp) then
                value = 0.0_dp
                return
            end if
        end do
        do i = 1, n
            value = value*u(i,i)
        end do
        if (.not. (ieee_is_finite(value))) then
            value = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine det_real

    subroutine inv_real(a, inverse, status, check_finite)
        real(dp), intent(in) :: a(:, :) !! square nonsingular matrix, possibly empty
        real(dp), intent(out) :: inverse(:, :) !! inverse A**(-1), shape of a; zero on error
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        real(dp), allocatable :: factors(:, :), identity(:, :)
        integer, allocatable :: pivots(:)
        integer :: n, i
        logical :: check

        inverse = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (size(a,2) /= n .or. size(inverse,1) /= n .or. size(inverse,2) /= n) return
        check = .true.
        if (present(check_finite)) check = check_finite
        if (check) then
            if (.not. (all(ieee_is_finite(a)))) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        allocate(factors(n,n),identity(n,n),pivots(n))
        identity = 0.0_dp
        do i = 1, n
            identity(i,i) = 1.0_dp
        end do
        call lu_factor(a,factors,pivots,status)
        if (status /= linalg_status_success) return
        call lu_solve(factors,pivots,identity,inverse,status)
        if (status /= linalg_status_success) inverse = 0.0_dp
    end subroutine inv_real
    subroutine det_complex(a, value, status, check_finite)
        complex(dp), intent(in) :: a(:, :) !! square coefficient matrix, possibly singular or empty
        complex(dp), intent(out) :: value !! determinant, exactly zero for singular matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        complex(dp), allocatable :: l(:, :), u(:, :)
        integer, allocatable :: perm(:)
        integer :: n, i, j, parity

        value = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (size(a,2) /= n) return
        if (n == 0) then
            value = cmplx(1.0_dp,0.0_dp,dp)
            status = linalg_status_success
            return
        end if
        allocate(l(n,n),u(n,n),perm(n))
        call lu_indices(a,perm,l,u,status,check_finite)
        if (status /= linalg_status_success) return
        ! The parity of the one-based inverse row permutation determines sign.
        parity = 0
        do i = 1, n-1
            do j = i+1, n
                if (perm(j) < perm(i)) parity = 1-parity
            end do
        end do
        value = cmplx(1.0_dp,0.0_dp,dp)
        if (parity == 1) value = -value
        do i = 1, n
            if (u(i,i) == cmplx(0.0_dp,0.0_dp,dp)) then
                value = cmplx(0.0_dp,0.0_dp,dp)
                return
            end if
        end do
        do i = 1, n
            value = value*u(i,i)
        end do
        if (.not. (ieee_is_finite(real(value,dp)) .and. ieee_is_finite(aimag(value)))) then
            value = cmplx(0.0_dp,0.0_dp,dp)
            status = linalg_status_no_convergence
        end if
    end subroutine det_complex

    subroutine inv_complex(a, inverse, status, check_finite)
        complex(dp), intent(in) :: a(:, :) !! square nonsingular matrix, possibly empty
        complex(dp), intent(out) :: inverse(:, :) !! inverse A**(-1), shape of a; zero on error
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        complex(dp), allocatable :: factors(:, :), identity(:, :)
        integer, allocatable :: pivots(:)
        integer :: n, i
        logical :: check

        inverse = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (size(a,2) /= n .or. size(inverse,1) /= n .or. size(inverse,2) /= n) return
        check = .true.
        if (present(check_finite)) check = check_finite
        if (check) then
            if (.not. (all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a))))) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        allocate(factors(n,n),identity(n,n),pivots(n))
        identity = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1, n
            identity(i,i) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        call lu_factor(a,factors,pivots,status)
        if (status /= linalg_status_success) return
        call lu_solve(factors,pivots,identity,inverse,status)
        if (status /= linalg_status_success) inverse = cmplx(0.0_dp,0.0_dp,dp)
    end subroutine inv_complex
end module scifort_linalg_det_inv
