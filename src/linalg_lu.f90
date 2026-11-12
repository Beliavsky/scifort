! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rectangular pivoted LU decomposition, including singular matrices.
! Uses independent partial-pivoting elimination. One-based permutation
! indices perm satisfy A = L(perm,:) * U, in matrix form A = P L U.
module scifort_linalg_lu
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: lu, lu_indices, lu_permuted

    interface lu
        module procedure lu_real
        module procedure lu_complex
    end interface
    interface lu_indices
        module procedure lu_indices_real
        module procedure lu_indices_complex
    end interface
    interface lu_permuted
        module procedure lu_permuted_real
        module procedure lu_permuted_complex
    end interface
contains
    subroutine lu_indices_real(a, perm, l, u, status, check_finite)
        real(dp), intent(in) :: a(:, :) !! m-by-n input, including rectangular and singular matrices
        integer, intent(out) :: perm(:) !! one-based row indices: A = L(perm,:) * U, length m
        real(dp), intent(out) :: l(:, :) !! unit lower trapezoid, shape (m,min(m,n))
        real(dp), intent(out) :: u(:, :) !! upper trapezoid, shape (min(m,n),n)
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        real(dp), allocatable :: work(:, :)
        real(dp) :: temp
        integer, allocatable :: rows(:)
        integer :: m, n, k, i, j, p, original
        logical :: check

        perm = 0
        l = 0.0_dp
        u = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (size(perm) /= m .or. size(l,1) /= m .or. size(l,2) /= k) return
        if (size(u,1) /= k .or. size(u,2) /= n) return
        check = .true.
        if (present(check_finite)) check = check_finite
        if (check) then
            if (.not. all(ieee_is_finite(a))) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        allocate(work(m,n),rows(m))
        work = a
        do i = 1, m
            rows(i) = i
        end do
        do j = 1, k
            p = j-1 + maxloc(abs(work(j:m,j)),dim=1)
            if (p /= j) then
                original = rows(j)
                rows(j) = rows(p)
                rows(p) = original
                do i = 1, n
                    temp = work(j,i)
                    work(j,i) = work(p,i)
                    work(p,i) = temp
                end do
            end if
            ! A zero pivot means that the column has no remaining multipliers.
            ! Singular input is accepted, as in scipy.linalg.lu.
            if (work(j,j) == 0.0_dp) cycle
            do i = j+1, m
                work(i,j) = work(i,j)/work(j,j)
            end do
            do i = j+1, m
                work(i,j+1:n) = work(i,j+1:n) - work(i,j)*work(j,j+1:n)
            end do
        end do
        if (.not. (all(ieee_is_finite(work)))) then
            status = linalg_status_no_convergence
            return
        end if
        do i = 1, m
            perm(rows(i)) = i
        end do
        do j = 1, k
            l(j,j) = 1.0_dp
            if (j < m) l(j+1:m,j) = work(j+1:m,j)
            u(j,j:n) = work(j,j:n)
        end do
        status = linalg_status_success
    end subroutine lu_indices_real

    subroutine lu_real(a, p, l, u, status, check_finite)
        real(dp), intent(in) :: a(:, :) !! m-by-n matrix
        real(dp), intent(out) :: p(:, :) !! m-by-m permutation matrix, A = P L U
        real(dp), intent(out) :: l(:, :) !! m-by-min(m,n) lower factor
        real(dp), intent(out) :: u(:, :) !! min(m,n)-by-n upper factor
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        integer, allocatable :: perm(:)
        integer :: m, i

        p = 0.0_dp
        l = 0.0_dp
        u = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        if (size(p,1) /= m .or. size(p,2) /= m) return
        allocate(perm(m))
        call lu_indices_real(a,perm,l,u,status,check_finite)
        if (status /= linalg_status_success) return
        do i = 1, m
            p(i,perm(i)) = 1.0_dp
        end do
    end subroutine lu_real

    subroutine lu_permuted_real(a, pl, u, status, check_finite)
        real(dp), intent(in) :: a(:, :) !! m-by-n matrix
        real(dp), intent(out) :: pl(:, :) !! m-by-min(m,n) matrix P*L
        real(dp), intent(out) :: u(:, :) !! min(m,n)-by-n upper factor
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        real(dp), allocatable :: l(:, :)
        integer, allocatable :: perm(:)
        integer :: m, n, k, i

        pl = 0.0_dp
        u = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (size(pl,1) /= m .or. size(pl,2) /= k) return
        if (size(u,1) /= k .or. size(u,2) /= n) return
        allocate(l(m,k),perm(m))
        call lu_indices_real(a,perm,l,u,status,check_finite)
        if (status /= linalg_status_success) return
        do i = 1, m
            pl(i,:) = l(perm(i),:)
        end do
    end subroutine lu_permuted_real
    subroutine lu_indices_complex(a, perm, l, u, status, check_finite)
        complex(dp), intent(in) :: a(:, :) !! m-by-n input, including rectangular and singular matrices
        integer, intent(out) :: perm(:) !! one-based row indices: A = L(perm,:) * U, length m
        complex(dp), intent(out) :: l(:, :) !! unit lower trapezoid, shape (m,min(m,n))
        complex(dp), intent(out) :: u(:, :) !! upper trapezoid, shape (min(m,n),n)
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        complex(dp), allocatable :: work(:, :)
        complex(dp) :: temp
        integer, allocatable :: rows(:)
        integer :: m, n, k, i, j, p, original
        logical :: check

        perm = 0
        l = cmplx(0.0_dp,0.0_dp,dp)
        u = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (size(perm) /= m .or. size(l,1) /= m .or. size(l,2) /= k) return
        if (size(u,1) /= k .or. size(u,2) /= n) return
        check = .true.
        if (present(check_finite)) check = check_finite
        if (check) then
            if (.not. all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        allocate(work(m,n),rows(m))
        work = a
        do i = 1, m
            rows(i) = i
        end do
        do j = 1, k
            ! LAPACK complex pivoting compares |real(z)| + |imag(z)|.
            p = j-1 + maxloc(abs(real(work(j:m,j),dp)) + &
                abs(aimag(work(j:m,j))),dim=1)
            if (p /= j) then
                original = rows(j)
                rows(j) = rows(p)
                rows(p) = original
                do i = 1, n
                    temp = work(j,i)
                    work(j,i) = work(p,i)
                    work(p,i) = temp
                end do
            end if
            ! A zero pivot means that the column has no remaining multipliers.
            ! Singular input is accepted, as in scipy.linalg.lu.
            if (work(j,j) == cmplx(0.0_dp,0.0_dp,dp)) cycle
            do i = j+1, m
                work(i,j) = work(i,j)/work(j,j)
            end do
            do i = j+1, m
                work(i,j+1:n) = work(i,j+1:n) - work(i,j)*work(j,j+1:n)
            end do
        end do
        if (.not. (all(ieee_is_finite(real(work,dp))) .and. all(ieee_is_finite(aimag(work))))) then
            status = linalg_status_no_convergence
            return
        end if
        do i = 1, m
            perm(rows(i)) = i
        end do
        do j = 1, k
            l(j,j) = cmplx(1.0_dp,0.0_dp,dp)
            if (j < m) l(j+1:m,j) = work(j+1:m,j)
            u(j,j:n) = work(j,j:n)
        end do
        status = linalg_status_success
    end subroutine lu_indices_complex

    subroutine lu_complex(a, p, l, u, status, check_finite)
        complex(dp), intent(in) :: a(:, :) !! m-by-n matrix
        complex(dp), intent(out) :: p(:, :) !! m-by-m permutation matrix, A = P L U
        complex(dp), intent(out) :: l(:, :) !! m-by-min(m,n) lower factor
        complex(dp), intent(out) :: u(:, :) !! min(m,n)-by-n upper factor
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        integer, allocatable :: perm(:)
        integer :: m, i

        p = cmplx(0.0_dp,0.0_dp,dp)
        l = cmplx(0.0_dp,0.0_dp,dp)
        u = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(a,1)
        if (size(p,1) /= m .or. size(p,2) /= m) return
        allocate(perm(m))
        call lu_indices_complex(a,perm,l,u,status,check_finite)
        if (status /= linalg_status_success) return
        do i = 1, m
            p(i,perm(i)) = cmplx(1.0_dp,0.0_dp,dp)
        end do
    end subroutine lu_complex

    subroutine lu_permuted_complex(a, pl, u, status, check_finite)
        complex(dp), intent(in) :: a(:, :) !! m-by-n matrix
        complex(dp), intent(out) :: pl(:, :) !! m-by-min(m,n) matrix P*L
        complex(dp), intent(out) :: u(:, :) !! min(m,n)-by-n upper factor
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! finite input validation, default true
        complex(dp), allocatable :: l(:, :)
        integer, allocatable :: perm(:)
        integer :: m, n, k, i

        pl = cmplx(0.0_dp,0.0_dp,dp)
        u = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (size(pl,1) /= m .or. size(pl,2) /= k) return
        if (size(u,1) /= k .or. size(u,2) /= n) return
        allocate(l(m,k),perm(m))
        call lu_indices_complex(a,perm,l,u,status,check_finite)
        if (status /= linalg_status_success) return
        do i = 1, m
            pl(i,:) = l(perm(i),:)
        end do
    end subroutine lu_permuted_complex
end module scifort_linalg_lu
