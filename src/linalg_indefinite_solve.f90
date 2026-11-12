! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Solves for symmetric/Hermitian indefinite matrices via SciFort's LDL factors.
! Independent triangular and 1x1/2x2 block substitutions; no external BLAS.
module scifort_linalg_indefinite_solve
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_ldl, only : ldl
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: ldl_solve, solve_indefinite

    ! Reuse (lu, d, perm) from ldl for any number of RHS matrices.
    interface ldl_solve
        module procedure ldl_solve_real_vector
        module procedure ldl_solve_real_matrix
        module procedure ldl_solve_complex_vector
        module procedure ldl_solve_complex_matrix
    end interface ldl_solve

    ! Factor and solve in one call. Equivalent to the sym/her specializations
    ! of scipy.linalg.solve, not a replacement for its general solve mode.
    interface solve_indefinite
        module procedure solve_indefinite_real_vector
        module procedure solve_indefinite_real_matrix
        module procedure solve_indefinite_complex_vector
        module procedure solve_indefinite_complex_matrix
    end interface solve_indefinite

contains

    subroutine ldl_solve_real_vector(lu, d, perm, b, x, status, lower)
        real(dp), intent(in) :: lu(:, :) !! unpermuted LDL outer factor, shape (n,n)
        real(dp), intent(in) :: d(:, :) !! 1x1/2x2-block diagonal factor, shape (n,n)
        integer, intent(in) :: perm(:) !! one-based row permutation from ldl, length n
        real(dp), intent(in) :: b(:) !! right-hand-side vector, length n
        real(dp), intent(out) :: x(:) !! solution vector; zero on error
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! true for lower LDL (default)
        real(dp), allocatable :: bb(:, :), xx(:, :)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(b)
        if (size(x) /= n) return
        allocate(bb(n,1),xx(n,1))
        bb(:,1) = b
        call ldl_solve_real_matrix(lu,d,perm,bb,xx,status,lower)
        if (status == linalg_status_success) x = xx(:,1)
    end subroutine ldl_solve_real_vector

    subroutine ldl_solve_real_matrix(lu, d, perm, b, x, status, lower)
        real(dp), intent(in) :: lu(:, :) !! unpermuted LDL outer factor, shape (n,n)
        real(dp), intent(in) :: d(:, :) !! 1x1/2x2-block diagonal factor, shape (n,n)
        integer, intent(in) :: perm(:) !! one-based row permutation from ldl, length n
        real(dp), intent(in) :: b(:, :) !! RHS matrix, shape (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solutions, shape (n,nrhs); zero on error
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! true for lower LDL (default)
        complex(dp), allocatable :: xc(:, :)
        integer :: n, nrhs

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(lu,1)
        nrhs = size(b,2)
        if (size(lu,2) /= n .or. size(d,1) /= n .or. size(d,2) /= n) return
        if (size(perm) /= n .or. size(b,1) /= n) return
        if (size(x,1) /= n .or. size(x,2) /= nrhs) return
        allocate(xc(n,nrhs))
        call ldl_solve_complex_matrix(cmplx(lu,kind=dp),cmplx(d,kind=dp),perm, &
            cmplx(b,kind=dp),xc,status,lower=lower,hermitian=.false.)
        if (status == linalg_status_success) x = real(xc,dp)
    end subroutine ldl_solve_real_matrix

    subroutine ldl_solve_complex_vector(lu, d, perm, b, x, status, lower, hermitian)
        complex(dp), intent(in) :: lu(:, :) !! unpermuted LDL outer factor, shape (n,n)
        complex(dp), intent(in) :: d(:, :) !! 1x1/2x2-block diagonal factor, shape (n,n)
        integer, intent(in) :: perm(:) !! one-based row permutation from ldl, length n
        complex(dp), intent(in) :: b(:) !! RHS vector, length n
        complex(dp), intent(out) :: x(:) !! solution vector; zero on error
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! true for lower LDL (default)
        logical, intent(in), optional :: hermitian !! conjugate transpose (default true)
        complex(dp), allocatable :: bb(:, :), xx(:, :)
        integer :: n

        x = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(b)
        if (size(x) /= n) return
        allocate(bb(n,1),xx(n,1))
        bb(:,1) = b
        call ldl_solve_complex_matrix(lu,d,perm,bb,xx,status,lower,hermitian)
        if (status == linalg_status_success) x = xx(:,1)
    end subroutine ldl_solve_complex_vector

    subroutine ldl_solve_complex_matrix(lu, d, perm, b, x, status, lower, hermitian)
        complex(dp), intent(in) :: lu(:, :) !! outer factor from ldl, shape (n,n)
        complex(dp), intent(in) :: d(:, :) !! block-diagonal matrix from ldl, shape (n,n)
        integer, intent(in) :: perm(:) !! one-based row indices making lu(perm,:) triangular
        complex(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solutions, shape (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! lower outer triangle (default true)
        logical, intent(in), optional :: hermitian !! use adjoint (default true), else transpose
        complex(dp), allocatable :: w(:, :)
        complex(dp) :: a11, a12, a21, a22, det, z1, z2, factor
        real(dp) :: scale, tolerance
        integer :: i, j, k, n, nrhs, rhs
        logical :: lo, conj_trans, block_pair

        x = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(lu,1)
        nrhs = size(b,2)
        if (size(lu,2) /= n .or. size(d,1) /= n .or. size(d,2) /= n) return
        if (size(perm) /= n .or. size(b,1) /= n) return
        if (size(x,1) /= n .or. size(x,2) /= nrhs) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        lo = .true.
        if (present(lower)) lo = lower
        conj_trans = .true.
        if (present(hermitian)) conj_trans = hermitian
        if (any(perm < 1) .or. any(perm > n)) return
        do i = 1, n
            if (count(perm == i) /= 1) return
        end do
        status = linalg_status_no_convergence
        if (.not. finite_complex_matrix(lu) .or. .not. finite_complex_matrix(d)) return
        if (.not. finite_complex_matrix(b)) return

        ! Factor structure is checked before substituting. In particular,
        ! no two 2x2 blocks may share a diagonal entry.
        status = linalg_status_invalid_shape
        do i = 1, n
            if (abs(lu(perm(i),i)-cmplx(1.0_dp,0.0_dp,dp)) > &
                64.0_dp*epsilon(1.0_dp)) return
            do j = 1, n
                if (abs(i-j) > 1 .and. abs(d(i,j)) /= 0.0_dp) return
                if (lo .and. j > i .and. abs(lu(perm(i),j)) /= 0.0_dp) return
                if (.not. lo .and. j < i .and. abs(lu(perm(i),j)) /= 0.0_dp) return
            end do
        end do
        do i = 1, n-1
            if ((abs(d(i,i+1)) == 0.0_dp) .neqv. (abs(d(i+1,i)) == 0.0_dp)) return
        end do

        ! P A P^T = T D T^H (or T D T^T), with T = lu(perm,:).
        ! First solve T*y = P*b, using the *selected* triangular factor.
        allocate(w(n,nrhs))
        w = b(perm,:)
        if (lo) then
            do i = 1, n
                do j = 1, i-1
                    w(i,:) = w(i,:) - lu(perm(i),j)*w(j,:)
                end do
            end do
        else
            do i = n, 1, -1
                do j = i+1, n
                    w(i,:) = w(i,:) - lu(perm(i),j)*w(j,:)
                end do
            end do
        end if

        ! Solve D*z = y. Unlike a diagonal solve, 2x2 pivots must be
        ! inverted as coupled systems, including zero diagonal elements.
        status = linalg_status_singular
        k = 1
        do while (k <= n)
            block_pair = .false.
            if (k < n) block_pair = abs(d(k+1,k)) > 0.0_dp
            if (block_pair) then
                if (k+1 < n) then
                    if (abs(d(k+1,k+2)) > 0.0_dp) then
                        status = linalg_status_invalid_shape
                        return
                    end if
                end if
                scale = max(abs(d(k,k)),abs(d(k,k+1)),abs(d(k+1,k)),abs(d(k+1,k+1)))
                if (scale == 0.0_dp) return
                a11 = d(k,k)/scale
                a12 = d(k,k+1)/scale
                a21 = d(k+1,k)/scale
                a22 = d(k+1,k+1)/scale
                det = a11*a22-a12*a21
                tolerance = 8.0_dp*epsilon(1.0_dp)* &
                    max(abs(a11*a22),abs(a12*a21))
                if (abs(det) <= tolerance) return
                do rhs = 1, nrhs
                    z1 = w(k,rhs)/scale
                    z2 = w(k+1,rhs)/scale
                    w(k,rhs) = (a22*z1-a12*z2)/det
                    w(k+1,rhs) = (a11*z2-a21*z1)/det
                end do
                k = k+2
            else
                if (abs(d(k,k)) == 0.0_dp) return
                w(k,:) = w(k,:)/d(k,k)
                k = k+1
            end if
        end do

        ! Final solve T^H*x_perm = z or T^T*x_perm = z.
        if (lo) then
            do i = n, 1, -1
                do j = i+1, n
                    factor = lu(perm(j),i)
                    if (conj_trans) factor = conjg(factor)
                    w(i,:) = w(i,:) - factor*w(j,:)
                end do
            end do
        else
            do i = 1, n
                do j = 1, i-1
                    factor = lu(perm(j),i)
                    if (conj_trans) factor = conjg(factor)
                    w(i,:) = w(i,:) - factor*w(j,:)
                end do
            end do
        end if
        status = linalg_status_no_convergence
        if (.not. finite_complex_matrix(w)) return
        do i = 1, n
            x(perm(i),:) = w(i,:)
        end do
        status = linalg_status_success
    end subroutine ldl_solve_complex_matrix

    subroutine solve_indefinite_real_vector(a, b, x, status, lower)
        real(dp), intent(in) :: a(:, :) !! symmetric matrix; selected triangle authoritative
        real(dp), intent(in) :: b(:) !! right-hand side, length n
        real(dp), intent(out) :: x(:) !! solution, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! true for lower input triangle (default)
        real(dp), allocatable :: bb(:, :), xx(:, :)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(b)
        if (size(x) /= n) return
        allocate(bb(n,1),xx(n,1))
        bb(:,1) = b
        call solve_indefinite_real_matrix(a,bb,xx,status,lower)
        if (status == linalg_status_success) x = xx(:,1)
    end subroutine solve_indefinite_real_vector

    subroutine solve_indefinite_real_matrix(a, b, x, status, lower)
        real(dp), intent(in) :: a(:, :) !! symmetric matrix; selected triangle authoritative
        real(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solutions, shape (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! true for lower input triangle (default)
        real(dp), allocatable :: lu(:, :), d(:, :)
        integer, allocatable :: perm(:)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (size(a,2) /= n .or. size(b,1) /= n) return
        if (size(x,1) /= n .or. size(x,2) /= size(b,2)) return
        allocate(lu(n,n),d(n,n),perm(n))
        call ldl(a,lu,d,perm,status,lower)
        if (status /= linalg_status_success) return
        call ldl_solve(lu,d,perm,b,x,status,lower)
    end subroutine solve_indefinite_real_matrix

    subroutine solve_indefinite_complex_vector(a, b, x, status, lower, hermitian)
        complex(dp), intent(in) :: a(:, :) !! Hermitian or complex symmetric matrix
        complex(dp), intent(in) :: b(:) !! RHS vector, length n
        complex(dp), intent(out) :: x(:) !! solution vector; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! true for lower input triangle (default)
        logical, intent(in), optional :: hermitian !! true for Hermitian (default)
        complex(dp), allocatable :: bb(:, :), xx(:, :)
        integer :: n

        x = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(b)
        if (size(x) /= n) return
        allocate(bb(n,1),xx(n,1))
        bb(:,1) = b
        call solve_indefinite_complex_matrix(a,bb,xx,status,lower,hermitian)
        if (status == linalg_status_success) x = xx(:,1)
    end subroutine solve_indefinite_complex_vector

    subroutine solve_indefinite_complex_matrix(a, b, x, status, lower, hermitian)
        complex(dp), intent(in) :: a(:, :) !! Hermitian or complex symmetric matrix
        complex(dp), intent(in) :: b(:, :) !! RHS matrix, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solution matrix, shape (n,nrhs); zero on error
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! true for lower input triangle (default)
        logical, intent(in), optional :: hermitian !! true for Hermitian (default)
        complex(dp), allocatable :: lu(:, :), d(:, :)
        integer, allocatable :: perm(:)
        integer :: n

        x = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (size(a,2) /= n .or. size(b,1) /= n) return
        if (size(x,1) /= n .or. size(x,2) /= size(b,2)) return
        allocate(lu(n,n),d(n,n),perm(n))
        call ldl(a,lu,d,perm,status,lower,hermitian)
        if (status /= linalg_status_success) return
        call ldl_solve(lu,d,perm,b,x,status,lower,hermitian)
    end subroutine solve_indefinite_complex_matrix

    pure logical function finite_complex_matrix(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! matrix to test for finite components
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex_matrix

end module scifort_linalg_indefinite_solve
