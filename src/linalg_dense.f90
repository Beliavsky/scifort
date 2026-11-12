! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style dense triangular solves and reusable pivoted LU factorization.
! Independent implementations of standard substitution and partial-pivoting
! algorithms; no external BLAS/LAPACK dependency.
module scifort_linalg_dense
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence, &
        linalg_status_not_positive_definite
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: solve_triangular, lu_factor, lu_solve
    public :: cholesky, cho_factor, cho_solve

    interface solve_triangular
        module procedure solve_triangular_real_vector
        module procedure solve_triangular_real_matrix
        module procedure solve_triangular_complex_vector
        module procedure solve_triangular_complex_matrix
    end interface solve_triangular

    interface lu_factor
        module procedure lu_factor_real
        module procedure lu_factor_complex
    end interface lu_factor

    interface lu_solve
        module procedure lu_solve_real_vector
        module procedure lu_solve_real_matrix
        module procedure lu_solve_complex_vector
        module procedure lu_solve_complex_matrix
    end interface lu_solve

    interface cholesky
        module procedure cholesky_real
        module procedure cholesky_complex
    end interface cholesky

    interface cho_factor
        module procedure cho_factor_real
        module procedure cho_factor_complex
    end interface cho_factor

    interface cho_solve
        module procedure cho_solve_real_vector
        module procedure cho_solve_real_matrix
        module procedure cho_solve_complex_vector
        module procedure cho_solve_complex_matrix
    end interface cho_solve

contains
    subroutine solve_triangular_real_vector(a, b, x, status, lower, trans, unit_diagonal)
        real(dp), intent(in) :: a(:, :) !! coefficient matrix, square; selected triangle used
        real(dp), intent(in) :: b(:) !! right-hand side of length n
        real(dp), intent(out) :: x(:) !! length-n solution; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! use lower triangle; default false
        character(len=1), intent(in), optional :: trans !! 'N', 'T', 'C'; default 'N'
        logical, intent(in), optional :: unit_diagonal !! replace diagonal with ones; default false
        real(dp), allocatable :: rhs(:, :), result(:, :)
        integer :: n

        x = 0.0_dp
        n = size(b)
        status = linalg_status_invalid_shape
        if (n < 1 .or. size(x) /= n) return
        allocate(rhs(n,1), result(n,1))
        rhs(:,1) = b
        call solve_triangular_real_matrix(a,rhs,result,status,lower,trans,unit_diagonal)
        if (status == linalg_status_success) x = result(:,1)
    end subroutine solve_triangular_real_vector

    subroutine solve_triangular_real_matrix(a, b, x, status, lower, trans, unit_diagonal)
        real(dp), intent(in) :: a(:, :) !! n-by-n matrix; only selected triangle used
        real(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solution array, same shape as b; zero on failure
        integer, intent(out) :: status !! success, invalid_shape, singular, no_convergence
        logical, intent(in), optional :: lower !! solve using lower triangle; default false
        character(len=1), intent(in), optional :: trans !! normal, transpose, adjoint; default 'N'
        logical, intent(in), optional :: unit_diagonal !! implicit unit diagonal; default false
        real(dp) :: work, factor, diagonal
        integer :: n, nrhs, i, j, k, istart, iend, istep, kstart, kend
        logical :: use_lower, unit, forward
        character(len=1) :: op

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        nrhs = size(b,2)
        if (n < 1 .or. size(a,2) /= n) return
        if (size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        op = 'N'
        if (present(trans)) op = trans
        if (op /= 'N' .and. op /= 'T' .and. op /= 'C') return
        use_lower = .false.
        if (present(lower)) use_lower = lower
        unit = .false.
        if (present(unit_diagonal)) unit = unit_diagonal
        status = linalg_status_no_convergence
        if (.not. (all(ieee_is_finite(a)))) return
        if (.not. (all(ieee_is_finite(b)))) return
        if (.not. unit) then
            do i = 1, n
                if (abs(a(i,i)) == 0.0_dp) then
                    status = linalg_status_singular
                    return
                end if
            end do
        end if
        x = b
        forward = (op == 'N' .and. use_lower) .or. (op /= 'N' .and. .not. use_lower)
        if (forward) then
            istart = 1
            iend = n
            istep = 1
        else
            istart = n
            iend = 1
            istep = -1
        end if
        do j = 1, nrhs
            do i = istart, iend, istep
                work = x(i,j)
                if (forward) then
                    kstart = 1
                    kend = i-1
                else
                    kstart = i+1
                    kend = n
                end if
                do k = kstart, kend
                    if (op == 'N') then
                        factor = a(i,k)
                    else
                        factor = a(k,i)
                    end if
                    work = work - factor*x(k,j)
                end do
                if (.not. unit) then
                    diagonal = a(i,i)
                    work = work/diagonal
                end if
                x(i,j) = work
            end do
        end do
        if (.not. (all(ieee_is_finite(x)))) then
            x = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine solve_triangular_real_matrix

    subroutine lu_factor_real(a, lu, ipiv, status)
        real(dp), intent(in) :: a(:, :) !! n-by-n nonsingular coefficient matrix
        real(dp), intent(out) :: lu(:, :) !! packed unit-lower and upper factors, n-by-n
        integer, intent(out) :: ipiv(:) !! one-based row swaps at each elimination step, length n
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp) :: tmp
        integer :: n, k, p, i, j

        lu = 0.0_dp
        ipiv = 0
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (n < 1 .or. size(a,2) /= n .or. size(lu,1) /= n .or. &
            size(lu,2) /= n .or. size(ipiv) /= n) return
        status = linalg_status_no_convergence
        if (.not. (all(ieee_is_finite(a)))) return
        lu = a
        do k = 1, n
            p = k-1 + maxloc(abs(lu(k:n,k)),dim=1)
            if (abs(lu(p,k)) == 0.0_dp) then
                lu = 0.0_dp
                ipiv = 0
                status = linalg_status_singular
                return
            end if
            ipiv(k) = p
            if (p /= k) then
                do j = 1, n
                    tmp = lu(k,j)
                    lu(k,j) = lu(p,j)
                    lu(p,j) = tmp
                end do
            end if
            do i = k+1, n
                lu(i,k) = lu(i,k)/lu(k,k)
            end do
            do j = k+1, n
                do i = k+1, n
                    lu(i,j) = lu(i,j) - lu(i,k)*lu(k,j)
                end do
            end do
        end do
        if (.not. (all(ieee_is_finite(lu)))) then
            lu = 0.0_dp
            ipiv = 0
            return
        end if
        status = linalg_status_success
    end subroutine lu_factor_real

    subroutine lu_solve_real_vector(lu, ipiv, b, x, status, trans)
        real(dp), intent(in) :: lu(:, :) !! factors returned by lu_factor, n-by-n
        integer, intent(in) :: ipiv(:) !! one-based pivots from lu_factor
        real(dp), intent(in) :: b(:) !! vector right-hand side of length n
        real(dp), intent(out) :: x(:) !! solution, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(len=1), intent(in), optional :: trans !! 'N', 'T', 'C'; default 'N'
        real(dp), allocatable :: rhs(:, :), result(:, :)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(b)
        if (n < 1 .or. size(x) /= n) return
        allocate(rhs(n,1), result(n,1))
        rhs(:,1) = b
        call lu_solve_real_matrix(lu,ipiv,rhs,result,status,trans)
        if (status == linalg_status_success) x = result(:,1)
    end subroutine lu_solve_real_vector

    subroutine lu_solve_real_matrix(lu, ipiv, b, x, status, trans)
        real(dp), intent(in) :: lu(:, :) !! compact LU factors from lu_factor, n-by-n
        integer, intent(in) :: ipiv(:) !! one-based row swaps, length n
        real(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solutions, same shape as b; zero on failure
        integer, intent(out) :: status !! success, invalid_shape, singular, no_convergence
        character(len=1), intent(in), optional :: trans !! 'N','T','C'; default 'N'
        real(dp) :: factor, diagonal, tmp
        integer :: n, nrhs, i, j, k, p
        character(len=1) :: op

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(lu,1)
        nrhs = size(b,2)
        if (n < 1 .or. size(lu,2) /= n .or. size(ipiv) /= n) return
        if (size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        op = 'N'
        if (present(trans)) op = trans
        if (op /= 'N' .and. op /= 'T' .and. op /= 'C') return
        do i = 1, n
            if (ipiv(i) < i .or. ipiv(i) > n) return
        end do
        status = linalg_status_no_convergence
        if (.not. (all(ieee_is_finite(lu)))) return
        if (.not. (all(ieee_is_finite(b)))) return
        do i = 1, n
            if (abs(lu(i,i)) == 0.0_dp) then
                status = linalg_status_singular
                return
            end if
        end do
        x = b
        if (op == 'N') then
            ! Apply the permutation P in the same order as the factorization.
            do i = 1, n
                p = ipiv(i)
                if (p /= i) then
                    do j = 1, nrhs
                        tmp = x(i,j)
                        x(i,j) = x(p,j)
                        x(p,j) = tmp
                    end do
                end if
            end do
            ! Unit-lower L followed by upper U.
            do j = 1, nrhs
                do i = 1, n
                    do k = 1, i-1
                        x(i,j) = x(i,j)-lu(i,k)*x(k,j)
                    end do
                end do
                do i = n, 1, -1
                    do k = i+1, n
                        x(i,j) = x(i,j)-lu(i,k)*x(k,j)
                    end do
                    x(i,j) = x(i,j)/lu(i,i)
                end do
            end do
        else
            ! Transpose or conjugate transpose: U^op, L^op, then P^T.
            do j = 1, nrhs
                do i = 1, n
                    do k = 1, i-1
                        factor = lu(k,i)
                        x(i,j) = x(i,j)-factor*x(k,j)
                    end do
                    diagonal = lu(i,i)
                    x(i,j) = x(i,j)/diagonal
                end do
                do i = n, 1, -1
                    do k = i+1, n
                        factor = lu(k,i)
                        x(i,j) = x(i,j)-factor*x(k,j)
                    end do
                end do
            end do
            do i = n, 1, -1
                p = ipiv(i)
                if (p /= i) then
                    do j = 1, nrhs
                        tmp = x(i,j)
                        x(i,j) = x(p,j)
                        x(p,j) = tmp
                    end do
                end if
            end do
        end if
        if (.not. (all(ieee_is_finite(x)))) then
            x = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine lu_solve_real_matrix

    subroutine solve_triangular_complex_vector(a, b, x, status, lower, trans, unit_diagonal)
        complex(dp), intent(in) :: a(:, :) !! coefficient matrix, square; selected triangle used
        complex(dp), intent(in) :: b(:) !! right-hand side of length n
        complex(dp), intent(out) :: x(:) !! length-n solution; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! use lower triangle; default false
        character(len=1), intent(in), optional :: trans !! 'N', 'T', 'C'; default 'N'
        logical, intent(in), optional :: unit_diagonal !! replace diagonal with ones; default false
        complex(dp), allocatable :: rhs(:, :), result(:, :)
        integer :: n

        x = 0.0_dp
        n = size(b)
        status = linalg_status_invalid_shape
        if (n < 1 .or. size(x) /= n) return
        allocate(rhs(n,1), result(n,1))
        rhs(:,1) = b
        call solve_triangular_complex_matrix(a,rhs,result,status,lower,trans,unit_diagonal)
        if (status == linalg_status_success) x = result(:,1)
    end subroutine solve_triangular_complex_vector

    subroutine solve_triangular_complex_matrix(a, b, x, status, lower, trans, unit_diagonal)
        complex(dp), intent(in) :: a(:, :) !! n-by-n matrix; only selected triangle used
        complex(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solution array, same shape as b; zero on failure
        integer, intent(out) :: status !! success, invalid_shape, singular, no_convergence
        logical, intent(in), optional :: lower !! solve using lower triangle; default false
        character(len=1), intent(in), optional :: trans !! normal, transpose, adjoint; default 'N'
        logical, intent(in), optional :: unit_diagonal !! implicit unit diagonal; default false
        complex(dp) :: work, factor, diagonal
        integer :: n, nrhs, i, j, k, istart, iend, istep, kstart, kend
        logical :: use_lower, unit, forward
        character(len=1) :: op

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        nrhs = size(b,2)
        if (n < 1 .or. size(a,2) /= n) return
        if (size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        op = 'N'
        if (present(trans)) op = trans
        if (op /= 'N' .and. op /= 'T' .and. op /= 'C') return
        use_lower = .false.
        if (present(lower)) use_lower = lower
        unit = .false.
        if (present(unit_diagonal)) unit = unit_diagonal
        status = linalg_status_no_convergence
        if (.not. (all(ieee_is_finite(real(a,dp))) .and. &
            all(ieee_is_finite(aimag(a))))) return
        if (.not. (all(ieee_is_finite(real(b,dp))) .and. &
            all(ieee_is_finite(aimag(b))))) return
        if (.not. unit) then
            do i = 1, n
                if (abs(a(i,i)) == 0.0_dp) then
                    status = linalg_status_singular
                    return
                end if
            end do
        end if
        x = b
        forward = (op == 'N' .and. use_lower) .or. (op /= 'N' .and. .not. use_lower)
        if (forward) then
            istart = 1
            iend = n
            istep = 1
        else
            istart = n
            iend = 1
            istep = -1
        end if
        do j = 1, nrhs
            do i = istart, iend, istep
                work = x(i,j)
                if (forward) then
                    kstart = 1
                    kend = i-1
                else
                    kstart = i+1
                    kend = n
                end if
                do k = kstart, kend
                    if (op == 'N') then
                        factor = a(i,k)
                    else
                        factor = a(k,i)
                        if (op == 'C') factor = conjg(factor)
                    end if
                    work = work - factor*x(k,j)
                end do
                if (.not. unit) then
                    diagonal = a(i,i)
                        if (op == 'C') diagonal = conjg(diagonal)
                    work = work/diagonal
                end if
                x(i,j) = work
            end do
        end do
        if (.not. (all(ieee_is_finite(real(x,dp))) .and. &
            all(ieee_is_finite(aimag(x))))) then
            x = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine solve_triangular_complex_matrix

    subroutine lu_factor_complex(a, lu, ipiv, status)
        complex(dp), intent(in) :: a(:, :) !! n-by-n nonsingular coefficient matrix
        complex(dp), intent(out) :: lu(:, :) !! packed unit-lower and upper factors, n-by-n
        integer, intent(out) :: ipiv(:) !! one-based row swaps at each elimination step, length n
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp) :: tmp
        integer :: n, k, p, i, j

        lu = 0.0_dp
        ipiv = 0
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (n < 1 .or. size(a,2) /= n .or. size(lu,1) /= n .or. &
            size(lu,2) /= n .or. size(ipiv) /= n) return
        status = linalg_status_no_convergence
        if (.not. (all(ieee_is_finite(real(a,dp))) .and. &
            all(ieee_is_finite(aimag(a))))) return
        lu = a
        do k = 1, n
            p = k-1 + maxloc(abs(lu(k:n,k)),dim=1)
            if (abs(lu(p,k)) == 0.0_dp) then
                lu = 0.0_dp
                ipiv = 0
                status = linalg_status_singular
                return
            end if
            ipiv(k) = p
            if (p /= k) then
                do j = 1, n
                    tmp = lu(k,j)
                    lu(k,j) = lu(p,j)
                    lu(p,j) = tmp
                end do
            end if
            do i = k+1, n
                lu(i,k) = lu(i,k)/lu(k,k)
            end do
            do j = k+1, n
                do i = k+1, n
                    lu(i,j) = lu(i,j) - lu(i,k)*lu(k,j)
                end do
            end do
        end do
        if (.not. (all(ieee_is_finite(real(lu,dp))) .and. &
            all(ieee_is_finite(aimag(lu))))) then
            lu = 0.0_dp
            ipiv = 0
            return
        end if
        status = linalg_status_success
    end subroutine lu_factor_complex

    subroutine lu_solve_complex_vector(lu, ipiv, b, x, status, trans)
        complex(dp), intent(in) :: lu(:, :) !! factors returned by lu_factor, n-by-n
        integer, intent(in) :: ipiv(:) !! one-based pivots from lu_factor
        complex(dp), intent(in) :: b(:) !! vector right-hand side of length n
        complex(dp), intent(out) :: x(:) !! solution, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(len=1), intent(in), optional :: trans !! 'N', 'T', 'C'; default 'N'
        complex(dp), allocatable :: rhs(:, :), result(:, :)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(b)
        if (n < 1 .or. size(x) /= n) return
        allocate(rhs(n,1), result(n,1))
        rhs(:,1) = b
        call lu_solve_complex_matrix(lu,ipiv,rhs,result,status,trans)
        if (status == linalg_status_success) x = result(:,1)
    end subroutine lu_solve_complex_vector

    subroutine lu_solve_complex_matrix(lu, ipiv, b, x, status, trans)
        complex(dp), intent(in) :: lu(:, :) !! compact LU factors from lu_factor, n-by-n
        integer, intent(in) :: ipiv(:) !! one-based row swaps, length n
        complex(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solutions, same shape as b; zero on failure
        integer, intent(out) :: status !! success, invalid_shape, singular, no_convergence
        character(len=1), intent(in), optional :: trans !! 'N','T','C'; default 'N'
        complex(dp) :: factor, diagonal, tmp
        integer :: n, nrhs, i, j, k, p
        character(len=1) :: op

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(lu,1)
        nrhs = size(b,2)
        if (n < 1 .or. size(lu,2) /= n .or. size(ipiv) /= n) return
        if (size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        op = 'N'
        if (present(trans)) op = trans
        if (op /= 'N' .and. op /= 'T' .and. op /= 'C') return
        do i = 1, n
            if (ipiv(i) < i .or. ipiv(i) > n) return
        end do
        status = linalg_status_no_convergence
        if (.not. (all(ieee_is_finite(real(lu,dp))) .and. &
            all(ieee_is_finite(aimag(lu))))) return
        if (.not. (all(ieee_is_finite(real(b,dp))) .and. &
            all(ieee_is_finite(aimag(b))))) return
        do i = 1, n
            if (abs(lu(i,i)) == 0.0_dp) then
                status = linalg_status_singular
                return
            end if
        end do
        x = b
        if (op == 'N') then
            ! Apply the permutation P in the same order as the factorization.
            do i = 1, n
                p = ipiv(i)
                if (p /= i) then
                    do j = 1, nrhs
                        tmp = x(i,j)
                        x(i,j) = x(p,j)
                        x(p,j) = tmp
                    end do
                end if
            end do
            ! Unit-lower L followed by upper U.
            do j = 1, nrhs
                do i = 1, n
                    do k = 1, i-1
                        x(i,j) = x(i,j)-lu(i,k)*x(k,j)
                    end do
                end do
                do i = n, 1, -1
                    do k = i+1, n
                        x(i,j) = x(i,j)-lu(i,k)*x(k,j)
                    end do
                    x(i,j) = x(i,j)/lu(i,i)
                end do
            end do
        else
            ! Transpose or conjugate transpose: U^op, L^op, then P^T.
            do j = 1, nrhs
                do i = 1, n
                    do k = 1, i-1
                        factor = lu(k,i)
                        if (op == 'C') factor = conjg(factor)
                        x(i,j) = x(i,j)-factor*x(k,j)
                    end do
                    diagonal = lu(i,i)
                        if (op == 'C') diagonal = conjg(diagonal)
                    x(i,j) = x(i,j)/diagonal
                end do
                do i = n, 1, -1
                    do k = i+1, n
                        factor = lu(k,i)
                        if (op == 'C') factor = conjg(factor)
                        x(i,j) = x(i,j)-factor*x(k,j)
                    end do
                end do
            end do
            do i = n, 1, -1
                p = ipiv(i)
                if (p /= i) then
                    do j = 1, nrhs
                        tmp = x(i,j)
                        x(i,j) = x(p,j)
                        x(p,j) = tmp
                    end do
                end if
            end do
        end if
        if (.not. (all(ieee_is_finite(real(x,dp))) .and. &
            all(ieee_is_finite(aimag(x))))) then
            x = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine lu_solve_complex_matrix


    ! Dense Cholesky factorization in the requested triangular storage.
    ! Only the selected triangle contributes to the factorization.  The
    ! opposite triangle of the output is deliberately initialized to zero.

    subroutine cho_factor_real(a, f, status, lower)
        real(dp), intent(in) :: a(:, :) !! square positive-definite matrix; selected triangle only
        real(dp), intent(out) :: f(:, :) !! triangular factor with unused triangle set to zero
        integer, intent(out) :: status !! success, invalid_shape, not_positive_definite, no_convergence
        logical, intent(in), optional :: lower !! lower L if true, upper U if false (default)
        integer :: n, i, j
        real(dp) :: diagonal
        real(dp) :: work
        logical :: use_lower

        f = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (n < 1 .or. size(a, 2) /= n) return
        if (size(f, 1) /= n .or. size(f, 2) /= n) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        use_lower = .false.
        if (present(lower)) use_lower = lower
        do j = 1, n
            diagonal = a(j, j)
            if (j > 1) then
                if (use_lower) then
                    diagonal = diagonal - sum(f(j, 1:j-1)**2)
                else
                    diagonal = diagonal - sum(f(1:j-1, j)**2)
                end if
            end if
            if (.not. ieee_is_finite(diagonal)) then
                f = 0.0_dp
                return
            end if
            if (diagonal <= 0.0_dp) then
                f = 0.0_dp
                status = linalg_status_not_positive_definite
                return
            end if
            f(j, j) = sqrt(diagonal)
            do i = j+1, n
                if (use_lower) then
                    work = a(i, j)
                    if (j > 1) work = work - sum(f(i, 1:j-1)*f(j, 1:j-1))
                    work = work / f(j, j)
                    f(i, j) = work
                else
                    work = a(j, i)
                    if (j > 1) work = work - sum(f(1:j-1, j)*f(1:j-1, i))
                    work = work / f(j, j)
                    f(j, i) = work
                end if
                if (.not. ieee_is_finite(work)) then
                    f = 0.0_dp
                    return
                end if
            end do
        end do
        status = linalg_status_success
    end subroutine cho_factor_real

    subroutine cholesky_real(a, f, status, lower)
        real(dp), intent(in) :: a(:, :) !! square positive-definite coefficient matrix
        real(dp), intent(out) :: f(:, :) !! clean lower or upper triangular Cholesky factor
        integer, intent(out) :: status !! success or a linalg_status_* failure
        logical, intent(in), optional :: lower !! choose lower L, default upper U

        call cho_factor_real(a, f, status, lower)
    end subroutine cholesky_real

    subroutine cho_solve_real_vector(f, b, x, status, lower)
        real(dp), intent(in) :: f(:, :) !! triangular Cholesky factor of the original matrix
        real(dp), intent(in) :: b(:) !! right-hand side, length n
        real(dp), intent(out) :: x(:) !! solution, zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in), optional :: lower !! true for L, default upper U
        real(dp), allocatable :: rhs(:, :), solution(:, :)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(b)
        if (n < 1 .or. size(x) /= n) return
        allocate(rhs(n, 1), solution(n, 1))
        rhs(:, 1) = b
        call cho_solve_real_matrix(f, rhs, solution, status, lower)
        if (status == linalg_status_success) x = solution(:, 1)
    end subroutine cho_solve_real_vector

    subroutine cho_solve_real_matrix(f, b, x, status, lower)
        real(dp), intent(in) :: f(:, :) !! lower L or upper U from cho_factor
        real(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solution, shape (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in), optional :: lower !! true for L, default upper U
        real(dp), allocatable :: work(:, :)
        integer :: n
        logical :: use_lower

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(f, 1)
        if (n < 1 .or. size(f, 2) /= n) return
        if (size(b, 1) /= n .or. size(x, 1) /= n .or. size(x, 2) /= size(b, 2)) return
        use_lower = .false.
        if (present(lower)) use_lower = lower
        allocate(work(n, size(b, 2)))
        if (use_lower) then
            call solve_triangular(f, b, work, status, lower=.true.)
            if (status /= linalg_status_success) return
            call solve_triangular(f, work, x, status, lower=.true., trans='C')
        else
            call solve_triangular(f, b, work, status, trans='C')
            if (status /= linalg_status_success) return
            call solve_triangular(f, work, x, status)
        end if
        if (status /= linalg_status_success) x = 0.0_dp
    end subroutine cho_solve_real_matrix

    subroutine cho_factor_complex(a, f, status, lower)
        complex(dp), intent(in) :: a(:, :) !! square positive-definite matrix; selected triangle only
        complex(dp), intent(out) :: f(:, :) !! triangular factor with unused triangle set to zero
        integer, intent(out) :: status !! success, invalid_shape, not_positive_definite, no_convergence
        logical, intent(in), optional :: lower !! lower L if true, upper U if false (default)
        integer :: n, i, j
        real(dp) :: diagonal
        complex(dp) :: work
        logical :: use_lower

        f = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (n < 1 .or. size(a, 2) /= n) return
        if (size(f, 1) /= n .or. size(f, 2) /= n) return
        status = linalg_status_no_convergence
        if (.not. (all(ieee_is_finite(real(a, dp))) .and. &
            all(ieee_is_finite(aimag(a))))) return
        use_lower = .false.
        if (present(lower)) use_lower = lower
        do j = 1, n
            diagonal = real(a(j, j), dp)
            if (j > 1) then
                if (use_lower) then
                    diagonal = diagonal - sum(abs(f(j, 1:j-1))**2)
                else
                    diagonal = diagonal - sum(abs(f(1:j-1, j))**2)
                end if
            end if
            if (.not. ieee_is_finite(diagonal)) then
                f = cmplx(0.0_dp, 0.0_dp, dp)
                return
            end if
            if (diagonal <= 0.0_dp) then
                f = cmplx(0.0_dp, 0.0_dp, dp)
                status = linalg_status_not_positive_definite
                return
            end if
            f(j, j) = sqrt(diagonal)
            do i = j+1, n
                if (use_lower) then
                    work = a(i, j)
                    if (j > 1) work = work - sum(f(i, 1:j-1)*conjg(f(j, 1:j-1)))
                    work = work / f(j, j)
                    f(i, j) = work
                else
                    work = a(j, i)
                    if (j > 1) work = work - sum(conjg(f(1:j-1, j))*f(1:j-1, i))
                    work = work / f(j, j)
                    f(j, i) = work
                end if
                if (.not. (ieee_is_finite(real(work, dp)) .and. ieee_is_finite(aimag(work)))) then
                    f = cmplx(0.0_dp, 0.0_dp, dp)
                    return
                end if
            end do
        end do
        status = linalg_status_success
    end subroutine cho_factor_complex

    subroutine cholesky_complex(a, f, status, lower)
        complex(dp), intent(in) :: a(:, :) !! square positive-definite coefficient matrix
        complex(dp), intent(out) :: f(:, :) !! clean lower or upper triangular Cholesky factor
        integer, intent(out) :: status !! success or a linalg_status_* failure
        logical, intent(in), optional :: lower !! choose lower L, default upper U

        call cho_factor_complex(a, f, status, lower)
    end subroutine cholesky_complex

    subroutine cho_solve_complex_vector(f, b, x, status, lower)
        complex(dp), intent(in) :: f(:, :) !! triangular Cholesky factor of the original matrix
        complex(dp), intent(in) :: b(:) !! right-hand side, length n
        complex(dp), intent(out) :: x(:) !! solution, zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in), optional :: lower !! true for L, default upper U
        complex(dp), allocatable :: rhs(:, :), solution(:, :)
        integer :: n

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(b)
        if (n < 1 .or. size(x) /= n) return
        allocate(rhs(n, 1), solution(n, 1))
        rhs(:, 1) = b
        call cho_solve_complex_matrix(f, rhs, solution, status, lower)
        if (status == linalg_status_success) x = solution(:, 1)
    end subroutine cho_solve_complex_vector

    subroutine cho_solve_complex_matrix(f, b, x, status, lower)
        complex(dp), intent(in) :: f(:, :) !! lower L or upper U from cho_factor
        complex(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solution, shape (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in), optional :: lower !! true for L, default upper U
        complex(dp), allocatable :: work(:, :)
        integer :: n
        logical :: use_lower

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(f, 1)
        if (n < 1 .or. size(f, 2) /= n) return
        if (size(b, 1) /= n .or. size(x, 1) /= n .or. size(x, 2) /= size(b, 2)) return
        use_lower = .false.
        if (present(lower)) use_lower = lower
        allocate(work(n, size(b, 2)))
        if (use_lower) then
            call solve_triangular(f, b, work, status, lower=.true.)
            if (status /= linalg_status_success) return
            call solve_triangular(f, work, x, status, lower=.true., trans='C')
        else
            call solve_triangular(f, b, work, status, trans='C')
            if (status /= linalg_status_success) return
            call solve_triangular(f, work, x, status)
        end if
        if (status /= linalg_status_success) x = cmplx(0.0_dp, 0.0_dp, dp)
    end subroutine cho_solve_complex_matrix

end module scifort_linalg_dense
