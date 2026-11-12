! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-compatible packed-band direct solvers. This independent Fortran
! implementation stores only the active band, including partial-pivot fill.
! No LAPACK or external package is needed.
module scifort_linalg_banded
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_not_positive_definite, linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: solve_banded, solveh_banded, cholesky_banded, cho_solve_banded
    public :: lu_factor_banded, lu_solve_banded

    interface solve_banded
        module procedure solve_banded_real_vector
        module procedure solve_banded_real_matrix
        module procedure solve_banded_complex_vector
        module procedure solve_banded_complex_matrix
    end interface solve_banded
    interface lu_factor_banded
        module procedure lu_factor_banded_real
        module procedure lu_factor_banded_complex
    end interface lu_factor_banded
    interface lu_solve_banded
        module procedure lu_solve_banded_real_vector
        module procedure lu_solve_banded_real_matrix
        module procedure lu_solve_banded_complex_vector
        module procedure lu_solve_banded_complex_matrix
    end interface lu_solve_banded
    interface solveh_banded
        module procedure solveh_banded_real_vector
        module procedure solveh_banded_real_matrix
        module procedure solveh_banded_complex_vector
        module procedure solveh_banded_complex_matrix
    end interface solveh_banded
    interface cholesky_banded
        module procedure cholesky_banded_real
        module procedure cholesky_banded_complex
    end interface cholesky_banded
    interface cho_solve_banded
        module procedure cho_solve_banded_real_vector
        module procedure cho_solve_banded_real_matrix
        module procedure cho_solve_banded_complex_vector
        module procedure cho_solve_banded_complex_matrix
    end interface cho_solve_banded

contains
    subroutine solve_banded_real_vector(ab, lower, upper, b, x, status)
        real(dp), intent(in) :: ab(:, :) !! packed diagonals, shape (lower+upper+1,n)
        integer, intent(in) :: lower !! number of lower diagonals, 0..n-1
        integer, intent(in) :: upper !! number of upper diagonals, 0..n-1
        real(dp), intent(in) :: b(:) !! right-hand side, length n
        real(dp), intent(out) :: x(:) !! solution, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), allocatable :: rhs(:, :), result(:, :)
        integer :: n

        x = 0.0_dp
        n = size(b)
        status = linalg_status_invalid_shape
        if (size(x) /= n .or. n < 1) return
        allocate(rhs(n, 1), result(n, 1))
        rhs(:, 1) = b
        call solve_banded_real_matrix(ab, lower, upper, rhs, result, status)
        if (status == linalg_status_success) x = result(:, 1)
    end subroutine solve_banded_real_vector

    subroutine solve_banded_real_matrix(ab, lower, upper, b, x, status)
        real(dp), intent(in) :: ab(:, :) !! compact band (lower+upper+1,n)
        integer, intent(in) :: lower !! number of subdiagonals
        integer, intent(in) :: upper !! number of superdiagonals
        real(dp), intent(in) :: b(:, :) !! right-hand sides (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solution (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), allocatable :: lu(:, :)
        integer, allocatable :: ipiv(:)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(ab, 2)
        if (n < 1 .or. lower < 0 .or. upper < 0) return
        if (lower >= n .or. upper >= n) return
        if (size(ab, 1) /= lower+upper+1) return
        if (size(b, 1) /= n .or. size(x, 1) /= n .or. size(x, 2) /= size(b, 2)) return
        allocate(lu(upper+2*lower+1,n), ipiv(n))
        call lu_factor_banded_real(ab, lower, upper, lu, ipiv, status)
        if (status /= linalg_status_success) return
        call lu_solve_banded_real_matrix(lu, lower, upper, ipiv, b, x, status)
    end subroutine solve_banded_real_matrix
    subroutine cholesky_banded_real(ab, cb, status, lower)
        real(dp), intent(in) :: ab(:, :) !! packed Hermitian/symmetric positive definite A
        real(dp), intent(out) :: cb(:, :) !! packed L if lower else packed U, A=LL^H=U^HU
        integer, intent(out) :: status !! success/invalid_shape/not_positive_definite
        logical, intent(in), optional :: lower !! .true. for lower storage; default .false.
        real(dp), allocatable :: factor(:, :)
        real(dp) :: entry, work
        real(dp) :: diagonal
        integer :: n, bw, i, j, k, kfirst, diag_index
        logical :: use_lower

        cb = 0.0_dp
        status = linalg_status_invalid_shape
        bw = size(ab, 1)-1
        n = size(ab, 2)
        if (n < 1 .or. bw < 0 .or. bw >= n) return
        if (any(shape(cb) /= shape(ab))) return
        status = linalg_status_not_positive_definite
        if (.not. finite_real(ab)) return
        use_lower = .false.
        if (present(lower)) use_lower = lower
        diag_index = bw+1
        if (use_lower) diag_index = 1
        allocate(factor(bw+1,n))
        factor = 0.0_dp
        do j = 1, n
            diagonal = ab(diag_index,j)
            do k = max(1,j-bw), j-1
                diagonal = diagonal - factor(1+j-k,k)**2
            end do
            if (.not. ieee_is_finite(diagonal) .or. diagonal <= 0.0_dp) return
            factor(1,j) = sqrt(diagonal)
            do i = j+1, min(n,j+bw)
                if (use_lower) then
                    entry = ab(1+i-j,j)
                else
                    entry = ab(bw+1+j-i, i)
                end if
                work = entry
                kfirst = max(1, j-bw, i-bw)
                do k = kfirst,j-1
                    work = work - factor(1+i-k,k)*factor(1+j-k,k)
                end do
                factor(1+i-j,j) = work/factor(1,j)
            end do
        end do
        if (.not. finite_real(factor)) return
        if (use_lower) then
            cb = factor
        else
            do j = 1, n
                do i = max(1, j-bw), j
                    cb(bw+1+i-j, j) = factor(1+j-i,i)
                end do
            end do
        end if
        status = linalg_status_success
    end subroutine cholesky_banded_real
    subroutine cho_solve_banded_real_vector(cb, b, x, status, lower)
        real(dp), intent(in) :: cb(:, :) !! packed Cholesky factor from cholesky_banded
        real(dp), intent(in) :: b(:) !! right-hand side n
        real(dp), intent(out) :: x(:) !! solution n
        integer, intent(out) :: status !! status code
        logical, intent(in), optional :: lower !! match factor storage; default .false.
        real(dp), allocatable :: rhs(:, :), result(:, :)
        integer :: n

        x = 0.0_dp
        n = size(b)
        status = linalg_status_invalid_shape
        if (n < 1 .or. size(x) /= n) return
        allocate(rhs(n,1),result(n,1))
        rhs(:,1) = b
        call cho_solve_banded_real_matrix(cb,rhs,result,status,lower)
        if (status == linalg_status_success) x = result(:,1)
    end subroutine cho_solve_banded_real_vector

    subroutine cho_solve_banded_real_matrix(cb, b, x, status, lower)
        real(dp), intent(in) :: cb(:, :) !! packed triangular Cholesky factor, shape (bw+1,n)
        real(dp), intent(in) :: b(:, :) !! right-hand sides (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solutions (n,nrhs); zero on failure
        integer, intent(out) :: status !! status code
        logical, intent(in), optional :: lower !! .true. for packed L; default packed U
        real(dp), allocatable :: factor(:, :)
        integer :: n, bw, i, j
        logical :: use_lower

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(cb,2)
        bw = size(cb,1)-1
        if (n < 1 .or. bw < 0 .or. bw >= n) return
        if (size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= size(b,2)) return
        status = linalg_status_no_convergence
        if (.not. finite_real(cb) .or. .not. finite_real(b)) return
        status = linalg_status_singular
        use_lower = .false.
        if (present(lower)) use_lower = lower
        allocate(factor(bw+1,n))
        factor = 0.0_dp
        if (use_lower) then
            factor = cb
        else
            do j = 1, n
                do i = j, min(n,j+bw)
                    factor(1+i-j,j) = cb(bw+1+j-i,i)
                end do
            end do
        end if
        if (any(abs(factor(1,:)) <= 0.0_dp)) return
        x = b
        ! Solve L*y = b, then L^H*x = y, reusing banded Cholesky.
        do j = 1, n
            x(j,:) = x(j,:)/factor(1,j)
            do i = j+1, min(n,j+bw)
                x(i,:) = x(i,:) - factor(1+i-j,j)*x(j,:)
            end do
        end do
        do j = n, 1, -1
            x(j,:) = x(j,:)/factor(1,j)
            do i = max(1,j-bw), j-1
                x(i,:) = x(i,:) - factor(1+j-i,i)*x(j,:)
            end do
        end do
        status = linalg_status_no_convergence
        if (.not. finite_real(x)) then
            x = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine cho_solve_banded_real_matrix
    subroutine solveh_banded_real_vector(ab, b, x, status, lower)
        real(dp), intent(in) :: ab(:, :) !! packed symmetric/Hermitian positive definite A
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(out) :: x(:) !! solution
        integer, intent(out) :: status !! success or diagnostic
        logical, intent(in), optional :: lower !! packed triangle used, default upper
        real(dp), allocatable :: cb(:, :)

        x = 0.0_dp
        allocate(cb(size(ab,1),size(ab,2)))
        call cholesky_banded_real(ab,cb,status,lower)
        if (status /= linalg_status_success) return
        call cho_solve_banded_real_vector(cb,b,x,status,lower)
    end subroutine solveh_banded_real_vector

    subroutine solveh_banded_real_matrix(ab, b, x, status, lower)
        real(dp), intent(in) :: ab(:, :) !! packed symmetric/Hermitian positive definite A
        real(dp), intent(in) :: b(:, :) !! right-hand sides (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solutions (n,nrhs)
        integer, intent(out) :: status !! success or diagnostic
        logical, intent(in), optional :: lower !! packed triangle used, default upper
        real(dp), allocatable :: cb(:, :)

        x = 0.0_dp
        allocate(cb(size(ab,1),size(ab,2)))
        call cholesky_banded_real(ab,cb,status,lower)
        if (status /= linalg_status_success) return
        call cho_solve_banded_real_matrix(cb,b,x,status,lower)
    end subroutine solveh_banded_real_matrix
    subroutine solve_banded_complex_vector(ab, lower, upper, b, x, status)
        complex(dp), intent(in) :: ab(:, :) !! packed diagonals, shape (lower+upper+1,n)
        integer, intent(in) :: lower !! number of lower diagonals, 0..n-1
        integer, intent(in) :: upper !! number of upper diagonals, 0..n-1
        complex(dp), intent(in) :: b(:) !! right-hand side, length n
        complex(dp), intent(out) :: x(:) !! solution, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        complex(dp), allocatable :: rhs(:, :), result(:, :)
        integer :: n

        x = cmplx(0.0_dp, 0.0_dp, dp)
        n = size(b)
        status = linalg_status_invalid_shape
        if (size(x) /= n .or. n < 1) return
        allocate(rhs(n, 1), result(n, 1))
        rhs(:, 1) = b
        call solve_banded_complex_matrix(ab, lower, upper, rhs, result, status)
        if (status == linalg_status_success) x = result(:, 1)
    end subroutine solve_banded_complex_vector

    subroutine solve_banded_complex_matrix(ab, lower, upper, b, x, status)
        complex(dp), intent(in) :: ab(:, :) !! compact band (lower+upper+1,n)
        integer, intent(in) :: lower !! number of subdiagonals
        integer, intent(in) :: upper !! number of superdiagonals
        complex(dp), intent(in) :: b(:, :) !! right-hand sides (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solution (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        complex(dp), allocatable :: lu(:, :)
        integer, allocatable :: ipiv(:)
        integer :: n

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(ab, 2)
        if (n < 1 .or. lower < 0 .or. upper < 0) return
        if (lower >= n .or. upper >= n) return
        if (size(ab, 1) /= lower+upper+1) return
        if (size(b, 1) /= n .or. size(x, 1) /= n .or. size(x, 2) /= size(b, 2)) return
        allocate(lu(upper+2*lower+1,n), ipiv(n))
        call lu_factor_banded_complex(ab, lower, upper, lu, ipiv, status)
        if (status /= linalg_status_success) return
        call lu_solve_banded_complex_matrix(lu, lower, upper, ipiv, b, x, status)
    end subroutine solve_banded_complex_matrix
    subroutine cholesky_banded_complex(ab, cb, status, lower)
        complex(dp), intent(in) :: ab(:, :) !! packed Hermitian/symmetric positive definite A
        complex(dp), intent(out) :: cb(:, :) !! packed L if lower else packed U, A=LL^H=U^HU
        integer, intent(out) :: status !! success/invalid_shape/not_positive_definite
        logical, intent(in), optional :: lower !! .true. for lower storage; default .false.
        complex(dp), allocatable :: factor(:, :)
        complex(dp) :: entry, work
        real(dp) :: diagonal
        integer :: n, bw, i, j, k, kfirst, diag_index
        logical :: use_lower

        cb = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        bw = size(ab, 1)-1
        n = size(ab, 2)
        if (n < 1 .or. bw < 0 .or. bw >= n) return
        if (any(shape(cb) /= shape(ab))) return
        status = linalg_status_not_positive_definite
        if (.not. finite_complex(ab)) return
        use_lower = .false.
        if (present(lower)) use_lower = lower
        diag_index = bw+1
        if (use_lower) diag_index = 1
        if (any(abs(aimag(ab(diag_index, :))) > &
                128.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(real(ab(diag_index, :), dp))))) return
        allocate(factor(bw+1,n))
        factor = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, n
            diagonal = real(ab(diag_index,j), dp)
            do k = max(1,j-bw), j-1
                diagonal = diagonal - abs(factor(1+j-k,k))**2
            end do
            if (.not. ieee_is_finite(diagonal) .or. diagonal <= 0.0_dp) return
            factor(1,j) = sqrt(diagonal)
            do i = j+1, min(n,j+bw)
                if (use_lower) then
                    entry = ab(1+i-j,j)
                else
                    entry = conjg(ab(bw+1+j-i, i))
                end if
                work = entry
                kfirst = max(1, j-bw, i-bw)
                do k = kfirst,j-1
                    work = work - factor(1+i-k,k)*conjg(factor(1+j-k,k))
                end do
                factor(1+i-j,j) = work/factor(1,j)
            end do
        end do
        if (.not. finite_complex(factor)) return
        if (use_lower) then
            cb = factor
        else
            do j = 1, n
                do i = max(1, j-bw), j
                    cb(bw+1+i-j, j) = conjg(factor(1+j-i,i))
                end do
            end do
        end if
        status = linalg_status_success
    end subroutine cholesky_banded_complex
    subroutine cho_solve_banded_complex_vector(cb, b, x, status, lower)
        complex(dp), intent(in) :: cb(:, :) !! packed Cholesky factor from cholesky_banded
        complex(dp), intent(in) :: b(:) !! right-hand side n
        complex(dp), intent(out) :: x(:) !! solution n
        integer, intent(out) :: status !! status code
        logical, intent(in), optional :: lower !! match factor storage; default .false.
        complex(dp), allocatable :: rhs(:, :), result(:, :)
        integer :: n

        x = cmplx(0.0_dp, 0.0_dp, dp)
        n = size(b)
        status = linalg_status_invalid_shape
        if (n < 1 .or. size(x) /= n) return
        allocate(rhs(n,1),result(n,1))
        rhs(:,1) = b
        call cho_solve_banded_complex_matrix(cb,rhs,result,status,lower)
        if (status == linalg_status_success) x = result(:,1)
    end subroutine cho_solve_banded_complex_vector

    subroutine cho_solve_banded_complex_matrix(cb, b, x, status, lower)
        complex(dp), intent(in) :: cb(:, :) !! packed triangular Cholesky factor, shape (bw+1,n)
        complex(dp), intent(in) :: b(:, :) !! right-hand sides (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solutions (n,nrhs); zero on failure
        integer, intent(out) :: status !! status code
        logical, intent(in), optional :: lower !! .true. for packed L; default packed U
        complex(dp), allocatable :: factor(:, :)
        integer :: n, bw, i, j
        logical :: use_lower

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(cb,2)
        bw = size(cb,1)-1
        if (n < 1 .or. bw < 0 .or. bw >= n) return
        if (size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= size(b,2)) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(cb) .or. .not. finite_complex(b)) return
        status = linalg_status_singular
        use_lower = .false.
        if (present(lower)) use_lower = lower
        allocate(factor(bw+1,n))
        factor = cmplx(0.0_dp, 0.0_dp, dp)
        if (use_lower) then
            factor = cb
        else
            do j = 1, n
                do i = j, min(n,j+bw)
                    factor(1+i-j,j) = conjg(cb(bw+1+j-i,i))
                end do
            end do
        end if
        if (any(abs(factor(1,:)) <= 0.0_dp)) return
        x = b
        ! Solve L*y = b, then L^H*x = y, reusing banded Cholesky.
        do j = 1, n
            x(j,:) = x(j,:)/factor(1,j)
            do i = j+1, min(n,j+bw)
                x(i,:) = x(i,:) - factor(1+i-j,j)*x(j,:)
            end do
        end do
        do j = n, 1, -1
            x(j,:) = x(j,:)/conjg(factor(1,j))
            do i = max(1,j-bw), j-1
                x(i,:) = x(i,:) - conjg(factor(1+j-i,i))*x(j,:)
            end do
        end do
        status = linalg_status_no_convergence
        if (.not. finite_complex(x)) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        status = linalg_status_success
    end subroutine cho_solve_banded_complex_matrix
    subroutine solveh_banded_complex_vector(ab, b, x, status, lower)
        complex(dp), intent(in) :: ab(:, :) !! packed symmetric/Hermitian positive definite A
        complex(dp), intent(in) :: b(:) !! right-hand side
        complex(dp), intent(out) :: x(:) !! solution
        integer, intent(out) :: status !! success or diagnostic
        logical, intent(in), optional :: lower !! packed triangle used, default upper
        complex(dp), allocatable :: cb(:, :)

        x = cmplx(0.0_dp, 0.0_dp, dp)
        allocate(cb(size(ab,1),size(ab,2)))
        call cholesky_banded_complex(ab,cb,status,lower)
        if (status /= linalg_status_success) return
        call cho_solve_banded_complex_vector(cb,b,x,status,lower)
    end subroutine solveh_banded_complex_vector

    subroutine solveh_banded_complex_matrix(ab, b, x, status, lower)
        complex(dp), intent(in) :: ab(:, :) !! packed symmetric/Hermitian positive definite A
        complex(dp), intent(in) :: b(:, :) !! right-hand sides (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solutions (n,nrhs)
        integer, intent(out) :: status !! success or diagnostic
        logical, intent(in), optional :: lower !! packed triangle used, default upper
        complex(dp), allocatable :: cb(:, :)

        x = cmplx(0.0_dp, 0.0_dp, dp)
        allocate(cb(size(ab,1),size(ab,2)))
        call cholesky_banded_complex(ab,cb,status,lower)
        if (status /= linalg_status_success) return
        call cho_solve_banded_complex_matrix(cb,b,x,status,lower)
    end subroutine solveh_banded_complex_matrix
    ! This is a compact, reusable analogue of LAPACK's GBTRF storage.
    ! The second padding band accounts for pivot fill; row indices inside
    ! LU are offset by upper+lower+1 relative to the true diagonal.
    subroutine lu_factor_banded_real(ab, lower, upper, lu, ipiv, status)
        real(dp), intent(in) :: ab(:, :) !! compact input band (lower+upper+1,n)
        integer, intent(in) :: lower !! number of subdiagonals, 0..n-1
        integer, intent(in) :: upper !! number of superdiagonals, 0..n-1
        real(dp), intent(out) :: lu(:, :) !! factored extended band (upper+2*lower+1,n)
        integer, intent(out) :: ipiv(:) !! one-based pivot row for each elimination step
        integer, intent(out) :: status !! success / invalid_shape / singular / no_convergence
        real(dp) :: value, multiplier, pivot_value
        real(dp) :: pivot_abs, candidate
        integer :: n, offset, i, j, k, last_row, last_column, pivot_row

        lu = 0.0_dp
        ipiv = 0
        status = linalg_status_invalid_shape
        n = size(ab,2)
        if (n < 1 .or. lower < 0 .or. upper < 0) return
        if (lower >= n .or. upper >= n) return
        if (size(ab,1) /= lower+upper+1) return
        if (size(lu,1) /= upper+2*lower+1 .or. size(lu,2) /= n) return
        if (size(ipiv) /= n) return
        status = linalg_status_no_convergence
        if (.not. finite_real(ab)) return
        offset = upper+lower+1
        do j=1,n
            do i=max(1,j-upper),min(n,j+lower)
                lu(offset+i-j,j)=ab(upper+1+i-j,j)
            end do
        end do
        do k=1,n
            last_row=min(n,k+lower)
            last_column=min(n,k+upper+lower)
            pivot_row=k
            pivot_abs=abs(lu(offset,k))
            do i=k+1,last_row
                candidate=abs(lu(offset+i-k,k))
                if (candidate > pivot_abs) then
                    pivot_row=i
                    pivot_abs=candidate
                end if
            end do
            if (.not. ieee_is_finite(pivot_abs)) exit
            if (pivot_abs == 0.0_dp) then
                status=linalg_status_singular
                exit
            end if
            ipiv(k)=pivot_row
            if (pivot_row /= k) then
                do j=k,last_column
                    value=lu(offset+k-j,j)
                    lu(offset+k-j,j)=lu(offset+pivot_row-j,j)
                    lu(offset+pivot_row-j,j)=value
                end do
            end if
            pivot_value=lu(offset,k)
            do i=k+1,last_row
                if (lu(offset+i-k,k) == 0.0_dp) cycle
                multiplier=lu(offset+i-k,k)/pivot_value
                lu(offset+i-k,k)=multiplier
                do j=k+1,last_column
                    lu(offset+i-j,j)=lu(offset+i-j,j)-multiplier*lu(offset+k-j,j)
                end do
            end do
        end do
        if (k > n) then
            if (finite_real(lu)) status=linalg_status_success
        end if
        if (status /= linalg_status_success) then
            lu = 0.0_dp
            ipiv = 0
        end if
    end subroutine lu_factor_banded_real

    subroutine lu_solve_banded_real_vector(lu, lower, upper, ipiv, b, x, status, trans)
        real(dp), intent(in) :: lu(:, :) !! factors produced by lu_factor_banded
        integer, intent(in) :: lower !! number of subdiagonals
        integer, intent(in) :: upper !! number of superdiagonals
        integer, intent(in) :: ipiv(:) !! one-based pivot rows returned by factorization
        real(dp), intent(in) :: b(:) !! right-hand side, length n
        real(dp), intent(out) :: x(:) !! solution, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(len=1), intent(in), optional :: trans !! 'N', 'T', or 'C'
        real(dp), allocatable :: rhs(:, :), sol(:, :)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(lu,2)
        if (n < 1 .or. size(b) /= n .or. size(x) /= n) return
        allocate(rhs(n,1),sol(n,1))
        rhs(:,1)=b
        call lu_solve_banded_real_matrix(lu,lower,upper,ipiv,rhs,sol,status,trans)
        if (status == linalg_status_success) x=sol(:,1)
    end subroutine lu_solve_banded_real_vector

    subroutine lu_solve_banded_real_matrix(lu, lower, upper, ipiv, b, x, status, trans)
        real(dp), intent(in) :: lu(:, :) !! factored extended band (upper+2*lower+1,n)
        integer, intent(in) :: lower !! number of subdiagonals
        integer, intent(in) :: upper !! number of superdiagonals
        integer, intent(in) :: ipiv(:) !! one-based pivot rows returned by factorization
        real(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solutions, shape (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(len=1), intent(in), optional :: trans !! N: A, T: transpose, C: conjugate transpose
        real(dp) :: temp(size(b,2)), multiplier
        integer :: n, offset, k, i, j, p
        character(len=1) :: op

        x = 0.0_dp
        status=linalg_status_invalid_shape
        n=size(lu,2)
        if (n < 1 .or. lower < 0 .or. upper < 0) return
        if (lower >= n .or. upper >= n) return
        if (size(lu,1) /= upper+2*lower+1) return
        if (size(ipiv) /= n .or. size(b,1) /= n) return
        if (size(x,1) /= n .or. size(x,2) /= size(b,2)) return
        op='N'
        if (present(trans)) op=trans
        if (op /= 'N' .and. op /= 'T' .and. op /= 'C') return
        do k=1,n
            if (ipiv(k) < k .or. ipiv(k) > min(n,k+lower)) return
        end do
        status=linalg_status_no_convergence
        if (.not. finite_real(lu) .or. .not. finite_real(b)) return
        offset=upper+lower+1
        status=linalg_status_singular
        if (any(abs(lu(offset,:)) == 0.0_dp)) return
        status=linalg_status_no_convergence
        x=b
        if (op == 'N') then
            ! Apply pivot k, then L_k to the RHS in the original factor order.
            do k=1,n
                p=ipiv(k)
                if (p /= k) then
                    temp=x(k,:)
                    x(k,:)=x(p,:)
                    x(p,:)=temp
                end if
                do i=k+1,min(n,k+lower)
                    x(i,:)=x(i,:)-lu(offset+i-k,k)*x(k,:)
                end do
            end do
            do k=n,1,-1
                do j=k+1,min(n,k+upper+lower)
                    x(k,:)=x(k,:)-lu(offset+k-j,j)*x(j,:)
                end do
                x(k,:)=x(k,:)/lu(offset,k)
            end do
        else
            ! A^T/A^H: U^T/U^H solve, followed by L_k^T/L_k^H
            ! and inverse pivot k, in reverse elimination order.
            do k=1,n
                do i=max(1,k-upper-lower),k-1
                    multiplier=lu(offset+i-k,k)
                    if (op == 'C') multiplier=multiplier
                    x(k,:)=x(k,:)-multiplier*x(i,:)
                end do
                multiplier=lu(offset,k)
                if (op == 'C') multiplier=multiplier
                x(k,:)=x(k,:)/multiplier
            end do
            do k=n,1,-1
                do i=k+1,min(n,k+lower)
                    multiplier=lu(offset+i-k,k)
                    if (op == 'C') multiplier=multiplier
                    x(k,:)=x(k,:)-multiplier*x(i,:)
                end do
                p=ipiv(k)
                if (p /= k) then
                    temp=x(k,:)
                    x(k,:)=x(p,:)
                    x(p,:)=temp
                end if
            end do
        end if
        if (.not. finite_real(x)) then
            x=0.0_dp
            return
        end if
        status=linalg_status_success
    end subroutine lu_solve_banded_real_matrix

    ! This is a compact, reusable analogue of LAPACK's GBTRF storage.
    ! The second padding band accounts for pivot fill; row indices inside
    ! LU are offset by upper+lower+1 relative to the true diagonal.
    subroutine lu_factor_banded_complex(ab, lower, upper, lu, ipiv, status)
        complex(dp), intent(in) :: ab(:, :) !! compact input band (lower+upper+1,n)
        integer, intent(in) :: lower !! number of subdiagonals, 0..n-1
        integer, intent(in) :: upper !! number of superdiagonals, 0..n-1
        complex(dp), intent(out) :: lu(:, :) !! factored extended band (upper+2*lower+1,n)
        integer, intent(out) :: ipiv(:) !! one-based pivot row for each elimination step
        integer, intent(out) :: status !! success / invalid_shape / singular / no_convergence
        complex(dp) :: value, multiplier, pivot_value
        real(dp) :: pivot_abs, candidate
        integer :: n, offset, i, j, k, last_row, last_column, pivot_row

        lu = cmplx(0.0_dp,0.0_dp,dp)
        ipiv = 0
        status = linalg_status_invalid_shape
        n = size(ab,2)
        if (n < 1 .or. lower < 0 .or. upper < 0) return
        if (lower >= n .or. upper >= n) return
        if (size(ab,1) /= lower+upper+1) return
        if (size(lu,1) /= upper+2*lower+1 .or. size(lu,2) /= n) return
        if (size(ipiv) /= n) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(ab)) return
        offset = upper+lower+1
        do j=1,n
            do i=max(1,j-upper),min(n,j+lower)
                lu(offset+i-j,j)=ab(upper+1+i-j,j)
            end do
        end do
        do k=1,n
            last_row=min(n,k+lower)
            last_column=min(n,k+upper+lower)
            pivot_row=k
            pivot_abs=abs(lu(offset,k))
            do i=k+1,last_row
                candidate=abs(lu(offset+i-k,k))
                if (candidate > pivot_abs) then
                    pivot_row=i
                    pivot_abs=candidate
                end if
            end do
            if (.not. ieee_is_finite(pivot_abs)) exit
            if (pivot_abs == 0.0_dp) then
                status=linalg_status_singular
                exit
            end if
            ipiv(k)=pivot_row
            if (pivot_row /= k) then
                do j=k,last_column
                    value=lu(offset+k-j,j)
                    lu(offset+k-j,j)=lu(offset+pivot_row-j,j)
                    lu(offset+pivot_row-j,j)=value
                end do
            end if
            pivot_value=lu(offset,k)
            do i=k+1,last_row
                if (lu(offset+i-k,k) == cmplx(0.0_dp,0.0_dp,dp)) cycle
                multiplier=lu(offset+i-k,k)/pivot_value
                lu(offset+i-k,k)=multiplier
                do j=k+1,last_column
                    lu(offset+i-j,j)=lu(offset+i-j,j)-multiplier*lu(offset+k-j,j)
                end do
            end do
        end do
        if (k > n) then
            if (finite_complex(lu)) status=linalg_status_success
        end if
        if (status /= linalg_status_success) then
            lu = cmplx(0.0_dp,0.0_dp,dp)
            ipiv = 0
        end if
    end subroutine lu_factor_banded_complex

    subroutine lu_solve_banded_complex_vector(lu, lower, upper, ipiv, b, x, status, trans)
        complex(dp), intent(in) :: lu(:, :) !! factors produced by lu_factor_banded
        integer, intent(in) :: lower !! number of subdiagonals
        integer, intent(in) :: upper !! number of superdiagonals
        integer, intent(in) :: ipiv(:) !! one-based pivot rows returned by factorization
        complex(dp), intent(in) :: b(:) !! right-hand side, length n
        complex(dp), intent(out) :: x(:) !! solution, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(len=1), intent(in), optional :: trans !! 'N', 'T', or 'C'
        complex(dp), allocatable :: rhs(:, :), sol(:, :)
        integer :: n

        x = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(lu,2)
        if (n < 1 .or. size(b) /= n .or. size(x) /= n) return
        allocate(rhs(n,1),sol(n,1))
        rhs(:,1)=b
        call lu_solve_banded_complex_matrix(lu,lower,upper,ipiv,rhs,sol,status,trans)
        if (status == linalg_status_success) x=sol(:,1)
    end subroutine lu_solve_banded_complex_vector

    subroutine lu_solve_banded_complex_matrix(lu, lower, upper, ipiv, b, x, status, trans)
        complex(dp), intent(in) :: lu(:, :) !! factored extended band (upper+2*lower+1,n)
        integer, intent(in) :: lower !! number of subdiagonals
        integer, intent(in) :: upper !! number of superdiagonals
        integer, intent(in) :: ipiv(:) !! one-based pivot rows returned by factorization
        complex(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solutions, shape (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(len=1), intent(in), optional :: trans !! N: A, T: transpose, C: conjugate transpose
        complex(dp) :: temp(size(b,2)), multiplier
        integer :: n, offset, k, i, j, p
        character(len=1) :: op

        x = cmplx(0.0_dp,0.0_dp,dp)
        status=linalg_status_invalid_shape
        n=size(lu,2)
        if (n < 1 .or. lower < 0 .or. upper < 0) return
        if (lower >= n .or. upper >= n) return
        if (size(lu,1) /= upper+2*lower+1) return
        if (size(ipiv) /= n .or. size(b,1) /= n) return
        if (size(x,1) /= n .or. size(x,2) /= size(b,2)) return
        op='N'
        if (present(trans)) op=trans
        if (op /= 'N' .and. op /= 'T' .and. op /= 'C') return
        do k=1,n
            if (ipiv(k) < k .or. ipiv(k) > min(n,k+lower)) return
        end do
        status=linalg_status_no_convergence
        if (.not. finite_complex(lu) .or. .not. finite_complex(b)) return
        offset=upper+lower+1
        status=linalg_status_singular
        if (any(abs(lu(offset,:)) == 0.0_dp)) return
        status=linalg_status_no_convergence
        x=b
        if (op == 'N') then
            ! Apply pivot k, then L_k to the RHS in the original factor order.
            do k=1,n
                p=ipiv(k)
                if (p /= k) then
                    temp=x(k,:)
                    x(k,:)=x(p,:)
                    x(p,:)=temp
                end if
                do i=k+1,min(n,k+lower)
                    x(i,:)=x(i,:)-lu(offset+i-k,k)*x(k,:)
                end do
            end do
            do k=n,1,-1
                do j=k+1,min(n,k+upper+lower)
                    x(k,:)=x(k,:)-lu(offset+k-j,j)*x(j,:)
                end do
                x(k,:)=x(k,:)/lu(offset,k)
            end do
        else
            ! A^T/A^H: U^T/U^H solve, followed by L_k^T/L_k^H
            ! and inverse pivot k, in reverse elimination order.
            do k=1,n
                do i=max(1,k-upper-lower),k-1
                    multiplier=lu(offset+i-k,k)
                    if (op == 'C') multiplier=conjg(multiplier)
                    x(k,:)=x(k,:)-multiplier*x(i,:)
                end do
                multiplier=lu(offset,k)
                if (op == 'C') multiplier=conjg(multiplier)
                x(k,:)=x(k,:)/multiplier
            end do
            do k=n,1,-1
                do i=k+1,min(n,k+lower)
                    multiplier=lu(offset+i-k,k)
                    if (op == 'C') multiplier=conjg(multiplier)
                    x(k,:)=x(k,:)-multiplier*x(i,:)
                end do
                p=ipiv(k)
                if (p /= k) then
                    temp=x(k,:)
                    x(k,:)=x(p,:)
                    x(p,:)=temp
                end if
            end do
        end if
        if (.not. finite_complex(x)) then
            x=cmplx(0.0_dp,0.0_dp,dp)
            return
        end if
        status=linalg_status_success
    end subroutine lu_solve_banded_complex_matrix

    pure logical function finite_real(a) result(ok)
        real(dp), intent(in) :: a(:, :) !! real matrix tested for finite entries
        ok = all(ieee_is_finite(a))
    end function finite_real

    pure logical function finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! complex matrix tested for finite entries
        ok = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex
end module scifort_linalg_banded
