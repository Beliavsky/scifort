! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Unblocked Householder QR (optionally column-pivoted) and
! QR/SVD minimum-norm least-squares solves.
module scifort_linalg_qr
    use scifort_kinds, only : dp
    use scifort_linalg_svd, only : svd_lstsq
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: qr, qr_factor, qr_unpack, qr_solve, lstsq

    interface qr
        module procedure qr_real
        module procedure qr_complex
    end interface qr

    interface qr_factor
        module procedure qr_factor_real
        module procedure qr_factor_complex
    end interface qr_factor

    interface qr_unpack
        module procedure qr_unpack_real
        module procedure qr_unpack_complex
    end interface qr_unpack

    interface qr_solve
        module procedure qr_solve_real_vector
        module procedure qr_solve_real_matrix
        module procedure qr_solve_complex_vector
        module procedure qr_solve_complex_matrix
    end interface qr_solve

    interface lstsq
        module procedure lstsq_real_vector
        module procedure lstsq_real_matrix
        module procedure lstsq_complex_vector
        module procedure lstsq_complex_matrix
    end interface lstsq

contains
    pure logical function finite_real(x) result(valid)
        real(dp), intent(in) :: x(:, :) !! values to inspect for NaNs and infinities
        valid = all(ieee_is_finite(x))
    end function finite_real

    pure real(dp) function safe_norm_real(v) result(norm)
        real(dp), intent(in) :: v(:) !! input vector, any finite values
        real(dp) :: scale
        scale = 0.0_dp
        if (size(v) > 0) scale = maxval(abs(v))
        norm = 0.0_dp
        if (scale > 0.0_dp) norm = scale * sqrt(sum((abs(v)/scale)**2))
    end function safe_norm_real

    subroutine qr_factor_real(a, af, tau, status, pivots)
        real(dp), intent(in) :: a(:, :) !! m-by-n matrix, both dimensions positive
        real(dp), intent(out) :: af(:, :) !! compact QR, same shape as a; zero on failure
        real(dp), intent(out) :: tau(:) !! min(m,n) Householder coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(out), optional :: pivots(:) !! final one-based column order; enables pivoting
        real(dp) :: alpha, beta, step
        real(dp), allocatable :: column(:)
        real(dp) :: length, trail, best, candidate
        integer :: m, n, kmax, k, j, p, col

        af = 0.0_dp
        tau = 0.0_dp
        if (present(pivots)) pivots = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        kmax = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(af,1) /= m .or. size(af,2) /= n .or. size(tau) /= kmax) return
        if (present(pivots)) then
            if (size(pivots) /= n) return
        end if
        status = linalg_status_no_convergence
        if (.not. finite_real(a)) return
        af = a
        if (present(pivots)) then
            do j = 1, n
                pivots(j) = j
            end do
        end if
        allocate(column(m))
        do k = 1, kmax
            if (present(pivots)) then
                p = k
                best = safe_norm_real(af(k:m,k))
                do j = k+1, n
                    candidate = safe_norm_real(af(k:m,j))
                    if (candidate > best) then
                        best = candidate
                        p = j
                    end if
                end do
                if (p /= k) then
                    column = af(:,k)
                    af(:,k) = af(:,p)
                    af(:,p) = column
                    col = pivots(k)
                    pivots(k) = pivots(p)
                    pivots(p) = col
                end if
            end if
            if (k == m) cycle
            trail = safe_norm_real(af(k+1:m,k))
            if (trail == 0.0_dp) cycle
            length = safe_norm_real(af(k:m,k))
            alpha = af(k,k)
            beta = -sign(length,alpha)
            tau(k) = 1.0_dp + abs(alpha)/length
            af(k+1:m,k) = af(k+1:m,k)/(alpha-beta)
            af(k,k) = beta
            do j = k+1, n
                step = (af(k,j) + dot_product(af(k+1:m,k),af(k+1:m,j)))*tau(k)
                af(k,j) = af(k,j)-step
                af(k+1:m,j) = af(k+1:m,j)-af(k+1:m,k)*step
            end do
        end do
        if (.not. finite_real(af)) then
            af = 0.0_dp
            tau = 0.0_dp
            if (present(pivots)) pivots = 0
            return
        end if
        status = linalg_status_success
    end subroutine qr_factor_real

    subroutine qr_unpack_real(af, tau, q, r, status)
        real(dp), intent(in) :: af(:, :) !! compact QR from qr_factor
        real(dp), intent(in) :: tau(:) !! corresponding Householder coefficients
        real(dp), intent(out) :: q(:, :) !! m-by-min(m,n) or m-by-m orthogonal/unitary Q
        real(dp), intent(out) :: r(:, :) !! upper trapezoidal R, same row count as q columns
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp) :: step
        integer :: m, n, kmax, qcols, k, j, i

        q = 0.0_dp
        r = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(af,1)
        n = size(af,2)
        kmax = min(m,n)
        qcols = size(q,2)
        if (m < 1 .or. n < 1) return
        if (size(tau) /= kmax) return
        if (size(q,1) /= m .or. (qcols /= m .and. qcols /= kmax)) return
        if (size(r,1) /= qcols .or. size(r,2) /= n) return
        status = linalg_status_no_convergence
        if (.not. finite_real(af)) return
        if (.not. all(ieee_is_finite(tau))) return
        do i = 1, qcols
            q(i,i) = 1.0_dp
            do j = i, n
                r(i,j) = af(i,j)
            end do
        end do
        ! Q = H(1) ... H(kmax), so build from the reflectors in reverse.
        do k = kmax, 1, -1
            if (tau(k) == 0.0_dp) cycle
            do j = 1, qcols
                step = (q(k,j)+dot_product(af(k+1:m,k),q(k+1:m,j)))*tau(k)
                q(k,j) = q(k,j)-step
                q(k+1:m,j) = q(k+1:m,j)-af(k+1:m,k)*step
            end do
        end do
        if (.not. finite_real(q) .or. .not. finite_real(r)) then
            q = 0.0_dp
            r = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine qr_unpack_real

    subroutine qr_real(a, q, r, status, pivots)
        real(dp), intent(in) :: a(:, :) !! m-by-n input matrix
        real(dp), intent(out) :: q(:, :) !! reduced or complete Q, size m-by-k or m-by-m
        real(dp), intent(out) :: r(:, :) !! corresponding k-by-n or m-by-n R
        integer, intent(out) :: status !! success, invalid_shape, or no_convergence
        integer, intent(out), optional :: pivots(:) !! final one-based column order; pivoted QR
        real(dp), allocatable :: af(:, :), tau(:)
        integer :: m, n, kmax

        q = 0.0_dp
        r = 0.0_dp
        if (present(pivots)) pivots = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        kmax = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(q,1) /= m .or. (size(q,2) /= m .and. size(q,2) /= kmax)) return
        if (size(r,1) /= size(q,2) .or. size(r,2) /= n) return
        allocate(af(m,n),tau(kmax))
        call qr_factor_real(a,af,tau,status,pivots)
        if (status /= linalg_status_success) return
        call qr_unpack_real(af,tau,q,r,status)
    end subroutine qr_real
    pure logical function finite_complex(x) result(valid)
        complex(dp), intent(in) :: x(:, :) !! values to inspect for NaNs and infinities
        valid = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end function finite_complex

    pure real(dp) function safe_norm_complex(v) result(norm)
        complex(dp), intent(in) :: v(:) !! input vector, any finite values
        real(dp) :: scale
        scale = 0.0_dp
        if (size(v) > 0) scale = maxval(abs(v))
        norm = 0.0_dp
        if (scale > 0.0_dp) norm = scale * sqrt(sum((abs(v)/scale)**2))
    end function safe_norm_complex

    subroutine qr_factor_complex(a, af, tau, status, pivots)
        complex(dp), intent(in) :: a(:, :) !! m-by-n matrix, both dimensions positive
        complex(dp), intent(out) :: af(:, :) !! compact QR, same shape as a; zero on failure
        complex(dp), intent(out) :: tau(:) !! min(m,n) Householder coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(out), optional :: pivots(:) !! final one-based column order; enables pivoting
        complex(dp) :: alpha, beta, step
        complex(dp), allocatable :: column(:)
        real(dp) :: length, trail, best, candidate
        integer :: m, n, kmax, k, j, p, col

        af = cmplx(0.0_dp,0.0_dp,dp)
        tau = cmplx(0.0_dp,0.0_dp,dp)
        if (present(pivots)) pivots = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        kmax = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(af,1) /= m .or. size(af,2) /= n .or. size(tau) /= kmax) return
        if (present(pivots)) then
            if (size(pivots) /= n) return
        end if
        status = linalg_status_no_convergence
        if (.not. finite_complex(a)) return
        af = a
        if (present(pivots)) then
            do j = 1, n
                pivots(j) = j
            end do
        end if
        allocate(column(m))
        do k = 1, kmax
            if (present(pivots)) then
                p = k
                best = safe_norm_complex(af(k:m,k))
                do j = k+1, n
                    candidate = safe_norm_complex(af(k:m,j))
                    if (candidate > best) then
                        best = candidate
                        p = j
                    end if
                end do
                if (p /= k) then
                    column = af(:,k)
                    af(:,k) = af(:,p)
                    af(:,p) = column
                    col = pivots(k)
                    pivots(k) = pivots(p)
                    pivots(p) = col
                end if
            end if
            if (k == m) cycle
            trail = safe_norm_complex(af(k+1:m,k))
            if (trail == 0.0_dp) cycle
            length = safe_norm_complex(af(k:m,k))
            alpha = af(k,k)
            beta = cmplx(-length,0.0_dp,dp)
            if (abs(alpha) > 0.0_dp) beta = -length*(alpha/abs(alpha))
            tau(k) = cmplx(1.0_dp + abs(alpha)/length,0.0_dp,dp)
            af(k+1:m,k) = af(k+1:m,k)/(alpha-beta)
            af(k,k) = beta
            do j = k+1, n
                step = (af(k,j) + dot_product(af(k+1:m,k),af(k+1:m,j)))*tau(k)
                af(k,j) = af(k,j)-step
                af(k+1:m,j) = af(k+1:m,j)-af(k+1:m,k)*step
            end do
        end do
        if (.not. finite_complex(af)) then
            af = cmplx(0.0_dp,0.0_dp,dp)
            tau = cmplx(0.0_dp,0.0_dp,dp)
            if (present(pivots)) pivots = 0
            return
        end if
        status = linalg_status_success
    end subroutine qr_factor_complex

    subroutine qr_unpack_complex(af, tau, q, r, status)
        complex(dp), intent(in) :: af(:, :) !! compact QR from qr_factor
        complex(dp), intent(in) :: tau(:) !! corresponding Householder coefficients
        complex(dp), intent(out) :: q(:, :) !! m-by-min(m,n) or m-by-m orthogonal/unitary Q
        complex(dp), intent(out) :: r(:, :) !! upper trapezoidal R, same row count as q columns
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp) :: step
        integer :: m, n, kmax, qcols, k, j, i

        q = cmplx(0.0_dp,0.0_dp,dp)
        r = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(af,1)
        n = size(af,2)
        kmax = min(m,n)
        qcols = size(q,2)
        if (m < 1 .or. n < 1) return
        if (size(tau) /= kmax) return
        if (size(q,1) /= m .or. (qcols /= m .and. qcols /= kmax)) return
        if (size(r,1) /= qcols .or. size(r,2) /= n) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(af)) return
        if (.not. all(ieee_is_finite(real(tau,dp)))) return
        if (.not. all(ieee_is_finite(aimag(tau)))) return
        do i = 1, qcols
            q(i,i) = 1.0_dp
            do j = i, n
                r(i,j) = af(i,j)
            end do
        end do
        ! Q = H(1) ... H(kmax), so build from the reflectors in reverse.
        do k = kmax, 1, -1
            if (tau(k) == cmplx(0.0_dp,0.0_dp,dp)) cycle
            do j = 1, qcols
                step = (q(k,j)+dot_product(af(k+1:m,k),q(k+1:m,j)))*tau(k)
                q(k,j) = q(k,j)-step
                q(k+1:m,j) = q(k+1:m,j)-af(k+1:m,k)*step
            end do
        end do
        if (.not. finite_complex(q) .or. .not. finite_complex(r)) then
            q = cmplx(0.0_dp,0.0_dp,dp)
            r = cmplx(0.0_dp,0.0_dp,dp)
            return
        end if
        status = linalg_status_success
    end subroutine qr_unpack_complex

    subroutine qr_complex(a, q, r, status, pivots)
        complex(dp), intent(in) :: a(:, :) !! m-by-n input matrix
        complex(dp), intent(out) :: q(:, :) !! reduced or complete Q, size m-by-k or m-by-m
        complex(dp), intent(out) :: r(:, :) !! corresponding k-by-n or m-by-n R
        integer, intent(out) :: status !! success, invalid_shape, or no_convergence
        integer, intent(out), optional :: pivots(:) !! final one-based column order; pivoted QR
        complex(dp), allocatable :: af(:, :), tau(:)
        integer :: m, n, kmax

        q = cmplx(0.0_dp,0.0_dp,dp)
        r = cmplx(0.0_dp,0.0_dp,dp)
        if (present(pivots)) pivots = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        kmax = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(q,1) /= m .or. (size(q,2) /= m .and. size(q,2) /= kmax)) return
        if (size(r,1) /= size(q,2) .or. size(r,2) /= n) return
        allocate(af(m,n),tau(kmax))
        call qr_factor_complex(a,af,tau,status,pivots)
        if (status /= linalg_status_success) return
        call qr_unpack_complex(af,tau,q,r,status)
    end subroutine qr_complex
    subroutine lstsq_real_vector(a, b, x, status, rank, rcond, method)
        real(dp), intent(in) :: a(:, :) !! m-by-n coefficient matrix
        real(dp), intent(in) :: b(:) !! m-element right-hand side
        real(dp), intent(out) :: x(:) !! n-element least-squares/minimum-norm solution
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! estimated rank (0 on invalid input)
        real(dp), intent(in), optional :: rcond !! relative threshold; default max(m,n)*epsilon
        character(*), intent(in), optional :: method !! 'auto' (default), 'qr', or 'svd'
        real(dp), allocatable :: rhs(:, :), solution(:, :)
        integer :: m, n

        x = 0.0_dp
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (m < 1 .or. n < 1 .or. size(b) /= m .or. size(x) /= n) return
        allocate(rhs(m,1),solution(n,1))
        rhs(:,1) = b
        call lstsq_real_matrix(a,rhs,solution,status,rank,rcond,method)
        if (status == linalg_status_success) x = solution(:,1)
    end subroutine lstsq_real_vector

    subroutine lstsq_real_matrix(a, b, x, status, rank, rcond, method)
        real(dp), intent(in) :: a(:, :) !! m-by-n coefficient matrix
        real(dp), intent(in) :: b(:, :) !! right-hand sides (m,nrhs)
        real(dp), intent(out) :: x(:, :) !! solutions (n,nrhs), zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! estimated QR or SVD rank
        real(dp), intent(in), optional :: rcond !! relative threshold; default max(m,n)*epsilon
        character(*), intent(in), optional :: method !! 'auto' (default), 'qr', or 'svd'
        real(dp), allocatable :: af(:, :), tau(:), y(:, :), z(:, :), work(:, :)
        real(dp) :: step
        integer, allocatable :: pivots(:)
        real(dp) :: cutoff, maxdiag, threshold
        integer :: m, n, nrhs, kmax, j, k, i, numerical_rank
        character(:), allocatable :: choice

        x = 0.0_dp
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        nrhs = size(b,2)
        kmax = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(b,1) /= m .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        choice = 'auto'
        if (present(method)) choice = trim(adjustl(method))
        if (choice /= 'auto' .and. choice /= 'qr' .and. choice /= 'svd') return
        status = linalg_status_no_convergence
        if (.not. finite_real(a) .or. .not. finite_real(b)) return
        if (choice == 'svd') then
            call svd_lstsq(a, b, x, status, rank, rcond)
            return
        end if
        if (m >= n) then
            allocate(af(m,n),tau(n),pivots(n))
            call qr_factor_real(a,af,tau,status,pivots)
        else
            allocate(af(n,m),tau(m),pivots(m))
            call qr_factor_real(transpose(a),af,tau,status,pivots)
        end if
        if (status /= linalg_status_success) return
        maxdiag = 0.0_dp
        do i = 1, kmax
            maxdiag = max(maxdiag,abs(af(i,i)))
        end do
        threshold = cutoff*maxdiag
        numerical_rank = 0
        do i = 1, kmax
            if (abs(af(i,i)) > threshold) numerical_rank = numerical_rank+1
        end do
        if (present(rank)) rank = numerical_rank
        if (numerical_rank < kmax) then
            if (choice == 'auto') then
                call svd_lstsq(a, b, x, status, rank, rcond)
            else
                status = linalg_status_singular
            end if
            return
        end if
        status = linalg_status_no_convergence
        if (m >= n) then
            allocate(y(m,nrhs),z(n,nrhs))
            ! Apply Q^H to B. The stored reflectors are Hermitian.
            y = b
            do k = 1, n
                if (tau(k) == 0.0_dp) cycle
                do j = 1, nrhs
                    step = (y(k,j)+dot_product(af(k+1:m,k),y(k+1:m,j)))*tau(k)
                    y(k,j) = y(k,j)-step
                    y(k+1:m,j) = y(k+1:m,j)-af(k+1:m,k)*step
                end do
            end do
            z = 0.0_dp
            do j = 1, nrhs
                do i = n, 1, -1
                    step = y(i,j)
                    do k = i+1, n
                        step = step-af(i,k)*z(k,j)
                    end do
                    z(i,j) = step/af(i,i)
                end do
            end do
            do i = 1, n
                x(pivots(i),:) = z(i,:)
            end do
        else
            allocate(y(m,nrhs),work(n,nrhs))
            ! A^H P = Q R: R^H y = P^T B, and X = Q [y; 0].
            y = 0.0_dp
            do j = 1, nrhs
                do i = 1, m
                    step = b(pivots(i),j)
                    do k = 1, i-1
                        step = step-af(k,i)*y(k,j)
                    end do
                    y(i,j) = step/af(i,i)
                end do
            end do
            work = 0.0_dp
            work(1:m,:) = y
            do k = m, 1, -1
                if (tau(k) == 0.0_dp) cycle
                do j = 1, nrhs
                    step = (work(k,j)+dot_product(af(k+1:n,k),work(k+1:n,j)))*tau(k)
                    work(k,j) = work(k,j)-step
                    work(k+1:n,j) = work(k+1:n,j)-af(k+1:n,k)*step
                end do
            end do
            x = work
        end if
        if (.not. finite_real(x)) then
            x = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine lstsq_real_matrix

    subroutine lstsq_complex_vector(a, b, x, status, rank, rcond, method)
        complex(dp), intent(in) :: a(:, :) !! m-by-n coefficient matrix
        complex(dp), intent(in) :: b(:) !! m-element right-hand side
        complex(dp), intent(out) :: x(:) !! n-element least-squares/minimum-norm solution
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! estimated rank (0 on invalid input)
        real(dp), intent(in), optional :: rcond !! relative threshold; default max(m,n)*epsilon
        character(*), intent(in), optional :: method !! 'auto' (default), 'qr', or 'svd'
        complex(dp), allocatable :: rhs(:, :), solution(:, :)
        integer :: m, n

        x = cmplx(0.0_dp,0.0_dp,dp)
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (m < 1 .or. n < 1 .or. size(b) /= m .or. size(x) /= n) return
        allocate(rhs(m,1),solution(n,1))
        rhs(:,1) = b
        call lstsq_complex_matrix(a,rhs,solution,status,rank,rcond,method)
        if (status == linalg_status_success) x = solution(:,1)
    end subroutine lstsq_complex_vector

    subroutine lstsq_complex_matrix(a, b, x, status, rank, rcond, method)
        complex(dp), intent(in) :: a(:, :) !! m-by-n coefficient matrix
        complex(dp), intent(in) :: b(:, :) !! right-hand sides (m,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solutions (n,nrhs), zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! estimated QR or SVD rank
        real(dp), intent(in), optional :: rcond !! relative threshold; default max(m,n)*epsilon
        character(*), intent(in), optional :: method !! 'auto' (default), 'qr', or 'svd'
        complex(dp), allocatable :: af(:, :), tau(:), y(:, :), z(:, :), work(:, :)
        complex(dp) :: step
        integer, allocatable :: pivots(:)
        real(dp) :: cutoff, maxdiag, threshold
        integer :: m, n, nrhs, kmax, j, k, i, numerical_rank
        character(:), allocatable :: choice

        x = cmplx(0.0_dp,0.0_dp,dp)
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        nrhs = size(b,2)
        kmax = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(b,1) /= m .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        choice = 'auto'
        if (present(method)) choice = trim(adjustl(method))
        if (choice /= 'auto' .and. choice /= 'qr' .and. choice /= 'svd') return
        status = linalg_status_no_convergence
        if (.not. finite_complex(a) .or. .not. finite_complex(b)) return
        if (choice == 'svd') then
            call svd_lstsq(a, b, x, status, rank, rcond)
            return
        end if
        if (m >= n) then
            allocate(af(m,n),tau(n),pivots(n))
            call qr_factor_complex(a,af,tau,status,pivots)
        else
            allocate(af(n,m),tau(m),pivots(m))
            call qr_factor_complex(transpose(conjg(a)),af,tau,status,pivots)
        end if
        if (status /= linalg_status_success) return
        maxdiag = 0.0_dp
        do i = 1, kmax
            maxdiag = max(maxdiag,abs(af(i,i)))
        end do
        threshold = cutoff*maxdiag
        numerical_rank = 0
        do i = 1, kmax
            if (abs(af(i,i)) > threshold) numerical_rank = numerical_rank+1
        end do
        if (present(rank)) rank = numerical_rank
        if (numerical_rank < kmax) then
            if (choice == 'auto') then
                call svd_lstsq(a, b, x, status, rank, rcond)
            else
                status = linalg_status_singular
            end if
            return
        end if
        status = linalg_status_no_convergence
        if (m >= n) then
            allocate(y(m,nrhs),z(n,nrhs))
            ! Apply Q^H to B. The stored reflectors are Hermitian.
            y = b
            do k = 1, n
                if (tau(k) == cmplx(0.0_dp,0.0_dp,dp)) cycle
                do j = 1, nrhs
                    step = (y(k,j)+dot_product(af(k+1:m,k),y(k+1:m,j)))*tau(k)
                    y(k,j) = y(k,j)-step
                    y(k+1:m,j) = y(k+1:m,j)-af(k+1:m,k)*step
                end do
            end do
            z = cmplx(0.0_dp,0.0_dp,dp)
            do j = 1, nrhs
                do i = n, 1, -1
                    step = y(i,j)
                    do k = i+1, n
                        step = step-af(i,k)*z(k,j)
                    end do
                    z(i,j) = step/af(i,i)
                end do
            end do
            do i = 1, n
                x(pivots(i),:) = z(i,:)
            end do
        else
            allocate(y(m,nrhs),work(n,nrhs))
            ! A^H P = Q R: R^H y = P^T B, and X = Q [y; 0].
            y = cmplx(0.0_dp,0.0_dp,dp)
            do j = 1, nrhs
                do i = 1, m
                    step = b(pivots(i),j)
                    do k = 1, i-1
                        step = step-conjg(af(k,i))*y(k,j)
                    end do
                    y(i,j) = step/conjg(af(i,i))
                end do
            end do
            work = cmplx(0.0_dp,0.0_dp,dp)
            work(1:m,:) = y
            do k = m, 1, -1
                if (tau(k) == cmplx(0.0_dp,0.0_dp,dp)) cycle
                do j = 1, nrhs
                    step = (work(k,j)+dot_product(af(k+1:n,k),work(k+1:n,j)))*tau(k)
                    work(k,j) = work(k,j)-step
                    work(k+1:n,j) = work(k+1:n,j)-af(k+1:n,k)*step
                end do
            end do
            x = work
        end if
        if (.not. finite_complex(x)) then
            x = cmplx(0.0_dp,0.0_dp,dp)
            return
        end if
        status = linalg_status_success
    end subroutine lstsq_complex_matrix

    subroutine qr_solve_real_vector(af, tau, pivots, b, x, status, rank, rcond)
        real(dp), intent(in) :: af(:, :) !! compact pivoted QR from qr_factor, m >= n
        real(dp), intent(in) :: tau(:) !! corresponding length-n Householder coefficients
        integer, intent(in) :: pivots(:) !! final one-based column permutation, length n
        real(dp), intent(in) :: b(:) !! m-element right-hand side
        real(dp), intent(out) :: x(:) !! n-element least-squares solution; zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! estimated numerical rank
        real(dp), intent(in), optional :: rcond !! relative diagonal threshold; default m*epsilon
        real(dp), allocatable :: rhs(:, :), solution(:, :)
        integer :: m, n

        x = 0.0_dp
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(af,1)
        n = size(af,2)
        if (m < n .or. n < 1 .or. size(b) /= m .or. size(x) /= n) return
        allocate(rhs(m,1),solution(n,1))
        rhs(:,1) = b
        call qr_solve_real_matrix(af,tau,pivots,rhs,solution,status,rank,rcond)
        if (status == linalg_status_success) x = solution(:,1)
    end subroutine qr_solve_real_vector

    subroutine qr_solve_real_matrix(af, tau, pivots, b, x, status, rank, rcond)
        real(dp), intent(in) :: af(:, :) !! compact pivoted QR from qr_factor, m >= n
        real(dp), intent(in) :: tau(:) !! corresponding length-n Householder coefficients
        integer, intent(in) :: pivots(:) !! final one-based column permutation, length n
        real(dp), intent(in) :: b(:, :) !! right-hand sides, (m,nrhs)
        real(dp), intent(out) :: x(:, :) !! least-squares solutions, (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! estimated numerical rank
        real(dp), intent(in), optional :: rcond !! relative diagonal threshold; default m*epsilon
        real(dp), allocatable :: y(:, :), z(:, :)
        real(dp) :: step
        real(dp) :: cutoff, maxdiag
        integer :: m, n, nrhs, i, j, k, numerical_rank

        x = 0.0_dp
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(af,1)
        n = size(af,2)
        nrhs = size(b,2)
        if (m < n .or. n < 1 .or. size(tau) /= n .or. size(pivots) /= n) return
        if (size(b,1) /= m .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        cutoff = real(m,dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        do i = 1, n
            if (pivots(i) < 1 .or. pivots(i) > n) return
            if (count(pivots == pivots(i)) /= 1) return
        end do
        status = linalg_status_no_convergence
        if (.not. finite_real(af) .or. .not. finite_real(b)) return
        if (.not. all(ieee_is_finite(tau))) return
        maxdiag = 0.0_dp
        do i = 1, n
            maxdiag = max(maxdiag,abs(af(i,i)))
        end do
        numerical_rank = 0
        do i = 1, n
            if (abs(af(i,i)) > cutoff*maxdiag) numerical_rank = numerical_rank+1
        end do
        if (present(rank)) rank = numerical_rank
        if (numerical_rank < n) then
            status = linalg_status_singular
            return
        end if
        status = linalg_status_no_convergence
        allocate(y(m,nrhs),z(n,nrhs))
        y = b
        ! Apply the Hermitian Householder reflectors to obtain Q^H B.
        do k = 1, n
            if (tau(k) == 0.0_dp) cycle
            do j = 1, nrhs
                step = (y(k,j)+dot_product(af(k+1:m,k),y(k+1:m,j)))*tau(k)
                y(k,j) = y(k,j)-step
                y(k+1:m,j) = y(k+1:m,j)-af(k+1:m,k)*step
            end do
        end do
        z = 0.0_dp
        do j = 1, nrhs
            do i = n, 1, -1
                step = y(i,j)
                do k = i+1, n
                    step = step-af(i,k)*z(k,j)
                end do
                z(i,j) = step/af(i,i)
            end do
        end do
        do i = 1, n
            x(pivots(i),:) = z(i,:)
        end do
        if (.not. finite_real(x)) then
            x = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine qr_solve_real_matrix

    subroutine qr_solve_complex_vector(af, tau, pivots, b, x, status, rank, rcond)
        complex(dp), intent(in) :: af(:, :) !! compact pivoted QR from qr_factor, m >= n
        complex(dp), intent(in) :: tau(:) !! corresponding length-n Householder coefficients
        integer, intent(in) :: pivots(:) !! final one-based column permutation, length n
        complex(dp), intent(in) :: b(:) !! m-element right-hand side
        complex(dp), intent(out) :: x(:) !! n-element least-squares solution; zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! estimated numerical rank
        real(dp), intent(in), optional :: rcond !! relative diagonal threshold; default m*epsilon
        complex(dp), allocatable :: rhs(:, :), solution(:, :)
        integer :: m, n

        x = cmplx(0.0_dp,0.0_dp,dp)
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(af,1)
        n = size(af,2)
        if (m < n .or. n < 1 .or. size(b) /= m .or. size(x) /= n) return
        allocate(rhs(m,1),solution(n,1))
        rhs(:,1) = b
        call qr_solve_complex_matrix(af,tau,pivots,rhs,solution,status,rank,rcond)
        if (status == linalg_status_success) x = solution(:,1)
    end subroutine qr_solve_complex_vector

    subroutine qr_solve_complex_matrix(af, tau, pivots, b, x, status, rank, rcond)
        complex(dp), intent(in) :: af(:, :) !! compact pivoted QR from qr_factor, m >= n
        complex(dp), intent(in) :: tau(:) !! corresponding length-n Householder coefficients
        integer, intent(in) :: pivots(:) !! final one-based column permutation, length n
        complex(dp), intent(in) :: b(:, :) !! right-hand sides, (m,nrhs)
        complex(dp), intent(out) :: x(:, :) !! least-squares solutions, (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! estimated numerical rank
        real(dp), intent(in), optional :: rcond !! relative diagonal threshold; default m*epsilon
        complex(dp), allocatable :: y(:, :), z(:, :)
        complex(dp) :: step
        real(dp) :: cutoff, maxdiag
        integer :: m, n, nrhs, i, j, k, numerical_rank

        x = cmplx(0.0_dp,0.0_dp,dp)
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(af,1)
        n = size(af,2)
        nrhs = size(b,2)
        if (m < n .or. n < 1 .or. size(tau) /= n .or. size(pivots) /= n) return
        if (size(b,1) /= m .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        cutoff = real(m,dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        do i = 1, n
            if (pivots(i) < 1 .or. pivots(i) > n) return
            if (count(pivots == pivots(i)) /= 1) return
        end do
        status = linalg_status_no_convergence
        if (.not. finite_complex(af) .or. .not. finite_complex(b)) return
        if (.not. all(ieee_is_finite(real(tau,dp)))) return
        if (.not. all(ieee_is_finite(aimag(tau)))) return
        maxdiag = 0.0_dp
        do i = 1, n
            maxdiag = max(maxdiag,abs(af(i,i)))
        end do
        numerical_rank = 0
        do i = 1, n
            if (abs(af(i,i)) > cutoff*maxdiag) numerical_rank = numerical_rank+1
        end do
        if (present(rank)) rank = numerical_rank
        if (numerical_rank < n) then
            status = linalg_status_singular
            return
        end if
        status = linalg_status_no_convergence
        allocate(y(m,nrhs),z(n,nrhs))
        y = b
        ! Apply the Hermitian Householder reflectors to obtain Q^H B.
        do k = 1, n
            if (tau(k) == cmplx(0.0_dp,0.0_dp,dp)) cycle
            do j = 1, nrhs
                step = (y(k,j)+dot_product(af(k+1:m,k),y(k+1:m,j)))*tau(k)
                y(k,j) = y(k,j)-step
                y(k+1:m,j) = y(k+1:m,j)-af(k+1:m,k)*step
            end do
        end do
        z = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1, nrhs
            do i = n, 1, -1
                step = y(i,j)
                do k = i+1, n
                    step = step-af(i,k)*z(k,j)
                end do
                z(i,j) = step/af(i,i)
            end do
        end do
        do i = 1, n
            x(pivots(i),:) = z(i,:)
        end do
        if (.not. finite_complex(x)) then
            x = cmplx(0.0_dp,0.0_dp,dp)
            return
        end if
        status = linalg_status_success
    end subroutine qr_solve_complex_matrix

end module scifort_linalg_qr
