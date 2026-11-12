! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Complex Hermitian Lanczos/Rayleigh-Ritz eigenpairs with two-pass full
! reorthogonalization.  The small real-symmetric Ritz problem is solved by a
! dependency-free Jacobi method.  This is one pass (no implicit restarting).
module scifort_sparse_linalg_complex_eigen
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, &
        csr_complex_matvec, csc_complex_to_csr
    implicit none
    private

    integer, parameter, public :: complex_eigen_success = 0
    integer, parameter, public :: complex_eigen_no_convergence = 1
    integer, parameter, public :: complex_eigen_invalid_input = -1
    integer, parameter, public :: complex_eigen_breakdown = -2
    public :: eigsh_complex

    interface eigsh_complex
        module procedure eigsh_complex_csr, eigsh_complex_csc
    end interface eigsh_complex

contains

    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Candidate complex scalar.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    pure real(dp) function complex_norm(z) result(value)
        complex(dp), intent(in) :: z(:) !! Vector to measure with scaled real norm.
        value = norm2(abs(z))
    end function complex_norm

    pure integer function find_col(a, row, column) result(pos)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(in) :: row !! One-based row.
        integer, intent(in) :: column !! One-based column to find.
        integer :: lo, hi, mid

        lo = a%indptr(row)
        hi = a%indptr(row+1)-1
        pos = 0
        do while (lo <= hi)
            mid = lo+(hi-lo)/2
            if (a%indices(mid) == column) then
                pos = mid
                return
            else if (a%indices(mid) < column) then
                lo = mid+1
            else
                hi = mid-1
            end if
        end do
    end function find_col

    pure logical function is_hermitian(a) result(ok)
        type(csr_complex_matrix), intent(in) :: a !! Valid square CSR matrix.
        integer :: i, j, p, q
        complex(dp) :: reflected

        ok = .false.
        if (.not. valid_csr_complex(a)) return
        if (a%n_rows /= a%n_cols) return
        do i = 1,a%n_rows
            do p = a%indptr(i),a%indptr(i+1)-1
                j = a%indices(p)
                q = find_col(a,j,i)
                reflected = cmplx(0.0_dp,0.0_dp,dp)
                if (q > 0) reflected = conjg(a%data(q))
                if (abs(a%data(p)-reflected) > &
                    1.0e-12_dp*max(1.0_dp,abs(a%data(p)),abs(reflected))) return
            end do
        end do
        ok = .true.
    end function is_hermitian

    subroutine jacobi_symmetric(matrix, eigenvalues, eigenvectors, info)
        real(dp), intent(in) :: matrix(:,:) !! Real symmetric projected matrix (m,m).
        real(dp), intent(out) :: eigenvalues(:) !! Eigenvalues (m).
        real(dp), intent(out) :: eigenvectors(:,:) !! Orthonormal columns (m,m).
        integer, intent(out) :: info !! 0 for convergence, -2 on numerical failure.
        real(dp), allocatable :: work(:,:)
        real(dp) :: app, aqq, apq, tau, t, c, s, arp, arq, wrp, wrq, threshold, largest
        integer :: m, i, j, p, q, sweep

        m = size(matrix,1)
        eigenvectors = 0.0_dp
        do i = 1,m
            eigenvectors(i,i) = 1.0_dp
        end do
        work = matrix
        info = complex_eigen_breakdown
        if (m == 0) then
            info = 0
            return
        end if
        do sweep = 1, max(30,12*m)
            largest = 0.0_dp
            do p = 1,m-1
                do q = p+1,m
                    largest = max(largest,abs(work(p,q)))
                end do
            end do
            threshold = 4.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(work)))
            if (largest <= threshold) then
                eigenvalues = [(work(i,i),i=1,m)]
                info = 0
                return
            end if
            do p = 1,m-1
                do q = p+1,m
                    apq = work(p,q)
                    if (abs(apq) <= threshold/max(1,m)) cycle
                    app = work(p,p)
                    aqq = work(q,q)
                    tau = (aqq-app)/(2.0_dp*apq)
                    if (tau >= 0.0_dp) then
                        t = 1.0_dp/(abs(tau)+sqrt(1.0_dp+tau*tau))
                    else
                        t = -1.0_dp/(abs(tau)+sqrt(1.0_dp+tau*tau))
                    end if
                    c = 1.0_dp/sqrt(1.0_dp+t*t)
                    s = t*c
                    do j = 1,m
                        if (j == p .or. j == q) cycle
                        arp = work(j,p)
                        arq = work(j,q)
                        work(j,p) = c*arp-s*arq
                        work(p,j) = work(j,p)
                        work(j,q) = s*arp+c*arq
                        work(q,j) = work(j,q)
                    end do
                    work(p,p) = app-t*apq
                    work(q,q) = aqq+t*apq
                    work(p,q) = 0.0_dp
                    work(q,p) = 0.0_dp
                    do j = 1,m
                        wrp = eigenvectors(j,p)
                        wrq = eigenvectors(j,q)
                        eigenvectors(j,p) = c*wrp-s*wrq
                        eigenvectors(j,q) = s*wrp+c*wrq
                    end do
                end do
            end do
        end do
    end subroutine jacobi_symmetric

    subroutine eigsh_complex_csc(a, k, eigenvalues, eigenvectors, info, which, tol, &
            ncv, v0, iterations, residual_norms, check_hermitian)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical complex CSC input.
        integer, intent(in) :: k !! Requested eigenpair count, 1 <= k < n.
        real(dp), allocatable, intent(out) :: eigenvalues(:) !! Real Hermitian Ritz values.
        complex(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Complex Ritz vectors (n,k).
        integer, intent(out) :: info !! 0 converged, 1 insufficient Lanczos basis, negative input/error.
        character(len=*), intent(in), optional :: which !! LA, SA, LM, SM (default LM).
        real(dp), intent(in), optional :: tol !! Relative true-residual tolerance, default 1e-8.
        integer, intent(in), optional :: ncv !! Maximum Lanczos basis, k+1..n; default min(n,max(2*k+8,24)).
        complex(dp), intent(in), optional :: v0(:) !! Starting vector of length n.
        integer, intent(out), optional :: iterations !! Sparse matvec count.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True Euclidean residuals.
        logical, intent(in), optional :: check_hermitian !! Verify matrix Hermitian, default true.
        type(csr_complex_matrix) :: csr
        integer :: status

        info = complex_eigen_invalid_input
        if (present(iterations)) iterations = 0
        allocate(eigenvalues(0), eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        deallocate(eigenvalues,eigenvectors)
        if (present(residual_norms)) deallocate(residual_norms)
        call eigsh_complex_csr(csr,k,eigenvalues,eigenvectors,info,which,tol, &
            ncv,v0,iterations,residual_norms,check_hermitian)
    end subroutine eigsh_complex_csc

    subroutine eigsh_complex_csr(a, k, eigenvalues, eigenvectors, info, which, tol, &
            ncv, v0, iterations, residual_norms, check_hermitian)
        type(csr_complex_matrix), intent(in) :: a !! Square canonical complex CSR input.
        integer, intent(in) :: k !! Requested eigenpair count, 1 <= k < n.
        real(dp), allocatable, intent(out) :: eigenvalues(:) !! Real Hermitian Ritz values.
        complex(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Complex Ritz vectors (n,k).
        integer, intent(out) :: info !! 0 converged, 1 insufficient Lanczos basis, negative input/error.
        character(len=*), intent(in), optional :: which !! LA, SA, LM, SM (default LM).
        real(dp), intent(in), optional :: tol !! Relative true-residual tolerance, default 1e-8.
        integer, intent(in), optional :: ncv !! Maximum Lanczos basis, k+1..n; default min(n,max(2*k+8,24)).
        complex(dp), intent(in), optional :: v0(:) !! Starting vector of length n.
        integer, intent(out), optional :: iterations !! Sparse matvec count.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True Euclidean residuals.
        logical, intent(in), optional :: check_hermitian !! Verify matrix Hermitian, default true.
        complex(dp), allocatable :: basis(:,:), z(:), y(:), ay(:), seed(:), seed_best(:)
        real(dp), allocatable :: alpha(:), beta(:), small(:,:), values(:), vecs(:,:), residuals(:)
        integer, allocatable :: order(:)
        integer :: n, mmax, j, t, status, idx, best, i, order_temp, sweeps, available
        real(dp) :: tolerance, beta_next, matrix_bound, &
            magnitude, candidate, max_score
        complex(dp) :: dot
        character(len=2) :: target
        logical :: validate_hermitian

        info = complex_eigen_invalid_input
        if (present(iterations)) iterations = 0
        allocate(eigenvalues(0), eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csr_complex(a)) return
        n = a%n_rows
        if (n /= a%n_cols .or. k < 1 .or. k >= n) return
        target = 'LM'
        if (present(which)) then
            if (len_trim(which) /= 2) return
            target = which
        end if
        if (target /= 'LA' .and. target /= 'SA' .and. &
            target /= 'LM' .and. target /= 'SM') return
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        if (tolerance == 0.0_dp) tolerance = 100.0_dp*epsilon(1.0_dp)
        mmax = min(n,max(2*k+8,24))
        if (present(ncv)) mmax = ncv
        if (mmax <= k .or. mmax > n) return
        if (present(v0)) then
            if (size(v0) /= n) return
            if (.not. all(finite_complex(v0))) return
            if (complex_norm(v0) == 0.0_dp) return
        end if
        validate_hermitian = .true.
        if (present(check_hermitian)) validate_hermitian = check_hermitian
        if (validate_hermitian) then
            if (.not. is_hermitian(a)) return
        end if
        matrix_bound = 0.0_dp
        do i = 1,n
            matrix_bound = max(matrix_bound,sum(abs(a%data(a%indptr(i):a%indptr(i+1)-1))))
        end do
        if (.not. ieee_is_finite(matrix_bound)) return

        allocate(basis(n,mmax),z(n),y(n),ay(n),seed(n),seed_best(n), &
            alpha(mmax),beta(mmax),small(mmax,mmax),values(mmax), &
            vecs(mmax,mmax),order(mmax),residuals(k))
        if (present(v0)) then
            basis(:,1) = v0/complex_norm(v0)
        else
            do i = 1,n
                basis(i,1) = cmplx(sin(real(i+2,dp)*0.719_dp)+cos(real(i+1,dp)*1.337_dp), &
                    cos(real(i+3,dp)*0.511_dp)-sin(real(i+1,dp)*0.239_dp),dp)
            end do
            basis(:,1) = basis(:,1)/complex_norm(basis(:,1))
        end if
        beta = 0.0_dp
        alpha = 0.0_dp
        available = mmax
        do j = 1,mmax
            z = csr_complex_matvec(a,basis(:,j),status)
            if (status /= 0) then
                info = complex_eigen_breakdown
                return
            end if
            if (j > 1) z = z - beta(j-1)*basis(:,j-1)
            alpha(j) = real(dot_product(basis(:,j),z),dp)
            z = z - alpha(j)*basis(:,j)
            ! Reorthogonalize twice to handle clustered eigenvalues.
            do sweeps = 1,2
                do t = 1,j
                    dot = dot_product(basis(:,t),z)
                    z = z - dot*basis(:,t)
                end do
            end do
            if (.not. all(finite_complex(z)) .or. .not. ieee_is_finite(alpha(j))) then
                info = complex_eigen_breakdown
                return
            end if
            if (j == mmax) exit
            beta_next = complex_norm(z)
            if (.not. ieee_is_finite(beta_next)) then
                info = complex_eigen_breakdown
                return
            end if
            if (beta_next <= 50.0_dp*epsilon(1.0_dp)*max(1.0_dp,matrix_bound)) then
                ! An invariant Krylov subspace may have too few distinct
                ! eigenvalues.  Start an independent block to capture repeats.
                ! Its link in the real projected tridiagonal is exactly zero.
                beta(j) = 0.0_dp
                best = 0
                max_score = 0.0_dp
                do idx = 1,n
                    seed = cmplx(0.0_dp,0.0_dp,dp)
                    seed(idx) = cmplx(1.0_dp,0.0_dp,dp)
                    do sweeps = 1,2
                        do t = 1,j
                            seed = seed-dot_product(basis(:,t),seed)*basis(:,t)
                        end do
                    end do
                    candidate = complex_norm(seed)
                    if (candidate > max_score) then
                        max_score = candidate
                        best = idx
                        seed_best = seed
                    end if
                    if (candidate > 0.5_dp) exit
                end do
                if (best == 0 .or. max_score < 1.0e-10_dp) then
                    available = j
                    exit
                end if
                basis(:,j+1) = seed_best/max_score
            else
                beta(j) = beta_next
                basis(:,j+1) = z/beta_next
            end if
        end do
        if (present(iterations)) iterations = available
        small(:available,:available) = 0.0_dp
        do j = 1,available
            small(j,j) = alpha(j)
            if (j < available) then
                small(j,j+1) = beta(j)
                small(j+1,j) = beta(j)
            end if
        end do
        call jacobi_symmetric(small(:available,:available), values(:available), &
            vecs(:available,:available),status)
        if (status /= 0) then
            info = complex_eigen_breakdown
            return
        end if
        do i = 1,available
            order(i) = i
        end do
        ! Stable insertion sort in requested selection order; output is in
        ! ascending numerical order for SA/SM, descending for LA/LM.
        do i = 2,available
            order_temp = order(i)
            j = i-1
            do while (j >= 1)
                magnitude = values(order(j))
                candidate = values(order_temp)
                if (target == 'LM' .or. target == 'SM') then
                    magnitude = abs(magnitude)
                    candidate = abs(candidate)
                end if
                if (target == 'LA' .or. target == 'LM') then
                    if (magnitude >= candidate) exit
                else
                    if (magnitude <= candidate) exit
                end if
                order(j+1) = order(j)
                j = j-1
            end do
            order(j+1) = order_temp
        end do
        deallocate(eigenvalues,eigenvectors)
        allocate(eigenvalues(k),eigenvectors(n,k))
        eigenvalues = 0.0_dp
        eigenvectors = cmplx(0.0_dp,0.0_dp,dp)
        residuals = huge(1.0_dp)
        info = complex_eigen_no_convergence
        if (available < k) then
            if (present(residual_norms)) residual_norms = residuals
            return
        end if
        do i = 1,k
            eigenvalues(i) = values(order(i))
            y = matmul(basis(:,:available), cmplx(vecs(:available,order(i)),0.0_dp,dp))
            y = y/complex_norm(y)
            ! Deterministic eigenvector phase (largest component real-positive).
            idx = maxloc(abs(y),dim=1)
            if (abs(y(idx)) > 0.0_dp) y = y * conjg(y(idx))/abs(y(idx))
            eigenvectors(:,i) = y
            ay = csr_complex_matvec(a,y,status)
            if (status /= 0) then
                info = complex_eigen_breakdown
                return
            end if
            residuals(i) = complex_norm(ay-eigenvalues(i)*y)
        end do
        if (.not. all(ieee_is_finite(residuals))) then
            info = complex_eigen_breakdown
        else if (all(residuals <= tolerance*max(1.0_dp,matrix_bound))) then
            info = complex_eigen_success
        end if
        if (present(residual_norms)) residual_norms = residuals
    end subroutine eigsh_complex_csr

end module scifort_sparse_linalg_complex_eigen
