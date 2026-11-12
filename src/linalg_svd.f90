! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Economy and full SVD, singular values, and minimum-norm solves,
! using one-sided Jacobi rotations.
! Adapted mathematically from SciFort's independently implemented polar
! Jacobi iteration; no external LAPACK/BLAS dependency.
module scifort_linalg_svd
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf
    implicit none
    private

    public :: svd, svdvals, diagsvd, pinv, svd_lstsq, cond2, matrix_rank

    interface svd
        module procedure svd_real
        module procedure svd_complex
    end interface

    interface svdvals
        module procedure svdvals_real
        module procedure svdvals_complex
    end interface

    interface cond2
        module procedure cond2_real
        module procedure cond2_complex
    end interface

    interface matrix_rank
        module procedure matrix_rank_real
        module procedure matrix_rank_complex
    end interface

    interface pinv
        module procedure pinv_real
        module procedure pinv_complex
    end interface

    interface svd_lstsq
        module procedure svd_lstsq_real
        module procedure svd_lstsq_complex
    end interface

contains

    pure logical function finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! input matrix to validate
        ok = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex

    subroutine svd_real(a, u, s, vh, status, max_sweeps, full_matrices)
        real(dp), intent(in) :: a(:, :) !! finite m-by-n real matrix, m,n >= 1
        real(dp), intent(out) :: u(:, :) !! m-by-k or m-by-m left singular vectors
        real(dp), intent(out) :: s(:) !! k singular values in descending order
        real(dp), intent(out) :: vh(:, :) !! k-by-n or n-by-n transposed right vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! Jacobi sweep limit, default 120
        logical, intent(in), optional :: full_matrices !! full square factors; default false
        complex(dp), allocatable :: uc(:, :), vhc(:, :)
        integer :: m, n, k, nu, nv
        logical :: full

        u = 0.0_dp
        s = 0.0_dp
        vh = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a, 1)
        n = size(a, 2)
        k = min(m, n)
        if (m < 1 .or. n < 1) return
        full = .false.
        if (present(full_matrices)) full = full_matrices
        nu = k
        nv = k
        if (full) then
            nu = m
            nv = n
        end if
        if (size(u,1) /= m .or. size(u,2) /= nu .or. size(s) /= k .or. &
            size(vh,1) /= nv .or. size(vh,2) /= n) return
        allocate(uc(m,nu),vhc(nv,n))
        call svd_complex(cmplx(a, 0.0_dp, dp), uc, s, vhc, status, &
            max_sweeps=max_sweeps, full_matrices=full)
        if (status /= linalg_status_success) return
        ! Jacobi rotations and orthogonal completion preserve real arithmetic
        ! for real input. Reject unexpected imaginary contamination.
        if (max(maxval(abs(aimag(uc))),maxval(abs(aimag(vhc)))) > 1.0e-11_dp) then
            status = linalg_status_no_convergence
            s = 0.0_dp
            return
        end if
        u = real(uc, dp)
        vh = real(vhc, dp)
    end subroutine svd_real

    subroutine svd_complex(a, u, s, vh, status, max_sweeps, full_matrices)
        complex(dp), intent(in) :: a(:, :) !! finite m-by-n complex matrix, m,n >= 1
        complex(dp), intent(out) :: u(:, :) !! m-by-k or m-by-m unitary-column factor
        real(dp), intent(out) :: s(:) !! k singular values in descending order
        complex(dp), intent(out) :: vh(:, :) !! k-by-n or n-by-n conjugate right vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! Jacobi sweep limit, default 120
        logical, intent(in), optional :: full_matrices !! full square factors; default false
        complex(dp), allocatable :: ut(:, :), vt(:, :), full_v(:, :)
        integer :: m, n, k, limit, nu, nv
        logical :: full

        u = cmplx(0.0_dp, 0.0_dp, dp)
        s = 0.0_dp
        vh = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        m = size(a, 1)
        n = size(a, 2)
        k = min(m, n)
        if (m < 1 .or. n < 1) return
        full = .false.
        if (present(full_matrices)) full = full_matrices
        nu = k
        nv = k
        if (full) then
            nu = m
            nv = n
        end if
        if (size(u,1) /= m .or. size(u,2) /= nu .or. size(s) /= k .or. &
            size(vh,1) /= nv .or. size(vh,2) /= n) return
        limit = 120
        if (present(max_sweeps)) limit = max_sweeps
        if (limit < 1) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(a)) return
        if (m >= n) then
            allocate(vt(n,n))
            call tall_svd(a, u(:,1:k), s, vt, status, limit)
            if (status == linalg_status_success) then
                vh(1:k,:) = transpose(conjg(vt))
                if (full .and. m > k) call complete_unitary(u, k, status)
            end if
        else
            ! A^H = Ut S Vt^H => A = Vt S Ut^H.
            allocate(ut(n,m))
            call tall_svd(transpose(conjg(a)), ut, s, u(:,1:k), status, limit)
            if (status == linalg_status_success) then
                if (full) then
                    allocate(full_v(n,n))
                    full_v = cmplx(0.0_dp,0.0_dp,dp)
                    full_v(:,1:k) = ut
                    call complete_unitary(full_v,k,status)
                    if (status == linalg_status_success) vh = transpose(conjg(full_v))
                else
                    vh(1:k,:) = transpose(conjg(ut))
                end if
            end if
        end if
        if (status /= linalg_status_success) then
            u = cmplx(0.0_dp, 0.0_dp, dp)
            s = 0.0_dp
            vh = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end subroutine svd_complex

    subroutine complete_unitary(q, k, status)
        complex(dp), intent(inout) :: q(:, :) !! n-by-n matrix with k orthonormal columns
        integer, intent(in) :: k !! number of known columns, 1 <= k <= n
        integer, intent(out) :: status !! success or no_convergence
        complex(dp), allocatable :: trial(:)
        real(dp), allocatable :: residual_diagonal(:)
        real(dp) :: norm_trial
        integer :: n, j, i, pass, pivot

        n = size(q,1)
        status = linalg_status_no_convergence
        if (size(q,2) /= n .or. k < 1 .or. k > n) return
        if (k == n) then
            status = linalg_status_success
            return
        end if
        allocate(trial(n), residual_diagonal(n))
        ! Greedily choose the coordinate direction farthest from the current
        ! span. Two modified Gram-Schmidt passes suppress loss of orthogonality.
        residual_diagonal = max(0.0_dp, 1.0_dp - sum(abs(q(:,1:k))**2,dim=2))
        do j = k+1, n
            pivot = maxloc(residual_diagonal,dim=1)
            trial = cmplx(0.0_dp,0.0_dp,dp)
            trial(pivot) = cmplx(1.0_dp,0.0_dp,dp)
            do pass = 1, 2
                do i = 1, j-1
                    trial = trial - dot_product(q(:,i),trial)*q(:,i)
                end do
            end do
            norm_trial = sqrt(sum(abs(trial)**2))
            if (norm_trial <= 32.0_dp*epsilon(1.0_dp)) return
            q(:,j) = trial/norm_trial
            residual_diagonal = max(0.0_dp,residual_diagonal - abs(q(:,j))**2)
        end do
        if (.not. finite_complex(q)) return
        status = linalg_status_success
    end subroutine complete_unitary

    subroutine svdvals_real(a, s, status, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! finite m-by-n real matrix, m,n >= 1
        real(dp), intent(out) :: s(:) !! min(m,n) singular values, descending
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! Jacobi sweep limit, default 120
        integer :: m, n, k, limit

        s = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (m < 1 .or. n < 1 .or. size(s) /= k) return
        limit = 120
        if (present(max_sweeps)) limit = max_sweeps
        if (limit < 1) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        if (m >= n) then
            call tall_svdvals_real(a, s, status, limit)
        else
            call tall_svdvals_real(transpose(a), s, status, limit)
        end if
    end subroutine svdvals_real

    subroutine svdvals_complex(a, s, status, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! finite m-by-n complex matrix, m,n >= 1
        real(dp), intent(out) :: s(:) !! min(m,n) singular values, descending
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! Jacobi sweep limit, default 120
        integer :: m, n, k, limit

        s = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (m < 1 .or. n < 1 .or. size(s) /= k) return
        limit = 120
        if (present(max_sweeps)) limit = max_sweeps
        if (limit < 1) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(a)) return
        if (m >= n) then
            call tall_svdvals_complex(a, s, status, limit)
        else
            call tall_svdvals_complex(transpose(conjg(a)), s, status, limit)
        end if
    end subroutine svdvals_complex

    ! Singular-value-only Jacobi sweeps. Unlike tall_svd, these routines
    ! neither accumulate right rotations nor construct left singular vectors.
    ! The real path also avoids promoting a real matrix to complex arithmetic.
    subroutine tall_svdvals_real(a, s, status, limit)
        real(dp), intent(in) :: a(:, :) !! tall real matrix, rows >= columns
        real(dp), intent(out) :: s(:) !! singular values in descending order
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in) :: limit !! maximum sweeps, >= 1
        real(dp), allocatable :: w(:, :), column(:)
        real(dp) :: scale, na, nb, gamma, g, phase, angle, c, sn
        real(dp) :: floor_sigma
        integer :: n, i, j, sweep
        logical :: changed, converged

        status = linalg_status_no_convergence
        n = size(a,2)
        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            s = 0.0_dp
            status = linalg_status_success
            return
        end if
        allocate(w(size(a,1),n),column(size(a,1)))
        w = a/scale
        converged = .false.
        do sweep = 1, limit
            changed = .false.
            do i = 1, n-1
                do j = i+1, n
                    na = sum(w(:,i)**2)
                    nb = sum(w(:,j)**2)
                    if (na == 0.0_dp .or. nb == 0.0_dp) cycle
                    gamma = dot_product(w(:,i),w(:,j))
                    g = abs(gamma)
                    if (g <= 32.0_dp*epsilon(1.0_dp)*sqrt(na)*sqrt(nb)) cycle
                    phase = sign(1.0_dp,gamma)
                    if (nb >= na) then
                        angle = 0.5_dp*atan2(2.0_dp*g,nb-na)
                    else
                        angle = -0.5_dp*atan2(2.0_dp*g,na-nb)
                    end if
                    if (abs(angle) <= epsilon(1.0_dp)) cycle
                    c = cos(angle)
                    sn = sin(angle)
                    column = c*w(:,i) - sn*phase*w(:,j)
                    w(:,j) = sn*w(:,i) + c*phase*w(:,j)
                    w(:,i) = column
                    changed = .true.
                end do
            end do
            if (.not. changed) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) return
        do j = 1, n
            s(j) = sqrt(sum(w(:,j)**2))
        end do
        call sort_singular_values(s)
        floor_sigma = 8.0_dp*epsilon(1.0_dp)*s(1)
        where (s <= floor_sigma) s = 0.0_dp
        s = s*scale
        if (.not. all(ieee_is_finite(s))) then
            s = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine tall_svdvals_real

    subroutine tall_svdvals_complex(a, s, status, limit)
        complex(dp), intent(in) :: a(:, :) !! tall complex matrix, rows >= columns
        real(dp), intent(out) :: s(:) !! singular values in descending order
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in) :: limit !! maximum sweeps, >= 1
        complex(dp), allocatable :: w(:, :), column(:)
        complex(dp) :: gamma, phase
        real(dp) :: scale, na, nb, g, angle, c, sn, floor_sigma
        integer :: n, i, j, sweep
        logical :: changed, converged

        status = linalg_status_no_convergence
        n = size(a,2)
        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            s = 0.0_dp
            status = linalg_status_success
            return
        end if
        allocate(w(size(a,1),n),column(size(a,1)))
        w = a/scale
        converged = .false.
        do sweep = 1, limit
            changed = .false.
            do i = 1, n-1
                do j = i+1, n
                    na = sum(abs(w(:,i))**2)
                    nb = sum(abs(w(:,j))**2)
                    if (na == 0.0_dp .or. nb == 0.0_dp) cycle
                    gamma = dot_product(w(:,i),w(:,j))
                    g = abs(gamma)
                    if (g <= 32.0_dp*epsilon(1.0_dp)*sqrt(na)*sqrt(nb)) cycle
                    phase = conjg(gamma)/g
                    if (nb >= na) then
                        angle = 0.5_dp*atan2(2.0_dp*g,nb-na)
                    else
                        angle = -0.5_dp*atan2(2.0_dp*g,na-nb)
                    end if
                    if (abs(angle) <= epsilon(1.0_dp)) cycle
                    c = cos(angle)
                    sn = sin(angle)
                    column = c*w(:,i) - sn*phase*w(:,j)
                    w(:,j) = sn*w(:,i) + c*phase*w(:,j)
                    w(:,i) = column
                    changed = .true.
                end do
            end do
            if (.not. changed) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) return
        do j = 1, n
            s(j) = sqrt(sum(abs(w(:,j))**2))
        end do
        call sort_singular_values(s)
        floor_sigma = 8.0_dp*epsilon(1.0_dp)*s(1)
        where (s <= floor_sigma) s = 0.0_dp
        s = s*scale
        if (.not. all(ieee_is_finite(s))) then
            s = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine tall_svdvals_complex

    subroutine sort_singular_values(s)
        real(dp), intent(inout) :: s(:) !! singular values, sorted in place
        real(dp) :: swap_sigma
        integer :: j, best

        do j = 1, size(s)-1
            best = j-1+maxloc(s(j:),dim=1)
            if (best == j) cycle
            swap_sigma = s(j)
            s(j) = s(best)
            s(best) = swap_sigma
        end do
    end subroutine sort_singular_values

    subroutine cond2_real(a, condition, status, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! finite rectangular matrix, m,n >= 1
        real(dp), intent(out) :: condition !! 2-norm spectral condition number, +Inf if singular
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! Jacobi sweep limit, default 120
        real(dp), allocatable :: s(:)
        integer :: m, n

        condition = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (m < 1 .or. n < 1) return
        allocate(s(min(m,n)))
        call svdvals_real(a,s,status,max_sweeps)
        if (status /= linalg_status_success) return
        call singular_values_condition(s,condition)
    end subroutine cond2_real

    subroutine cond2_complex(a, condition, status, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! finite rectangular matrix, m,n >= 1
        real(dp), intent(out) :: condition !! 2-norm spectral condition number, +Inf if singular
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! Jacobi sweep limit, default 120
        real(dp), allocatable :: s(:)
        integer :: m, n

        condition = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (m < 1 .or. n < 1) return
        allocate(s(min(m,n)))
        call svdvals_complex(a,s,status,max_sweeps)
        if (status /= linalg_status_success) return
        call singular_values_condition(s,condition)
    end subroutine cond2_complex

    subroutine singular_values_condition(s, condition)
        real(dp), intent(in) :: s(:) !! descending nonnegative singular values, nonempty
        real(dp), intent(out) :: condition !! finite condition or positive IEEE infinity

        if (s(size(s)) == 0.0_dp) then
            condition = ieee_value(1.0_dp,ieee_positive_inf)
        else if (s(size(s)) <= s(1)/huge(1.0_dp)) then
            condition = ieee_value(1.0_dp,ieee_positive_inf)
        else
            condition = s(1)/s(size(s))
        end if
    end subroutine singular_values_condition

    subroutine matrix_rank_real(a, rank, status, tol, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! finite rectangular matrix, m,n >= 1
        integer, intent(out) :: rank !! numerical rank, 0 to min(m,n)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: tol !! absolute singular-value threshold, >= 0
        integer, intent(in), optional :: max_sweeps !! Jacobi sweep limit, default 120
        real(dp), allocatable :: s(:)
        real(dp) :: cutoff
        integer :: m, n

        rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (m < 1 .or. n < 1) return
        if (present(tol)) then
            if (.not. ieee_is_finite(tol) .or. tol < 0.0_dp) return
        end if
        allocate(s(min(m,n)))
        call svdvals_real(a,s,status,max_sweeps)
        if (status /= linalg_status_success) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)*s(1)
        if (present(tol)) cutoff = tol
        rank = count(s > cutoff)
    end subroutine matrix_rank_real

    subroutine matrix_rank_complex(a, rank, status, tol, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! finite rectangular matrix, m,n >= 1
        integer, intent(out) :: rank !! numerical rank, 0 to min(m,n)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: tol !! absolute singular-value threshold, >= 0
        integer, intent(in), optional :: max_sweeps !! Jacobi sweep limit, default 120
        real(dp), allocatable :: s(:)
        real(dp) :: cutoff
        integer :: m, n

        rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (m < 1 .or. n < 1) return
        if (present(tol)) then
            if (.not. ieee_is_finite(tol) .or. tol < 0.0_dp) return
        end if
        allocate(s(min(m,n)))
        call svdvals_complex(a,s,status,max_sweeps)
        if (status /= linalg_status_success) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)*s(1)
        if (present(tol)) cutoff = tol
        rank = count(s > cutoff)
    end subroutine matrix_rank_complex

    subroutine diagsvd(s, m, n, d, status)
        real(dp), intent(in) :: s(:) !! min(m,n) values to place on the diagonal
        integer, intent(in) :: m !! output row count, >= 0
        integer, intent(in) :: n !! output column count, >= 0
        real(dp), intent(out) :: d(:, :) !! m-by-n rectangular diagonal matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        integer :: i

        d = 0.0_dp
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0) return
        if (size(s) /= min(m,n)) return
        if (size(d,1) /= m .or. size(d,2) /= n) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(s))) return
        do i = 1, size(s)
            d(i,i) = s(i)
        end do
        status = linalg_status_success
    end subroutine diagsvd

    subroutine tall_svd(a, u, s, v, status, limit)
        complex(dp), intent(in) :: a(:, :) !! finite input matrix with m >= n
        complex(dp), intent(out) :: u(:, :) !! m-by-n orthonormal columns
        real(dp), intent(out) :: s(:) !! descending singular values
        complex(dp), intent(out) :: v(:, :) !! n-by-n right singular vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in) :: limit !! maximum number of sweeps, >= 1
        complex(dp), allocatable :: w(:, :), column(:), rcolumn(:), best_column(:)
        complex(dp) :: gamma, phase
        real(dp) :: scale, na, nb, g, angle, c, sn, floor_sigma, biggest
        real(dp) :: norm_u, candidate, best_norm, swap_sigma
        integer :: m, n, i, j, sweep, best, pivot, pass, idx
        logical :: changed, converged

        m = size(a, 1)
        n = size(a, 2)
        u = cmplx(0.0_dp, 0.0_dp, dp)
        v = cmplx(0.0_dp, 0.0_dp, dp)
        s = 0.0_dp
        status = linalg_status_no_convergence
        do j = 1, n
            v(j,j) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            do j = 1, n
                u(j,j) = cmplx(1.0_dp, 0.0_dp, dp)
            end do
            status = linalg_status_success
            return
        end if
        allocate(w(m,n),column(m),rcolumn(n),best_column(m))
        w = a / scale
        converged = .false.
        do sweep = 1, limit
            changed = .false.
            do i = 1, n-1
                do j = i+1, n
                    na = sum(abs(w(:,i))**2)
                    nb = sum(abs(w(:,j))**2)
                    if (na == 0.0_dp .or. nb == 0.0_dp) cycle
                    gamma = dot_product(w(:,i),w(:,j))
                    g = abs(gamma)
                    if (g <= 32.0_dp*epsilon(1.0_dp)*sqrt(na)*sqrt(nb)) cycle
                    phase = conjg(gamma)/g
                    if (nb >= na) then
                        angle = 0.5_dp*atan2(2.0_dp*g,nb-na)
                    else
                        angle = -0.5_dp*atan2(2.0_dp*g,na-nb)
                    end if
                    if (abs(angle) <= epsilon(1.0_dp)) cycle
                    c = cos(angle)
                    sn = sin(angle)
                    column = c*w(:,i) - sn*phase*w(:,j)
                    w(:,j) = sn*w(:,i) + c*phase*w(:,j)
                    w(:,i) = column
                    rcolumn = c*v(:,i) - sn*phase*v(:,j)
                    v(:,j) = sn*v(:,i) + c*phase*v(:,j)
                    v(:,i) = rcolumn
                    changed = .true.
                end do
            end do
            if (.not. changed) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) return
        do j = 1, n
            s(j) = sqrt(sum(abs(w(:,j))**2))
        end do
        do j = 1, n-1
            best = j - 1 + maxloc(s(j:n),dim=1)
            if (best == j) cycle
            swap_sigma = s(j)
            s(j) = s(best)
            s(best) = swap_sigma
            column = w(:,j)
            w(:,j) = w(:,best)
            w(:,best) = column
            rcolumn = v(:,j)
            v(:,j) = v(:,best)
            v(:,best) = rcolumn
        end do
        biggest = maxval(s)
        floor_sigma = 8.0_dp*epsilon(1.0_dp)*biggest
        do j = 1, n
            if (s(j) > floor_sigma) then
                column = w(:,j)/s(j)
            else
                ! Complete an orthonormal basis without dividing by a
                ! near-zero singular value.
                s(j) = 0.0_dp
                best_norm = -1.0_dp
                do pivot = 1, m
                    column = cmplx(0.0_dp, 0.0_dp, dp)
                    column(pivot) = cmplx(1.0_dp, 0.0_dp, dp)
                    do pass = 1, 2
                        do idx = 1, j-1
                            column = column-dot_product(u(:,idx),column)*u(:,idx)
                        end do
                    end do
                    candidate = sqrt(sum(abs(column)**2))
                    if (candidate > best_norm) then
                        best_norm = candidate
                        best_column = column
                    end if
                end do
                column = best_column
            end if
            do pass = 1, 2
                do idx = 1, j-1
                    column = column-dot_product(u(:,idx),column)*u(:,idx)
                end do
            end do
            norm_u = sqrt(sum(abs(column)**2))
            if (norm_u <= 32.0_dp*epsilon(1.0_dp)) return
            u(:,j) = column/norm_u
        end do
        s = s*scale
        if (.not. all(ieee_is_finite(s)) .or. .not. finite_complex(u) .or. &
            .not. finite_complex(v)) return
        status = linalg_status_success
    end subroutine tall_svd

    subroutine svd_lstsq_real(a, b, x, status, rank, rcond)
        real(dp), intent(in) :: a(:, :) !! m-by-n coefficient matrix, m,n >= 1
        real(dp), intent(in) :: b(:, :) !! m-by-nrhs right-hand sides
        real(dp), intent(out) :: x(:, :) !! n-by-nrhs minimum-norm solutions
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(out), optional :: rank !! estimated numerical rank from singular values
        real(dp), intent(in), optional :: rcond !! nonnegative relative cutoff; default max(m,n)*epsilon
        real(dp), allocatable :: u(:, :), vh(:, :), s(:), projection(:, :)
        real(dp) :: cutoff, threshold
        integer :: m, n, k, j, numerical_rank

        x = 0.0_dp
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(b,1) /= m .or. size(x,1) /= n .or. size(x,2) /= size(b,2)) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        allocate(u(m,k),vh(k,n),s(k),projection(k,size(b,2)))
        call svd_real(a,u,s,vh,status)
        if (status /= linalg_status_success) return
        threshold = cutoff*s(1)
        numerical_rank = count(s > threshold)
        if (present(rank)) rank = numerical_rank
        projection = matmul(transpose(u),b)
        do j = 1, k
            if (s(j) > threshold) then
                projection(j,:) = (projection(j,:)/s(j))
            else
                projection(j,:) = 0.0_dp
            end if
        end do
        x = matmul(transpose(vh),projection)
        if (.not. all(ieee_is_finite(x))) then
            x = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine svd_lstsq_real

    subroutine svd_lstsq_complex(a, b, x, status, rank, rcond)
        complex(dp), intent(in) :: a(:, :) !! m-by-n coefficient matrix, m,n >= 1
        complex(dp), intent(in) :: b(:, :) !! m-by-nrhs right-hand sides
        complex(dp), intent(out) :: x(:, :) !! n-by-nrhs minimum-norm solutions
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(out), optional :: rank !! estimated numerical rank from singular values
        real(dp), intent(in), optional :: rcond !! nonnegative relative cutoff; default max(m,n)*epsilon
        complex(dp), allocatable :: u(:, :), vh(:, :), projection(:, :)
        real(dp), allocatable :: s(:)
        real(dp) :: cutoff, threshold
        integer :: m, n, k, j, numerical_rank

        x = cmplx(0.0_dp,0.0_dp,dp)
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(b,1) /= m .or. size(x,1) /= n .or. size(x,2) /= size(b,2)) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(a) .or. .not. finite_complex(b)) return
        allocate(u(m,k),vh(k,n),s(k),projection(k,size(b,2)))
        call svd_complex(a,u,s,vh,status)
        if (status /= linalg_status_success) return
        threshold = cutoff*s(1)
        numerical_rank = count(s > threshold)
        if (present(rank)) rank = numerical_rank
        projection = matmul(transpose(conjg(u)),b)
        do j = 1, k
            if (s(j) > threshold) then
                projection(j,:) = projection(j,:)/s(j)
            else
                projection(j,:) = cmplx(0.0_dp,0.0_dp,dp)
            end if
        end do
        x = matmul(transpose(conjg(vh)),projection)
        if (.not. finite_complex(x)) then
            x = cmplx(0.0_dp,0.0_dp,dp)
            status = linalg_status_no_convergence
        end if
    end subroutine svd_lstsq_complex

    subroutine pinv_real(a, ap, status, rank, rcond)
        real(dp), intent(in) :: a(:, :) !! m-by-n finite input matrix, m,n >= 1
        real(dp), intent(out) :: ap(:, :) !! n-by-m Moore-Penrose pseudoinverse
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! rank after singular-value cutoff
        real(dp), intent(in), optional :: rcond !! relative cutoff, default max(m,n)*epsilon
        real(dp), allocatable :: u(:, :), vh(:, :), s(:), vr(:, :)
        real(dp) :: cutoff
        integer :: m, n, k, j

        ap = 0.0_dp
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(ap,1) /= n .or. size(ap,2) /= m) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        allocate(u(m,k),vh(k,n),s(k),vr(n,k))
        call svd_real(a,u,s,vh,status)
        if (status /= linalg_status_success) return
        vr = transpose(vh)
        do j = 1, k
            if (s(j) > cutoff*s(1)) then
                vr(:,j) = vr(:,j)/s(j)
                if (present(rank)) rank = rank + 1
            else
                vr(:,j) = 0.0_dp
            end if
        end do
        ap = matmul(vr,transpose(u))
        if (.not. all(ieee_is_finite(ap))) then
            ap = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine pinv_real

    subroutine pinv_complex(a, ap, status, rank, rcond)
        complex(dp), intent(in) :: a(:, :) !! m-by-n finite input matrix, m,n >= 1
        complex(dp), intent(out) :: ap(:, :) !! n-by-m Moore-Penrose pseudoinverse
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(out), optional :: rank !! rank after singular-value cutoff
        real(dp), intent(in), optional :: rcond !! relative cutoff, default max(m,n)*epsilon
        complex(dp), allocatable :: u(:, :), vh(:, :), vr(:, :)
        real(dp), allocatable :: s(:)
        real(dp) :: cutoff
        integer :: m, n, k, j

        ap = cmplx(0.0_dp,0.0_dp,dp)
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(ap,1) /= n .or. size(ap,2) /= m) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        allocate(u(m,k),vh(k,n),s(k),vr(n,k))
        call svd_complex(a,u,s,vh,status)
        if (status /= linalg_status_success) return
        vr = transpose(conjg(vh))
        do j = 1, k
            if (s(j) > cutoff*s(1)) then
                vr(:,j) = vr(:,j)/s(j)
                if (present(rank)) rank = rank + 1
            else
                vr(:,j) = cmplx(0.0_dp,0.0_dp,dp)
            end if
        end do
        ap = matmul(vr,transpose(conjg(u)))
        if (.not. finite_complex(ap)) then
            ap = cmplx(0.0_dp,0.0_dp,dp)
            status = linalg_status_no_convergence
        end if
    end subroutine pinv_complex

end module scifort_linalg_svd
