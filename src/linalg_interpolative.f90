! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Deterministic interpolative decomposition (ID) from column-pivoted QR.
! SciPy analogues: scipy.linalg.interpolative.interp_decomp,
! reconstruct_skel_matrix, reconstruct_interp_matrix,
! reconstruct_matrix_from_id, and id_to_svd.
!
! Matrix A approximately equals B*P where B consists of columns of A
! and P has a k-by-k identity submatrix. Column indices are one-based.
module scifort_linalg_interpolative
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_qr, only : qr_factor, qr
    use scifort_linalg_svd, only : svd_lstsq, svd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: interp_decomp, reconstruct_skel_matrix
    public :: reconstruct_interp_matrix, reconstruct_matrix_from_id, id_to_svd

    interface interp_decomp
        module procedure interp_real_rank
        module procedure interp_complex_rank
        module procedure interp_real_tol
        module procedure interp_complex_tol
    end interface interp_decomp
    interface reconstruct_skel_matrix
        module procedure skeleton_real, skeleton_complex
    end interface reconstruct_skel_matrix
    interface reconstruct_interp_matrix
        module procedure interp_matrix_real, interp_matrix_complex
    end interface reconstruct_interp_matrix
    interface reconstruct_matrix_from_id
        module procedure reconstruct_real, reconstruct_complex
    end interface reconstruct_matrix_from_id
    interface id_to_svd
        module procedure id_to_svd_real, id_to_svd_complex
    end interface id_to_svd

contains

    pure logical function finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! complex input values
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex

    pure logical function valid_indices(idx) result(ok)
        integer, intent(in) :: idx(:) !! one-based permutation of columns 1:n
        logical :: seen(size(idx))
        integer :: i, n
        n = size(idx)
        ok = .false.
        if (any(idx < 1) .or. any(idx > n)) return
        seen = .false.
        do i = 1, n
            if (seen(idx(i))) return
            seen(idx(i)) = .true.
        end do
        ok = .true.
    end function valid_indices

    ! Frobenius-norm tail of the triangular R, scaled before squaring.
    ! Selecting a rank from ||R(k+1:,:)||_F gives a deterministic
    ! Frobenius error certificate for the interpolative approximation.
    integer function rank_from_qr_real(af, eps, scale) result(k)
        real(dp), intent(in) :: af(:, :) !! Householder QR with upper-triangular R
        real(dp), intent(in) :: eps !! requested relative Frobenius error, 0 <= eps < 1
        real(dp), intent(in) :: scale !! positive scale of original input
        real(dp), allocatable :: row_energy(:)
        real(dp) :: total, tail
        integer :: i, r, n
        r = min(size(af,1),size(af,2))
        n = size(af,2)
        allocate(row_energy(r))
        do i = 1, r
            row_energy(i) = sum((abs(af(i,i:n))/scale)**2)
        end do
        total = sum(row_energy)
        tail = total
        k = 0
        if (total == 0.0_dp) return
        do i = 1, r
            if (tail <= eps**2 * total) return
            k = i
            tail = sum(row_energy(i+1:r))
        end do
    end function rank_from_qr_real

    integer function rank_from_qr_complex(af, eps, scale) result(k)
        complex(dp), intent(in) :: af(:, :) !! complex Householder QR
        real(dp), intent(in) :: eps !! requested relative Frobenius error
        real(dp), intent(in) :: scale !! positive input scale
        real(dp), allocatable :: row_energy(:)
        real(dp) :: total, tail
        integer :: i, r, n
        r = min(size(af,1),size(af,2))
        n = size(af,2)
        allocate(row_energy(r))
        do i = 1, r
            row_energy(i) = sum((abs(af(i,i:n))/scale)**2)
        end do
        total = sum(row_energy)
        tail = total
        k = 0
        if (total == 0.0_dp) return
        do i = 1, r
            if (tail <= eps**2 * total) return
            k = i
            tail = sum(row_energy(i+1:r))
        end do
    end function rank_from_qr_complex

    subroutine interp_real_rank(a, k, idx, proj, status)
        real(dp), intent(in) :: a(:, :) !! m-by-n real matrix
        integer, intent(in) :: k !! requested rank, 0 <= k <= min(m,n)
        integer, allocatable, intent(out) :: idx(:) !! one-based pivoted column indices
        real(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), allocatable :: af(:, :), tau(:)
        integer :: m, n, i
        allocate(idx(0),proj(0,0))
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (k < 0 .or. k > min(m,n)) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        deallocate(idx,proj)
        allocate(idx(n),proj(k,n-k))
        idx = [(i, i=1,n)]
        proj = 0.0_dp
        if (m == 0 .or. n == 0 .or. k == 0) then
            status = linalg_status_success
            return
        end if
        allocate(af(m,n),tau(min(m,n)))
        call qr_factor(a,af,tau,status,idx)
        if (status /= linalg_status_success) return
        call projection_real(a,af,k,idx,proj,status)
        if (status /= linalg_status_success) proj = 0.0_dp
    end subroutine interp_real_rank

    subroutine interp_complex_rank(a, k, idx, proj, status)
        complex(dp), intent(in) :: a(:, :) !! m-by-n complex matrix
        integer, intent(in) :: k !! requested rank, 0 <= k <= min(m,n)
        integer, allocatable, intent(out) :: idx(:) !! one-based pivoted column indices
        complex(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), allocatable :: af(:, :), tau(:)
        integer :: m, n, i
        allocate(idx(0),proj(0,0))
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (k < 0 .or. k > min(m,n)) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(a)) return
        deallocate(idx,proj)
        allocate(idx(n),proj(k,n-k))
        idx = [(i, i=1,n)]
        proj = cmplx(0.0_dp,0.0_dp,dp)
        if (m == 0 .or. n == 0 .or. k == 0) then
            status = linalg_status_success
            return
        end if
        allocate(af(m,n),tau(min(m,n)))
        call qr_factor(a,af,tau,status,idx)
        if (status /= linalg_status_success) return
        call projection_complex(a,af,k,idx,proj,status)
        if (status /= linalg_status_success) proj = cmplx(0.0_dp,0.0_dp,dp)
    end subroutine interp_complex_rank

    subroutine interp_real_tol(a, eps, k, idx, proj, status)
        real(dp), intent(in) :: a(:, :) !! m-by-n real matrix
        real(dp), intent(in) :: eps !! relative Frobenius error, 0 <= eps < 1
        integer, intent(out) :: k !! selected rank
        integer, allocatable, intent(out) :: idx(:) !! one-based pivoted column indices
        real(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), allocatable :: af(:, :), tau(:)
        real(dp) :: scale
        integer :: m, n, i
        k = 0
        allocate(idx(0),proj(0,0))
        status = linalg_status_invalid_shape
        if (.not. ieee_is_finite(eps)) return
        if (eps < 0.0_dp .or. eps >= 1.0_dp) return
        m = size(a,1)
        n = size(a,2)
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        deallocate(idx,proj)
        allocate(idx(n))
        idx = [(i, i=1,n)]
        if (m == 0 .or. n == 0) then
            allocate(proj(0,n))
            status = linalg_status_success
            return
        end if
        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            allocate(proj(0,n))
            status = linalg_status_success
            return
        end if
        allocate(af(m,n),tau(min(m,n)))
        call qr_factor(a,af,tau,status,idx)
        if (status /= linalg_status_success) then
            allocate(proj(0,0))
            return
        end if
        k = rank_from_qr_real(af,eps,scale)
        allocate(proj(k,n-k))
        proj = 0.0_dp
        if (k == 0) return
        call projection_real(a,af,k,idx,proj,status)
        if (status /= linalg_status_success) then
            k = 0
            proj = 0.0_dp
        end if
    end subroutine interp_real_tol

    subroutine interp_complex_tol(a, eps, k, idx, proj, status)
        complex(dp), intent(in) :: a(:, :) !! m-by-n complex matrix
        real(dp), intent(in) :: eps !! relative Frobenius error, 0 <= eps < 1
        integer, intent(out) :: k !! selected rank
        integer, allocatable, intent(out) :: idx(:) !! one-based pivoted column indices
        complex(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), allocatable :: af(:, :), tau(:)
        real(dp) :: scale
        integer :: m, n, i
        k = 0
        allocate(idx(0),proj(0,0))
        status = linalg_status_invalid_shape
        if (.not. ieee_is_finite(eps)) return
        if (eps < 0.0_dp .or. eps >= 1.0_dp) return
        m = size(a,1)
        n = size(a,2)
        status = linalg_status_no_convergence
        if (.not. finite_complex(a)) return
        deallocate(idx,proj)
        allocate(idx(n))
        idx = [(i, i=1,n)]
        if (m == 0 .or. n == 0) then
            allocate(proj(0,n))
            status = linalg_status_success
            return
        end if
        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            allocate(proj(0,n))
            status = linalg_status_success
            return
        end if
        allocate(af(m,n),tau(min(m,n)))
        call qr_factor(a,af,tau,status,idx)
        if (status /= linalg_status_success) then
            allocate(proj(0,0))
            return
        end if
        k = rank_from_qr_complex(af,eps,scale)
        allocate(proj(k,n-k))
        proj = cmplx(0.0_dp,0.0_dp,dp)
        if (k == 0) return
        call projection_complex(a,af,k,idx,proj,status)
        if (status /= linalg_status_success) then
            k = 0
            proj = cmplx(0.0_dp,0.0_dp,dp)
        end if
    end subroutine interp_complex_tol

    subroutine projection_real(a, af, k, idx, proj, status)
        real(dp), intent(in) :: a(:, :) !! original m-by-n matrix
        real(dp), intent(in) :: af(:, :) !! pivoted QR factors of a
        integer, intent(in) :: k !! selected rank
        integer, intent(in) :: idx(:) !! pivoted columns
        real(dp), intent(out) :: proj(:, :) !! interpolation coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp) :: cutoff, rhs, scale
        integer :: n, i, j
        n = size(a,2)
        status = linalg_status_success
        proj = 0.0_dp
        if (k == 0 .or. k == n) return
        scale = maxval(abs(af(1:k,1:k)))
        cutoff = real(max(size(a,1),n),dp) * epsilon(1.0_dp) * scale
        if (any(abs([(af(i,i),i=1,k)]) <= cutoff)) then
            ! A rank-deficient selected skeleton needs a minimum-norm solve.
            call svd_lstsq(a(:,idx(1:k)),a(:,idx(k+1:n)),proj,status)
            return
        end if
        do j = 1, n-k
            do i = k, 1, -1
                rhs = af(i,k+j)
                if (i < k) rhs = rhs-sum(af(i,i+1:k)*proj(i+1:k,j))
                proj(i,j) = rhs/af(i,i)
            end do
        end do
        if (.not. all(ieee_is_finite(proj))) status = linalg_status_no_convergence
    end subroutine projection_real

    subroutine projection_complex(a, af, k, idx, proj, status)
        complex(dp), intent(in) :: a(:, :) !! original m-by-n matrix
        complex(dp), intent(in) :: af(:, :) !! pivoted QR factors of a
        integer, intent(in) :: k !! selected rank
        integer, intent(in) :: idx(:) !! pivoted columns
        complex(dp), intent(out) :: proj(:, :) !! interpolation coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp) :: rhs
        real(dp) :: cutoff, scale
        integer :: n, i, j
        n = size(a,2)
        status = linalg_status_success
        proj = cmplx(0.0_dp,0.0_dp,dp)
        if (k == 0 .or. k == n) return
        scale = maxval(abs(af(1:k,1:k)))
        cutoff = real(max(size(a,1),n),dp) * epsilon(1.0_dp) * scale
        if (any(abs([(af(i,i),i=1,k)]) <= cutoff)) then
            call svd_lstsq(a(:,idx(1:k)),a(:,idx(k+1:n)),proj,status)
            return
        end if
        do j = 1, n-k
            do i = k, 1, -1
                rhs = af(i,k+j)
                ! Do not conjugate: T = inv(R11) R12, not inv(R11^H) R12.
                if (i < k) rhs = rhs-sum(af(i,i+1:k)*proj(i+1:k,j))
                proj(i,j) = rhs/af(i,i)
            end do
        end do
        if (.not. finite_complex(proj)) status = linalg_status_no_convergence
    end subroutine projection_complex

    subroutine skeleton_real(a, k, idx, b, status)
        real(dp), intent(in) :: a(:, :) !! original m-by-n real matrix
        integer, intent(in) :: k !! number of selected columns
        integer, intent(in) :: idx(:) !! one-based permutation, length n
        real(dp), intent(out) :: b(:, :) !! m-by-k skeleton column matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        b = 0.0_dp
        status = linalg_status_invalid_shape
        if (k < 0 .or. k > min(size(a,1),size(a,2))) return
        if (size(idx) /= size(a,2) .or. size(b,1) /= size(a,1) .or. size(b,2) /= k) return
        if (.not. valid_indices(idx)) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        b = a(:,idx(1:k))
        status = linalg_status_success
    end subroutine skeleton_real

    subroutine skeleton_complex(a, k, idx, b, status)
        complex(dp), intent(in) :: a(:, :) !! original m-by-n complex matrix
        integer, intent(in) :: k !! number of selected columns
        integer, intent(in) :: idx(:) !! one-based permutation, length n
        complex(dp), intent(out) :: b(:, :) !! m-by-k skeleton column matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        b = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (k < 0 .or. k > min(size(a,1),size(a,2))) return
        if (size(idx) /= size(a,2) .or. size(b,1) /= size(a,1) .or. size(b,2) /= k) return
        if (.not. valid_indices(idx)) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(a)) return
        b = a(:,idx(1:k))
        status = linalg_status_success
    end subroutine skeleton_complex

    subroutine interp_matrix_real(idx, proj, p, status)
        integer, intent(in) :: idx(:) !! one-based permutation of n columns
        real(dp), intent(in) :: proj(:, :) !! k-by-(n-k) coefficients
        real(dp), intent(out) :: p(:, :) !! k-by-n interpolation matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        integer :: i, k, n
        p = 0.0_dp
        status = linalg_status_invalid_shape
        k = size(proj,1)
        n = size(idx)
        if (k > n .or. size(proj,2) /= n-k) return
        if (size(p,1) /= k .or. size(p,2) /= n) return
        if (.not. valid_indices(idx)) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(proj))) return
        do i = 1, k
            p(i,idx(i)) = 1.0_dp
        end do
        do i = k+1, n
            p(:,idx(i)) = proj(:,i-k)
        end do
        status = linalg_status_success
    end subroutine interp_matrix_real

    subroutine interp_matrix_complex(idx, proj, p, status)
        integer, intent(in) :: idx(:) !! one-based permutation of n columns
        complex(dp), intent(in) :: proj(:, :) !! k-by-(n-k) coefficients
        complex(dp), intent(out) :: p(:, :) !! k-by-n interpolation matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        integer :: i, k, n
        p = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        k = size(proj,1)
        n = size(idx)
        if (k > n .or. size(proj,2) /= n-k) return
        if (size(p,1) /= k .or. size(p,2) /= n) return
        if (.not. valid_indices(idx)) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(proj)) return
        do i = 1, k
            p(i,idx(i)) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        do i = k+1, n
            p(:,idx(i)) = proj(:,i-k)
        end do
        status = linalg_status_success
    end subroutine interp_matrix_complex

    subroutine reconstruct_real(b, idx, proj, a, status)
        real(dp), intent(in) :: b(:, :) !! m-by-k skeleton matrix
        integer, intent(in) :: idx(:) !! one-based column permutation
        real(dp), intent(in) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        real(dp), intent(out) :: a(:, :) !! reconstructed m-by-n matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        integer :: i, k, n
        a = 0.0_dp
        status = linalg_status_invalid_shape
        k = size(b,2)
        n = size(idx)
        if (k > n .or. size(proj,1) /= k .or. size(proj,2) /= n-k) return
        if (size(a,1) /= size(b,1) .or. size(a,2) /= n) return
        if (.not. valid_indices(idx)) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(b)) .or. .not. all(ieee_is_finite(proj))) return
        do i = 1, k
            a(:,idx(i)) = b(:,i)
        end do
        do i = k+1, n
            a(:,idx(i)) = matmul(b,proj(:,i-k))
        end do
        if (all(ieee_is_finite(a))) status = linalg_status_success
    end subroutine reconstruct_real

    subroutine reconstruct_complex(b, idx, proj, a, status)
        complex(dp), intent(in) :: b(:, :) !! m-by-k skeleton matrix
        integer, intent(in) :: idx(:) !! one-based column permutation
        complex(dp), intent(in) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        complex(dp), intent(out) :: a(:, :) !! reconstructed m-by-n matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        integer :: i, k, n
        a = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        k = size(b,2)
        n = size(idx)
        if (k > n .or. size(proj,1) /= k .or. size(proj,2) /= n-k) return
        if (size(a,1) /= size(b,1) .or. size(a,2) /= n) return
        if (.not. valid_indices(idx)) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(b) .or. .not. finite_complex(proj)) return
        do i = 1, k
            a(:,idx(i)) = b(:,i)
        end do
        do i = k+1, n
            a(:,idx(i)) = matmul(b,proj(:,i-k))
        end do
        if (finite_complex(a)) status = linalg_status_success
    end subroutine reconstruct_complex

    ! The ID-to-SVD conversion uses QR on the small skeleton and SVD on
    ! the k-by-n core. It never forms the original m-by-n matrix.
    subroutine id_to_svd_real(b, idx, proj, u, s, v, status)
        real(dp), intent(in) :: b(:, :) !! m-by-k skeleton matrix
        integer, intent(in) :: idx(:) !! one-based column permutation
        real(dp), intent(in) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        real(dp), intent(out) :: u(:, :) !! m-by-k left singular vectors
        real(dp), intent(out) :: s(:) !! k singular values, descending
        real(dp), intent(out) :: v(:, :) !! n-by-k right singular vectors (not V^T)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), allocatable :: q(:, :), r(:, :), p(:, :), core(:, :), uu(:, :), vh(:, :)
        integer :: k, m, n
        u = 0.0_dp
        s = 0.0_dp
        v = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(b,1)
        k = size(b,2)
        n = size(idx)
        if (k > min(m,n) .or. size(proj,1) /= k .or. size(proj,2) /= n-k) return
        if (size(u,1) /= m .or. size(u,2) /= k .or. size(s) /= k) return
        if (size(v,1) /= n .or. size(v,2) /= k .or. .not. valid_indices(idx)) return
        status = linalg_status_success
        if (k == 0) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(b)) .or. .not. all(ieee_is_finite(proj))) return
        allocate(q(m,k),r(k,k),p(k,n),core(k,n),uu(k,k),vh(k,n))
        call qr(b,q,r,status)
        if (status /= linalg_status_success) return
        call interp_matrix_real(idx,proj,p,status)
        if (status /= linalg_status_success) return
        core = matmul(r,p)
        call svd(core,uu,s,vh,status)
        if (status /= linalg_status_success) return
        u = matmul(q,uu)
        v = transpose(vh)
    end subroutine id_to_svd_real

    subroutine id_to_svd_complex(b, idx, proj, u, s, v, status)
        complex(dp), intent(in) :: b(:, :) !! m-by-k skeleton matrix
        integer, intent(in) :: idx(:) !! one-based column permutation
        complex(dp), intent(in) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        complex(dp), intent(out) :: u(:, :) !! m-by-k left singular vectors
        real(dp), intent(out) :: s(:) !! k singular values, descending
        complex(dp), intent(out) :: v(:, :) !! n-by-k right singular vectors (not V^H)
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), allocatable :: q(:, :), r(:, :), p(:, :), core(:, :), uu(:, :), vh(:, :)
        integer :: k, m, n
        u = cmplx(0.0_dp,0.0_dp,dp)
        s = 0.0_dp
        v = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(b,1)
        k = size(b,2)
        n = size(idx)
        if (k > min(m,n) .or. size(proj,1) /= k .or. size(proj,2) /= n-k) return
        if (size(u,1) /= m .or. size(u,2) /= k .or. size(s) /= k) return
        if (size(v,1) /= n .or. size(v,2) /= k .or. .not. valid_indices(idx)) return
        status = linalg_status_success
        if (k == 0) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(b) .or. .not. finite_complex(proj)) return
        allocate(q(m,k),r(k,k),p(k,n),core(k,n),uu(k,k),vh(k,n))
        call qr(b,q,r,status)
        if (status /= linalg_status_success) return
        call interp_matrix_complex(idx,proj,p,status)
        if (status /= linalg_status_success) return
        core = matmul(r,p)
        call svd(core,uu,s,vh,status)
        if (status /= linalg_status_success) return
        u = matmul(q,uu)
        v = transpose(conjg(vh))
    end subroutine id_to_svd_complex
end module scifort_linalg_interpolative
