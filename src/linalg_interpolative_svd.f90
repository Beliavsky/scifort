! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Deterministic interpolative SVD, analogous to
! scipy.linalg.interpolative.svd(..., rand=False).
!
! The ID chooses actual columns A(:,idx(1:k)); its small-core SVD yields
! A ~ U * diag(s) * V^H with U(m,k), s(k), V(n,k). No full m-by-n
! reconstruction is required. Integer rank and real relative-tolerance
! arguments select the two generic interfaces.
module scifort_linalg_interpolative_svd
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_interpolative, only : interp_decomp, id_to_svd
    implicit none
    private

    public :: interpolative_svd

    interface interpolative_svd
        module procedure interpolative_svd_real_rank
        module procedure interpolative_svd_complex_rank
        module procedure interpolative_svd_real_tol
        module procedure interpolative_svd_complex_tol
    end interface interpolative_svd

contains

    subroutine interpolative_svd_real_rank(a, k, u, s, v, status)
        real(dp), intent(in) :: a(:, :) !! m-by-n real matrix
        integer, intent(in) :: k !! target rank, 0 <= k <= min(m,n)
        real(dp), intent(out) :: u(:, :) !! m-by-k orthonormal left vectors
        real(dp), intent(out) :: s(:) !! k singular values in decreasing order
        real(dp), intent(out) :: v(:, :) !! n-by-k orthonormal right vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), allocatable :: b(:, :), proj(:, :)
        integer, allocatable :: idx(:)
        integer :: m, n

        m = size(a,1)
        n = size(a,2)
        u = 0.0_dp
        s = 0.0_dp
        v = 0.0_dp
        status = linalg_status_invalid_shape
        if (k < 0 .or. k > min(m,n)) return
        if (size(u,1) /= m .or. size(u,2) /= k) return
        if (size(s) /= k .or. size(v,1) /= n .or. size(v,2) /= k) return

        call interp_decomp(a,k,idx,proj,status)
        if (status /= linalg_status_success .or. k == 0) return
        allocate(b(m,k))
        b = a(:,idx(1:k))
        call id_to_svd(b,idx,proj,u,s,v,status)
        if (status /= linalg_status_success) then
            u = 0.0_dp
            s = 0.0_dp
            v = 0.0_dp
        end if
    end subroutine interpolative_svd_real_rank

    subroutine interpolative_svd_complex_rank(a, k, u, s, v, status)
        complex(dp), intent(in) :: a(:, :) !! m-by-n complex matrix
        integer, intent(in) :: k !! target rank, 0 <= k <= min(m,n)
        complex(dp), intent(out) :: u(:, :) !! m-by-k orthonormal left vectors
        real(dp), intent(out) :: s(:) !! k nonnegative singular values
        complex(dp), intent(out) :: v(:, :) !! n-by-k orthonormal right vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), allocatable :: b(:, :), proj(:, :)
        integer, allocatable :: idx(:)
        integer :: m, n

        m = size(a,1)
        n = size(a,2)
        u = cmplx(0.0_dp,0.0_dp,dp)
        s = 0.0_dp
        v = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (k < 0 .or. k > min(m,n)) return
        if (size(u,1) /= m .or. size(u,2) /= k) return
        if (size(s) /= k .or. size(v,1) /= n .or. size(v,2) /= k) return

        call interp_decomp(a,k,idx,proj,status)
        if (status /= linalg_status_success .or. k == 0) return
        allocate(b(m,k))
        b = a(:,idx(1:k))
        call id_to_svd(b,idx,proj,u,s,v,status)
        if (status /= linalg_status_success) then
            u = cmplx(0.0_dp,0.0_dp,dp)
            s = 0.0_dp
            v = cmplx(0.0_dp,0.0_dp,dp)
        end if
    end subroutine interpolative_svd_complex_rank

    subroutine interpolative_svd_real_tol(a, eps, k, u, s, v, status)
        real(dp), intent(in) :: a(:, :) !! m-by-n real matrix
        real(dp), intent(in) :: eps !! relative Frobenius tolerance, 0 <= eps < 1
        integer, intent(out) :: k !! selected rank, 0..min(m,n)
        real(dp), allocatable, intent(out) :: u(:, :) !! allocated m-by-k factor
        real(dp), allocatable, intent(out) :: s(:) !! allocated k singular values
        real(dp), allocatable, intent(out) :: v(:, :) !! allocated n-by-k factor
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), allocatable :: b(:, :), proj(:, :)
        integer, allocatable :: idx(:)
        integer :: m, n

        m = size(a,1)
        n = size(a,2)
        k = 0
        allocate(u(0,0),s(0),v(0,0))
        call interp_decomp(a,eps,k,idx,proj,status)
        if (status /= linalg_status_success) then
            k = 0
            return
        end if
        deallocate(u,s,v)
        allocate(u(m,k),s(k),v(n,k))
        u = 0.0_dp
        s = 0.0_dp
        v = 0.0_dp
        if (k == 0) return
        allocate(b(m,k))
        b = a(:,idx(1:k))
        call id_to_svd(b,idx,proj,u,s,v,status)
        if (status /= linalg_status_success) then
            k = 0
            deallocate(u,s,v)
            allocate(u(0,0),s(0),v(0,0))
        end if
    end subroutine interpolative_svd_real_tol

    subroutine interpolative_svd_complex_tol(a, eps, k, u, s, v, status)
        complex(dp), intent(in) :: a(:, :) !! m-by-n complex matrix
        real(dp), intent(in) :: eps !! relative Frobenius tolerance, 0 <= eps < 1
        integer, intent(out) :: k !! selected rank, 0..min(m,n)
        complex(dp), allocatable, intent(out) :: u(:, :) !! allocated m-by-k factor
        real(dp), allocatable, intent(out) :: s(:) !! allocated k singular values
        complex(dp), allocatable, intent(out) :: v(:, :) !! allocated n-by-k factor
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), allocatable :: b(:, :), proj(:, :)
        integer, allocatable :: idx(:)
        integer :: m, n

        m = size(a,1)
        n = size(a,2)
        k = 0
        allocate(u(0,0),s(0),v(0,0))
        call interp_decomp(a,eps,k,idx,proj,status)
        if (status /= linalg_status_success) then
            k = 0
            return
        end if
        deallocate(u,s,v)
        allocate(u(m,k),s(k),v(n,k))
        u = cmplx(0.0_dp,0.0_dp,dp)
        s = 0.0_dp
        v = cmplx(0.0_dp,0.0_dp,dp)
        if (k == 0) return
        allocate(b(m,k))
        b = a(:,idx(1:k))
        call id_to_svd(b,idx,proj,u,s,v,status)
        if (status /= linalg_status_success) then
            k = 0
            deallocate(u,s,v)
            allocate(u(0,0),s(0),v(0,0))
        end if
    end subroutine interpolative_svd_complex_tol

end module scifort_linalg_interpolative_svd
