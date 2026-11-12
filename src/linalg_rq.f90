! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! RQ factorization for real and complex rectangular matrices.
! The orthogonal/unitary factor is obtained by a QR factorization of
! J_n * A**H * J_m, where J is the order-reversal permutation matrix.
! No LAPACK dependency or explicit reversal matrix is required.
module scifort_linalg_rq
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_qr, only : qr, qr_factor
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: rq, rq_r

    interface rq
        module procedure rq_real
        module procedure rq_complex
        module procedure rq_real_batch
        module procedure rq_complex_batch
    end interface rq

    interface rq_r
        module procedure rq_r_real
        module procedure rq_r_complex
        module procedure rq_r_real_batch
        module procedure rq_r_complex_batch
    end interface rq_r

contains

    subroutine rq_real(a, r, q, status, mode, check_finite)
        real(dp), intent(in) :: a(:, :) !! m-by-n input matrix; zero dimensions allowed
        real(dp), intent(out) :: r(:, :) !! m-by-n (full) or m-by-k (economic)
        real(dp), intent(out) :: q(:, :) !! n-by-n (full) or k-by-n (economic)
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default) or 'economic'
        logical, intent(in), optional :: check_finite !! whether to check input; default true
        real(dp), allocatable :: b(:, :), qb(:, :), rb(:, :)
        integer :: m, n, k, nq, i, j
        character(:), allocatable :: selected

        r = 0.0_dp
        q = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        selected = 'full'
        if (present(mode)) selected = trim(adjustl(mode))
        select case (selected)
        case ('full')
            nq = n
        case ('economic')
            nq = k
        case default
            return
        end select
        if (size(r,1) /= m .or. size(r,2) /= nq) return
        if (size(q,1) /= nq .or. size(q,2) /= n) return
        ! The QR backend always rejects non-finite matrices, even if the
        ! optional up-front check is disabled.
        if (present(check_finite)) then
            if (check_finite .and. .not. all(ieee_is_finite(a))) then
                status = linalg_status_no_convergence
                return
            end if
        else
            if (.not. all(ieee_is_finite(a))) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        status = linalg_status_success
        if (m == 0 .or. n == 0) then
            if (selected == 'full') then
                do i = 1, n
                    q(i,i) = 1.0_dp
                end do
            end if
            return
        end if

        allocate(b(n,m),qb(n,nq),rb(nq,m))
        do j = 1, m
            do i = 1, n
                b(i,j) = a(m+1-j,n+1-i)
            end do
        end do
        call qr(b,qb,rb,status)
        if (status /= linalg_status_success) return
        do j = 1, nq
            do i = 1, m
                r(i,j) = rb(nq+1-j,m+1-i)
            end do
        end do
        do j = 1, n
            do i = 1, nq
                q(i,j) = qb(n+1-j,nq+1-i)
            end do
        end do
    end subroutine rq_real

    subroutine rq_complex(a, r, q, status, mode, check_finite)
        complex(dp), intent(in) :: a(:, :) !! m-by-n complex input; zero dimensions allowed
        complex(dp), intent(out) :: r(:, :) !! m-by-n or m-by-k R factor
        complex(dp), intent(out) :: q(:, :) !! n-by-n or k-by-n unitary-row Q factor
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default) or 'economic'
        logical, intent(in), optional :: check_finite !! whether to check input; default true
        complex(dp), allocatable :: b(:, :), qb(:, :), rb(:, :)
        integer :: m, n, k, nq, i, j
        character(:), allocatable :: selected

        r = cmplx(0.0_dp,0.0_dp,dp)
        q = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        selected = 'full'
        if (present(mode)) selected = trim(adjustl(mode))
        select case (selected)
        case ('full')
            nq = n
        case ('economic')
            nq = k
        case default
            return
        end select
        if (size(r,1) /= m .or. size(r,2) /= nq) return
        if (size(q,1) /= nq .or. size(q,2) /= n) return
        if (present(check_finite)) then
            if (check_finite) then
                if (.not. finite_complex(a)) then
                    status = linalg_status_no_convergence
                    return
                end if
            end if
        else
            if (.not. finite_complex(a)) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        status = linalg_status_success
        if (m == 0 .or. n == 0) then
            if (selected == 'full') then
                do i = 1, n
                    q(i,i) = cmplx(1.0_dp,0.0_dp,dp)
                end do
            end if
            return
        end if

        allocate(b(n,m),qb(n,nq),rb(nq,m))
        do j = 1, m
            do i = 1, n
                b(i,j) = conjg(a(m+1-j,n+1-i))
            end do
        end do
        call qr(b,qb,rb,status)
        if (status /= linalg_status_success) return
        do j = 1, nq
            do i = 1, m
                r(i,j) = conjg(rb(nq+1-j,m+1-i))
            end do
        end do
        do j = 1, n
            do i = 1, nq
                q(i,j) = conjg(qb(n+1-j,nq+1-i))
            end do
        end do
    end subroutine rq_complex

    ! The R-only path avoids allocating/building the unitary factor.
    subroutine rq_r_real(a, r, status, check_finite)
        real(dp), intent(in) :: a(:, :) !! m-by-n input matrix
        real(dp), intent(out) :: r(:, :) !! m-by-n full-mode R factor
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! whether to check input
        real(dp), allocatable :: b(:, :), af(:, :), tau(:)
        integer :: m, n, k, i, j, row, col

        r = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (size(r,1) /= m .or. size(r,2) /= n) return
        if (present(check_finite)) then
            if (check_finite .and. .not. all(ieee_is_finite(a))) then
                status = linalg_status_no_convergence
                return
            end if
        else
            if (.not. all(ieee_is_finite(a))) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        status = linalg_status_success
        if (m == 0 .or. n == 0) return
        k = min(m,n)
        allocate(b(n,m),af(n,m),tau(k))
        do j = 1, m
            do i = 1, n
                b(i,j) = a(m+1-j,n+1-i)
            end do
        end do
        call qr_factor(b,af,tau,status)
        if (status /= linalg_status_success) return
        do j = 1, n
            row = n+1-j
            if (row > k) cycle
            do i = 1, m
                col = m+1-i
                if (row <= col) r(i,j) = af(row,col)
            end do
        end do
    end subroutine rq_r_real

    subroutine rq_r_complex(a, r, status, check_finite)
        complex(dp), intent(in) :: a(:, :) !! m-by-n complex input matrix
        complex(dp), intent(out) :: r(:, :) !! m-by-n R factor
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! whether to check input
        complex(dp), allocatable :: b(:, :), af(:, :), tau(:)
        integer :: m, n, k, i, j, row, col

        r = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        if (size(r,1) /= m .or. size(r,2) /= n) return
        if (present(check_finite)) then
            if (check_finite) then
                if (.not. finite_complex(a)) then
                    status = linalg_status_no_convergence
                    return
                end if
            end if
        else
            if (.not. finite_complex(a)) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        status = linalg_status_success
        if (m == 0 .or. n == 0) return
        k = min(m,n)
        allocate(b(n,m),af(n,m),tau(k))
        do j = 1, m
            do i = 1, n
                b(i,j) = conjg(a(m+1-j,n+1-i))
            end do
        end do
        call qr_factor(b,af,tau,status)
        if (status /= linalg_status_success) return
        do j = 1, n
            row = n+1-j
            if (row > k) cycle
            do i = 1, m
                col = m+1-i
                if (row <= col) r(i,j) = conjg(af(row,col))
            end do
        end do
    end subroutine rq_r_complex

    pure logical function finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! complex matrix
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex

    subroutine rq_real_batch(a, r, q, status, mode, check_finite)
        real(dp), intent(in) :: a(:, :, :) !! batch-by-m-by-n real matrices
        real(dp), intent(out) :: r(:, :, :) !! batch-by-m-by-n or batch-by-m-by-k
        real(dp), intent(out) :: q(:, :, :) !! batch-by-n-by-n or batch-by-k-by-n
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' or 'economic'
        logical, intent(in), optional :: check_finite !! whether to check input
        integer :: nb, m, n, nq, i
        character(:), allocatable :: selected

        r = 0.0_dp
        q = 0.0_dp
        status = linalg_status_invalid_shape
        nb = size(a,1)
        m = size(a,2)
        n = size(a,3)
        selected = 'full'
        if (present(mode)) selected = trim(adjustl(mode))
        select case (selected)
        case ('full')
            nq = n
        case ('economic')
            nq = min(m,n)
        case default
            return
        end select
        if (nb < 1) return
        if (any(shape(r) /= [nb,m,nq]) .or. any(shape(q) /= [nb,nq,n])) return
        do i = 1, nb
            call rq_real(a(i,:,:),r(i,:,:),q(i,:,:),status,selected,check_finite)
            if (status /= linalg_status_success) then
                r = 0.0_dp
                q = 0.0_dp
                return
            end if
        end do
    end subroutine rq_real_batch

    subroutine rq_complex_batch(a, r, q, status, mode, check_finite)
        complex(dp), intent(in) :: a(:, :, :) !! batch-by-m-by-n complex matrices
        complex(dp), intent(out) :: r(:, :, :) !! batch-by-m-by-n or batch-by-m-by-k
        complex(dp), intent(out) :: q(:, :, :) !! batch-by-n-by-n or batch-by-k-by-n
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' or 'economic'
        logical, intent(in), optional :: check_finite !! whether to check input
        integer :: nb, m, n, nq, i
        character(:), allocatable :: selected

        r = cmplx(0.0_dp,0.0_dp,dp)
        q = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        nb = size(a,1)
        m = size(a,2)
        n = size(a,3)
        selected = 'full'
        if (present(mode)) selected = trim(adjustl(mode))
        select case (selected)
        case ('full')
            nq = n
        case ('economic')
            nq = min(m,n)
        case default
            return
        end select
        if (nb < 1) return
        if (any(shape(r) /= [nb,m,nq]) .or. any(shape(q) /= [nb,nq,n])) return
        do i = 1, nb
            call rq_complex(a(i,:,:),r(i,:,:),q(i,:,:),status,selected,check_finite)
            if (status /= linalg_status_success) then
                r = cmplx(0.0_dp,0.0_dp,dp)
                q = cmplx(0.0_dp,0.0_dp,dp)
                return
            end if
        end do
    end subroutine rq_complex_batch

    subroutine rq_r_real_batch(a, r, status, check_finite)
        real(dp), intent(in) :: a(:, :, :) !! batch-by-m-by-n matrices
        real(dp), intent(out) :: r(:, :, :) !! same shape as a
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! whether to check input
        integer :: nb, i

        r = 0.0_dp
        status = linalg_status_invalid_shape
        nb = size(a,1)
        if (nb < 1 .or. any(shape(r) /= shape(a))) return
        do i = 1, nb
            call rq_r_real(a(i,:,:),r(i,:,:),status,check_finite)
            if (status /= linalg_status_success) then
                r = 0.0_dp
                return
            end if
        end do
    end subroutine rq_r_real_batch

    subroutine rq_r_complex_batch(a, r, status, check_finite)
        complex(dp), intent(in) :: a(:, :, :) !! batch-by-m-by-n matrices
        complex(dp), intent(out) :: r(:, :, :) !! same shape as a
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: check_finite !! whether to check input
        integer :: nb, i

        r = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        nb = size(a,1)
        if (nb < 1 .or. any(shape(r) /= shape(a))) return
        do i = 1, nb
            call rq_r_complex(a(i,:,:),r(i,:,:),status,check_finite)
            if (status /= linalg_status_success) then
                r = cmplx(0.0_dp,0.0_dp,dp)
                return
            end if
        end do
    end subroutine rq_r_complex_batch

end module scifort_linalg_rq
