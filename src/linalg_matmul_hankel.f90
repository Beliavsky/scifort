! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! FFT-based Hankel matrix products via column reversal: H = T * J,
! where T is Toeplitz and J reverses the column order.  This module
! reuses scifort_linalg_matmul_toeplitz, including its FFT implementation.
module scifort_linalg_matmul_hankel
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_matmul_toeplitz, only : matmul_toeplitz
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: matmul_hankel

    interface matmul_hankel
        module procedure matmul_hankel_real_real_vector
        module procedure matmul_hankel_real_real_matrix
        module procedure matmul_hankel_complex_complex_vector
        module procedure matmul_hankel_complex_complex_matrix
        module procedure matmul_hankel_real_complex_vector
        module procedure matmul_hankel_real_complex_matrix
        module procedure matmul_hankel_complex_real_vector
        module procedure matmul_hankel_complex_real_matrix
    end interface matmul_hankel

contains

    ! Map H(i,j) = s(i+j-1) to T(i,j) = H(i,m-j+1).
    ! tc(i) = s(m+i-1); tr(j) = s(m-j+1), including tr(1) = tc(1).
    ! The caller must pass tr explicitly because Toeplitz's default row
    ! is the CONJUGATE of the column, not the required Hankel anti-diagonal.
    pure subroutine hankel_to_toeplitz_real(c, tc, tr, r)
        real(dp), intent(in) :: c(:) !! Hankel first column, length n
        real(dp), intent(out) :: tc(:) !! Toeplitz first column, length n
        real(dp), intent(out) :: tr(:) !! Toeplitz first row, length m
        real(dp), intent(in), optional :: r(:) !! Hankel last row, length m
        integer :: i, n, m, k

        n = size(c)
        m = size(tr)
        tc = 0.0_dp
        tr = 0.0_dp
        if (n == 0 .or. m == 0) return
        do i = 1, n
            k = m+i-1
            if (k <= n) then
                tc(i) = c(k)
            else if (present(r)) then
                tc(i) = r(k-n+1)
            end if
        end do
        do i = 1, m
            k = m-i+1
            if (k <= n) then
                tr(i) = c(k)
            else if (present(r)) then
                tr(i) = r(k-n+1)
            end if
        end do
    end subroutine hankel_to_toeplitz_real

    pure subroutine hankel_to_toeplitz_complex(c, tc, tr, r)
        complex(dp), intent(in) :: c(:) !! Hankel first column, length n
        complex(dp), intent(out) :: tc(:) !! Toeplitz first column, length n
        complex(dp), intent(out) :: tr(:) !! Toeplitz first row, length m
        complex(dp), intent(in), optional :: r(:) !! Hankel last row, length m
        integer :: i, n, m, k

        n = size(c)
        m = size(tr)
        tc = (0.0_dp,0.0_dp)
        tr = (0.0_dp,0.0_dp)
        if (n == 0 .or. m == 0) return
        do i = 1, n
            k = m+i-1
            if (k <= n) then
                tc(i) = c(k)
            else if (present(r)) then
                tc(i) = r(k-n+1)
            end if
        end do
        do i = 1, m
            k = m-i+1
            if (k <= n) then
                tr(i) = c(k)
            else if (present(r)) then
                tr(i) = r(k-n+1)
            end if
        end do
    end subroutine hankel_to_toeplitz_complex

    subroutine matmul_hankel_real_real_vector(c, x, y, status, r, check_finite)
        real(dp), intent(in) :: c(:) !! Hankel first column, length n
        real(dp), intent(in) :: x(:) !! input vector, length m
        real(dp), intent(out) :: y(:) !! output vector, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! Hankel last row; r(1) ignored
        logical, intent(in), optional :: check_finite !! check inputs, default true
        real(dp), allocatable :: work(:,:)

        y = 0.0_dp
        status = linalg_status_invalid_shape
        if (size(y) /= size(c)) return
        allocate(work(size(c),1))
        call matmul_hankel_real_real_matrix(c,reshape(x,[size(x),1]),work, &
            status,r,check_finite)
        if (status == linalg_status_success) y = work(:,1)
    end subroutine matmul_hankel_real_real_vector

    subroutine matmul_hankel_real_real_matrix(c, x, y, status, r, check_finite)
        real(dp), intent(in) :: c(:) !! Hankel first column, length n
        real(dp), intent(in) :: x(:,:) !! input matrix, shape (m,nrhs)
        real(dp), intent(out) :: y(:,:) !! output matrix, shape (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! Hankel last row; r(1) ignored
        logical, intent(in), optional :: check_finite !! check inputs, default true
        real(dp), allocatable :: tc(:), tr(:)
        integer :: n, m
        logical :: require_finite

        y = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(c)
        m = n
        if (present(r)) m = size(r)
        if (size(x,1) /= m .or. size(y,1) /= n .or. &
            size(y,2) /= size(x,2)) return
        if (n > 0 .and. m > 0) then
            if (n > huge(n) - m + 1) return
        end if
        require_finite = .true.
        if (present(check_finite)) require_finite = check_finite
        if (require_finite) then
            status = linalg_status_no_convergence
            if (.not. all(ieee_is_finite(c))) return
            if (present(r)) then
                if (.not. all(ieee_is_finite(r))) return
            end if
        end if
        allocate(tc(n),tr(m))
        call hankel_to_toeplitz_real(c,tc,tr,r)
        ! Reversing the RHS reverses the columns of the implied Hankel matrix.
        call matmul_toeplitz(tc,x(m:1:-1,:),y,status,r=tr,check_finite=check_finite)
    end subroutine matmul_hankel_real_real_matrix

    subroutine matmul_hankel_complex_complex_vector(c, x, y, status, r, check_finite)
        complex(dp), intent(in) :: c(:) !! Hankel first column, length n
        complex(dp), intent(in) :: x(:) !! input vector, length m
        complex(dp), intent(out) :: y(:) !! output vector, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! Hankel last row; r(1) ignored
        logical, intent(in), optional :: check_finite !! check inputs, default true
        complex(dp), allocatable :: work(:,:)

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        if (size(y) /= size(c)) return
        allocate(work(size(c),1))
        call matmul_hankel_complex_complex_matrix(c,reshape(x,[size(x),1]),work, &
            status,r,check_finite)
        if (status == linalg_status_success) y = work(:,1)
    end subroutine matmul_hankel_complex_complex_vector

    subroutine matmul_hankel_complex_complex_matrix(c, x, y, status, r, check_finite)
        complex(dp), intent(in) :: c(:) !! Hankel first column, length n
        complex(dp), intent(in) :: x(:,:) !! input matrix, shape (m,nrhs)
        complex(dp), intent(out) :: y(:,:) !! output matrix, shape (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! Hankel last row; r(1) ignored
        logical, intent(in), optional :: check_finite !! check inputs, default true
        complex(dp), allocatable :: tc(:), tr(:)
        integer :: n, m
        logical :: require_finite

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        m = n
        if (present(r)) m = size(r)
        if (size(x,1) /= m .or. size(y,1) /= n .or. &
            size(y,2) /= size(x,2)) return
        if (n > 0 .and. m > 0) then
            if (n > huge(n) - m + 1) return
        end if
        require_finite = .true.
        if (present(check_finite)) require_finite = check_finite
        if (require_finite) then
            status = linalg_status_no_convergence
            if (.not. finite_complex_vector(c)) return
            if (present(r)) then
                if (.not. finite_complex_vector(r)) return
            end if
        end if
        allocate(tc(n),tr(m))
        call hankel_to_toeplitz_complex(c,tc,tr,r)
        ! Reversing the RHS reverses the columns of the implied Hankel matrix.
        call matmul_toeplitz(tc,x(m:1:-1,:),y,status,r=tr,check_finite=check_finite)
    end subroutine matmul_hankel_complex_complex_matrix

    subroutine matmul_hankel_real_complex_vector(c, x, y, status, r, check_finite)
        real(dp), intent(in) :: c(:) !! Hankel first column, length n
        complex(dp), intent(in) :: x(:) !! input vector, length m
        complex(dp), intent(out) :: y(:) !! output vector, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! Hankel last row; r(1) ignored
        logical, intent(in), optional :: check_finite !! check inputs, default true
        complex(dp), allocatable :: work(:,:)

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        if (size(y) /= size(c)) return
        allocate(work(size(c),1))
        call matmul_hankel_real_complex_matrix(c,reshape(x,[size(x),1]),work, &
            status,r,check_finite)
        if (status == linalg_status_success) y = work(:,1)
    end subroutine matmul_hankel_real_complex_vector

    subroutine matmul_hankel_real_complex_matrix(c, x, y, status, r, check_finite)
        real(dp), intent(in) :: c(:) !! Hankel first column, length n
        complex(dp), intent(in) :: x(:,:) !! input matrix, shape (m,nrhs)
        complex(dp), intent(out) :: y(:,:) !! output matrix, shape (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! Hankel last row; r(1) ignored
        logical, intent(in), optional :: check_finite !! check inputs, default true
        real(dp), allocatable :: tc(:), tr(:)
        integer :: n, m
        logical :: require_finite

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        m = n
        if (present(r)) m = size(r)
        if (size(x,1) /= m .or. size(y,1) /= n .or. &
            size(y,2) /= size(x,2)) return
        if (n > 0 .and. m > 0) then
            if (n > huge(n) - m + 1) return
        end if
        require_finite = .true.
        if (present(check_finite)) require_finite = check_finite
        if (require_finite) then
            status = linalg_status_no_convergence
            if (.not. all(ieee_is_finite(c))) return
            if (present(r)) then
                if (.not. all(ieee_is_finite(r))) return
            end if
        end if
        allocate(tc(n),tr(m))
        call hankel_to_toeplitz_real(c,tc,tr,r)
        ! Reversing the RHS reverses the columns of the implied Hankel matrix.
        call matmul_toeplitz(tc,x(m:1:-1,:),y,status,r=tr,check_finite=check_finite)
    end subroutine matmul_hankel_real_complex_matrix

    subroutine matmul_hankel_complex_real_vector(c, x, y, status, r, check_finite)
        complex(dp), intent(in) :: c(:) !! Hankel first column, length n
        real(dp), intent(in) :: x(:) !! input vector, length m
        complex(dp), intent(out) :: y(:) !! output vector, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! Hankel last row; r(1) ignored
        logical, intent(in), optional :: check_finite !! check inputs, default true
        complex(dp), allocatable :: work(:,:)

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        if (size(y) /= size(c)) return
        allocate(work(size(c),1))
        call matmul_hankel_complex_real_matrix(c,reshape(x,[size(x),1]),work, &
            status,r,check_finite)
        if (status == linalg_status_success) y = work(:,1)
    end subroutine matmul_hankel_complex_real_vector

    subroutine matmul_hankel_complex_real_matrix(c, x, y, status, r, check_finite)
        complex(dp), intent(in) :: c(:) !! Hankel first column, length n
        real(dp), intent(in) :: x(:,:) !! input matrix, shape (m,nrhs)
        complex(dp), intent(out) :: y(:,:) !! output matrix, shape (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! Hankel last row; r(1) ignored
        logical, intent(in), optional :: check_finite !! check inputs, default true
        complex(dp), allocatable :: tc(:), tr(:)
        integer :: n, m
        logical :: require_finite

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        m = n
        if (present(r)) m = size(r)
        if (size(x,1) /= m .or. size(y,1) /= n .or. &
            size(y,2) /= size(x,2)) return
        if (n > 0 .and. m > 0) then
            if (n > huge(n) - m + 1) return
        end if
        require_finite = .true.
        if (present(check_finite)) require_finite = check_finite
        if (require_finite) then
            status = linalg_status_no_convergence
            if (.not. finite_complex_vector(c)) return
            if (present(r)) then
                if (.not. finite_complex_vector(r)) return
            end if
        end if
        allocate(tc(n),tr(m))
        call hankel_to_toeplitz_complex(c,tc,tr,r)
        ! Reversing the RHS reverses the columns of the implied Hankel matrix.
        call matmul_toeplitz(tc,x(m:1:-1,:),y,status,r=tr,check_finite=check_finite)
    end subroutine matmul_hankel_complex_real_matrix

    pure logical function finite_complex_vector(z) result(ok)
        complex(dp), intent(in) :: z(:) !! vector to check
        ok = all(ieee_is_finite(real(z,dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_complex_vector

end module scifort_linalg_matmul_hankel
