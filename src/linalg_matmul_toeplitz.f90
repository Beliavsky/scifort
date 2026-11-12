! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! FFT-based multiplication by a rectangular Toeplitz matrix.
! Independent implementation of the circulant embedding identity.
! See docs/linalg_matmul_toeplitz.md.
module scifort_linalg_matmul_toeplitz
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_signal_spectral, only : fft_plan, spectral_status_success
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: matmul_toeplitz

    interface matmul_toeplitz
        module procedure matmul_toeplitz_real_vector
        module procedure matmul_toeplitz_real_matrix
        module procedure matmul_toeplitz_complex_vector
        module procedure matmul_toeplitz_complex_matrix
        module procedure matmul_toeplitz_real_complex_vector
        module procedure matmul_toeplitz_real_complex_matrix
        module procedure matmul_toeplitz_complex_real_vector
        module procedure matmul_toeplitz_complex_real_matrix
    end interface matmul_toeplitz

contains

    subroutine matmul_toeplitz_real_vector(c, x, y, status, r, check_finite)
        real(dp), intent(in) :: c(:) !! first column of Toeplitz matrix; length n
        real(dp), intent(in) :: x(:) !! vector to multiply; length m
        real(dp), intent(out) :: y(:) !! product; length n, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! first row, length m; r(1) ignored
        logical, intent(in), optional :: check_finite !! check all inputs, default true
        real(dp), allocatable :: product(:,:)

        y = 0.0_dp
        status = linalg_status_invalid_shape
        if (size(y) /= size(c)) return
        allocate(product(size(c),1))
        call matmul_toeplitz_real_matrix(c,reshape(x,[size(x),1]),product, &
            status,r,check_finite)
        if (status == linalg_status_success) y = product(:,1)
    end subroutine matmul_toeplitz_real_vector

    subroutine matmul_toeplitz_complex_vector(c, x, y, status, r, check_finite)
        complex(dp), intent(in) :: c(:) !! first column of Toeplitz matrix; length n
        complex(dp), intent(in) :: x(:) !! vector to multiply; length m
        complex(dp), intent(out) :: y(:) !! product; length n, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! first row, length m; r(1) ignored
        logical, intent(in), optional :: check_finite !! check all inputs, default true
        complex(dp), allocatable :: product(:,:)

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        if (size(y) /= size(c)) return
        allocate(product(size(c),1))
        call matmul_toeplitz_complex_matrix(c,reshape(x,[size(x),1]),product, &
            status,r,check_finite)
        if (status == linalg_status_success) y = product(:,1)
    end subroutine matmul_toeplitz_complex_vector

    subroutine matmul_toeplitz_real_complex_vector(c, x, y, status, r, check_finite)
        real(dp), intent(in) :: c(:) !! real first column; length n
        complex(dp), intent(in) :: x(:) !! complex multiplier; length m
        complex(dp), intent(out) :: y(:) !! complex product; length n, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! real first row, length m; r(1) ignored
        logical, intent(in), optional :: check_finite !! check all inputs, default true
        complex(dp), allocatable :: product(:,:)

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        if (size(y) /= size(c)) return
        allocate(product(size(c),1))
        call matmul_toeplitz_real_complex_matrix(c,reshape(x,[size(x),1]),product, &
            status,r,check_finite)
        if (status == linalg_status_success) y = product(:,1)
    end subroutine matmul_toeplitz_real_complex_vector

    subroutine matmul_toeplitz_complex_real_vector(c, x, y, status, r, check_finite)
        complex(dp), intent(in) :: c(:) !! complex first column; length n
        real(dp), intent(in) :: x(:) !! real multiplier; length m
        complex(dp), intent(out) :: y(:) !! complex product; length n, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! complex first row, length m; r(1) ignored
        logical, intent(in), optional :: check_finite !! check all inputs, default true
        complex(dp), allocatable :: product(:,:)

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        if (size(y) /= size(c)) return
        allocate(product(size(c),1))
        call matmul_toeplitz_complex_real_matrix(c,reshape(x,[size(x),1]),product, &
            status,r,check_finite)
        if (status == linalg_status_success) y = product(:,1)
    end subroutine matmul_toeplitz_complex_real_vector

    subroutine matmul_toeplitz_real_matrix(c, x, y, status, r, check_finite)
        real(dp), intent(in) :: c(:) !! real first column, length n
        real(dp), intent(in) :: x(:,:) !! multiplier, shape (m,nrhs)
        real(dp), intent(out) :: y(:,:) !! product, shape (n,nrhs), zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! real first row, length m; r(1) ignored
        logical, intent(in), optional :: check_finite !! check all inputs, default true
        complex(dp), allocatable :: product(:,:)

        y = 0.0_dp
        allocate(product(size(y,1),size(y,2)))
        if (present(r)) then
            call fft_toeplitz_product(cmplx(c,kind=dp),cmplx(x,kind=dp), &
                product,status,cmplx(r,kind=dp),check_finite)
        else
            call fft_toeplitz_product(cmplx(c,kind=dp),cmplx(x,kind=dp), &
                product,status,check_finite=check_finite)
        end if
        if (status == linalg_status_success) y = real(product,dp)
    end subroutine matmul_toeplitz_real_matrix

    subroutine matmul_toeplitz_complex_matrix(c, x, y, status, r, check_finite)
        complex(dp), intent(in) :: c(:) !! complex first column, length n
        complex(dp), intent(in) :: x(:,:) !! complex multiplier, shape (m,nrhs)
        complex(dp), intent(out) :: y(:,:) !! complex product, shape (n,nrhs), zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! complex first row, length m; r(1) ignored
        logical, intent(in), optional :: check_finite !! check all inputs, default true

        call fft_toeplitz_product(c,x,y,status,r,check_finite)
    end subroutine matmul_toeplitz_complex_matrix

    subroutine matmul_toeplitz_real_complex_matrix(c, x, y, status, r, check_finite)
        real(dp), intent(in) :: c(:) !! real first column, length n
        complex(dp), intent(in) :: x(:,:) !! complex multiplier, shape (m,nrhs)
        complex(dp), intent(out) :: y(:,:) !! complex product, shape (n,nrhs), zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! real first row, length m; r(1) ignored
        logical, intent(in), optional :: check_finite !! check all inputs, default true

        if (present(r)) then
            call fft_toeplitz_product(cmplx(c,kind=dp),x,y,status, &
                cmplx(r,kind=dp),check_finite)
        else
            call fft_toeplitz_product(cmplx(c,kind=dp),x,y,status, &
                check_finite=check_finite)
        end if
    end subroutine matmul_toeplitz_real_complex_matrix

    subroutine matmul_toeplitz_complex_real_matrix(c, x, y, status, r, check_finite)
        complex(dp), intent(in) :: c(:) !! complex first column, length n
        real(dp), intent(in) :: x(:,:) !! real multiplier, shape (m,nrhs)
        complex(dp), intent(out) :: y(:,:) !! complex product, shape (n,nrhs), zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! complex first row, length m; r(1) ignored
        logical, intent(in), optional :: check_finite !! check all inputs, default true

        call fft_toeplitz_product(c,cmplx(x,kind=dp),y,status,r,check_finite)
    end subroutine matmul_toeplitz_complex_real_matrix

    subroutine fft_toeplitz_product(c, x, y, status, r, check_finite)
        complex(dp), intent(in) :: c(:) !! Toeplitz first column, length n
        complex(dp), intent(in) :: x(:,:) !! matrix to multiply, shape (m,nrhs)
        complex(dp), intent(out) :: y(:,:) !! product, shape (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! first row, length m; r(1) ignored
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs, default true
        type(fft_plan) :: plan
        complex(dp), allocatable :: frequency(:)
        complex(dp), allocatable :: scratch(:)
        complex(dp), allocatable :: product(:,:)
        integer :: n
        integer :: m
        integer :: nrhs
        integer :: p
        integer :: j
        integer :: fft_status
        logical :: require_finite

        y = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        m = n
        if (present(r)) m = size(r)
        nrhs = size(x,2)
        if (size(x,1) /= m .or. size(y,1) /= n .or. size(y,2) /= nrhs) return
        require_finite = .true.
        if (present(check_finite)) require_finite = check_finite
        status = linalg_status_no_convergence
        if (require_finite) then
            if (.not. finite_complex_vector(c)) return
            if (.not. finite_complex_matrix(x)) return
            if (present(r)) then
                if (.not. finite_complex_vector(r)) return
            end if
        end if

        ! The empty matrix product is always a correctly shaped zero array.
        if (n == 0 .or. m == 0 .or. nrhs == 0) then
            status = linalg_status_success
            return
        end if
        if (n > huge(n) - m + 1) then
            status = linalg_status_invalid_shape
            return
        end if
        p = n + m - 1
        call plan%initialize(p,fft_status)
        if (fft_status /= spectral_status_success) return
        allocate(frequency(p),scratch(p),product(n,nrhs))
        frequency = (0.0_dp,0.0_dp)
        frequency(1:n) = c
        if (present(r)) then
            frequency(n+1:p) = r(m:2:-1)
        else
            frequency(n+1:p) = conjg(c(n:2:-1))
        end if
        call plan%execute(frequency,fft_status)
        if (fft_status /= spectral_status_success) return
        if (.not. finite_complex_vector(frequency)) return
        do j = 1, nrhs
            scratch = (0.0_dp,0.0_dp)
            scratch(1:m) = x(:,j)
            call plan%execute(scratch,fft_status)
            if (fft_status /= spectral_status_success) return
            scratch = scratch * frequency
            if (.not. finite_complex_vector(scratch)) return
            call plan%inverse(scratch,fft_status)
            if (fft_status /= spectral_status_success) return
            if (.not. finite_complex_vector(scratch(1:n))) return
            product(:,j) = scratch(1:n)
        end do
        y = product
        status = linalg_status_success
    end subroutine fft_toeplitz_product

    pure logical function finite_complex_vector(z) result(ok)
        complex(dp), intent(in) :: z(:) !! vector to test
        ok = all(ieee_is_finite(real(z,dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_complex_vector

    pure logical function finite_complex_matrix(z) result(ok)
        complex(dp), intent(in) :: z(:,:) !! matrix to test
        ok = all(ieee_is_finite(real(z,dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_complex_matrix

end module scifort_linalg_matmul_toeplitz
