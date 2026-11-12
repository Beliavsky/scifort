! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! FFT-based multiplication by the Toeplitz convolution matrix without
! constructing its dense representation. SciFort extension to scipy.linalg.
module scifort_linalg_matmul_convolution_matrix
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_convolution_matrix, only : convolution_matrix_rows
    use scifort_linalg_matmul_toeplitz, only : matmul_toeplitz
    implicit none
    private
    public :: matmul_convolution_matrix

    interface matmul_convolution_matrix
        module procedure convmul_real_real_vector
        module procedure convmul_real_real_matrix
        module procedure convmul_real_complex_vector
        module procedure convmul_real_complex_matrix
        module procedure convmul_complex_real_vector
        module procedure convmul_complex_real_matrix
        module procedure convmul_complex_complex_vector
        module procedure convmul_complex_complex_matrix
    end interface matmul_convolution_matrix

contains

    pure subroutine real_to_toeplitz(a, col, row, shift)
        real(dp), intent(in) :: a(:) !! convolution coefficients, size m
        real(dp), intent(out) :: col(:) !! first column, length k
        real(dp), intent(out) :: row(:) !! first row, length n
        integer, intent(in) :: shift !! offset into full convolution
        integer :: i, coefficient

        col = 0.0_dp
        row = 0.0_dp
        do i = 1,size(col)
            coefficient = i + shift
            if (coefficient >= 1 .and. coefficient <= size(a)) col(i) = a(coefficient)
        end do
        do i = 1,size(row)
            coefficient = 2 - i + shift
            if (coefficient >= 1 .and. coefficient <= size(a)) row(i) = a(coefficient)
        end do
    end subroutine real_to_toeplitz

    pure subroutine complex_to_toeplitz(a, col, row, shift)
        complex(dp), intent(in) :: a(:) !! convolution coefficients, size m
        complex(dp), intent(out) :: col(:) !! first column, length k
        complex(dp), intent(out) :: row(:) !! first row, length n
        integer, intent(in) :: shift !! offset into full convolution
        integer :: i, coefficient

        col = (0.0_dp,0.0_dp)
        row = (0.0_dp,0.0_dp)
        do i = 1,size(col)
            coefficient = i + shift
            if (coefficient >= 1 .and. coefficient <= size(a)) col(i) = a(coefficient)
        end do
        do i = 1,size(row)
            coefficient = 2 - i + shift
            if (coefficient >= 1 .and. coefficient <= size(a)) row(i) = a(coefficient)
        end do
    end subroutine complex_to_toeplitz

    subroutine get_geometry(m, n, mode, rows, shift, status)
        integer, intent(in) :: m !! convolution coefficient count
        integer, intent(in) :: n !! RHS length
        character(*), intent(in), optional :: mode !! convolution mode
        integer, intent(out) :: rows !! output row count
        integer, intent(out) :: shift !! shift in full convolution
        integer, intent(out) :: status !! linalg status

        rows = convolution_matrix_rows(m,n,mode,status)
        shift = 0
        if (status /= linalg_status_success) return
        if (present(mode)) then
            select case (mode)
            case ('same')
                shift = (min(m,n)-1)/2
            case ('valid')
                shift = min(m,n)-1
            end select
        end if
        if (rows > huge(rows) - n + 1) status = linalg_status_invalid_shape
    end subroutine get_geometry

    subroutine convmul_real_real_vector(a, x, y, status, mode, check_finite)
        real(dp), intent(in) :: a(:) !! convolution kernel, length m
        real(dp), intent(in) :: x(:) !! input vector, length n
        real(dp), intent(out) :: y(:) !! output, length k, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        logical, intent(in), optional :: check_finite !! reject nonfinite, default true
        real(dp), allocatable :: result(:,:)

        y = 0.0_dp
        allocate(result(size(y),1))
        call convmul_real_real_matrix(a,reshape(x,[size(x),1]),result, &
            status,mode,check_finite)
        if (status == linalg_status_success) y = result(:,1)
    end subroutine convmul_real_real_vector

    subroutine convmul_real_real_matrix(a, x, y, status, mode, check_finite)
        real(dp), intent(in) :: a(:) !! convolution kernel, length m
        real(dp), intent(in) :: x(:,:) !! inputs, shape (n,nrhs)
        real(dp), intent(out) :: y(:,:) !! outputs, shape (k,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        logical, intent(in), optional :: check_finite !! reject nonfinite, default true
        real(dp), allocatable :: col(:), row(:)
        integer :: k, shift, n

        y = 0.0_dp
        n = size(x,1)
        call get_geometry(size(a),n,mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(y,1) /= k .or. size(y,2) /= size(x,2)) then
            status = linalg_status_invalid_shape
            return
        end if
        allocate(col(k),row(n))
        call real_to_toeplitz(a,col,row,shift)
        call matmul_toeplitz(col,x,y,status,r=row,check_finite=check_finite)
    end subroutine convmul_real_real_matrix

    subroutine convmul_real_complex_vector(a, x, y, status, mode, check_finite)
        real(dp), intent(in) :: a(:) !! convolution kernel, length m
        complex(dp), intent(in) :: x(:) !! input vector, length n
        complex(dp), intent(out) :: y(:) !! output, length k, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        logical, intent(in), optional :: check_finite !! reject nonfinite, default true
        complex(dp), allocatable :: result(:,:)

        y = (0.0_dp,0.0_dp)
        allocate(result(size(y),1))
        call convmul_real_complex_matrix(a,reshape(x,[size(x),1]),result, &
            status,mode,check_finite)
        if (status == linalg_status_success) y = result(:,1)
    end subroutine convmul_real_complex_vector

    subroutine convmul_real_complex_matrix(a, x, y, status, mode, check_finite)
        real(dp), intent(in) :: a(:) !! convolution kernel, length m
        complex(dp), intent(in) :: x(:,:) !! inputs, shape (n,nrhs)
        complex(dp), intent(out) :: y(:,:) !! outputs, shape (k,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        logical, intent(in), optional :: check_finite !! reject nonfinite, default true
        real(dp), allocatable :: col(:), row(:)
        integer :: k, shift, n

        y = (0.0_dp,0.0_dp)
        n = size(x,1)
        call get_geometry(size(a),n,mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(y,1) /= k .or. size(y,2) /= size(x,2)) then
            status = linalg_status_invalid_shape
            return
        end if
        allocate(col(k),row(n))
        call real_to_toeplitz(a,col,row,shift)
        call matmul_toeplitz(col,x,y,status,r=row,check_finite=check_finite)
    end subroutine convmul_real_complex_matrix

    subroutine convmul_complex_real_vector(a, x, y, status, mode, check_finite)
        complex(dp), intent(in) :: a(:) !! convolution kernel, length m
        real(dp), intent(in) :: x(:) !! input vector, length n
        complex(dp), intent(out) :: y(:) !! output, length k, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        logical, intent(in), optional :: check_finite !! reject nonfinite, default true
        complex(dp), allocatable :: result(:,:)

        y = (0.0_dp,0.0_dp)
        allocate(result(size(y),1))
        call convmul_complex_real_matrix(a,reshape(x,[size(x),1]),result, &
            status,mode,check_finite)
        if (status == linalg_status_success) y = result(:,1)
    end subroutine convmul_complex_real_vector

    subroutine convmul_complex_real_matrix(a, x, y, status, mode, check_finite)
        complex(dp), intent(in) :: a(:) !! convolution kernel, length m
        real(dp), intent(in) :: x(:,:) !! inputs, shape (n,nrhs)
        complex(dp), intent(out) :: y(:,:) !! outputs, shape (k,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        logical, intent(in), optional :: check_finite !! reject nonfinite, default true
        complex(dp), allocatable :: col(:), row(:)
        integer :: k, shift, n

        y = (0.0_dp,0.0_dp)
        n = size(x,1)
        call get_geometry(size(a),n,mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(y,1) /= k .or. size(y,2) /= size(x,2)) then
            status = linalg_status_invalid_shape
            return
        end if
        allocate(col(k),row(n))
        call complex_to_toeplitz(a,col,row,shift)
        call matmul_toeplitz(col,x,y,status,r=row,check_finite=check_finite)
    end subroutine convmul_complex_real_matrix

    subroutine convmul_complex_complex_vector(a, x, y, status, mode, check_finite)
        complex(dp), intent(in) :: a(:) !! convolution kernel, length m
        complex(dp), intent(in) :: x(:) !! input vector, length n
        complex(dp), intent(out) :: y(:) !! output, length k, zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        logical, intent(in), optional :: check_finite !! reject nonfinite, default true
        complex(dp), allocatable :: result(:,:)

        y = (0.0_dp,0.0_dp)
        allocate(result(size(y),1))
        call convmul_complex_complex_matrix(a,reshape(x,[size(x),1]),result, &
            status,mode,check_finite)
        if (status == linalg_status_success) y = result(:,1)
    end subroutine convmul_complex_complex_vector

    subroutine convmul_complex_complex_matrix(a, x, y, status, mode, check_finite)
        complex(dp), intent(in) :: a(:) !! convolution kernel, length m
        complex(dp), intent(in) :: x(:,:) !! inputs, shape (n,nrhs)
        complex(dp), intent(out) :: y(:,:) !! outputs, shape (k,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        logical, intent(in), optional :: check_finite !! reject nonfinite, default true
        complex(dp), allocatable :: col(:), row(:)
        integer :: k, shift, n

        y = (0.0_dp,0.0_dp)
        n = size(x,1)
        call get_geometry(size(a),n,mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(y,1) /= k .or. size(y,2) /= size(x,2)) then
            status = linalg_status_invalid_shape
            return
        end if
        allocate(col(k),row(n))
        call complex_to_toeplitz(a,col,row,shift)
        call matmul_toeplitz(col,x,y,status,r=row,check_finite=check_finite)
    end subroutine convmul_complex_complex_matrix

end module scifort_linalg_matmul_convolution_matrix
