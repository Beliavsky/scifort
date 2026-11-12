! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! One-dimensional convolution matrices, following scipy.linalg.convolution_matrix.
! A(i,j) = a(i-j+shift+1), with out-of-range coefficients set to zero.
! The offset shift is 0 (full), floor((min(m,n)-1)/2) (same),
! or min(m,n)-1 (valid). This also handles kernels longer than the input.
module scifort_linalg_convolution_matrix
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private
    public :: convolution_matrix, convolution_matrix_rows

    interface convolution_matrix
        module procedure convolution_matrix_real
        module procedure convolution_matrix_real_batch
        module procedure convolution_matrix_complex
        module procedure convolution_matrix_complex_batch
        module procedure convolution_matrix_integer
        module procedure convolution_matrix_integer_batch
    end interface convolution_matrix

contains

    pure subroutine convolution_geometry(m, n, mode, rows, shift, status)
        integer, intent(in) :: m !! coefficient vector length, positive
        integer, intent(in) :: n !! input vector length / matrix column count, positive
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        integer, intent(out) :: rows !! number of output rows; zero on failure
        integer, intent(out) :: shift !! offset in full convolution; zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(:), allocatable :: which

        rows = 0
        shift = 0
        status = linalg_status_invalid_shape
        if (m <= 0 .or. n <= 0) return
        which = 'full'
        if (present(mode)) which = mode
        select case (which)
        case ('full')
            if (m > huge(m) - n + 1) return
            rows = m + n - 1
        case ('same')
            rows = max(m,n)
            shift = (min(m,n) - 1) / 2
        case ('valid')
            rows = max(m,n) - min(m,n) + 1
            shift = min(m,n) - 1
        case default
            return
        end select
        status = linalg_status_success
    end subroutine convolution_geometry

    integer function convolution_matrix_rows(m, n, mode, status) result(rows)
        integer, intent(in) :: m !! kernel length (positive)
        integer, intent(in) :: n !! input length / number of matrix columns (positive)
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', or 'valid'
        integer, intent(out), optional :: status !! linalg_status_*; invalid_shape on error
        integer :: shift, result_status

        call convolution_geometry(m,n,mode,rows,shift,result_status)
        if (present(status)) status = result_status
    end function convolution_matrix_rows

    pure subroutine convolution_matrix_real(a, matrix, status, mode)
        real(dp), intent(in) :: a(:) !! coefficients, length m >= 1
        real(dp), intent(out) :: matrix(:,:) !! matrix, shape (k,n); k from convolution_matrix_rows
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', or 'valid'
        integer :: k, shift, i, j, coefficient

        matrix = 0.0_dp
        call convolution_geometry(size(a),size(matrix,2),mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(matrix,1) /= k) then
            status = linalg_status_invalid_shape
            return
        end if
        do j=1,size(matrix,2)
            do i=1,k
                coefficient = i-j+shift+1
                if (coefficient >= 1 .and. coefficient <= size(a)) &
                    matrix(i,j) = a(coefficient)
            end do
        end do
    end subroutine convolution_matrix_real

    pure subroutine convolution_matrix_real_batch(a, matrix, status, mode)
        real(dp), intent(in) :: a(:,:) !! batched coefficients, shape (batch,m)
        real(dp), intent(out) :: matrix(:,:,:) !! output, shape (batch,k,n)
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        integer :: k, shift, batch, s

        matrix = 0.0_dp
        call convolution_geometry(size(a,2),size(matrix,3),mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(matrix,1) /= size(a,1) .or. size(matrix,2) /= k) then
            status = linalg_status_invalid_shape
            return
        end if
        do batch=1,size(a,1)
            call convolution_matrix_real(a(batch,:),matrix(batch,:,:),s,mode)
            if (s /= linalg_status_success) then
                status = s
                matrix = 0.0_dp
                return
            end if
        end do
    end subroutine convolution_matrix_real_batch

    pure subroutine convolution_matrix_complex(a, matrix, status, mode)
        complex(dp), intent(in) :: a(:) !! coefficients, length m >= 1
        complex(dp), intent(out) :: matrix(:,:) !! matrix, shape (k,n); k from convolution_matrix_rows
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', or 'valid'
        integer :: k, shift, i, j, coefficient

        matrix = (0.0_dp,0.0_dp)
        call convolution_geometry(size(a),size(matrix,2),mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(matrix,1) /= k) then
            status = linalg_status_invalid_shape
            return
        end if
        do j=1,size(matrix,2)
            do i=1,k
                coefficient = i-j+shift+1
                if (coefficient >= 1 .and. coefficient <= size(a)) &
                    matrix(i,j) = a(coefficient)
            end do
        end do
    end subroutine convolution_matrix_complex

    pure subroutine convolution_matrix_complex_batch(a, matrix, status, mode)
        complex(dp), intent(in) :: a(:,:) !! batched coefficients, shape (batch,m)
        complex(dp), intent(out) :: matrix(:,:,:) !! output, shape (batch,k,n)
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        integer :: k, shift, batch, s

        matrix = (0.0_dp,0.0_dp)
        call convolution_geometry(size(a,2),size(matrix,3),mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(matrix,1) /= size(a,1) .or. size(matrix,2) /= k) then
            status = linalg_status_invalid_shape
            return
        end if
        do batch=1,size(a,1)
            call convolution_matrix_complex(a(batch,:),matrix(batch,:,:),s,mode)
            if (s /= linalg_status_success) then
                status = s
                matrix = (0.0_dp,0.0_dp)
                return
            end if
        end do
    end subroutine convolution_matrix_complex_batch

    pure subroutine convolution_matrix_integer(a, matrix, status, mode)
        integer, intent(in) :: a(:) !! coefficients, length m >= 1
        integer, intent(out) :: matrix(:,:) !! matrix, shape (k,n); k from convolution_matrix_rows
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', or 'valid'
        integer :: k, shift, i, j, coefficient

        matrix = 0
        call convolution_geometry(size(a),size(matrix,2),mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(matrix,1) /= k) then
            status = linalg_status_invalid_shape
            return
        end if
        do j=1,size(matrix,2)
            do i=1,k
                coefficient = i-j+shift+1
                if (coefficient >= 1 .and. coefficient <= size(a)) &
                    matrix(i,j) = a(coefficient)
            end do
        end do
    end subroutine convolution_matrix_integer

    pure subroutine convolution_matrix_integer_batch(a, matrix, status, mode)
        integer, intent(in) :: a(:,:) !! batched coefficients, shape (batch,m)
        integer, intent(out) :: matrix(:,:,:) !! output, shape (batch,k,n)
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        integer :: k, shift, batch, s

        matrix = 0
        call convolution_geometry(size(a,2),size(matrix,3),mode,k,shift,status)
        if (status /= linalg_status_success) return
        if (size(matrix,1) /= size(a,1) .or. size(matrix,2) /= k) then
            status = linalg_status_invalid_shape
            return
        end if
        do batch=1,size(a,1)
            call convolution_matrix_integer(a(batch,:),matrix(batch,:,:),s,mode)
            if (s /= linalg_status_success) then
                status = s
                matrix = 0
                return
            end if
        end do
    end subroutine convolution_matrix_integer_batch

end module scifort_linalg_convolution_matrix
