! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Construct rectangular SVD diagonal matrices without SVD recomputation.
module scifort_linalg_diagsvd
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private
    public :: diagsvd

    interface diagsvd
        module procedure diagsvd_i32_1
        module procedure diagsvd_i32_2
        module procedure diagsvd_i64_1
        module procedure diagsvd_i64_2
        module procedure diagsvd_r_1
        module procedure diagsvd_r_2
        module procedure diagsvd_c_1
        module procedure diagsvd_c_2
    end interface diagsvd

contains

    pure subroutine diagsvd_i32_1(singular_values, m, n, result_matrix, status)
        integer(int32), intent(in) :: singular_values(:) !! singular values; length min(m,n)
        integer, intent(in) :: m !! number of output rows, >=0
        integer, intent(in) :: n !! number of output columns, >=0
        integer(int32), intent(out) :: result_matrix(:, :) !! diagonal SVD sigma matrix
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i

        result_matrix = 0_int32
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0) return
        if (size(singular_values) /= min(m, n)) return
        if (size(result_matrix, 1) /= m .or. size(result_matrix, 2) /= n) return
        do i = 1, min(m, n)
            result_matrix(i, i) = singular_values(i)
        end do
        status = linalg_status_success
    end subroutine diagsvd_i32_1

    pure subroutine diagsvd_i32_2(singular_values, m, n, result_matrix, status)
        integer(int32), intent(in) :: singular_values(:, :) !! singular values; length min(m,n)
        integer, intent(in) :: m !! number of output rows, >=0
        integer, intent(in) :: n !! number of output columns, >=0
        integer(int32), intent(out) :: result_matrix(:, :, :) !! diagonal SVD sigma matrix
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i
        integer :: b, batch

        result_matrix = 0_int32
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0) return
        batch = size(singular_values, 1)
        if (batch == 0) return ! zero-size batches are unsupported in SciPy
        if (size(singular_values, 2) /= min(m, n)) return
        if (size(result_matrix, 1) /= batch) return
        if (size(result_matrix, 2) /= m .or. size(result_matrix, 3) /= n) return
        do b = 1, batch
            do i = 1, min(m, n)
                result_matrix(b, i, i) = singular_values(b, i)
            end do
        end do
        status = linalg_status_success
    end subroutine diagsvd_i32_2

    pure subroutine diagsvd_i64_1(singular_values, m, n, result_matrix, status)
        integer(int64), intent(in) :: singular_values(:) !! singular values; length min(m,n)
        integer, intent(in) :: m !! number of output rows, >=0
        integer, intent(in) :: n !! number of output columns, >=0
        integer(int64), intent(out) :: result_matrix(:, :) !! diagonal SVD sigma matrix
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i

        result_matrix = 0_int64
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0) return
        if (size(singular_values) /= min(m, n)) return
        if (size(result_matrix, 1) /= m .or. size(result_matrix, 2) /= n) return
        do i = 1, min(m, n)
            result_matrix(i, i) = singular_values(i)
        end do
        status = linalg_status_success
    end subroutine diagsvd_i64_1

    pure subroutine diagsvd_i64_2(singular_values, m, n, result_matrix, status)
        integer(int64), intent(in) :: singular_values(:, :) !! singular values; length min(m,n)
        integer, intent(in) :: m !! number of output rows, >=0
        integer, intent(in) :: n !! number of output columns, >=0
        integer(int64), intent(out) :: result_matrix(:, :, :) !! diagonal SVD sigma matrix
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i
        integer :: b, batch

        result_matrix = 0_int64
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0) return
        batch = size(singular_values, 1)
        if (batch == 0) return ! zero-size batches are unsupported in SciPy
        if (size(singular_values, 2) /= min(m, n)) return
        if (size(result_matrix, 1) /= batch) return
        if (size(result_matrix, 2) /= m .or. size(result_matrix, 3) /= n) return
        do b = 1, batch
            do i = 1, min(m, n)
                result_matrix(b, i, i) = singular_values(b, i)
            end do
        end do
        status = linalg_status_success
    end subroutine diagsvd_i64_2

    pure subroutine diagsvd_r_1(singular_values, m, n, result_matrix, status)
        real(dp), intent(in) :: singular_values(:) !! singular values; length min(m,n)
        integer, intent(in) :: m !! number of output rows, >=0
        integer, intent(in) :: n !! number of output columns, >=0
        real(dp), intent(out) :: result_matrix(:, :) !! diagonal SVD sigma matrix
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i

        result_matrix = 0.0_dp
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0) return
        if (size(singular_values) /= min(m, n)) return
        if (size(result_matrix, 1) /= m .or. size(result_matrix, 2) /= n) return
        do i = 1, min(m, n)
            result_matrix(i, i) = singular_values(i)
        end do
        status = linalg_status_success
    end subroutine diagsvd_r_1

    pure subroutine diagsvd_r_2(singular_values, m, n, result_matrix, status)
        real(dp), intent(in) :: singular_values(:, :) !! singular values; length min(m,n)
        integer, intent(in) :: m !! number of output rows, >=0
        integer, intent(in) :: n !! number of output columns, >=0
        real(dp), intent(out) :: result_matrix(:, :, :) !! diagonal SVD sigma matrix
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i
        integer :: b, batch

        result_matrix = 0.0_dp
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0) return
        batch = size(singular_values, 1)
        if (batch == 0) return ! zero-size batches are unsupported in SciPy
        if (size(singular_values, 2) /= min(m, n)) return
        if (size(result_matrix, 1) /= batch) return
        if (size(result_matrix, 2) /= m .or. size(result_matrix, 3) /= n) return
        do b = 1, batch
            do i = 1, min(m, n)
                result_matrix(b, i, i) = singular_values(b, i)
            end do
        end do
        status = linalg_status_success
    end subroutine diagsvd_r_2

    pure subroutine diagsvd_c_1(singular_values, m, n, result_matrix, status)
        complex(dp), intent(in) :: singular_values(:) !! singular values; length min(m,n)
        integer, intent(in) :: m !! number of output rows, >=0
        integer, intent(in) :: n !! number of output columns, >=0
        complex(dp), intent(out) :: result_matrix(:, :) !! diagonal SVD sigma matrix
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i

        result_matrix = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0) return
        if (size(singular_values) /= min(m, n)) return
        if (size(result_matrix, 1) /= m .or. size(result_matrix, 2) /= n) return
        do i = 1, min(m, n)
            result_matrix(i, i) = singular_values(i)
        end do
        status = linalg_status_success
    end subroutine diagsvd_c_1

    pure subroutine diagsvd_c_2(singular_values, m, n, result_matrix, status)
        complex(dp), intent(in) :: singular_values(:, :) !! singular values; length min(m,n)
        integer, intent(in) :: m !! number of output rows, >=0
        integer, intent(in) :: n !! number of output columns, >=0
        complex(dp), intent(out) :: result_matrix(:, :, :) !! diagonal SVD sigma matrix
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i
        integer :: b, batch

        result_matrix = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0) return
        batch = size(singular_values, 1)
        if (batch == 0) return ! zero-size batches are unsupported in SciPy
        if (size(singular_values, 2) /= min(m, n)) return
        if (size(result_matrix, 1) /= batch) return
        if (size(result_matrix, 2) /= m .or. size(result_matrix, 3) /= n) return
        do b = 1, batch
            do i = 1, min(m, n)
                result_matrix(b, i, i) = singular_values(b, i)
            end do
        end do
        status = linalg_status_success
    end subroutine diagsvd_c_2

end module scifort_linalg_diagsvd
