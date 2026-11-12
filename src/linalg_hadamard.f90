! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Sylvester Hadamard matrices corresponding to scipy.linalg.hadamard.
! The parity of the binary overlap of row and column indices gives the
! same ordering as Sylvester's recursive construction without temporary arrays.
module scifort_linalg_hadamard
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private

    public :: hadamard

    interface hadamard
        module procedure hadamard_i32
        module procedure hadamard_i64
        module procedure hadamard_real
        module procedure hadamard_complex
    end interface hadamard

contains

    pure logical function valid_order(n) result(valid)
        integer, intent(in) :: n !! requested positive power-of-two order
        valid = .false.
        if (n <= 0) return
        valid = iand(n, n-1) == 0
    end function valid_order

    pure integer function entry_sign(i, j) result(value)
        integer, intent(in) :: i !! one-based row index
        integer, intent(in) :: j !! one-based column index
        value = 1 - 2*mod(popcnt(iand(i-1,j-1)),2)
    end function entry_sign

    pure subroutine hadamard_i32(n, h, status)
        integer, intent(in) :: n !! order, positive power of two
        integer(int32), intent(out) :: h(:,:) !! output, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, j

        h = 0_int32
        status = linalg_status_invalid_shape
        if (.not. valid_order(n)) return
        if (size(h,1) /= n .or. size(h,2) /= n) return
        do j=1,n
            do i=1,n
                h(i,j) = int(entry_sign(i,j),int32)
            end do
        end do
        status = linalg_status_success
    end subroutine hadamard_i32

    pure subroutine hadamard_i64(n, h, status)
        integer, intent(in) :: n !! order, positive power of two
        integer(int64), intent(out) :: h(:,:) !! output, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, j

        h = 0_int64
        status = linalg_status_invalid_shape
        if (.not. valid_order(n)) return
        if (size(h,1) /= n .or. size(h,2) /= n) return
        do j=1,n
            do i=1,n
                h(i,j) = int(entry_sign(i,j),int64)
            end do
        end do
        status = linalg_status_success
    end subroutine hadamard_i64

    pure subroutine hadamard_real(n, h, status)
        integer, intent(in) :: n !! order, positive power of two
        real(dp), intent(out) :: h(:,:) !! output, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, j

        h = 0.0_dp
        status = linalg_status_invalid_shape
        if (.not. valid_order(n)) return
        if (size(h,1) /= n .or. size(h,2) /= n) return
        do j=1,n
            do i=1,n
                h(i,j) = real(entry_sign(i,j),dp)
            end do
        end do
        status = linalg_status_success
    end subroutine hadamard_real

    pure subroutine hadamard_complex(n, h, status)
        integer, intent(in) :: n !! order, positive power of two
        complex(dp), intent(out) :: h(:,:) !! output, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, j

        h = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (.not. valid_order(n)) return
        if (size(h,1) /= n .or. size(h,2) /= n) return
        do j=1,n
            do i=1,n
                h(i,j) = cmplx(real(entry_sign(i,j),dp),0.0_dp,dp)
            end do
        end do
        status = linalg_status_success
    end subroutine hadamard_complex

end module scifort_linalg_hadamard
