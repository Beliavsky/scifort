! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Leslie age-structured population matrices matching scipy.linalg.leslie.
! Rank-1 inputs produce one matrix; rank-2 inputs produce a leading batch.
! Rank-2 input batch dimensions may broadcast when their size is one.
module scifort_linalg_leslie
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private

    public :: leslie

    interface leslie
        module procedure leslie_rr
        module procedure leslie_cc
        module procedure leslie_rc
        module procedure leslie_cr
        module procedure leslie_i32
        module procedure leslie_i64
        module procedure leslie_batch_rr
        module procedure leslie_batch_cc
        module procedure leslie_batch_rc
        module procedure leslie_batch_cr
    end interface leslie

contains

    pure subroutine leslie_rr(f, s, a, status)
        real(dp), intent(in) :: f(:) !! fecundities, length n >= 2
        real(dp), intent(in) :: s(:) !! survival coefficients, length n-1
        real(dp), intent(out) :: a(:,:) !! Leslie matrix, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, n

        a = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(f)
        if (n < 2 .or. size(s) /= n-1) return
        if (size(a,1) /= n .or. size(a,2) /= n) return
        a(1,:) = f
        do i=1,n-1
            a(i+1,i) = s(i)
        end do
        status = linalg_status_success
    end subroutine leslie_rr

    pure subroutine leslie_cc(f, s, a, status)
        complex(dp), intent(in) :: f(:) !! fecundities, length n >= 2
        complex(dp), intent(in) :: s(:) !! survival coefficients, length n-1
        complex(dp), intent(out) :: a(:,:) !! Leslie matrix, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, n

        a = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(f)
        if (n < 2 .or. size(s) /= n-1) return
        if (size(a,1) /= n .or. size(a,2) /= n) return
        a(1,:) = f
        do i=1,n-1
            a(i+1,i) = s(i)
        end do
        status = linalg_status_success
    end subroutine leslie_cc

    pure subroutine leslie_rc(f, s, a, status)
        real(dp), intent(in) :: f(:) !! fecundities, length n >= 2
        complex(dp), intent(in) :: s(:) !! survival coefficients, length n-1
        complex(dp), intent(out) :: a(:,:) !! Leslie matrix, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, n

        a = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(f)
        if (n < 2 .or. size(s) /= n-1) return
        if (size(a,1) /= n .or. size(a,2) /= n) return
        a(1,:) = f
        do i=1,n-1
            a(i+1,i) = s(i)
        end do
        status = linalg_status_success
    end subroutine leslie_rc

    pure subroutine leslie_cr(f, s, a, status)
        complex(dp), intent(in) :: f(:) !! fecundities, length n >= 2
        real(dp), intent(in) :: s(:) !! survival coefficients, length n-1
        complex(dp), intent(out) :: a(:,:) !! Leslie matrix, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, n

        a = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(f)
        if (n < 2 .or. size(s) /= n-1) return
        if (size(a,1) /= n .or. size(a,2) /= n) return
        a(1,:) = f
        do i=1,n-1
            a(i+1,i) = s(i)
        end do
        status = linalg_status_success
    end subroutine leslie_cr

    pure subroutine leslie_i32(f, s, a, status)
        integer(int32), intent(in) :: f(:) !! fecundities, length n >= 2
        integer(int32), intent(in) :: s(:) !! survival coefficients, length n-1
        integer(int32), intent(out) :: a(:,:) !! Leslie matrix, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, n

        a = 0_int32
        status = linalg_status_invalid_shape
        n = size(f)
        if (n < 2 .or. size(s) /= n-1) return
        if (size(a,1) /= n .or. size(a,2) /= n) return
        a(1,:) = f
        do i=1,n-1
            a(i+1,i) = s(i)
        end do
        status = linalg_status_success
    end subroutine leslie_i32

    pure subroutine leslie_i64(f, s, a, status)
        integer(int64), intent(in) :: f(:) !! fecundities, length n >= 2
        integer(int64), intent(in) :: s(:) !! survival coefficients, length n-1
        integer(int64), intent(out) :: a(:,:) !! Leslie matrix, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, n

        a = 0_int64
        status = linalg_status_invalid_shape
        n = size(f)
        if (n < 2 .or. size(s) /= n-1) return
        if (size(a,1) /= n .or. size(a,2) /= n) return
        a(1,:) = f
        do i=1,n-1
            a(i+1,i) = s(i)
        end do
        status = linalg_status_success
    end subroutine leslie_i64

    pure subroutine leslie_batch_rr(f, s, a, status)
        real(dp), intent(in) :: f(:,:) !! fecundities, shape (bf,n)
        real(dp), intent(in) :: s(:,:) !! survival coefficients, shape (bs,n-1)
        real(dp), intent(out) :: a(:,:,:) !! Leslie matrices, shape (b,n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: b, bf, bs, i, k, n, f_index, s_index

        a = 0.0_dp
        status = linalg_status_invalid_shape
        bf = size(f,1)
        bs = size(s,1)
        n = size(f,2)
        b = max(bf,bs)
        if (n < 2 .or. b < 1 .or. size(s,2) /= n-1) return
        if ((bf /= b .and. bf /= 1) .or. (bs /= b .and. bs /= 1)) return
        if (size(a,1) /= b .or. size(a,2) /= n .or. size(a,3) /= n) return
        do k=1,b
            f_index = min(k,bf)
            s_index = min(k,bs)
            a(k,1,:) = f(f_index,:)
            do i=1,n-1
                a(k,i+1,i) = s(s_index,i)
            end do
        end do
        status = linalg_status_success
    end subroutine leslie_batch_rr

    pure subroutine leslie_batch_cc(f, s, a, status)
        complex(dp), intent(in) :: f(:,:) !! fecundities, shape (bf,n)
        complex(dp), intent(in) :: s(:,:) !! survival coefficients, shape (bs,n-1)
        complex(dp), intent(out) :: a(:,:,:) !! Leslie matrices, shape (b,n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: b, bf, bs, i, k, n, f_index, s_index

        a = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        bf = size(f,1)
        bs = size(s,1)
        n = size(f,2)
        b = max(bf,bs)
        if (n < 2 .or. b < 1 .or. size(s,2) /= n-1) return
        if ((bf /= b .and. bf /= 1) .or. (bs /= b .and. bs /= 1)) return
        if (size(a,1) /= b .or. size(a,2) /= n .or. size(a,3) /= n) return
        do k=1,b
            f_index = min(k,bf)
            s_index = min(k,bs)
            a(k,1,:) = f(f_index,:)
            do i=1,n-1
                a(k,i+1,i) = s(s_index,i)
            end do
        end do
        status = linalg_status_success
    end subroutine leslie_batch_cc

    pure subroutine leslie_batch_rc(f, s, a, status)
        real(dp), intent(in) :: f(:,:) !! fecundities, shape (bf,n)
        complex(dp), intent(in) :: s(:,:) !! survival coefficients, shape (bs,n-1)
        complex(dp), intent(out) :: a(:,:,:) !! Leslie matrices, shape (b,n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: b, bf, bs, i, k, n, f_index, s_index

        a = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        bf = size(f,1)
        bs = size(s,1)
        n = size(f,2)
        b = max(bf,bs)
        if (n < 2 .or. b < 1 .or. size(s,2) /= n-1) return
        if ((bf /= b .and. bf /= 1) .or. (bs /= b .and. bs /= 1)) return
        if (size(a,1) /= b .or. size(a,2) /= n .or. size(a,3) /= n) return
        do k=1,b
            f_index = min(k,bf)
            s_index = min(k,bs)
            a(k,1,:) = f(f_index,:)
            do i=1,n-1
                a(k,i+1,i) = s(s_index,i)
            end do
        end do
        status = linalg_status_success
    end subroutine leslie_batch_rc

    pure subroutine leslie_batch_cr(f, s, a, status)
        complex(dp), intent(in) :: f(:,:) !! fecundities, shape (bf,n)
        real(dp), intent(in) :: s(:,:) !! survival coefficients, shape (bs,n-1)
        complex(dp), intent(out) :: a(:,:,:) !! Leslie matrices, shape (b,n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: b, bf, bs, i, k, n, f_index, s_index

        a = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        bf = size(f,1)
        bs = size(s,1)
        n = size(f,2)
        b = max(bf,bs)
        if (n < 2 .or. b < 1 .or. size(s,2) /= n-1) return
        if ((bf /= b .and. bf /= 1) .or. (bs /= b .and. bs /= 1)) return
        if (size(a,1) /= b .or. size(a,2) /= n .or. size(a,3) /= n) return
        do k=1,b
            f_index = min(k,bf)
            s_index = min(k,bs)
            a(k,1,:) = f(f_index,:)
            do i=1,n-1
                a(k,i+1,i) = s(s_index,i)
            end do
        end do
        status = linalg_status_success
    end subroutine leslie_batch_cr

end module scifort_linalg_leslie
