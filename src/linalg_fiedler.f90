! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Symmetric distance matrices following scipy.linalg.fiedler.
! For a complex input the output is real, F(i,j)=abs(a(i)-a(j)).
module scifort_linalg_fiedler
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private

    public :: fiedler

    interface fiedler
        module procedure fiedler_real
        module procedure fiedler_complex
        module procedure fiedler_integer
        module procedure fiedler_real_batch
        module procedure fiedler_complex_batch
        module procedure fiedler_integer_batch
    end interface fiedler

contains

    pure subroutine fiedler_real(a, f, status)
        real(dp), intent(in) :: a(:) !! input numbers, length n >= 0
        real(dp), intent(out) :: f(:,:) !! distance matrix, shape (n,n)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, i, j

        f = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a)
        if (size(f,1) /= n .or. size(f,2) /= n) return
        do j=1,n
            do i=1,j-1
                f(i,j) = abs(a(i)-a(j))
                f(j,i) = f(i,j)
            end do
        end do
        status = linalg_status_success
    end subroutine fiedler_real

    pure subroutine fiedler_complex(a, f, status)
        complex(dp), intent(in) :: a(:) !! complex input numbers, length n >= 0
        real(dp), intent(out) :: f(:,:) !! real nonnegative distance matrix, shape (n,n)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, i, j

        f = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a)
        if (size(f,1) /= n .or. size(f,2) /= n) return
        do j=1,n
            do i=1,j-1
                f(i,j) = abs(a(i)-a(j))
                f(j,i) = f(i,j)
            end do
        end do
        status = linalg_status_success
    end subroutine fiedler_complex

    pure subroutine fiedler_integer(a, f, status)
        integer, intent(in) :: a(:) !! integer input numbers, length n >= 0
        integer, intent(out) :: f(:,:) !! integer distance matrix, shape (n,n)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, i, j

        f = 0
        status = linalg_status_invalid_shape
        n = size(a)
        if (size(f,1) /= n .or. size(f,2) /= n) return
        do j=1,n
            do i=1,j-1
                f(i,j) = abs(a(i)-a(j))
                f(j,i) = f(i,j)
            end do
        end do
        status = linalg_status_success
    end subroutine fiedler_integer

    pure subroutine fiedler_real_batch(a, f, status)
        real(dp), intent(in) :: a(:,:) !! real batched vectors, shape (batch,n)
        real(dp), intent(out) :: f(:,:,:) !! distance matrices, shape (batch,n,n)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, k, s

        f = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        if (size(f,1) /= size(a,1) .or. &
            size(f,2) /= n .or. size(f,3) /= n) return
        do k=1,size(a,1)
            call fiedler_real(a(k,:),f(k,:,:),s)
            if (s /= linalg_status_success) then
                f = 0.0_dp
                status = s
                return
            end if
        end do
        status = linalg_status_success
    end subroutine fiedler_real_batch

    pure subroutine fiedler_complex_batch(a, f, status)
        complex(dp), intent(in) :: a(:,:) !! complex batched vectors, shape (batch,n)
        real(dp), intent(out) :: f(:,:,:) !! real distance matrices, shape (batch,n,n)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, k, s

        f = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        if (size(f,1) /= size(a,1) .or. &
            size(f,2) /= n .or. size(f,3) /= n) return
        do k=1,size(a,1)
            call fiedler_complex(a(k,:),f(k,:,:),s)
            if (s /= linalg_status_success) then
                f = 0.0_dp
                status = s
                return
            end if
        end do
        status = linalg_status_success
    end subroutine fiedler_complex_batch

    pure subroutine fiedler_integer_batch(a, f, status)
        integer, intent(in) :: a(:,:) !! integer batched vectors, shape (batch,n)
        integer, intent(out) :: f(:,:,:) !! integer distance matrices, shape (batch,n,n)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, k, s

        f = 0
        status = linalg_status_invalid_shape
        n = size(a,2)
        if (size(f,1) /= size(a,1) .or. &
            size(f,2) /= n .or. size(f,3) /= n) return
        do k=1,size(a,1)
            call fiedler_integer(a(k,:),f(k,:,:),s)
            if (s /= linalg_status_success) then
                f = 0
                status = s
                return
            end if
        end do
        status = linalg_status_success
    end subroutine fiedler_integer_batch

end module scifort_linalg_fiedler
