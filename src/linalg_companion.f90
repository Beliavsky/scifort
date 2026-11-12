! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Polynomial companion and Fiedler companion matrices, using the conventions
! of scipy.linalg.companion and scipy.linalg.fiedler_companion.
! Coefficients are in descending order. Outputs are caller allocated; batches
! have shape (batch, degree, degree), with one-based Fortran indexing.
module scifort_linalg_companion
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private

    public :: companion, fiedler_companion

    interface companion
        module procedure companion_real
        module procedure companion_complex
        module procedure companion_integer
        module procedure companion_real_batch
        module procedure companion_complex_batch
        module procedure companion_integer_batch
    end interface companion

    interface fiedler_companion
        module procedure fiedler_companion_real
        module procedure fiedler_companion_complex
        module procedure fiedler_companion_integer
        module procedure fiedler_companion_real_batch
        module procedure fiedler_companion_complex_batch
        module procedure fiedler_companion_integer_batch
    end interface fiedler_companion

contains

    pure subroutine companion_real(a, c, status)
        real(dp), intent(in) :: a(:) !! polynomial coefficients, descending powers; length >= 2
        real(dp), intent(out) :: c(:,:) !! companion matrix, shape (size(a)-1,size(a)-1)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, i

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a)-1
        if (n < 1 .or. size(c,1) /= n .or. size(c,2) /= n) return
        if (a(1) == 0.0_dp) return
        c(1,:) = -a(2:)/a(1)
        do i=2,n
            c(i,i-1) = 1.0_dp
        end do
        status = linalg_status_success
    end subroutine companion_real

    pure subroutine companion_complex(a, c, status)
        complex(dp), intent(in) :: a(:) !! polynomial coefficients, descending powers; length >= 2
        complex(dp), intent(out) :: c(:,:) !! companion matrix, shape (size(a)-1,size(a)-1)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, i

        c = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(a)-1
        if (n < 1 .or. size(c,1) /= n .or. size(c,2) /= n) return
        if (a(1) == (0.0_dp,0.0_dp)) return
        c(1,:) = -a(2:)/a(1)
        do i=2,n
            c(i,i-1) = (1.0_dp,0.0_dp)
        end do
        status = linalg_status_success
    end subroutine companion_complex

    pure subroutine companion_integer(a, c, status)
        integer, intent(in) :: a(:) !! integer polynomial coefficients, descending powers
        real(dp), intent(out) :: c(:,:) !! real companion matrix, shape (degree,degree)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)

        call companion_real(real(a,dp),c,status)
    end subroutine companion_integer

    pure subroutine companion_real_batch(a, c, status)
        real(dp), intent(in) :: a(:,:) !! polynomial batch, shape (batch,number of coefficients)
        real(dp), intent(out) :: c(:,:,:) !! companion matrices, shape (batch,degree,degree)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, k, s

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)-1
        if (n < 1 .or. size(c,1) /= size(a,1) .or. &
            size(c,2) /= n .or. size(c,3) /= n) return
        if (any(a(:,1) == 0.0_dp)) return
        do k=1,size(a,1)
            call companion_real(a(k,:),c(k,:,:),s)
            if (s /= linalg_status_success) then
                c = 0.0_dp
                status = s
                return
            end if
        end do
        status = linalg_status_success
    end subroutine companion_real_batch

    pure subroutine companion_complex_batch(a, c, status)
        complex(dp), intent(in) :: a(:,:) !! polynomial batch, shape (batch,number of coefficients)
        complex(dp), intent(out) :: c(:,:,:) !! companion matrices, shape (batch,degree,degree)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, k, s

        c = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(a,2)-1
        if (n < 1 .or. size(c,1) /= size(a,1) .or. &
            size(c,2) /= n .or. size(c,3) /= n) return
        if (any(a(:,1) == (0.0_dp,0.0_dp))) return
        do k=1,size(a,1)
            call companion_complex(a(k,:),c(k,:,:),s)
            if (s /= linalg_status_success) then
                c = (0.0_dp,0.0_dp)
                status = s
                return
            end if
        end do
        status = linalg_status_success
    end subroutine companion_complex_batch

    pure subroutine companion_integer_batch(a, c, status)
        integer, intent(in) :: a(:,:) !! integer polynomial batch, shape (batch,coefficients)
        real(dp), intent(out) :: c(:,:,:) !! companion matrices, shape (batch,degree,degree)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)

        call companion_real_batch(real(a,dp),c,status)
    end subroutine companion_integer_batch

    ! The Fiedler companion has <= 5 nonzero diagonals; unlike the ordinary
    ! companion, its coefficients are interleaved with identity entries.
    ! For degree one, the two companion forms are identical.
    pure subroutine fiedler_companion_real(a, c, status)
        real(dp), intent(in) :: a(:) !! coefficients in descending order; nonzero leading term
        real(dp), intent(out) :: c(:,:) !! Fiedler companion, shape (max(0,size(a)-1),same)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, k, j

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = max(0,size(a)-1)
        if (size(c,1) /= n .or. size(c,2) /= n) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        if (a(1) == 0.0_dp) return
        c(1,1) = -a(2)/a(1)
        if (n >= 2) c(2,1) = 1.0_dp
        do k=3,size(a)
            if (mod(k,2) == 1) then
                c(k-2,k-1) = -a(k)/a(1)
            else
                c(k-1,k-2) = -a(k)/a(1)
            end if
        end do
        do j=3,n,2
            c(j-2,j) = 1.0_dp
        end do
        do j=4,n,2
            c(j,j-2) = 1.0_dp
        end do
        status = linalg_status_success
    end subroutine fiedler_companion_real

    pure subroutine fiedler_companion_complex(a, c, status)
        complex(dp), intent(in) :: a(:) !! coefficients in descending order; nonzero leading term
        complex(dp), intent(out) :: c(:,:) !! Fiedler companion, shape (max(0,size(a)-1),same)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, k, j

        c = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = max(0,size(a)-1)
        if (size(c,1) /= n .or. size(c,2) /= n) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        if (a(1) == (0.0_dp,0.0_dp)) return
        c(1,1) = -a(2)/a(1)
        if (n >= 2) c(2,1) = (1.0_dp,0.0_dp)
        do k=3,size(a)
            if (mod(k,2) == 1) then
                c(k-2,k-1) = -a(k)/a(1)
            else
                c(k-1,k-2) = -a(k)/a(1)
            end if
        end do
        do j=3,n,2
            c(j-2,j) = (1.0_dp,0.0_dp)
        end do
        do j=4,n,2
            c(j,j-2) = (1.0_dp,0.0_dp)
        end do
        status = linalg_status_success
    end subroutine fiedler_companion_complex

    pure subroutine fiedler_companion_integer(a, c, status)
        integer, intent(in) :: a(:) !! integer polynomial coefficients in descending order
        real(dp), intent(out) :: c(:,:) !! real Fiedler companion matrix
        integer, intent(out) :: status !! success (0) or invalid_shape (1)

        call fiedler_companion_real(real(a,dp),c,status)
    end subroutine fiedler_companion_integer

    pure subroutine fiedler_companion_real_batch(a, c, status)
        real(dp), intent(in) :: a(:,:) !! polynomial batch, shape (batch,coefficients)
        real(dp), intent(out) :: c(:,:,:) !! Fiedler companions, shape (batch,degree,degree)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, k, s

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = max(0,size(a,2)-1)
        if (size(c,1) /= size(a,1) .or. &
            size(c,2) /= n .or. size(c,3) /= n) return
        if (n >= 1) then
            if (any(a(:,1) == 0.0_dp)) return
        end if
        do k=1,size(a,1)
            call fiedler_companion_real(a(k,:),c(k,:,:),s)
            if (s /= linalg_status_success) then
                c = 0.0_dp
                status = s
                return
            end if
        end do
        status = linalg_status_success
    end subroutine fiedler_companion_real_batch

    pure subroutine fiedler_companion_complex_batch(a, c, status)
        complex(dp), intent(in) :: a(:,:) !! polynomial batch, shape (batch,coefficients)
        complex(dp), intent(out) :: c(:,:,:) !! Fiedler companions, shape (batch,degree,degree)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)
        integer :: n, k, s

        c = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = max(0,size(a,2)-1)
        if (size(c,1) /= size(a,1) .or. &
            size(c,2) /= n .or. size(c,3) /= n) return
        if (n >= 1) then
            if (any(a(:,1) == (0.0_dp,0.0_dp))) return
        end if
        do k=1,size(a,1)
            call fiedler_companion_complex(a(k,:),c(k,:,:),s)
            if (s /= linalg_status_success) then
                c = (0.0_dp,0.0_dp)
                status = s
                return
            end if
        end do
        status = linalg_status_success
    end subroutine fiedler_companion_complex_batch

    pure subroutine fiedler_companion_integer_batch(a, c, status)
        integer, intent(in) :: a(:,:) !! integer polynomial batch, shape (batch,coefficients)
        real(dp), intent(out) :: c(:,:,:) !! Fiedler companions, shape (batch,degree,degree)
        integer, intent(out) :: status !! success (0) or invalid_shape (1)

        call fiedler_companion_real_batch(real(a,dp),c,status)
    end subroutine fiedler_companion_integer_batch

end module scifort_linalg_companion
