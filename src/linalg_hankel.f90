! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-compatible Hankel matrix constructors.  A Hankel matrix has
! constant anti-diagonals; its final row's first entry is always c(n).
module scifort_linalg_hankel
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: hankel

    interface hankel
        module procedure hankel_real
        module procedure hankel_complex
        module procedure hankel_real_complex
        module procedure hankel_complex_real
    end interface hankel

contains

    subroutine hankel_real(c, h, status, r)
        real(dp), intent(in) :: c(:) !! first column, length n
        real(dp), intent(out) :: h(:,:) !! Hankel matrix, shape (n,m)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! last row, length m; r(1) ignored
        integer :: n, m, i, j, k

        h = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(c)
        m = n
        if (present(r)) m = size(r)
        if (size(h,1) /= n .or. size(h,2) /= m) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(c))) return
        if (present(r)) then
            if (.not. all(ieee_is_finite(r))) return
        end if
        do j = 1, m
            do i = 1, n
                k = i+j-1
                if (k <= n) then
                    h(i,j) = c(k)
                else if (present(r)) then
                    h(i,j) = r(k-n+1)
                end if
            end do
        end do
        status = linalg_status_success
    end subroutine hankel_real

    subroutine hankel_complex(c, h, status, r)
        complex(dp), intent(in) :: c(:) !! first column, length n
        complex(dp), intent(out) :: h(:,:) !! Hankel matrix, shape (n,m)
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! last row, length m; r(1) ignored
        integer :: n, m, i, j, k

        h = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        m = n
        if (present(r)) m = size(r)
        if (size(h,1) /= n .or. size(h,2) /= m) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(c,dp))) .or. &
            .not. all(ieee_is_finite(aimag(c)))) return
        if (present(r)) then
            if (.not. all(ieee_is_finite(real(r,dp))) .or. &
                .not. all(ieee_is_finite(aimag(r)))) return
        end if
        do j = 1, m
            do i = 1, n
                k = i+j-1
                if (k <= n) then
                    h(i,j) = c(k)
                else if (present(r)) then
                    h(i,j) = r(k-n+1)
                end if
            end do
        end do
        status = linalg_status_success
    end subroutine hankel_complex

    subroutine hankel_real_complex(c, h, status, r)
        real(dp), intent(in) :: c(:) !! real first column, length n
        complex(dp), intent(out) :: h(:,:) !! complex Hankel matrix, shape (n,m)
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in) :: r(:) !! complex last row, length m; r(1) ignored

        call hankel_complex(cmplx(c,kind=dp),h,status,r)
    end subroutine hankel_real_complex

    subroutine hankel_complex_real(c, h, status, r)
        complex(dp), intent(in) :: c(:) !! complex first column, length n
        complex(dp), intent(out) :: h(:,:) !! complex Hankel matrix, shape (n,m)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in) :: r(:) !! real last row, length m; r(1) ignored

        call hankel_complex(c,h,status,cmplx(r,kind=dp))
    end subroutine hankel_complex_real

end module scifort_linalg_hankel
