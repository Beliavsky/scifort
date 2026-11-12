! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Explicit discrete Fourier transform matrix, matching scipy.linalg.dft.
! This O(n**2) constructor is intended for matrix-based calculations;
! for transforming vectors, use SciFort's FFT instead.
module scifort_linalg_dft
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private

    public :: dft

contains

    pure subroutine dft(n, a, status, scale)
        integer, intent(in) :: n !! DFT order, nonnegative
        complex(dp), intent(out) :: a(:,:) !! complex output matrix, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape/scale
        character(*), intent(in), optional :: scale !! 'none' (default), 'sqrtn', 'n'
        integer(int64) :: frequency
        integer :: i, j
        real(dp) :: factor, theta, two_pi

        a = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (n < 0) return
        if (size(a,1) /= n .or. size(a,2) /= n) return
        factor = 1.0_dp
        if (present(scale)) then
            select case (scale)
            case ('none')
                continue
            case ('sqrtn')
                if (n > 0) factor = 1.0_dp/sqrt(real(n,dp))
            case ('n')
                if (n > 0) factor = 1.0_dp/real(n,dp)
            case default
                return
            end select
        end if
        status = linalg_status_success
        if (n == 0) return

        two_pi = 2.0_dp*acos(-1.0_dp)
        do j=1,n
            do i=1,n
                ! Modular arithmetic keeps phases bounded and is safe for
                ! any practically allocatable square array with default ints.
                frequency = modulo(int(i-1,int64)*int(j-1,int64),int(n,int64))
                theta = -two_pi*real(frequency,dp)/real(n,dp)
                a(i,j) = factor*cmplx(cos(theta),sin(theta),dp)
            end do
        end do
        ! The constant column and row are exact unity before scaling.
        a(1,:) = cmplx(factor,0.0_dp,dp)
        a(:,1) = cmplx(factor,0.0_dp,dp)
    end subroutine dft

end module scifort_linalg_dft
