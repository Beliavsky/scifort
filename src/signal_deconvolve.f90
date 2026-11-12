! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent implementation of scipy.signal.deconvolve via polynomial long
! division. The coefficients are in convolution (descending-power) order.
module scifort_signal_deconvolve
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: deconvolution_status_success = 0
    integer, parameter, public :: deconvolution_status_invalid_input = 1

    interface deconvolve
        module procedure deconvolve_rr
        module procedure deconvolve_rc
        module procedure deconvolve_cr
        module procedure deconvolve_cc
    end interface deconvolve
    public :: deconvolve

contains

    subroutine deconvolve_rr(signal, divisor, quotient, remainder, status)
        real(dp), intent(in) :: signal(:) !! measured sequence
        real(dp), intent(in) :: divisor(:) !! known impulse response
        real(dp), allocatable, intent(out) :: quotient(:) !! recovered sequence
        real(dp), allocatable, intent(out) :: remainder(:) !! full-length residual
        integer, intent(out), optional :: status !! zero on success
        integer :: n, d, i, j

        allocate(quotient(0), remainder(0))
        if (present(status)) status = deconvolution_status_invalid_input
        n = size(signal)
        d = size(divisor)
        if (n == 0 .or. d == 0) return
        if (.not. all(ieee_is_finite(signal)) .or. .not. all(ieee_is_finite(divisor))) return
        if (divisor(1) == 0.0_dp) return
        remainder = signal
        if (d <= n) then
            deallocate(quotient)
            allocate(quotient(n-d+1))
            do i = 1, n-d+1
                quotient(i) = remainder(i)/divisor(1)
                do j = 1, d
                    remainder(i+j-1) = remainder(i+j-1) - quotient(i)*divisor(j)
                end do
            end do
        end if
        if (present(status)) status = deconvolution_status_success
    end subroutine deconvolve_rr

    subroutine deconvolve_cc(signal, divisor, quotient, remainder, status)
        complex(dp), intent(in) :: signal(:) !! measured sequence
        complex(dp), intent(in) :: divisor(:) !! known impulse response
        complex(dp), allocatable, intent(out) :: quotient(:) !! recovered sequence
        complex(dp), allocatable, intent(out) :: remainder(:) !! full-length residual
        integer, intent(out), optional :: status !! zero on success

        call deconvolve_complex_kernel(signal, divisor, quotient, remainder, status)
    end subroutine deconvolve_cc

    subroutine deconvolve_rc(signal, divisor, quotient, remainder, status)
        real(dp), intent(in) :: signal(:) !! measured real sequence
        complex(dp), intent(in) :: divisor(:) !! complex impulse response
        complex(dp), allocatable, intent(out) :: quotient(:) !! recovered sequence
        complex(dp), allocatable, intent(out) :: remainder(:) !! full-length residual
        integer, intent(out), optional :: status !! zero on success

        call deconvolve_complex_kernel(cmplx(signal,0.0_dp,dp), divisor, quotient, remainder, status)
    end subroutine deconvolve_rc

    subroutine deconvolve_cr(signal, divisor, quotient, remainder, status)
        complex(dp), intent(in) :: signal(:) !! measured complex sequence
        real(dp), intent(in) :: divisor(:) !! real impulse response
        complex(dp), allocatable, intent(out) :: quotient(:) !! recovered sequence
        complex(dp), allocatable, intent(out) :: remainder(:) !! full-length residual
        integer, intent(out), optional :: status !! zero on success

        call deconvolve_complex_kernel(signal, cmplx(divisor,0.0_dp,dp), quotient, remainder, status)
    end subroutine deconvolve_cr

    subroutine deconvolve_complex_kernel(signal, divisor, quotient, remainder, status)
        complex(dp), intent(in) :: signal(:) !! measured sequence
        complex(dp), intent(in) :: divisor(:) !! known impulse response
        complex(dp), allocatable, intent(out) :: quotient(:) !! recovered sequence
        complex(dp), allocatable, intent(out) :: remainder(:) !! full-length residual
        integer, intent(out), optional :: status !! zero on success
        integer :: n, d, i, j

        allocate(quotient(0), remainder(0))
        if (present(status)) status = deconvolution_status_invalid_input
        n = size(signal)
        d = size(divisor)
        if (n == 0 .or. d == 0) return
        if (.not. all(ieee_is_finite(real(signal,dp))) .or. &
            .not. all(ieee_is_finite(aimag(signal)))) return
        if (.not. all(ieee_is_finite(real(divisor,dp))) .or. &
            .not. all(ieee_is_finite(aimag(divisor)))) return
        if (abs(divisor(1)) == 0.0_dp) return
        remainder = signal
        if (d <= n) then
            deallocate(quotient)
            allocate(quotient(n-d+1))
            do i = 1, n-d+1
                quotient(i) = remainder(i)/divisor(1)
                do j = 1, d
                    remainder(i+j-1) = remainder(i+j-1) - quotient(i)*divisor(j)
                end do
            end do
        end if
        if (present(status)) status = deconvolution_status_success
    end subroutine deconvolve_complex_kernel

end module scifort_signal_deconvolve
