! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! First-quadrant complex zeros of erf(z). The implementation is original
! SciFort code based on the asymptotic zero expansion in NIST DLMF 7.13.1-2,
! followed by Newton refinement with the existing complex error function.

module scifort_erf_zeros
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_constants, only : scifort_pi
    use scifort_faddeeva, only : erf_complex
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: erf_zeros_status_success = 0
    integer, parameter, public :: erf_zeros_status_invalid_argument = 1
    integer, parameter, public :: erf_zeros_status_numerical_failure = 2

    public :: erf_zeros

contains

    pure subroutine erf_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! number of first-quadrant nonzero erf zeros requested, >= 1
        complex(dp), allocatable, intent(out) :: zeros(:) !! first nt zeros ordered by magnitude
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        integer :: info
        integer :: k

        if (nt < 1) then
            allocate(zeros(0))
            if (present(status)) status = erf_zeros_status_invalid_argument
            return
        end if

        allocate(zeros(nt))
        zeros = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        do k = 1, nt
            call refine_erf_zero(k, zeros(k), info)
            if (info /= erf_zeros_status_success) then
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = erf_zeros_status_success
    end subroutine erf_zeros

    pure subroutine refine_erf_zero(index, zero, status)
        integer, intent(in) :: index !! one-based zero index in increasing magnitude order
        complex(dp), intent(out) :: zero !! refined first-quadrant zero
        integer, intent(out) :: status !! status: 0 success; nonzero otherwise

        integer, parameter :: max_iterations = 20
        integer :: iteration
        real(dp) :: common
        real(dp) :: lambda
        real(dp) :: mu
        real(dp) :: x0
        real(dp) :: y0
        complex(dp) :: correction
        complex(dp) :: derivative
        complex(dp) :: value
        complex(dp) :: z

        lambda = sqrt((real(index, dp) - 0.125_dp) * scifort_pi)
        mu = log(lambda * sqrt(2.0_dp * scifort_pi))
        common = (1.0_dp - mu + 0.5_dp * mu * mu) / (16.0_dp * lambda**3)
        x0 = lambda - 0.25_dp * mu / lambda + common
        y0 = lambda + 0.25_dp * mu / lambda + common
        z = cmplx(x0, y0, kind=dp)

        do iteration = 1, max_iterations
            value = erf_complex(z)
            derivative = 2.0_dp / sqrt(scifort_pi) * exp(-z * z)
            if (abs(derivative) == 0.0_dp .or. .not. finite_complex(derivative)) exit
            correction = value / derivative
            z = z - correction
            if (.not. finite_complex(z)) exit
            if (abs(correction) <= 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(z))) then
                zero = z
                status = erf_zeros_status_success
                return
            end if
        end do

        zero = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        status = erf_zeros_status_numerical_failure
    end subroutine refine_erf_zero

    pure elemental function finite_complex(z) result(is_finite)
        complex(dp), intent(in) :: z !! complex value to classify
        logical :: is_finite

        is_finite = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

end module scifort_erf_zeros
