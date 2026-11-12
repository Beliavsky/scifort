! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! First-quadrant complex zeros of the Fresnel S and C integrals. The
! implementation is original SciFort code based on the asymptotic zero
! expansions in NIST DLMF 7.13.5-7.13.7 followed by Newton refinement with
! the existing complex Fresnel kernel.

module scifort_fresnel_zeros
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_constants, only : scifort_pi
    use scifort_faddeeva, only : fresnel_complex
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: fresnel_zeros_status_success = 0
    integer, parameter, public :: fresnel_zeros_status_invalid_argument = 1
    integer, parameter, public :: fresnel_zeros_status_numerical_failure = 2

    integer, parameter :: root_s = 1
    integer, parameter :: root_c = 2

    public :: fresnel_zeros
    public :: fresnelc_zeros
    public :: fresnels_zeros

contains

    pure subroutine fresnel_zeros(nt, zeros_sine, zeros_cosine, status)
        integer, intent(in) :: nt !! number of first-quadrant zeros requested, >= 1
        complex(dp), allocatable, intent(out) :: zeros_sine(:) !! first nt zeros of Fresnel S
        complex(dp), allocatable, intent(out) :: zeros_cosine(:) !! first nt zeros of Fresnel C
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        integer :: info

        call compute_fresnel_zeros(nt, root_s, zeros_sine, info)
        if (info /= fresnel_zeros_status_success) then
            allocate(zeros_cosine(0))
            if (present(status)) status = info
            return
        end if
        call compute_fresnel_zeros(nt, root_c, zeros_cosine, info)
        if (present(status)) status = info
    end subroutine fresnel_zeros

    pure subroutine fresnels_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! number of first-quadrant S zeros requested, >= 1
        complex(dp), allocatable, intent(out) :: zeros(:) !! first nt zeros of Fresnel S
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_fresnel_zeros(nt, root_s, zeros, status)
    end subroutine fresnels_zeros

    pure subroutine fresnelc_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! number of first-quadrant C zeros requested, >= 1
        complex(dp), allocatable, intent(out) :: zeros(:) !! first nt zeros of Fresnel C
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_fresnel_zeros(nt, root_c, zeros, status)
    end subroutine fresnelc_zeros

    pure subroutine compute_fresnel_zeros(nt, root_kind, zeros, status)
        integer, intent(in) :: nt !! number of first-quadrant zeros requested, >= 1
        integer, intent(in) :: root_kind !! root_s or root_c
        complex(dp), allocatable, intent(out) :: zeros(:) !! refined first-quadrant zeros
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        integer :: k
        integer :: info

        if (nt < 1 .or. (root_kind /= root_s .and. root_kind /= root_c)) then
            allocate(zeros(0))
            if (present(status)) status = fresnel_zeros_status_invalid_argument
            return
        end if

        allocate(zeros(nt))
        zeros = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        do k = 1, nt
            call refine_fresnel_zero(k, root_kind, zeros(k), info)
            if (info /= fresnel_zeros_status_success) then
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = fresnel_zeros_status_success
    end subroutine compute_fresnel_zeros

    pure subroutine refine_fresnel_zero(index, root_kind, zero, status)
        integer, intent(in) :: index !! one-based zero index in increasing magnitude order
        integer, intent(in) :: root_kind !! root_s or root_c
        complex(dp), intent(out) :: zero !! refined first-quadrant zero
        integer, intent(out) :: status !! status: 0 success; nonzero otherwise

        integer, parameter :: max_iterations = 24
        integer :: iteration
        real(dp) :: alpha
        real(dp) :: lambda
        real(dp) :: x0
        real(dp) :: y0
        complex(dp) :: c_value
        complex(dp) :: derivative
        complex(dp) :: correction
        complex(dp) :: phase
        complex(dp) :: s_value
        complex(dp) :: z

        if (root_kind == root_c) then
            lambda = sqrt(real(4 * index - 1, dp))
        else
            lambda = 2.0_dp * sqrt(real(index, dp))
        end if
        alpha = (2.0_dp / scifort_pi) * log(scifort_pi * lambda)
        x0 = lambda + alpha * (alpha * scifort_pi - 4.0_dp) / &
            (8.0_dp * scifort_pi * lambda**3)
        y0 = alpha / (2.0_dp * lambda)
        z = cmplx(x0, y0, kind=dp)

        do iteration = 1, max_iterations
            call fresnel_complex(z, s_value, c_value)
            phase = 0.5_dp * scifort_pi * z * z
            if (root_kind == root_c) then
                derivative = cos(phase)
                if (abs(derivative) == 0.0_dp) exit
                correction = c_value / derivative
            else
                derivative = sin(phase)
                if (abs(derivative) == 0.0_dp) exit
                correction = s_value / derivative
            end if
            z = z - correction
            if (.not. ieee_is_finite(real(z, dp)) .or. .not. ieee_is_finite(aimag(z))) exit
            if (abs(correction) <= 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(z))) then
                zero = z
                status = fresnel_zeros_status_success
                return
            end if
        end do

        zero = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        status = fresnel_zeros_status_numerical_failure
    end subroutine refine_fresnel_zero

end module scifort_fresnel_zeros
