! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Real Hurwitz zeta for s > 1 and q > 0, plus its derivative with respect
! to s.  The implementation uses the Euler-Maclaurin representation from
! DLMF 25.11.43 and is independent of SciPy source code.

module scifort_zeta
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    public :: hurwitz_zeta
    public :: hurwitz_zeta_derivative

    real(dp), parameter :: em_coeff(6) = [ &
        1.0_dp / 12.0_dp, &
        -1.0_dp / 720.0_dp, &
        1.0_dp / 30240.0_dp, &
        -1.0_dp / 1209600.0_dp, &
        1.0_dp / 47900160.0_dp, &
        -691.0_dp / 1307674368000.0_dp]

contains

    pure elemental function hurwitz_zeta(s, q) result(value)
        real(dp), intent(in) :: s !! exponent, strictly greater than one
        real(dp), intent(in) :: q !! positive Hurwitz shift
        real(dp) :: value
        real(dp) :: derivative

        call zeta_pair(s, q, value, derivative)
    end function hurwitz_zeta

    pure elemental function hurwitz_zeta_derivative(s, q) result(value)
        real(dp), intent(in) :: s !! exponent, strictly greater than one
        real(dp), intent(in) :: q !! positive Hurwitz shift
        real(dp) :: value
        real(dp) :: zeta_value

        call zeta_pair(s, q, zeta_value, value)
    end function hurwitz_zeta_derivative

    pure subroutine zeta_pair(s, q, value, derivative)
        real(dp), intent(in) :: s !! exponent, strictly greater than one
        real(dp), intent(in) :: q !! positive Hurwitz shift
        real(dp), intent(out) :: value !! Hurwitz zeta zeta(s,q)
        real(dp), intent(out) :: derivative !! derivative with respect to s
        integer, parameter :: n_shift = 24
        integer :: i, j, r
        real(dp) :: x, term, rising, harmonic, correction

        if (.not. ieee_is_finite(s) .or. .not. ieee_is_finite(q) .or. &
            s <= 1.0_dp .or. q <= 0.0_dp) then
            value = quiet_nan(s + q)
            derivative = value
            return
        end if

        value = 0.0_dp
        derivative = 0.0_dp
        do i = 0, n_shift - 1
            x = q + real(i, dp)
            term = exp(-s * log(x))
            value = value + term
            derivative = derivative - log(x) * term
        end do

        x = q + real(n_shift, dp)
        term = exp((1.0_dp - s) * log(x)) / (s - 1.0_dp)
        value = value + term
        derivative = derivative + term * (-log(x) - 1.0_dp / (s - 1.0_dp))

        term = 0.5_dp * exp(-s * log(x))
        value = value + term
        derivative = derivative - log(x) * term

        do i = 1, size(em_coeff)
            r = 2 * i - 1
            rising = 1.0_dp
            harmonic = 0.0_dp
            do j = 0, r - 1
                rising = rising * (s + real(j, dp))
                harmonic = harmonic + 1.0_dp / (s + real(j, dp))
            end do
            correction = em_coeff(i) * rising * exp(-(s + real(r, dp)) * log(x))
            value = value + correction
            derivative = derivative + correction * (harmonic - log(x))
        end do
    end subroutine zeta_pair

end module scifort_zeta
