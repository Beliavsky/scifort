! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Positive-argument digamma helper for likelihood derivatives.
!
! The recurrence psi(x + 1) = psi(x) + 1/x moves x to the asymptotic
! region. The asymptotic expansion then uses DLMF 5.11.2 through the
! Bernoulli-number term of order x**(-14). This module is internal for now;
! scifort_special does not yet promise a public full-real-axis digamma API.

module scifort_digamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    public :: digamma_positive

contains

    pure elemental function digamma_positive(x) result(y)
        real(dp), intent(in) :: x !! positive argument; +infinity is permitted
        real(dp) :: y

        real(dp) :: r
        real(dp) :: r2
        real(dp) :: z

        if (ieee_is_nan(x) .or. x <= 0.0_dp) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = x
            return
        end if

        z = x
        y = 0.0_dp
        do while (z < 8.0_dp)
            y = y - 1.0_dp / z
            z = z + 1.0_dp
        end do

        r = 1.0_dp / z
        r2 = r * r
        y = y + log(z) - 0.5_dp * r + r2 * ( &
            -1.0_dp / 12.0_dp + r2 * ( &
            1.0_dp / 120.0_dp + r2 * ( &
            -1.0_dp / 252.0_dp + r2 * ( &
            1.0_dp / 240.0_dp + r2 * ( &
            -1.0_dp / 132.0_dp + r2 * ( &
            691.0_dp / 32760.0_dp - r2 / 12.0_dp))))))
    end function digamma_positive

end module scifort_digamma
