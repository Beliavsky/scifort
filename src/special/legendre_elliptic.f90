! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real Legendre elliptic integrals implemented from the Carlson connection
! formulas in NIST DLMF 19.25 using SciFort's Carlson symmetric kernels.

module scifort_legendre_elliptic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_carlson_elliptic, only : elliprd, elliprf, elliprg
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    real(dp), parameter :: half_pi = 0.5_dp * scifort_pi

    public :: ellipe
    public :: ellipeinc
    public :: ellipk
    public :: ellipkinc
    public :: ellipkm1

contains

    pure elemental function ellipk(m) result(value)
        real(dp), intent(in) :: m !! Legendre parameter m; real branch requires m <= 1
        real(dp) :: value

        if (ieee_is_nan(m)) then
            value = quiet_nan(m)
        else if (m > 1.0_dp) then
            value = quiet_nan(m)
        else if (.not. ieee_is_finite(m)) then
            value = 0.0_dp
        else if (m == 0.0_dp) then
            value = half_pi
        else
            value = elliprf(0.0_dp, 1.0_dp - m, 1.0_dp)
        end if
    end function ellipk

    pure elemental function ellipkm1(p) result(value)
        real(dp), intent(in) :: p !! complementary parameter p = 1 - m; must be nonnegative
        real(dp) :: value

        if (ieee_is_nan(p)) then
            value = quiet_nan(p)
        else if (p < 0.0_dp) then
            value = quiet_nan(p)
        else if (.not. ieee_is_finite(p)) then
            value = 0.0_dp
        else if (p == 1.0_dp) then
            value = half_pi
        else
            value = elliprf(0.0_dp, p, 1.0_dp)
        end if
    end function ellipkm1

    pure elemental function ellipe(m) result(value)
        real(dp), intent(in) :: m !! Legendre parameter m; real branch requires m <= 1
        real(dp) :: value

        real(dp) :: complementary
        real(dp) :: transformed_e

        if (ieee_is_nan(m)) then
            value = quiet_nan(m)
        else if (m > 1.0_dp) then
            value = quiet_nan(m)
        else if (.not. ieee_is_finite(m)) then
            value = positive_infinity(m)
        else if (m == 0.0_dp) then
            value = half_pi
        else if (m == 1.0_dp) then
            value = 1.0_dp
        else if (m < 0.0_dp) then
            complementary = 1.0_dp - m
            if (1.0_dp / complementary < 0.125_dp) then
                transformed_e = ellipe_near_one(1.0_dp / complementary)
            else
                transformed_e = 2.0_dp * &
                    elliprg(0.0_dp, 1.0_dp / complementary, 1.0_dp)
            end if
            value = sqrt(complementary) * transformed_e
        else if (1.0_dp - m < 0.125_dp) then
            value = ellipe_near_one(1.0_dp - m)
        else
            value = 2.0_dp * elliprg(0.0_dp, 1.0_dp - m, 1.0_dp)
        end if
    end function ellipe

    pure elemental function ellipkinc(phi, m) result(value)
        real(dp), intent(in) :: phi !! real amplitude in radians
        real(dp), intent(in) :: m !! Legendre parameter m; real branch requires m <= 1
        real(dp) :: value

        real(dp) :: complete
        real(dp) :: cycles
        real(dp) :: local
        real(dp) :: reduced

        if (ieee_is_nan(phi) .or. ieee_is_nan(m)) then
            value = quiet_nan(phi)
            return
        end if
        if (m > 1.0_dp) then
            value = quiet_nan(phi)
            return
        end if
        if (.not. ieee_is_finite(phi)) then
            if (.not. ieee_is_finite(m)) then
                value = quiet_nan(phi)
            else if (phi > 0.0_dp) then
                value = positive_infinity(phi)
            else
                value = negative_infinity(phi)
            end if
            return
        end if
        if (.not. ieee_is_finite(m)) then
            value = sign(0.0_dp, phi)
            return
        end if
        if (phi == 0.0_dp .or. m == 0.0_dp) then
            value = phi
            return
        end if

        call reduce_amplitude(phi, reduced, cycles)
        local = ellipkinc_principal(reduced, m)
        complete = ellipk(m)
        value = add_periods(cycles, complete, local, phi)
    end function ellipkinc

    pure elemental function ellipeinc(phi, m) result(value)
        real(dp), intent(in) :: phi !! real amplitude in radians
        real(dp), intent(in) :: m !! Legendre parameter m; real branch requires m <= 1
        real(dp) :: value

        real(dp) :: complete
        real(dp) :: cycles
        real(dp) :: local
        real(dp) :: reduced

        if (ieee_is_nan(phi) .or. ieee_is_nan(m)) then
            value = quiet_nan(phi)
            return
        end if
        if (m > 1.0_dp) then
            value = quiet_nan(phi)
            return
        end if
        if (.not. ieee_is_finite(phi)) then
            if (phi > 0.0_dp) then
                value = positive_infinity(phi)
            else
                value = negative_infinity(phi)
            end if
            return
        end if
        if (phi == 0.0_dp .or. m == 0.0_dp) then
            value = phi
            return
        end if
        if (.not. ieee_is_finite(m)) then
            if (phi > 0.0_dp) then
                value = positive_infinity(phi)
            else
                value = negative_infinity(phi)
            end if
            return
        end if

        call reduce_amplitude(phi, reduced, cycles)
        local = ellipeinc_principal(reduced, m)
        complete = ellipe(m)
        value = add_periods(cycles, complete, local, phi)
    end function ellipeinc

    pure elemental function ellipkinc_principal(phi, m) result(value)
        real(dp), intent(in) :: phi !! amplitude reduced to [-pi/2, pi/2]
        real(dp), intent(in) :: m !! finite Legendre parameter, <= 1
        real(dp) :: value

        real(dp) :: c
        real(dp) :: c2
        real(dp) :: s
        real(dp) :: s2
        real(dp) :: y

        if (abs(phi) == half_pi) then
            value = sign(ellipk(m), phi)
            return
        end if

        s = sin(phi)
        c = cos(phi)
        s2 = s * s
        c2 = c * c
        y = 1.0_dp - m * s2
        if (y < 0.0_dp .and. y > -8.0_dp * epsilon(1.0_dp)) y = 0.0_dp
        value = s * elliprf(c2, y, 1.0_dp)
    end function ellipkinc_principal

    pure elemental function ellipeinc_principal(phi, m) result(value)
        real(dp), intent(in) :: phi !! amplitude reduced to [-pi/2, pi/2]
        real(dp), intent(in) :: m !! finite Legendre parameter, <= 1
        real(dp) :: value

        real(dp) :: c
        real(dp) :: c2
        real(dp) :: rf_value
        real(dp) :: rd_value
        real(dp) :: s
        real(dp) :: s2
        real(dp) :: y

        if (abs(phi) == half_pi) then
            value = sign(ellipe(m), phi)
            return
        end if

        s = sin(phi)
        if (m == 1.0_dp) then
            value = s
            return
        end if

        c = cos(phi)
        s2 = s * s
        c2 = c * c
        y = 1.0_dp - m * s2
        if (y < 0.0_dp .and. y > -8.0_dp * epsilon(1.0_dp)) y = 0.0_dp
        rf_value = elliprf(c2, y, 1.0_dp)
        if (m == 0.0_dp .or. s == 0.0_dp) then
            value = s * rf_value
            return
        end if
        rd_value = elliprd(c2, y, 1.0_dp)
        value = s * (rf_value - (m * s2 / 3.0_dp) * rd_value)
    end function ellipeinc_principal

    pure elemental function ellipe_near_one(p) result(value)
        real(dp), intent(in) :: p !! complementary parameter 1 - m, in [0, 1/8)
        real(dp) :: value

        integer :: j
        real(dp) :: coefficient
        real(dp) :: d_value
        real(dp) :: bracket
        real(dp) :: log_inv_kprime
        real(dp) :: power
        real(dp) :: term

        if (p == 0.0_dp) then
            value = 1.0_dp
            return
        end if

        log_inv_kprime = -0.5_dp * log(p)
        d_value = 2.0_dp * log(2.0_dp)
        coefficient = 1.0_dp
        power = p
        value = 1.0_dp
        do j = 0, 60
            bracket = log_inv_kprime + d_value - &
                1.0_dp / (real(2 * j + 1, dp) * real(2 * j + 2, dp))
            term = 0.5_dp * coefficient * power * bracket
            value = value + term
            if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * abs(value)) exit
            d_value = d_value - 2.0_dp / &
                (real(2 * j + 1, dp) * real(2 * j + 2, dp))
            coefficient = coefficient * &
                (real(j, dp) + 0.5_dp) * (real(j, dp) + 1.5_dp) / &
                ((real(j, dp) + 2.0_dp) * (real(j, dp) + 1.0_dp))
            power = power * p
        end do
    end function ellipe_near_one

    pure elemental subroutine reduce_amplitude(phi, reduced, cycles)
        real(dp), intent(in) :: phi !! finite amplitude in radians
        real(dp), intent(out) :: reduced !! equivalent amplitude in [-pi/2, pi/2]
        real(dp), intent(out) :: cycles !! integral number of pi periods removed

        reduced = atan2(sin(phi), cos(phi))
        if (reduced > half_pi) reduced = reduced - scifort_pi
        if (reduced < -half_pi) reduced = reduced + scifort_pi
        cycles = anint((phi - reduced) / scifort_pi)
    end subroutine reduce_amplitude

    pure elemental function add_periods(cycles, complete, local, selector) result(value)
        real(dp), intent(in) :: cycles !! integral number of pi periods
        real(dp), intent(in) :: complete !! corresponding complete elliptic integral
        real(dp), intent(in) :: local !! reduced-amplitude contribution
        real(dp), intent(in) :: selector !! argument whose sign/kind selects overflow result
        real(dp) :: value

        real(dp) :: base
        real(dp) :: period

        if (cycles == 0.0_dp) then
            value = local
            return
        end if
        if (.not. ieee_is_finite(complete)) then
            if (cycles > 0.0_dp) then
                value = positive_infinity(selector)
            else
                value = negative_infinity(selector)
            end if
            return
        end if

        period = 2.0_dp * complete
        if (abs(cycles) > huge(1.0_dp) / period) then
            if (cycles > 0.0_dp) then
                value = positive_infinity(selector)
            else
                value = negative_infinity(selector)
            end if
            return
        end if
        base = cycles * period
        if (abs(base) > huge(1.0_dp) - abs(local)) then
            if (base > 0.0_dp) then
                value = positive_infinity(selector)
            else
                value = negative_infinity(selector)
            end if
            return
        end if
        value = base + local
    end function add_periods

end module scifort_legendre_elliptic
