! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real Jacobi elliptic functions from the arithmetic-geometric-mean inversion
! described in NIST DLMF 22.20. The real SciPy-compatible domain is 0 <= m <= 1.

module scifort_jacobi_elliptic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_legendre_elliptic, only : ellipk
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter :: max_agm_iterations = 32
    real(dp), parameter :: half_pi = 0.5_dp * scifort_pi

    public :: ellipj

contains

    pure elemental subroutine ellipj(u, m, sn, cn, dn, ph)
        real(dp), intent(in) :: u !! real argument
        real(dp), intent(in) :: m !! Legendre parameter, in the closed interval [0, 1]
        real(dp), intent(out) :: sn !! Jacobi sn(u|m)
        real(dp), intent(out) :: cn !! Jacobi cn(u|m)
        real(dp), intent(out) :: dn !! Jacobi dn(u|m)
        real(dp), intent(out) :: ph !! Jacobi amplitude am(u|m), in radians

        real(dp) :: a(0:max_agm_iterations)
        real(dp) :: b
        real(dp) :: c(0:max_agm_iterations)
        real(dp) :: complete_k
        real(dp) :: e
        real(dp) :: four_k
        real(dp) :: old_a
        real(dp) :: phi
        real(dp) :: phase_shift
        real(dp) :: reduced_u
        real(dp) :: scale
        real(dp) :: t
        integer :: i
        integer :: n

        if (ieee_is_nan(u) .or. ieee_is_nan(m)) then
            call set_nan_outputs(u, sn, cn, dn, ph)
            return
        end if
        if (.not. ieee_is_finite(u) .or. .not. ieee_is_finite(m)) then
            call set_nan_outputs(u, sn, cn, dn, ph)
            return
        end if
        if (m < 0.0_dp .or. m > 1.0_dp) then
            call set_nan_outputs(u, sn, cn, dn, ph)
            return
        end if

        if (m == 0.0_dp) then
            sn = sin(u)
            cn = cos(u)
            dn = 1.0_dp
            ph = u
            return
        end if

        if (m == 1.0_dp) then
            sn = tanh(u)
            e = exp(-abs(u))
            cn = 2.0_dp * e / (1.0_dp + e * e)
            dn = cn
            if (u >= 0.0_dp) then
                ph = half_pi - 2.0_dp * atan(e)
            else
                ph = -half_pi + 2.0_dp * atan(e)
            end if
            return
        end if

        complete_k = ellipk(m)
        four_k = 4.0_dp * complete_k
        reduced_u = modulo(u, four_k)
        if (reduced_u > 2.0_dp * complete_k) then
            reduced_u = reduced_u - four_k
        end if

        ! am(u + 4 K | m) = am(u | m) + 2 pi. Compute the amplitude shift
        ! without first forming a potentially overflowing period count.
        phase_shift = (u - reduced_u) * (half_pi / complete_k)

        a(0) = 1.0_dp
        b = sqrt(1.0_dp - m)
        c(0) = sqrt(m)
        scale = 1.0_dp
        n = 0

        do while (abs(c(n) / a(n)) > epsilon(1.0_dp))
            if (n >= max_agm_iterations) then
                call set_nan_outputs(u, sn, cn, dn, ph)
                return
            end if
            old_a = a(n)
            n = n + 1
            c(n) = 0.5_dp * (old_a - b)
            a(n) = 0.5_dp * (old_a + b)
            b = sqrt(old_a * b)
            scale = 2.0_dp * scale
        end do

        phi = scale * a(n) * reduced_u
        do i = n, 1, -1
            t = c(i) * sin(phi) / a(i)
            t = max(-1.0_dp, min(1.0_dp, t))
            phi = 0.5_dp * (asin(t) + phi)
        end do

        sn = sin(phi)
        cn = cos(phi)
        dn = sqrt(max(0.0_dp, cn * cn + (1.0_dp - m) * sn * sn))
        ph = phi + phase_shift
    end subroutine ellipj

    pure elemental subroutine set_nan_outputs(source, sn, cn, dn, ph)
        real(dp), intent(in) :: source !! value used only to select binary64 NaN kind
        real(dp), intent(out) :: sn !! NaN output
        real(dp), intent(out) :: cn !! NaN output
        real(dp), intent(out) :: dn !! NaN output
        real(dp), intent(out) :: ph !! NaN output

        sn = quiet_nan(source)
        cn = sn
        dn = sn
        ph = sn
    end subroutine set_nan_outputs

end module scifort_jacobi_elliptic
