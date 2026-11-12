! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Complex sine/cosine and hyperbolic sine/cosine integrals.
!
! This is an independent implementation based on the defining series and
! interrelations in NIST DLMF 6.2. Small arguments use the convergent defining
! series directly. Larger sine/cosine-integral arguments reuse SciFort's
! principal E1 implementation; larger hyperbolic arguments reuse principal Ei.

module scifort_complex_sine_integrals
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_complex_exponential_integrals, only : exp1_complex, expi_complex
    use scifort_exponential_integrals, only : shichi_real => shichi, sici_real => sici
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: shichi_complex
    public :: sici_complex

    real(dp), parameter :: euler_gamma = &
        0.577215664901532860606512090082402431_dp
    real(dp), parameter :: series_radius_squared = 64.0_dp
    integer, parameter :: max_iterations = 10000

contains

    pure elemental subroutine sici_complex(z, si, ci)
        complex(dp), intent(in) :: z !! complex argument of the sine and cosine integrals
        complex(dp), intent(out) :: si !! principal sine integral Si(z)
        complex(dp), intent(out) :: ci !! principal cosine integral Ci(z)

        complex(dp) :: ci_work
        complex(dp) :: si_work
        complex(dp) :: w
        real(dp) :: ci_real
        real(dp) :: si_real
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            si = complex_nan()
            ci = complex_nan()
            return
        end if

        if (y == 0.0_dp) then
            if (.not. ieee_is_finite(x)) then
                if (x > 0.0_dp) then
                    si = cmplx(0.5_dp * scifort_pi, 0.0_dp, kind=dp)
                    ci = cmplx(0.0_dp, 0.0_dp, kind=dp)
                else
                    si = cmplx(-0.5_dp * scifort_pi, 0.0_dp, kind=dp)
                    ci = cmplx(0.0_dp, scifort_pi, kind=dp)
                end if
                return
            end if

            if (x == 0.0_dp) then
                si = z
                ci = cmplx(negative_infinity(0.0_dp), quiet_nan(0.0_dp), kind=dp)
                return
            end if

            call sici_real(x, si_real, ci_real)
            si = cmplx(si_real, 0.0_dp, kind=dp)
            if (x < 0.0_dp) then
                ci = cmplx(ci_real, scifort_pi, kind=dp)
            else
                ci = cmplx(ci_real, 0.0_dp, kind=dp)
            end if
            return
        end if

        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            si = complex_nan()
            ci = complex_nan()
            return
        end if

        if (abs(z) * abs(z) <= series_radius_squared) then
            call sici_series(z, si, ci)
            return
        end if

        if (x > 0.0_dp .or. (x == 0.0_dp .and. y > 0.0_dp)) then
            call sici_right_half_plane(z, si, ci)
        else
            w = -z
            call sici_right_half_plane(w, si_work, ci_work)
            si = -si_work
            if (y > 0.0_dp) then
                ci = ci_work + cmplx(0.0_dp, scifort_pi, kind=dp)
            else
                ci = ci_work - cmplx(0.0_dp, scifort_pi, kind=dp)
            end if
        end if
    end subroutine sici_complex

    pure elemental subroutine shichi_complex(z, shi, chi)
        complex(dp), intent(in) :: z !! complex argument of the hyperbolic integrals
        complex(dp), intent(out) :: shi !! principal hyperbolic sine integral Shi(z)
        complex(dp), intent(out) :: chi !! principal hyperbolic cosine integral Chi(z)

        complex(dp) :: ei_minus
        complex(dp) :: ei_plus
        real(dp) :: chi_real
        real(dp) :: correction
        real(dp) :: shi_real
        real(dp) :: x
        real(dp) :: y

        x = real(z, dp)
        y = aimag(z)

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            shi = complex_nan()
            chi = complex_nan()
            return
        end if

        if (y == 0.0_dp) then
            if (.not. ieee_is_finite(x)) then
                shi = cmplx(sign(positive_infinity(0.0_dp), x), 0.0_dp, kind=dp)
                chi = cmplx(positive_infinity(0.0_dp), 0.0_dp, kind=dp)
                return
            end if

            if (x == 0.0_dp) then
                shi = z
                chi = cmplx(negative_infinity(0.0_dp), quiet_nan(0.0_dp), kind=dp)
                return
            end if

            call shichi_real(x, shi_real, chi_real)
            shi = cmplx(shi_real, 0.0_dp, kind=dp)
            if (x < 0.0_dp) then
                chi = cmplx(chi_real, scifort_pi, kind=dp)
            else
                chi = cmplx(chi_real, 0.0_dp, kind=dp)
            end if
            return
        end if

        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            shi = complex_nan()
            chi = complex_nan()
            return
        end if

        if (abs(z) * abs(z) <= series_radius_squared) then
            call shichi_series(z, shi, chi)
            return
        end if

        ei_plus = expi_complex(z)
        ei_minus = expi_complex(-z)
        correction = sign(0.5_dp * scifort_pi, y)
        shi = 0.5_dp * (ei_plus - ei_minus) - cmplx(0.0_dp, correction, kind=dp)
        chi = 0.5_dp * (ei_plus + ei_minus) + cmplx(0.0_dp, correction, kind=dp)
    end subroutine shichi_complex

    pure subroutine sici_series(z, si, ci)
        complex(dp), intent(in) :: z !! finite nonzero argument inside the series disk
        complex(dp), intent(out) :: si !! sine-integral series sum
        complex(dp), intent(out) :: ci !! cosine-integral series sum on the principal branch

        complex(dp) :: ci_sum
        complex(dp) :: ci_term
        complex(dp) :: si_sum
        complex(dp) :: si_term
        complex(dp) :: z_squared
        integer :: k

        z_squared = z * z

        si_term = z
        si_sum = si_term
        do k = 1, max_iterations
            si_term = si_term * (-z_squared) * real(2 * k - 1, dp) / &
                (real(2 * k, dp) * real(2 * k + 1, dp) ** 2)
            si_sum = si_sum + si_term
            if (abs(si_term) <= 2.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(si_sum))) exit
        end do

        ci_term = -0.25_dp * z_squared
        ci_sum = ci_term
        do k = 2, max_iterations
            ci_term = ci_term * (-z_squared) * real(k - 1, dp) / &
                (2.0_dp * real(k * k * (2 * k - 1), dp))
            ci_sum = ci_sum + ci_term
            if (abs(ci_term) <= 2.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(ci_sum))) exit
        end do

        si = si_sum
        ci = euler_gamma + log(z) + ci_sum
    end subroutine sici_series

    pure subroutine shichi_series(z, shi, chi)
        complex(dp), intent(in) :: z !! finite nonzero argument inside the series disk
        complex(dp), intent(out) :: shi !! hyperbolic-sine-integral series sum
        complex(dp), intent(out) :: chi !! hyperbolic-cosine-integral principal series sum

        complex(dp) :: chi_sum
        complex(dp) :: chi_term
        complex(dp) :: shi_sum
        complex(dp) :: shi_term
        complex(dp) :: z_squared
        integer :: k

        z_squared = z * z

        shi_term = z
        shi_sum = shi_term
        do k = 1, max_iterations
            shi_term = shi_term * z_squared * real(2 * k - 1, dp) / &
                (real(2 * k, dp) * real(2 * k + 1, dp) ** 2)
            shi_sum = shi_sum + shi_term
            if (abs(shi_term) <= 2.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(shi_sum))) exit
        end do

        chi_term = 0.25_dp * z_squared
        chi_sum = chi_term
        do k = 2, max_iterations
            chi_term = chi_term * z_squared * real(k - 1, dp) / &
                (2.0_dp * real(k * k * (2 * k - 1), dp))
            chi_sum = chi_sum + chi_term
            if (abs(chi_term) <= 2.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(chi_sum))) exit
        end do

        shi = shi_sum
        chi = euler_gamma + log(z) + chi_sum
    end subroutine shichi_series

    pure subroutine sici_right_half_plane(z, si, ci)
        complex(dp), intent(in) :: z !! argument in the closed upper/right canonical half-plane
        complex(dp), intent(out) :: si !! sine integral at z
        complex(dp), intent(out) :: ci !! cosine integral at z

        complex(dp) :: e1_minus
        complex(dp) :: e1_plus
        complex(dp) :: iz

        iz = cmplx(-aimag(z), abs(real(z, dp)), kind=dp)
        e1_plus = exp1_complex(iz)
        e1_minus = exp1_complex(-iz)
        si = cmplx(0.5_dp * scifort_pi, 0.0_dp, kind=dp) + &
            (e1_plus - e1_minus) / cmplx(0.0_dp, 2.0_dp, kind=dp)
        ci = -0.5_dp * (e1_plus + e1_minus)
    end subroutine sici_right_half_plane

    pure function complex_nan() result(value)
        complex(dp) :: value

        real(dp) :: nan_value

        nan_value = quiet_nan(0.0_dp)
        value = cmplx(nan_value, nan_value, kind=dp)
    end function complex_nan

end module scifort_complex_sine_integrals
