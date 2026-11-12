! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real and complex Airy functions derived from DLMF 9.2, 9.4, 9.7, and the Airy ODE.

module scifort_airy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_jy, only : jv
    use scifort_bessel_order, only : ive, kve
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    public :: ai_zeros
    public :: airy
    public :: airye
    public :: airy_status_invalid_argument
    public :: airy_status_numerical_failure
    public :: airy_status_success
    public :: bi_zeros
    public :: itairy

    real(dp), parameter :: ai_zero = &
        0.355028053887817239260063186004183176398_dp
    real(dp), parameter :: aip_zero = &
        -0.258819403792806798405183560189203963479_dp
    real(dp), parameter :: bi_zero = &
        0.614926627446000735150922369093613553595_dp
    real(dp), parameter :: bip_zero = &
        0.448288357353826357914823710398828390866_dp
    real(dp), parameter :: one_third = 1.0_dp / 3.0_dp
    real(dp), parameter :: two_thirds = 2.0_dp / 3.0_dp
    real(dp), parameter :: inv_sqrt_three = 0.577350269189625764509148780501957456_dp
    real(dp), parameter :: inv_sqrt_pi = 0.564189583547756286948079451560772586_dp
    integer, parameter :: airy_status_success = 0
    integer, parameter :: airy_status_invalid_argument = 1
    integer, parameter :: airy_status_numerical_failure = 2
    real(dp), parameter :: oscillatory_limit_x = 1.0e22_dp
    real(dp), parameter :: gl_x(8) = [ &
        0.0950125098376374401853193354250_dp, &
        0.281603550779258913230460501460_dp, &
        0.458016777657227386342419442984_dp, &
        0.617876244402643748446671764049_dp, &
        0.755404408355003033895101194847_dp, &
        0.865631202387831743880467897712_dp, &
        0.944575023073232576077988415535_dp, &
        0.989400934991649932596154173450_dp]
    real(dp), parameter :: gl_w(8) = [ &
        0.189450610455068496285396723208_dp, &
        0.182603415044923588866763667969_dp, &
        0.169156519395002538189312079030_dp, &
        0.149595988816576732081501730547_dp, &
        0.124628971255533872052476282192_dp, &
        0.0951585116824927848099251076022_dp, &
        0.0622535239386478928628438369944_dp, &
        0.0271524594117540948517805724560_dp]
    real(dp), parameter :: complex_series_limit = 6.8_dp
    real(dp), parameter :: complex_scaled_taylor_zeta_limit = 1.0_dp
    real(dp), parameter :: complex_scaled_taylor_anchor = 8.0_dp
    integer, parameter :: complex_scaled_taylor_order = 24

    interface airy
        module procedure airy_real
        module procedure airy_complex
    end interface airy

    interface airye
        module procedure airye_real
        module procedure airye_complex
    end interface airye

contains

    pure elemental subroutine airy_real(x, ai, aip, bi, bip)
        real(dp), intent(in) :: x !! real Airy argument
        real(dp), intent(out) :: ai !! Airy Ai(x)
        real(dp), intent(out) :: aip !! derivative Ai'(x)
        real(dp), intent(out) :: bi !! Airy Bi(x)
        real(dp), intent(out) :: bip !! derivative Bi'(x)

        real(dp) :: ax
        real(dp) :: scaled_ai
        real(dp) :: scaled_aip
        real(dp) :: scaled_bi
        real(dp) :: scaled_bip
        real(dp) :: sqrt_ax
        real(dp) :: zeta

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            call set_all_nan(x, ai, aip, bi, bip)
            return
        end if

        ax = abs(x)
        if (ax <= 1.0_dp) then
            call airy_series(x, ai, aip, bi, bip)
            return
        end if

        sqrt_ax = sqrt(ax)
        if (zeta_would_overflow(ax, sqrt_ax)) then
            if (x > 0.0_dp) then
                ai = 0.0_dp
                aip = -0.0_dp
                bi = positive_infinity(x)
                bip = positive_infinity(x)
            else
                call set_all_nan(x, ai, aip, bi, bip)
            end if
            return
        end if
        zeta = two_thirds * ax * sqrt_ax

        if (x > 0.0_dp) then
            call airy_positive_scaled(x, zeta, scaled_ai, scaled_aip, scaled_bi, scaled_bip)
            ai = positive_scaled_exp(scaled_ai, -zeta, x)
            aip = -positive_scaled_exp(-scaled_aip, -zeta, x)
            bi = positive_scaled_exp(scaled_bi, zeta, x)
            bip = positive_scaled_exp(scaled_bip, zeta, x)
        else
            call airy_negative_bessel(ax, zeta, ai, aip, bi, bip)
        end if
    end subroutine airy_real

    pure elemental subroutine airy_complex(z, ai, aip, bi, bip)
        complex(dp), intent(in) :: z !! complex Airy argument
        complex(dp), intent(out) :: ai !! Airy Ai(z)
        complex(dp), intent(out) :: aip !! derivative Ai'(z)
        complex(dp), intent(out) :: bi !! Airy Bi(z)
        complex(dp), intent(out) :: bip !! derivative Bi'(z)

        complex(dp), parameter :: omega_minus = &
            cmplx(-0.5_dp, -0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: omega_plus = &
            cmplx(-0.5_dp, 0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: phase_minus = &
            cmplx(0.866025403784438646763723170753_dp, -0.5_dp, kind=dp)
        complex(dp), parameter :: phase_plus = &
            cmplx(0.866025403784438646763723170753_dp, 0.5_dp, kind=dp)
        complex(dp) :: ai_minus
        complex(dp) :: ai_plus
        complex(dp) :: aip_minus
        complex(dp) :: aip_plus
        real(dp) :: ai_real_value
        real(dp) :: aip_real_value
        real(dp) :: bi_real_value
        real(dp) :: bip_real_value

        if (.not. complex_is_finite(z)) then
            call set_all_nan_complex(z, ai, aip, bi, bip)
            return
        end if

        if (aimag(z) == 0.0_dp) then
            call airy_real(real(z, dp), ai_real_value, aip_real_value, &
                bi_real_value, bip_real_value)
            ai = cmplx(ai_real_value, 0.0_dp, kind=dp)
            aip = cmplx(aip_real_value, 0.0_dp, kind=dp)
            bi = cmplx(bi_real_value, 0.0_dp, kind=dp)
            bip = cmplx(bip_real_value, 0.0_dp, kind=dp)
            return
        end if

        if (abs(z) <= complex_series_limit) then
            call airy_series_complex(z, ai, aip, bi, bip)
            return
        end if

        call airy_ai_complex(z, ai, aip)
        call airy_ai_complex(z * omega_minus, ai_minus, aip_minus)
        call airy_ai_complex(z * omega_plus, ai_plus, aip_plus)

        ! DLMF 9.2.10 and its derivative.
        bi = phase_minus * ai_minus + phase_plus * ai_plus
        bip = phase_minus * omega_minus * aip_minus + &
            phase_plus * omega_plus * aip_plus
    end subroutine airy_complex

    pure elemental subroutine airye_real(x, eai, eaip, ebi, ebip)
        real(dp), intent(in) :: x !! real Airy argument
        real(dp), intent(out) :: eai !! exponentially scaled Ai(x)
        real(dp), intent(out) :: eaip !! exponentially scaled Ai'(x)
        real(dp), intent(out) :: ebi !! exponentially scaled Bi(x)
        real(dp), intent(out) :: ebip !! exponentially scaled Bi'(x)

        real(dp) :: ai
        real(dp) :: aip
        real(dp) :: bi
        real(dp) :: bip
        real(dp) :: sqrt_x
        real(dp) :: zeta

        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            call set_all_nan(x, eai, eaip, ebi, ebip)
            return
        end if

        if (x < 0.0_dp) then
            call airy(x, ai, aip, bi, bip)
            eai = quiet_nan(x)
            eaip = quiet_nan(x)
            ebi = bi
            ebip = bip
            return
        end if

        if (x <= 1.0_dp) then
            call airy_series(x, ai, aip, bi, bip)
            zeta = two_thirds * x * sqrt(x)
            eai = ai * exp(zeta)
            eaip = aip * exp(zeta)
            ebi = bi * exp(-zeta)
            ebip = bip * exp(-zeta)
            return
        end if

        sqrt_x = sqrt(x)
        if (zeta_would_overflow(x, sqrt_x)) then
            call airy_positive_scaled_asymptotic(x, eai, eaip, ebi, ebip)
            return
        end if
        zeta = two_thirds * x * sqrt_x
        call airy_positive_scaled(x, zeta, eai, eaip, ebi, ebip)
    end subroutine airye_real

    pure elemental subroutine airye_complex(z, eai, eaip, ebi, ebip)
        complex(dp), intent(in) :: z !! complex Airy argument
        complex(dp), intent(out) :: eai !! exponentially scaled Ai(z)
        complex(dp), intent(out) :: eaip !! exponentially scaled Ai'(z)
        complex(dp), intent(out) :: ebi !! exponentially scaled Bi(z)
        complex(dp), intent(out) :: ebip !! exponentially scaled Bi'(z)

        complex(dp) :: ai
        complex(dp) :: aip
        complex(dp) :: bi
        complex(dp) :: bip
        complex(dp) :: ai_scale
        complex(dp) :: zeta
        real(dp) :: bi_scale
        real(dp) :: eai_real
        real(dp) :: eaip_real
        real(dp) :: ebi_real
        real(dp) :: ebip_real

        if (.not. complex_is_finite(z)) then
            call set_all_nan_complex(z, eai, eaip, ebi, ebip)
            return
        end if

        if (aimag(z) == 0.0_dp .and. real(z, dp) >= 0.0_dp) then
            call airye_real(real(z, dp), eai_real, eaip_real, ebi_real, ebip_real)
            eai = cmplx(eai_real, 0.0_dp, kind=dp)
            eaip = cmplx(eaip_real, 0.0_dp, kind=dp)
            ebi = cmplx(ebi_real, 0.0_dp, kind=dp)
            ebip = cmplx(ebip_real, 0.0_dp, kind=dp)
            return
        end if

        if (abs(z) <= complex_series_limit) then
            call airy_series_complex(z, ai, aip, bi, bip)
            zeta = airy_zeta_complex(z)
            if (real(zeta, dp) > complex_scaled_taylor_zeta_limit) then
                call airy_ai_scaled_taylor(z, eai, eaip)
            else
                ai_scale = exp(zeta)
                eai = ai * ai_scale
                eaip = aip * ai_scale
            end if
            bi_scale = exp(-abs(real(zeta, dp)))
            ebi = bi * bi_scale
            ebip = bip * bi_scale
            return
        end if

        call airy_ai_scaled_complex(z, eai, eaip)
        call airy_bi_scaled_complex(z, ebi, ebip)
    end subroutine airye_complex

    pure elemental subroutine airy_series_complex(z, ai, aip, bi, bip)
        complex(dp), intent(in) :: z !! finite argument in the power-series region
        complex(dp), intent(out) :: ai !! Airy Ai(z)
        complex(dp), intent(out) :: aip !! derivative Ai'(z)
        complex(dp), intent(out) :: bi !! Airy Bi(z)
        complex(dp), intent(out) :: bip !! derivative Bi'(z)

        complex(dp) :: f
        complex(dp) :: fp
        complex(dp) :: g
        complex(dp) :: gp
        complex(dp) :: term_f
        complex(dp) :: term_g
        complex(dp) :: z3
        integer :: k
        integer :: n_f
        integer :: n_g

        f = cmplx(1.0_dp, 0.0_dp, kind=dp)
        fp = cmplx(0.0_dp, 0.0_dp, kind=dp)
        g = z
        gp = cmplx(1.0_dp, 0.0_dp, kind=dp)
        term_f = cmplx(1.0_dp, 0.0_dp, kind=dp)
        term_g = z
        z3 = z * z * z

        do k = 1, 120
            n_f = 3 * k
            n_g = 3 * k + 1
            term_f = term_f * z3 / real(n_f * (n_f - 1), dp)
            term_g = term_g * z3 / real(n_g * (n_g - 1), dp)
            f = f + term_f
            g = g + term_g
            if (z /= cmplx(0.0_dp, 0.0_dp, kind=dp)) then
                fp = fp + real(n_f, dp) * term_f / z
                gp = gp + real(n_g, dp) * term_g / z
            end if
            if (abs(term_f) <= epsilon(1.0_dp) * max(1.0_dp, abs(f)) .and. &
                abs(term_g) <= epsilon(1.0_dp) * max(1.0_dp, abs(g))) exit
        end do

        ai = ai_zero * f + aip_zero * g
        aip = ai_zero * fp + aip_zero * gp
        bi = bi_zero * f + bip_zero * g
        bip = bi_zero * fp + bip_zero * gp
    end subroutine airy_series_complex

    pure elemental subroutine airy_ai_complex(z, ai, aip)
        complex(dp), intent(in) :: z !! finite complex argument outside the series region
        complex(dp), intent(out) :: ai !! Airy Ai(z)
        complex(dp), intent(out) :: aip !! derivative Ai'(z)

        complex(dp), parameter :: rot_minus = &
            cmplx(0.5_dp, -0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: rot_plus = &
            cmplx(0.5_dp, 0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: rot2_minus = &
            cmplx(-0.5_dp, -0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: rot2_plus = &
            cmplx(-0.5_dp, 0.866025403784438646763723170753_dp, kind=dp)
        complex(dp) :: ai_minus
        complex(dp) :: ai_plus
        complex(dp) :: aip_minus
        complex(dp) :: aip_plus
        complex(dp) :: w
        real(dp) :: argument

        if (abs(z) <= complex_series_limit) then
            call airy_ai_series_complex(z, ai, aip)
            return
        end if

        argument = atan2(aimag(z), real(z, dp))
        if (abs(argument) < two_thirds * scifort_pi) then
            call airy_ai_asymptotic_direct(z, ai, aip)
            return
        end if

        ! DLMF 9.2.14, differentiated for Ai'.  Rotating -z by +/- pi/3
        ! keeps both asymptotic evaluations away from the negative-real sector.
        w = -z
        call airy_ai_asymptotic_direct(w * rot_plus, ai_plus, aip_plus)
        call airy_ai_asymptotic_direct(w * rot_minus, ai_minus, aip_minus)
        ai = rot_plus * ai_plus + rot_minus * ai_minus
        aip = -(rot2_plus * aip_plus + rot2_minus * aip_minus)
    end subroutine airy_ai_complex

    pure elemental subroutine airy_ai_series_complex(z, ai, aip)
        complex(dp), intent(in) :: z !! finite argument in the power-series region
        complex(dp), intent(out) :: ai !! Airy Ai(z)
        complex(dp), intent(out) :: aip !! derivative Ai'(z)
        complex(dp) :: bi_dummy
        complex(dp) :: bip_dummy

        call airy_series_complex(z, ai, aip, bi_dummy, bip_dummy)
    end subroutine airy_ai_series_complex

    pure elemental subroutine airy_ai_scaled_taylor(z, eai, eaip)
        complex(dp), intent(in) :: z !! finite moderate argument with positive real Airy phase
        complex(dp), intent(out) :: eai !! exp(zeta)*Ai(z)
        complex(dp), intent(out) :: eaip !! exp(zeta)*Ai'(z)

        complex(dp) :: ai
        complex(dp) :: ai_next
        complex(dp) :: aip
        complex(dp) :: aip_next
        complex(dp) :: anchor
        complex(dp) :: anchor_scale
        complex(dp) :: current
        complex(dp) :: step
        complex(dp) :: scaled_ai
        complex(dp) :: scaled_aip
        complex(dp) :: target_scale
        real(dp) :: radius
        integer :: k
        integer :: nstep

        radius = abs(z)
        anchor = complex_scaled_taylor_anchor * z / radius
        call airy_ai_asymptotic_scaled_direct(anchor, scaled_ai, scaled_aip)
        anchor_scale = exp(-airy_zeta_complex(anchor))
        ai = scaled_ai * anchor_scale
        aip = scaled_aip * anchor_scale

        nstep = max(1, ceiling(abs(anchor - z)))
        step = (z - anchor) / real(nstep, dp)
        current = anchor
        do k = 1, nstep
            call airy_taylor_step_complex(current, ai, aip, step, ai_next, aip_next)
            current = current + step
            ai = ai_next
            aip = aip_next
        end do

        target_scale = exp(airy_zeta_complex(z))
        eai = ai * target_scale
        eaip = aip * target_scale
    end subroutine airy_ai_scaled_taylor

    pure elemental subroutine airy_taylor_step_complex(z0, ai, aip, step, ai_next, aip_next)
        complex(dp), intent(in) :: z0 !! expansion center
        complex(dp), intent(in) :: ai !! Ai(z0)
        complex(dp), intent(in) :: aip !! Ai'(z0)
        complex(dp), intent(in) :: step !! complex continuation step from z0
        complex(dp), intent(out) :: ai_next !! Ai(z0+step)
        complex(dp), intent(out) :: aip_next !! Ai'(z0+step)

        complex(dp) :: coefficient(0:complex_scaled_taylor_order)
        integer :: n

        coefficient(0) = ai
        coefficient(1) = aip
        coefficient(2) = 0.5_dp * z0 * coefficient(0)
        do n = 1, complex_scaled_taylor_order - 2
            coefficient(n + 2) = (z0 * coefficient(n) + coefficient(n - 1)) / &
                real((n + 2) * (n + 1), dp)
        end do

        ai_next = coefficient(complex_scaled_taylor_order)
        aip_next = real(complex_scaled_taylor_order, dp) * &
            coefficient(complex_scaled_taylor_order)
        do n = complex_scaled_taylor_order - 1, 1, -1
            ai_next = coefficient(n) + step * ai_next
            aip_next = real(n, dp) * coefficient(n) + step * aip_next
        end do
        ai_next = coefficient(0) + step * ai_next
    end subroutine airy_taylor_step_complex

    pure elemental subroutine airy_ai_scaled_complex(z, eai, eaip)
        complex(dp), intent(in) :: z !! finite complex argument outside the series region
        complex(dp), intent(out) :: eai !! exp(zeta)*Ai(z), zeta=2*z*sqrt(z)/3
        complex(dp), intent(out) :: eaip !! exp(zeta)*Ai'(z), zeta=2*z*sqrt(z)/3

        complex(dp), parameter :: rot_minus = &
            cmplx(0.5_dp, -0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: rot_plus = &
            cmplx(0.5_dp, 0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: rot2_minus = &
            cmplx(-0.5_dp, -0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: rot2_plus = &
            cmplx(-0.5_dp, 0.866025403784438646763723170753_dp, kind=dp)
        complex(dp) :: eai_minus
        complex(dp) :: eai_plus
        complex(dp) :: eaip_minus
        complex(dp) :: eaip_plus
        complex(dp) :: factor_minus
        complex(dp) :: factor_plus
        complex(dp) :: u_minus
        complex(dp) :: u_plus
        complex(dp) :: w
        complex(dp) :: zeta
        complex(dp) :: zeta_minus
        complex(dp) :: zeta_plus
        real(dp) :: argument

        argument = atan2(aimag(z), real(z, dp))
        if (abs(argument) < two_thirds * scifort_pi) then
            call airy_ai_asymptotic_scaled_direct(z, eai, eaip)
            return
        end if

        ! Scale DLMF 9.2.14 before combining the rotated Ai branches.
        ! The real parts of both exponent differences are nonpositive in this sector.
        w = -z
        u_plus = w * rot_plus
        u_minus = w * rot_minus
        call airy_ai_asymptotic_scaled_direct(u_plus, eai_plus, eaip_plus)
        call airy_ai_asymptotic_scaled_direct(u_minus, eai_minus, eaip_minus)
        zeta = airy_zeta_complex(z)
        zeta_plus = airy_zeta_complex(u_plus)
        zeta_minus = airy_zeta_complex(u_minus)
        factor_plus = exp(zeta - zeta_plus)
        factor_minus = exp(zeta - zeta_minus)
        eai = rot_plus * eai_plus * factor_plus + rot_minus * eai_minus * factor_minus
        eaip = -(rot2_plus * eaip_plus * factor_plus + &
            rot2_minus * eaip_minus * factor_minus)
    end subroutine airy_ai_scaled_complex

    pure elemental subroutine airy_bi_scaled_complex(z, ebi, ebip)
        complex(dp), intent(in) :: z !! finite complex argument outside the series region
        complex(dp), intent(out) :: ebi !! exp(-abs(real(zeta)))*Bi(z)
        complex(dp), intent(out) :: ebip !! exp(-abs(real(zeta)))*Bi'(z)

        complex(dp), parameter :: omega_minus = &
            cmplx(-0.5_dp, -0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: omega_plus = &
            cmplx(-0.5_dp, 0.866025403784438646763723170753_dp, kind=dp)
        complex(dp), parameter :: phase_minus = &
            cmplx(0.866025403784438646763723170753_dp, -0.5_dp, kind=dp)
        complex(dp), parameter :: phase_plus = &
            cmplx(0.866025403784438646763723170753_dp, 0.5_dp, kind=dp)
        complex(dp) :: eai_minus
        complex(dp) :: eai_plus
        complex(dp) :: eaip_minus
        complex(dp) :: eaip_plus
        complex(dp) :: factor_minus
        complex(dp) :: factor_plus
        complex(dp) :: u_minus
        complex(dp) :: u_plus
        complex(dp) :: zeta
        complex(dp) :: zeta_minus
        complex(dp) :: zeta_plus
        real(dp) :: log_scale

        ! DLMF 9.2.10 reconstructed from exponentially scaled rotated Ai values.
        ! Combining -abs(real(zeta)) with each rotated exponent avoids overflow.
        u_minus = z * omega_minus
        u_plus = z * omega_plus
        call airy_ai_scaled_complex(u_minus, eai_minus, eaip_minus)
        call airy_ai_scaled_complex(u_plus, eai_plus, eaip_plus)
        zeta = airy_zeta_complex(z)
        zeta_minus = airy_zeta_complex(u_minus)
        zeta_plus = airy_zeta_complex(u_plus)
        log_scale = -abs(real(zeta, dp))
        factor_minus = exp(-zeta_minus + cmplx(log_scale, 0.0_dp, kind=dp))
        factor_plus = exp(-zeta_plus + cmplx(log_scale, 0.0_dp, kind=dp))
        ebi = phase_minus * eai_minus * factor_minus + &
            phase_plus * eai_plus * factor_plus
        ebip = phase_minus * omega_minus * eaip_minus * factor_minus + &
            phase_plus * omega_plus * eaip_plus * factor_plus
    end subroutine airy_bi_scaled_complex

    pure elemental subroutine airy_ai_asymptotic_scaled_direct(z, eai, eaip)
        complex(dp), intent(in) :: z !! finite large argument away from the negative-real sector
        complex(dp), intent(out) :: eai !! exp(zeta)*Ai(z) asymptotic value
        complex(dp), intent(out) :: eaip !! exp(zeta)*Ai'(z) asymptotic value

        complex(dp) :: sum_u
        complex(dp) :: sum_v
        complex(dp) :: sqrt_z
        complex(dp) :: z_quarter
        complex(dp) :: zeta

        sqrt_z = sqrt(z)
        z_quarter = sqrt(sqrt_z)
        zeta = two_thirds * z * sqrt_z
        call airy_asymptotic_sums(zeta, sum_u, sum_v)
        eai = 0.5_dp * inv_sqrt_pi * sum_u / z_quarter
        eaip = -0.5_dp * inv_sqrt_pi * z_quarter * sum_v
    end subroutine airy_ai_asymptotic_scaled_direct

    pure elemental function airy_zeta_complex(z) result(zeta)
        complex(dp), intent(in) :: z !! finite Airy argument
        complex(dp) :: zeta

        zeta = two_thirds * z * sqrt(z)
    end function airy_zeta_complex

    pure elemental subroutine airy_ai_asymptotic_direct(z, ai, aip)
        complex(dp), intent(in) :: z !! finite large argument away from the negative-real sector
        complex(dp), intent(out) :: ai !! asymptotic Airy Ai(z)
        complex(dp), intent(out) :: aip !! asymptotic derivative Ai'(z)

        complex(dp) :: exp_minus_zeta
        complex(dp) :: sum_u
        complex(dp) :: sum_v
        complex(dp) :: sqrt_z
        complex(dp) :: z_quarter
        complex(dp) :: zeta

        sqrt_z = sqrt(z)
        z_quarter = sqrt(sqrt_z)
        zeta = two_thirds * z * sqrt_z
        call airy_asymptotic_sums(zeta, sum_u, sum_v)
        exp_minus_zeta = exp(-zeta)
        ai = 0.5_dp * inv_sqrt_pi * exp_minus_zeta * sum_u / z_quarter
        aip = -0.5_dp * inv_sqrt_pi * z_quarter * exp_minus_zeta * sum_v
    end subroutine airy_ai_asymptotic_direct

    pure elemental subroutine airy_asymptotic_sums(zeta, sum_u, sum_v)
        complex(dp), intent(in) :: zeta !! nonzero Airy phase 2*z**(3/2)/3
        complex(dp), intent(out) :: sum_u !! DLMF 9.7.5 asymptotic multiplier
        complex(dp), intent(out) :: sum_v !! DLMF 9.7.6 asymptotic multiplier

        complex(dp) :: inv_power
        complex(dp) :: term_u
        complex(dp) :: term_v
        real(dp) :: coefficient_u
        real(dp) :: coefficient_v
        real(dp) :: magnitude
        real(dp) :: previous_u
        real(dp) :: previous_v
        integer :: k

        sum_u = cmplx(1.0_dp, 0.0_dp, kind=dp)
        sum_v = cmplx(1.0_dp, 0.0_dp, kind=dp)
        inv_power = cmplx(1.0_dp, 0.0_dp, kind=dp)
        coefficient_u = 1.0_dp
        previous_u = 1.0_dp
        previous_v = 1.0_dp

        do k = 1, 60
            coefficient_u = coefficient_u * &
                real((6 * k - 5) * (6 * k - 3) * (6 * k - 1), dp) / &
                real((2 * k - 1) * 216 * k, dp)
            coefficient_v = real(6 * k + 1, dp) / real(1 - 6 * k, dp) * coefficient_u
            inv_power = -inv_power / zeta
            term_u = coefficient_u * inv_power
            term_v = coefficient_v * inv_power

            magnitude = abs(term_u)
            if (k > 2 .and. magnitude > previous_u) exit
            sum_u = sum_u + term_u
            previous_u = magnitude

            magnitude = abs(term_v)
            if (k > 2 .and. magnitude > previous_v) exit
            sum_v = sum_v + term_v
            previous_v = magnitude

            if (abs(term_u) <= epsilon(1.0_dp) * max(1.0_dp, abs(sum_u)) .and. &
                abs(term_v) <= epsilon(1.0_dp) * max(1.0_dp, abs(sum_v))) exit
        end do
    end subroutine airy_asymptotic_sums

    pure elemental logical function complex_is_finite(z) result(is_finite)
        complex(dp), intent(in) :: z !! complex value to classify

        is_finite = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function complex_is_finite

    pure elemental subroutine set_all_nan_complex(kind_arg, a, b, c, d)
        complex(dp), intent(in) :: kind_arg !! complex-kind selector and NaN propagation source
        complex(dp), intent(out) :: a !! first NaN result
        complex(dp), intent(out) :: b !! second NaN result
        complex(dp), intent(out) :: c !! third NaN result
        complex(dp), intent(out) :: d !! fourth NaN result
        real(dp) :: nan_value

        nan_value = quiet_nan(real(kind_arg, dp))
        a = cmplx(nan_value, nan_value, kind=dp)
        b = a
        c = a
        d = a
    end subroutine set_all_nan_complex

    pure elemental subroutine airy_series(x, ai, aip, bi, bip)
        real(dp), intent(in) :: x !! argument with abs(x) <= 1
        real(dp), intent(out) :: ai !! Airy Ai(x)
        real(dp), intent(out) :: aip !! derivative Ai'(x)
        real(dp), intent(out) :: bi !! Airy Bi(x)
        real(dp), intent(out) :: bip !! derivative Bi'(x)

        integer :: k
        integer :: n_f
        integer :: n_g
        real(dp) :: f
        real(dp) :: fp
        real(dp) :: g
        real(dp) :: gp
        real(dp) :: term_f
        real(dp) :: term_g
        real(dp) :: x3

        f = 1.0_dp
        fp = 0.0_dp
        g = x
        gp = 1.0_dp
        term_f = 1.0_dp
        term_g = x
        x3 = x * x * x

        do k = 1, 80
            n_f = 3 * k
            n_g = 3 * k + 1
            term_f = term_f * x3 / real(n_f * (n_f - 1), dp)
            term_g = term_g * x3 / real(n_g * (n_g - 1), dp)
            f = f + term_f
            g = g + term_g
            if (x /= 0.0_dp) then
                fp = fp + real(n_f, dp) * term_f / x
                gp = gp + real(n_g, dp) * term_g / x
            end if
            if (abs(term_f) <= epsilon(1.0_dp) * max(1.0_dp, abs(f)) .and. &
                abs(term_g) <= epsilon(1.0_dp) * max(1.0_dp, abs(g))) exit
        end do

        ai = ai_zero * f + aip_zero * g
        aip = ai_zero * fp + aip_zero * gp
        bi = bi_zero * f + bip_zero * g
        bip = bi_zero * fp + bip_zero * gp
    end subroutine airy_series

    pure elemental subroutine airy_positive_scaled(x, zeta, eai, eaip, ebi, ebip)
        real(dp), intent(in) :: x !! positive real Airy argument
        real(dp), intent(in) :: zeta !! 2*x**(3/2)/3, finite and positive
        real(dp), intent(out) :: eai !! exp(zeta)*Ai(x)
        real(dp), intent(out) :: eaip !! exp(zeta)*Ai'(x)
        real(dp), intent(out) :: ebi !! exp(-zeta)*Bi(x)
        real(dp), intent(out) :: ebip !! exp(-zeta)*Bi'(x)

        real(dp) :: sqrt_x_over_three

        sqrt_x_over_three = sqrt(x * one_third)
        eai = sqrt_x_over_three * kve(one_third, zeta) / scifort_pi
        eaip = -x * inv_sqrt_three * kve(two_thirds, zeta) / scifort_pi
        ebi = sqrt_x_over_three * (ive(-one_third, zeta) + ive(one_third, zeta))
        ebip = x * inv_sqrt_three * (ive(-two_thirds, zeta) + ive(two_thirds, zeta))
    end subroutine airy_positive_scaled

    pure elemental subroutine airy_positive_scaled_asymptotic(x, eai, eaip, ebi, ebip)
        real(dp), intent(in) :: x !! positive argument whose Airy zeta overflows binary64
        real(dp), intent(out) :: eai !! asymptotic exp(zeta)*Ai(x)
        real(dp), intent(out) :: eaip !! asymptotic exp(zeta)*Ai'(x)
        real(dp), intent(out) :: ebi !! asymptotic exp(-zeta)*Bi(x)
        real(dp), intent(out) :: ebip !! asymptotic exp(-zeta)*Bi'(x)

        real(dp) :: x_quarter

        x_quarter = sqrt(sqrt(x))
        eai = 0.5_dp * inv_sqrt_pi / x_quarter
        eaip = -0.5_dp * inv_sqrt_pi * x_quarter
        ebi = inv_sqrt_pi / x_quarter
        ebip = inv_sqrt_pi * x_quarter
    end subroutine airy_positive_scaled_asymptotic

    pure elemental subroutine airy_negative_bessel(ax, zeta, ai, aip, bi, bip)
        real(dp), intent(in) :: ax !! positive magnitude of the negative Airy argument
        real(dp), intent(in) :: zeta !! 2*ax**(3/2)/3, finite and positive
        real(dp), intent(out) :: ai !! Airy Ai(-ax)
        real(dp), intent(out) :: aip !! derivative Ai'(-ax)
        real(dp), intent(out) :: bi !! Airy Bi(-ax)
        real(dp), intent(out) :: bip !! derivative Bi'(-ax)

        real(dp) :: j_m_one_third
        real(dp) :: j_m_two_thirds
        real(dp) :: j_one_third
        real(dp) :: j_two_thirds
        real(dp) :: sqrt_ax

        sqrt_ax = sqrt(ax)
        j_m_one_third = jv(-one_third, zeta)
        j_one_third = jv(one_third, zeta)
        j_m_two_thirds = jv(-two_thirds, zeta)
        j_two_thirds = jv(two_thirds, zeta)

        ai = (sqrt_ax / 3.0_dp) * (j_one_third + j_m_one_third)
        aip = (ax / 3.0_dp) * (j_two_thirds - j_m_two_thirds)
        bi = sqrt(ax * one_third) * (j_m_one_third - j_one_third)
        bip = ax * inv_sqrt_three * (j_m_two_thirds + j_two_thirds)
    end subroutine airy_negative_bessel

    pure elemental logical function zeta_would_overflow(ax, sqrt_ax) result(overflows)
        real(dp), intent(in) :: ax !! positive magnitude of the Airy argument
        real(dp), intent(in) :: sqrt_ax !! square root of ax

        if (ax <= 2.0_dp) then
            overflows = .false.
        else
            overflows = sqrt_ax > 1.5_dp * (huge(ax) / ax)
        end if
    end function zeta_would_overflow

    pure elemental function positive_scaled_exp(scale, exponent, kind_arg) result(y)
        real(dp), intent(in) :: scale !! nonnegative finite scaled function value
        real(dp), intent(in) :: exponent !! exponential correction to apply
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y
        real(dp) :: log_y
        real(dp) :: log_min_subnormal

        if (scale == 0.0_dp) then
            y = 0.0_dp
            return
        end if
        if (.not. ieee_is_finite(scale) .or. scale < 0.0_dp) then
            y = quiet_nan(kind_arg)
            return
        end if

        log_y = log(scale) + exponent
        if (log_y > log(huge(kind_arg))) then
            y = positive_infinity(kind_arg)
        else
            log_min_subnormal = log(tiny(kind_arg)) + log(epsilon(kind_arg))
            if (log_y < log_min_subnormal) then
                y = 0.0_dp
            else
                y = exp(log_y)
            end if
        end if
    end function positive_scaled_exp

    pure elemental subroutine itairy(x, apt, bpt, ant, bnt)
        real(dp), intent(in) :: x !! real upper integration limit
        real(dp), intent(out) :: apt !! integral of Ai(t) from 0 to x
        real(dp), intent(out) :: bpt !! integral of Bi(t) from 0 to x
        real(dp), intent(out) :: ant !! integral of Ai(-t) from 0 to x
        real(dp), intent(out) :: bnt !! integral of Bi(-t) from 0 to x

        real(dp) :: a
        real(dp) :: pa
        real(dp) :: pb
        real(dp) :: pn_a
        real(dp) :: pn_b

        if (ieee_is_nan(x)) then
            call set_all_nan(x, apt, bpt, ant, bnt)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                apt = one_third
                bpt = positive_infinity(x)
                ant = two_thirds
                bnt = 0.0_dp
            else
                apt = -two_thirds
                bpt = -0.0_dp
                ant = -one_third
                bnt = -positive_infinity(x)
            end if
            return
        end if

        a = abs(x)
        call itairy_positive(a, pa, pb, pn_a, pn_b)
        if (x >= 0.0_dp) then
            apt = pa
            bpt = pb
            ant = pn_a
            bnt = pn_b
        else
            apt = -pn_a
            bpt = -pn_b
            ant = -pa
            bnt = -pb
        end if
    end subroutine itairy

    pure subroutine ai_zeros(nt, a, ap, ai, aip, status)
        integer, intent(in) :: nt !! number of zeros requested, >= 1
        real(dp), allocatable, intent(out) :: a(:) !! first nt zeros of Ai(x)
        real(dp), allocatable, intent(out) :: ap(:) !! first nt zeros of Ai'(x)
        real(dp), allocatable, intent(out) :: ai(:) !! Ai evaluated at the Ai' zeros
        real(dp), allocatable, intent(out) :: aip(:) !! Ai' evaluated at the Ai zeros
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input
                                           !! or numerical failure

        integer :: info
        integer :: k
        real(dp) :: bi
        real(dp) :: bip
        real(dp) :: unused

        if (nt < 1) then
            allocate(a(0), ap(0), ai(0), aip(0))
            if (present(status)) status = airy_status_invalid_argument
            return
        end if

        allocate(a(nt), ap(nt), ai(nt), aip(nt))
        do k = 1, nt
            a(k) = airy_zero_seed(k, 1)
            call refine_airy_zero(1, a(k), info)
            if (info /= airy_status_success) then
                call set_zero_failure(a, ap, ai, aip, status)
                return
            end if
            call airy(a(k), unused, aip(k), bi, bip)

            ap(k) = airy_zero_seed(k, 2)
            call refine_airy_zero(2, ap(k), info)
            if (info /= airy_status_success) then
                call set_zero_failure(a, ap, ai, aip, status)
                return
            end if
            call airy(ap(k), ai(k), unused, bi, bip)
        end do
        if (present(status)) status = airy_status_success
    end subroutine ai_zeros

    pure subroutine bi_zeros(nt, b, bp, bi, bip, status)
        integer, intent(in) :: nt !! number of zeros requested, >= 1
        real(dp), allocatable, intent(out) :: b(:) !! first nt zeros of Bi(x)
        real(dp), allocatable, intent(out) :: bp(:) !! first nt zeros of Bi'(x)
        real(dp), allocatable, intent(out) :: bi(:) !! Bi evaluated at the Bi' zeros
        real(dp), allocatable, intent(out) :: bip(:) !! Bi' evaluated at the Bi zeros
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input
                                           !! or numerical failure

        integer :: info
        integer :: k
        real(dp) :: ai
        real(dp) :: aip
        real(dp) :: unused

        if (nt < 1) then
            allocate(b(0), bp(0), bi(0), bip(0))
            if (present(status)) status = airy_status_invalid_argument
            return
        end if

        allocate(b(nt), bp(nt), bi(nt), bip(nt))
        do k = 1, nt
            b(k) = airy_zero_seed(k, 3)
            call refine_airy_zero(3, b(k), info)
            if (info /= airy_status_success) then
                call set_zero_failure(b, bp, bi, bip, status)
                return
            end if
            call airy(b(k), ai, aip, unused, bip(k))

            bp(k) = airy_zero_seed(k, 4)
            call refine_airy_zero(4, bp(k), info)
            if (info /= airy_status_success) then
                call set_zero_failure(b, bp, bi, bip, status)
                return
            end if
            call airy(bp(k), ai, aip, bi(k), unused)
        end do
        if (present(status)) status = airy_status_success
    end subroutine bi_zeros

    pure subroutine itairy_positive(x, apt, bpt, ant, bnt)
        real(dp), intent(in) :: x !! nonnegative upper integration limit
        real(dp), intent(out) :: apt !! integral of Ai(t) from 0 to x
        real(dp), intent(out) :: bpt !! integral of Bi(t) from 0 to x
        real(dp), intent(out) :: ant !! integral of Ai(-t) from 0 to x
        real(dp), intent(out) :: bnt !! integral of Bi(-t) from 0 to x

        real(dp) :: series_apt
        real(dp) :: series_bpt
        real(dp) :: series_ant
        real(dp) :: series_bnt

        if (x <= 6.0_dp) then
            call itairy_series(x, apt, bpt, ant, bnt)
            return
        end if

        if (x <= 15.0_dp) then
            call itairy_series(x, series_apt, series_bpt, series_ant, series_bnt)
            bpt = series_bpt
        else
            call itairy_positive_bi_asymptotic(x, bpt)
        end if

        if (x <= 8.0_dp) then
            call itairy_moderate(x, apt, ant, bnt)
        else
            call itairy_positive_ai_asymptotic(x, apt)
            if (x <= 12.0_dp) then
                call itairy_moderate(x, series_apt, ant, bnt)
            else
                call itairy_oscillatory_asymptotic(x, ant, bnt)
            end if
        end if
    end subroutine itairy_positive

    pure subroutine itairy_series(x, apt, bpt, ant, bnt)
        real(dp), intent(in) :: x !! nonnegative argument in the convergent-series regime
        real(dp), intent(out) :: apt !! integral of Ai(t) from 0 to x
        real(dp), intent(out) :: bpt !! integral of Bi(t) from 0 to x
        real(dp), intent(out) :: ant !! integral of Ai(-t) from 0 to x
        real(dp), intent(out) :: bnt !! integral of Bi(-t) from 0 to x

        real(dp) :: f_int
        real(dp) :: fm_int
        real(dp) :: g_int
        real(dp) :: gm_int

        call airy_basis_integrals(x, f_int, g_int)
        call airy_basis_integrals(-x, fm_int, gm_int)
        apt = ai_zero * f_int + aip_zero * g_int
        bpt = bi_zero * f_int + bip_zero * g_int
        ant = -(ai_zero * fm_int + aip_zero * gm_int)
        bnt = -(bi_zero * fm_int + bip_zero * gm_int)
    end subroutine itairy_series

    pure subroutine airy_basis_integrals(x, f_int, g_int)
        real(dp), intent(in) :: x !! real series endpoint
        real(dp), intent(out) :: f_int !! integral of the unit-value Airy basis solution
        real(dp), intent(out) :: g_int !! integral of the unit-slope Airy basis solution

        integer :: k
        integer :: n_f
        integer :: n_g
        real(dp) :: add_f
        real(dp) :: add_g
        real(dp) :: term_f
        real(dp) :: term_g
        real(dp) :: x3

        term_f = 1.0_dp
        term_g = x
        f_int = x
        g_int = 0.5_dp * x * x
        x3 = x * x * x
        do k = 1, 160
            n_f = 3 * k
            n_g = 3 * k + 1
            term_f = term_f * x3 / real(n_f * (n_f - 1), dp)
            term_g = term_g * x3 / real(n_g * (n_g - 1), dp)
            add_f = term_f * x / real(n_f + 1, dp)
            add_g = term_g * x / real(n_g + 1, dp)
            f_int = f_int + add_f
            g_int = g_int + add_g
            if (abs(add_f) <= epsilon(1.0_dp) * max(1.0_dp, abs(f_int)) .and. &
                abs(add_g) <= epsilon(1.0_dp) * max(1.0_dp, abs(g_int))) exit
        end do
    end subroutine airy_basis_integrals

    pure subroutine itairy_moderate(x, apt, ant, bnt)
        real(dp), intent(in) :: x !! endpoint in [6, 12]
        real(dp), intent(out) :: apt !! Ai integral obtained by quadrature from the x=6 anchor
        real(dp), intent(out) :: ant !! Ai(-t) integral obtained by quadrature from x=6
        real(dp), intent(out) :: bnt !! Bi(-t) integral obtained by quadrature from x=6

        integer :: i
        integer :: j
        integer :: nseg
        real(dp) :: a
        real(dp) :: ai1
        real(dp) :: ai2
        real(dp) :: aim1
        real(dp) :: aim2
        real(dp) :: aip
        real(dp) :: b
        real(dp) :: bi
        real(dp) :: bim1
        real(dp) :: bim2
        real(dp) :: bip
        real(dp) :: bpt_anchor
        real(dp) :: half
        real(dp) :: mid
        real(dp) :: t1
        real(dp) :: t2

        call itairy_series(6.0_dp, apt, bpt_anchor, ant, bnt)
        if (x <= 6.0_dp) return
        nseg = max(1, ceiling(2.0_dp * (x - 6.0_dp)))
        do j = 1, nseg
            a = 6.0_dp + (x - 6.0_dp) * real(j - 1, dp) / real(nseg, dp)
            b = 6.0_dp + (x - 6.0_dp) * real(j, dp) / real(nseg, dp)
            mid = 0.5_dp * (a + b)
            half = 0.5_dp * (b - a)
            do i = 1, 8
                t1 = mid - half * gl_x(i)
                t2 = mid + half * gl_x(i)
                call airy(t1, ai1, aip, bi, bip)
                call airy(t2, ai2, aip, bi, bip)
                call airy(-t1, aim1, aip, bim1, bip)
                call airy(-t2, aim2, aip, bim2, bip)
                apt = apt + half * gl_w(i) * (ai1 + ai2)
                ant = ant + half * gl_w(i) * (aim1 + aim2)
                bnt = bnt + half * gl_w(i) * (bim1 + bim2)
            end do
        end do
    end subroutine itairy_moderate

    pure subroutine itairy_positive_ai_asymptotic(x, apt)
        real(dp), intent(in) :: x !! positive endpoint in the asymptotic regime
        real(dp), intent(out) :: apt !! integral of Ai(t) from 0 to x

        real(dp) :: a
        real(dp) :: ai
        real(dp) :: aip
        real(dp) :: b
        real(dp) :: bi
        real(dp) :: bip

        call airy_integral_coefficients(x, .false., a, b)
        call airy(x, ai, aip, bi, bip)
        apt = one_third + a * ai + b * aip
    end subroutine itairy_positive_ai_asymptotic

    pure subroutine itairy_positive_bi_asymptotic(x, bpt)
        real(dp), intent(in) :: x !! positive endpoint in the asymptotic regime
        real(dp), intent(out) :: bpt !! integral of Bi(t) from 0 to x

        real(dp) :: a
        real(dp) :: b
        real(dp) :: eai
        real(dp) :: eaip
        real(dp) :: ebi
        real(dp) :: ebip
        real(dp) :: scale
        real(dp) :: sqrt_x
        real(dp) :: zeta

        sqrt_x = sqrt(x)
        if (zeta_would_overflow(x, sqrt_x)) then
            bpt = positive_infinity(x)
            return
        end if
        zeta = two_thirds * x * sqrt_x
        call airy_integral_coefficients(x, .false., a, b)
        call airye(x, eai, eaip, ebi, ebip)
        scale = a * ebi + b * ebip
        bpt = positive_scaled_exp(scale, zeta, x)
    end subroutine itairy_positive_bi_asymptotic

    pure subroutine itairy_oscillatory_asymptotic(x, ant, bnt)
        real(dp), intent(in) :: x !! positive endpoint in the oscillatory asymptotic regime
        real(dp), intent(out) :: ant !! integral of Ai(-t) from 0 to x
        real(dp), intent(out) :: bnt !! integral of Bi(-t) from 0 to x

        real(dp) :: a
        real(dp) :: ai
        real(dp) :: aip
        real(dp) :: b
        real(dp) :: bi
        real(dp) :: bip

        if (x >= oscillatory_limit_x) then
            ant = two_thirds
            bnt = 0.0_dp
            return
        end if
        call airy_integral_coefficients(x, .true., a, b)
        call airy(-x, ai, aip, bi, bip)
        ant = two_thirds + a * ai - b * aip
        bnt = a * bi - b * bip
    end subroutine itairy_oscillatory_asymptotic

    pure subroutine airy_integral_coefficients(x, oscillatory, a, b)
        real(dp), intent(in) :: x !! positive asymptotic endpoint
        logical, intent(in) :: oscillatory !! true for y''=-x*y; false for y''=x*y
        real(dp), intent(out) :: a !! coefficient multiplying the function value
        real(dp), intent(out) :: b !! coefficient multiplying the derivative

        integer :: k
        real(dp) :: next_term
        real(dp) :: ratio
        real(dp) :: sign_factor
        real(dp) :: term_b

        if (oscillatory) then
            sign_factor = -1.0_dp
            term_b = -1.0_dp / x
        else
            sign_factor = 1.0_dp
            term_b = 1.0_dp / x
        end if
        b = term_b
        a = term_b / x
        do k = 1, 80
            ratio = real(3 * k - 2, dp) / x
            ratio = ratio * real(3 * k - 1, dp) / x
            ratio = ratio / x
            next_term = sign_factor * term_b * ratio
            if (k > 1 .and. abs(next_term) >= abs(term_b)) exit
            term_b = next_term
            b = b + term_b
            a = a + real(3 * k + 1, dp) * term_b / x
            if (abs(term_b) <= epsilon(1.0_dp) * max(1.0_dp, abs(b))) exit
        end do
    end subroutine airy_integral_coefficients

    pure elemental function airy_zero_seed(k, kind) result(x)
        integer, intent(in) :: k !! one-based Airy zero index
        integer, intent(in) :: kind !! 1=Ai, 2=Ai', 3=Bi, 4=Bi'
        real(dp) :: x
        real(dp) :: shift

        select case (kind)
        case (1, 4)
            shift = real(k, dp) - 0.25_dp
        case default
            shift = real(k, dp) - 0.75_dp
        end select
        x = -(1.5_dp * scifort_pi * shift)**two_thirds
    end function airy_zero_seed

    pure subroutine refine_airy_zero(kind, x, status)
        integer, intent(in) :: kind !! 1=Ai, 2=Ai', 3=Bi, 4=Bi'
        real(dp), intent(inout) :: x !! initial root estimate replaced by refined root
        integer, intent(out) :: status !! zero on convergence; nonzero on numerical failure

        integer :: iter
        real(dp) :: ai
        real(dp) :: aip
        real(dp) :: bi
        real(dp) :: bip
        real(dp) :: df
        real(dp) :: f
        real(dp) :: step

        do iter = 1, 24
            call airy(x, ai, aip, bi, bip)
            select case (kind)
            case (1)
                f = ai
                df = aip
            case (2)
                f = aip
                df = x * ai
            case (3)
                f = bi
                df = bip
            case default
                f = bip
                df = x * bi
            end select
            if (.not. ieee_is_finite(f) .or. .not. ieee_is_finite(df) .or. df == 0.0_dp) then
                status = airy_status_numerical_failure
                return
            end if
            step = f / df
            if (.not. ieee_is_finite(step)) then
                status = airy_status_numerical_failure
                return
            end if
            x = x - step
            if (abs(step) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x))) then
                status = airy_status_success
                return
            end if
        end do
        status = airy_status_numerical_failure
    end subroutine refine_airy_zero

    pure subroutine set_zero_failure(a, ap, value, deriv, status)
        real(dp), intent(inout) :: a(:) !! primary zero array to invalidate
        real(dp), intent(inout) :: ap(:) !! derivative-zero array to invalidate
        real(dp), intent(inout) :: value(:) !! function-value array to invalidate
        real(dp), intent(inout) :: deriv(:) !! derivative-value array to invalidate
        integer, intent(out), optional :: status !! receives numerical-failure status when present

        a = quiet_nan(0.0_dp)
        ap = quiet_nan(0.0_dp)
        value = quiet_nan(0.0_dp)
        deriv = quiet_nan(0.0_dp)
        if (present(status)) status = airy_status_numerical_failure
    end subroutine set_zero_failure

    pure elemental subroutine set_all_nan(kind_arg, a, b, c, d)
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp), intent(out) :: a !! first NaN result
        real(dp), intent(out) :: b !! second NaN result
        real(dp), intent(out) :: c !! third NaN result
        real(dp), intent(out) :: d !! fourth NaN result

        a = quiet_nan(kind_arg)
        b = quiet_nan(kind_arg)
        c = quiet_nan(kind_arg)
        d = quiet_nan(kind_arg)
    end subroutine set_all_nan

end module scifort_airy
