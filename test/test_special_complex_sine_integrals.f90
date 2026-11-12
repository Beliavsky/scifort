! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_complex_sine_integrals
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_special, only : exp1, expi, shichi, sici
    implicit none

    complex(dp) :: chi
    complex(dp) :: chi_array(2)
    complex(dp) :: ci
    complex(dp) :: ci_array(2)
    complex(dp) :: input_values(2)
    complex(dp) :: shi
    complex(dp) :: shi_array(2)
    complex(dp) :: si
    complex(dp) :: si_array(2)
    complex(dp) :: value
    complex(dp) :: value2
    real(dp) :: chi_real
    real(dp) :: ci_real
    real(dp) :: nan_value
    real(dp) :: shi_real
    real(dp) :: si_real

    call sici(cmplx(2.5_dp, 3.0_dp, kind=dp), si, ci)
    call check_complex_rel('sici ordinary si', si, &
        cmplx(4.5057358745639515_dp, 0.06863305018999655_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('sici ordinary ci', ci, &
        cmplx(0.07936442069069671_dp, -2.9355102629375414_dp, kind=dp), 3.0e-14_dp)

    call sici(cmplx(-2.0_dp, 3.0_dp, kind=dp), si, ci)
    call check_complex_rel('sici left upper si', si, &
        cmplx(-4.547513889562289_dp, 1.399196580646057_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('sici left upper ci', ci, &
        cmplx(1.4082925015208514_dp, 6.125210395619398_dp, kind=dp), 3.0e-14_dp)

    call sici(cmplx(1.0e-12_dp, 2.0e-12_dp, kind=dp), si, ci)
    call check_complex_rel('sici tiny si', si, &
        cmplx(1.0e-12_dp, 2.0e-12_dp, kind=dp), 2.0e-15_dp)
    call check_complex_rel('sici tiny ci', ci, &
        cmplx(-26.249086494809966_dp, 1.1071487177940904_dp, kind=dp), 3.0e-14_dp)

    call sici(cmplx(20.0_dp, 15.0_dp, kind=dp), si, ci)
    call check_complex_rel('sici large si', si, &
        cmplx(12556.904034618623_dp, 65715.00783089481_dp, kind=dp), 5.0e-14_dp)
    call check_complex_rel('sici large ci', ci, &
        cmplx(65715.00783090029_dp, -12555.333238302434_dp, kind=dp), 5.0e-14_dp)

    call shichi(cmplx(0.5_dp, 2.5_dp, kind=dp), shi, chi)
    call check_complex_rel('shichi ordinary shi', shi, &
        cmplx(0.11772029666668204_dp, 1.831091777729851_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('shichi ordinary chi', chi, &
        cmplx(0.29912435887648814_dp, 1.7395351121166562_dp, kind=dp), 3.0e-14_dp)

    call shichi(cmplx(-2.0_dp, -3.0_dp, kind=dp), shi, chi)
    call check_complex_rel('shichi left lower shi', shi, &
        cmplx(0.19318907627191972_dp, -2.645432555362369_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('shichi left lower chi', chi, &
        cmplx(-0.16836286832772052_dp, -0.5164767731384683_dp, kind=dp), 3.0e-14_dp)

    call shichi(cmplx(20.0_dp, 15.0_dp, kind=dp), shi, chi)
    call check_complex_rel('shichi large shi', shi, &
        cmplx(-1920173.2355022295_dp, 9838464.485123662_dp, kind=dp), 6.0e-14_dp)
    call check_complex_rel('shichi large chi', chi, &
        cmplx(-1920173.2355022295_dp, 9838464.485123662_dp, kind=dp), 6.0e-14_dp)

    call sici(cmplx(-2.0_dp, 0.0_dp, kind=dp), si, ci)
    call check_close('sici cut si real', real(si, dp), -1.6054129768026946_dp, 2.0e-15_dp)
    call check_close('sici cut ci real', real(ci, dp), 0.42298082877486476_dp, 2.0e-15_dp)
    call check_close('sici cut ci imag', aimag(ci), scifort_pi, 2.0e-15_dp)
    call sici(cmplx(-2.0_dp, -0.0_dp, kind=dp), si, ci)
    call check_close('sici cut minus-zero ci imag', aimag(ci), scifort_pi, 2.0e-15_dp)

    call shichi(cmplx(-2.0_dp, 0.0_dp, kind=dp), shi, chi)
    call check_close('shichi cut shi real', real(shi, dp), -2.5015674333549756_dp, 3.0e-15_dp)
    call check_close('shichi cut chi real', real(chi, dp), 2.452666922646915_dp, 3.0e-15_dp)
    call check_close('shichi cut chi imag', aimag(chi), scifort_pi, 2.0e-15_dp)

    call sici(cmplx(-2.0_dp, 1.0e-12_dp, kind=dp), si, ci)
    call check_close('sici upper near cut', aimag(ci), 3.1415926535900014_dp, 3.0e-14_dp)
    call sici(cmplx(-2.0_dp, -1.0e-12_dp, kind=dp), si, ci)
    call check_close('sici lower near cut', aimag(ci), -3.1415926535900014_dp, 3.0e-14_dp)

    value = cmplx(2.0_dp, 3.0_dp, kind=dp)
    call sici(value, si, ci)
    value2 = (exp1(cmplx(-aimag(value), real(value, dp), kind=dp)) - &
        exp1(cmplx(aimag(value), -real(value, dp), kind=dp))) / &
        cmplx(0.0_dp, 2.0_dp, kind=dp) + 0.5_dp * scifort_pi
    call check_complex_rel('sici E1 identity si', si, value2, 3.0e-14_dp)

    call shichi(value, shi, chi)
    value2 = 0.5_dp * (expi(value) - expi(-value)) - &
        cmplx(0.0_dp, 0.5_dp * scifort_pi, kind=dp)
    call check_complex_rel('shichi Ei identity shi', shi, value2, 3.0e-14_dp)

    call sici(2.5_dp, si_real, ci_real)
    call sici(cmplx(2.5_dp, 0.0_dp, kind=dp), si, ci)
    call check_close('sici exact positive real si parity', real(si, dp), si_real, 0.0_dp)
    call check_close('sici exact positive real ci parity', real(ci, dp), ci_real, 0.0_dp)
    call check_close('sici exact positive real imag', aimag(ci), 0.0_dp, 0.0_dp)

    call shichi(2.5_dp, shi_real, chi_real)
    call shichi(cmplx(2.5_dp, 0.0_dp, kind=dp), shi, chi)
    call check_close('shichi exact positive real shi parity', real(shi, dp), shi_real, 0.0_dp)
    call check_close('shichi exact positive real chi parity', real(chi, dp), chi_real, 0.0_dp)

    input_values = [cmplx(2.5_dp, 3.0_dp, kind=dp), cmplx(-2.0_dp, -3.0_dp, kind=dp)]
    call sici(input_values, si_array, ci_array)
    call check_complex_rel('sici elemental first', si_array(1), &
        cmplx(4.5057358745639515_dp, 0.06863305018999655_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('sici elemental second', ci_array(2), &
        cmplx(1.4082925015208514_dp, -6.125210395619398_dp, kind=dp), 3.0e-14_dp)

    call shichi(input_values, shi_array, chi_array)
    call check_complex_rel('shichi elemental first', shi_array(1), &
        cmplx(-0.49620226984975657_dp, 3.197655178639338_dp, kind=dp), 4.0e-14_dp)
    call check_complex_rel('shichi elemental second', chi_array(2), &
        cmplx(-0.16836286832772052_dp, -0.5164767731384683_dp, kind=dp), 3.0e-14_dp)

    call sici(cmplx(0.0_dp, 2.0_dp, kind=dp), si, ci)
    call check_complex_rel('sici imaginary axis si', si, &
        cmplx(0.0_dp, 2.5015674333549756_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('sici imaginary axis ci', ci, &
        cmplx(2.452666922646915_dp, 0.5_dp * scifort_pi, kind=dp), 3.0e-14_dp)

    call shichi(cmplx(0.0_dp, -2.0_dp, kind=dp), shi, chi)
    call check_complex_rel('shichi imaginary axis shi', shi, &
        cmplx(0.0_dp, -1.6054129768026946_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('shichi imaginary axis chi', chi, &
        cmplx(0.42298082877486476_dp, -0.5_dp * scifort_pi, kind=dp), 3.0e-14_dp)

    call sici(cmplx(0.0_dp, 0.0_dp, kind=dp), si, ci)
    call check_true('sici zero ci imag nan', ieee_is_nan(aimag(ci)))
    call shichi(cmplx(0.0_dp, 0.0_dp, kind=dp), shi, chi)
    call check_true('shichi zero chi imag nan', ieee_is_nan(aimag(chi)))

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    call sici(cmplx(nan_value, 1.0_dp, kind=dp), si, ci)
    call check_true('sici nan si real', ieee_is_nan(real(si, dp)))
    call check_true('sici nan ci imag', ieee_is_nan(aimag(ci)))
    call shichi(cmplx(1.0_dp, nan_value, kind=dp), shi, chi)
    call check_true('shichi nan shi real', ieee_is_nan(real(shi, dp)))
    call check_true('shichi nan chi imag', ieee_is_nan(aimag(chi)))

    print *, 'special complex sine/hyperbolic-integral tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed scalar value
        real(dp), intent(in) :: expected !! reference scalar value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_complex_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum scaled complex error

        if (actual == expected) return
        if (ieee_is_nan(real(actual, dp)) .or. ieee_is_nan(aimag(actual)) .or. &
                abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' scaled tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_rel

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_complex_sine_integrals
