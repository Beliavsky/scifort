! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_complex_exponential_integrals
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_special, only : exp1, expi
    implicit none

    complex(dp) :: actual
    complex(dp) :: array_values(3)
    complex(dp) :: input_values(3)
    complex(dp) :: lower_cut
    complex(dp) :: upper_cut
    complex(dp) :: value
    real(dp) :: nan_value

    call check_complex_rel('exp1 ordinary', exp1(cmplx(0.5_dp, 0.25_dp, kind=dp)), &
        cmplx(0.45944998817799976_dp, -0.26751135083437470_dp, kind=dp), 2.0e-14_dp)
    call check_complex_rel('exp1 left half-plane', exp1(cmplx(-1.0_dp, 0.5_dp, kind=dp)), &
        cmplx(-1.8831483186261848_dp, -1.8336478949049033_dp, kind=dp), 2.0e-14_dp)
    call check_complex_rel('exp1 continued fraction', exp1(cmplx(10.0_dp, 10.0_dp, kind=dp)), &
        cmplx(-7.669897841555350e-7_dp, 2.9581458278235625e-6_dp, kind=dp), 2.0e-14_dp)
    call check_complex_rel('exp1 near cut series', exp1(cmplx(-20.0_dp, 1.0e-8_dp, kind=dp)), &
        cmplx(-2.5615652664056595e7_dp, -2.899010055884898_dp, kind=dp), 3.0e-14_dp)

    call check_complex_rel('expi ordinary', expi(cmplx(0.5_dp, 0.25_dp, kind=dp)), &
        cmplx(0.5438975725612261_dp, 0.7867419334490524_dp, kind=dp), 2.0e-14_dp)
    call check_complex_rel('expi lower half-plane', expi(cmplx(3.0_dp, -2.0_dp, kind=dp)), &
        cmplx(2.8074890821669083_dp, -9.101927958386790_dp, kind=dp), 2.0e-14_dp)
    call check_complex_rel('expi near positive axis', expi(cmplx(10.0_dp, 1.0e-12_dp, kind=dp)), &
        cmplx(2492.2289762418786_dp, 2.2026465096303127e-9_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('expi near negative cut', expi(cmplx(-10.0_dp, 0.1_dp, kind=dp)), &
        cmplx(-4.132024797264409e-6_dp, 3.1415922005130423_dp, kind=dp), 2.0e-14_dp)

    upper_cut = exp1(cmplx(-1.0_dp, 0.0_dp, kind=dp))
    lower_cut = exp1(cmplx(-1.0_dp, -0.0_dp, kind=dp))
    call check_close('exp1 cut real +0', real(upper_cut, dp), -1.895117816355937_dp, 2.0e-15_dp)
    call check_close('exp1 cut imag +0', aimag(upper_cut), -scifort_pi, 2.0e-15_dp)
    call check_close('exp1 cut real -0', real(lower_cut, dp), -1.895117816355937_dp, 2.0e-15_dp)
    call check_close('exp1 cut imag -0', aimag(lower_cut), scifort_pi, 2.0e-15_dp)

    upper_cut = expi(cmplx(-1.0_dp, 0.0_dp, kind=dp))
    lower_cut = expi(cmplx(-1.0_dp, -0.0_dp, kind=dp))
    call check_close('expi cut real +0', real(upper_cut, dp), -0.2193839343955205_dp, 2.0e-15_dp)
    call check_close('expi cut imag +0', aimag(upper_cut), -0.0_dp, 0.0_dp)
    call check_close('expi cut real -0', real(lower_cut, dp), -0.2193839343955205_dp, 2.0e-15_dp)
    call check_close('expi cut imag -0', aimag(lower_cut), -0.0_dp, 0.0_dp)
    call check_true('expi cut signed zero +0', sign(1.0_dp, aimag(upper_cut)) < 0.0_dp)
    call check_true('expi cut signed zero -0', sign(1.0_dp, aimag(lower_cut)) < 0.0_dp)

    value = exp1(cmplx(2.0_dp, 3.0_dp, kind=dp))
    call check_complex_rel('exp1 conjugation', exp1(cmplx(2.0_dp, -3.0_dp, kind=dp)), &
        conjg(value), 3.0e-14_dp)
    value = expi(cmplx(-2.0_dp, 3.0_dp, kind=dp))
    call check_complex_rel('expi conjugation', expi(cmplx(-2.0_dp, -3.0_dp, kind=dp)), &
        conjg(value), 3.0e-14_dp)

    value = cmplx(2.0_dp, 3.0_dp, kind=dp)
    actual = expi(value) + exp1(-value)
    call check_complex_close('Ei-E1 relation upper', actual, &
        cmplx(0.0_dp, scifort_pi, kind=dp), 2.0e-14_dp)
    value = cmplx(2.0_dp, -3.0_dp, kind=dp)
    actual = expi(value) + exp1(-value)
    call check_complex_close('Ei-E1 relation lower', actual, &
        cmplx(0.0_dp, -scifort_pi, kind=dp), 2.0e-14_dp)

    call check_complex_rel('exp1 exact real dispatch', exp1(cmplx(0.5_dp, 0.0_dp, kind=dp)), &
        cmplx(exp1(0.5_dp), 0.0_dp, kind=dp), 0.0_dp)
    call check_complex_rel('expi exact real dispatch', expi(cmplx(10.0_dp, 0.0_dp, kind=dp)), &
        cmplx(expi(10.0_dp), 0.0_dp, kind=dp), 0.0_dp)

    input_values = [cmplx(0.5_dp, 0.25_dp, kind=dp), &
        cmplx(-1.0_dp, 0.5_dp, kind=dp), cmplx(3.0_dp, -2.0_dp, kind=dp)]
    array_values = exp1(input_values)
    call check_complex_rel('exp1 elemental array first', array_values(1), &
        cmplx(0.45944998817799976_dp, -0.26751135083437470_dp, kind=dp), 2.0e-14_dp)
    array_values = expi(input_values)
    call check_complex_rel('expi elemental array last', array_values(3), &
        cmplx(2.8074890821669083_dp, -9.101927958386790_dp, kind=dp), 2.0e-14_dp)

    value = exp1(cmplx(-700.0_dp, 1.0_dp, kind=dp))
    call check_true('exp1 large finite real part', abs(real(value, dp)) < huge(1.0_dp))
    call check_true('exp1 large finite imag part', abs(aimag(value)) < huge(1.0_dp))

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    value = exp1(cmplx(nan_value, 1.0_dp, kind=dp))
    call check_true('exp1 complex nan real', ieee_is_nan(real(value, dp)))
    call check_true('exp1 complex nan imag', ieee_is_nan(aimag(value)))
    value = expi(cmplx(1.0_dp, nan_value, kind=dp))
    call check_true('expi complex nan real', ieee_is_nan(real(value, dp)))
    call check_true('expi complex nan imag', ieee_is_nan(aimag(value)))

    print *, 'special complex exponential-integral tests passed'

contains

    subroutine check_close(name, actual_value, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual_value !! computed scalar value
        real(dp), intent(in) :: expected !! reference scalar value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual_value == expected) return
        if (ieee_is_nan(actual_value) .or. abs(actual_value - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual_value
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_complex_close(name, actual_value, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual_value !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum absolute complex error

        if (actual_value == expected) return
        if (ieee_is_nan(real(actual_value, dp)) .or. ieee_is_nan(aimag(actual_value)) .or. &
                abs(actual_value - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual_value
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_close

    subroutine check_complex_rel(name, actual_value, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual_value !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum scaled error

        if (actual_value == expected) return
        if (ieee_is_nan(real(actual_value, dp)) .or. ieee_is_nan(aimag(actual_value)) .or. &
                abs(actual_value - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual_value
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

end program test_special_complex_exponential_integrals
