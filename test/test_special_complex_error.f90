! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_complex_error
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_special, only : dawsn, dp, erf, erfc, erfi
    implicit none

    complex(dp) :: array_actual(3)
    complex(dp) :: array_input(3)
    complex(dp) :: z
    real(dp) :: nan_value

    call check_complex('erf ordinary', erf(cmplx(0.25_dp, 0.5_dp, kind=dp)), &
        cmplx(0.3511866593813804_dp, 0.5712866223943341_dp, kind=dp), 3.0e-14_dp)
    call check_complex('erfc ordinary', erfc(cmplx(-1.0_dp, 0.75_dp, kind=dp)), &
        cmplx(2.0967142233706118_dp, -0.23437075851891276_dp, kind=dp), 3.0e-14_dp)
    call check_complex('erfi ordinary', erfi(cmplx(2.0_dp, -1.0_dp, kind=dp)), &
        cmplx(-5.0491437034470374_dp, 0.5366435657785664_dp, kind=dp), 5.0e-14_dp)
    call check_complex('dawsn ordinary', dawsn(cmplx(0.25_dp, 0.5_dp, kind=dp)), &
        cmplx(0.37504331414825715_dp, 0.509600529121811_dp, kind=dp), 3.0e-14_dp)

    call check_complex('erf imaginary axis', erf(cmplx(0.0_dp, 1.5_dp, kind=dp)), &
        cmplx(0.0_dp, 4.5847332572844275_dp, kind=dp), 3.0e-14_dp)
    call check_complex('erfi imaginary axis', erfi(cmplx(0.0_dp, 1.5_dp, kind=dp)), &
        cmplx(0.0_dp, 0.9661051464753108_dp, kind=dp), 3.0e-14_dp)
    call check_complex('dawsn imaginary axis', dawsn(cmplx(0.0_dp, 1.5_dp, kind=dp)), &
        cmplx(0.0_dp, 8.123289304869214_dp, kind=dp), 5.0e-14_dp)

    z = cmplx(1.25_dp, -0.75_dp, kind=dp)
    call check_complex('erf+erfc identity', erf(z) + erfc(z), &
        cmplx(1.0_dp, 0.0_dp, kind=dp), 4.0e-14_dp)
    call check_complex('erf odd symmetry', erf(-z), -erf(z), 4.0e-14_dp)
    call check_complex('erfi rotation identity', erfi(z), &
        cmplx(0.0_dp, -1.0_dp, kind=dp) * erf(cmplx(0.0_dp, 1.0_dp, kind=dp) * z), &
        4.0e-14_dp)
    call check_complex('dawsn odd symmetry', dawsn(-z), -dawsn(z), 4.0e-14_dp)

    call check_close('real erf preserved', erf(0.5_dp), 0.5204998778130465_dp, 4.0e-16_dp)
    call check_close('real erfc preserved', erfc(0.5_dp), 0.4795001221869535_dp, 4.0e-16_dp)
    call check_close('real erfi preserved', erfi(0.5_dp), 0.6149520946965109_dp, 2.0e-15_dp)
    call check_close('real dawsn preserved', dawsn(0.5_dp), 0.4244363835020223_dp, 2.0e-15_dp)

    array_input = [cmplx(0.1_dp, 0.2_dp, kind=dp), &
        cmplx(-1.0_dp, 0.4_dp, kind=dp), cmplx(2.0_dp, -0.5_dp, kind=dp)]
    array_actual = erf(array_input)
    call check_complex('erf elemental 1', array_actual(1), erf(array_input(1)), 0.0_dp)
    call check_complex('erf elemental 2', array_actual(2), erf(array_input(2)), 0.0_dp)
    call check_complex('erf elemental 3', array_actual(3), erf(array_input(3)), 0.0_dp)

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    z = erf(cmplx(nan_value, 0.0_dp, kind=dp))
    call check_true('erf nan real', ieee_is_nan(real(z, dp)))
    call check_true('erf nan imaginary', ieee_is_nan(aimag(z)))

    print *, 'special complex error tests passed'

contains

    subroutine check_close(name, actual_value, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual_value !! computed value
        real(dp), intent(in) :: expected !! reference value
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

    subroutine check_complex(name, actual_value, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual_value !! computed complex value
        complex(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum scaled absolute error

        if (actual_value == expected) return
        if (ieee_is_nan(real(actual_value, dp)) .or. ieee_is_nan(aimag(actual_value)) .or. &
                abs(actual_value - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual_value
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_complex_error
