! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_complex_spence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_special, only : spence
    implicit none

    complex(dp) :: values(4)
    complex(dp) :: z
    complex(dp) :: lhs
    complex(dp) :: rhs
    real(dp) :: nan_value
    real(dp) :: pos_inf

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_dp, ieee_positive_inf)

    call check_complex_rel('spence 1-i', spence(cmplx(1.0_dp, -1.0_dp, kind=dp)), &
        cmplx(-0.20561675835602825_dp, 0.91596559417721902_dp, kind=dp), 4.0e-15_dp)
    call check_complex_rel('spence 1.5+2i', spence(cmplx(1.5_dp, 2.0_dp, kind=dp)), &
        cmplx(-0.86766908589075903_dp, -1.3956013375813165_dp, kind=dp), 4.0e-15_dp)
    call check_complex_rel('spence 3i', spence(cmplx(0.0_dp, 3.0_dp, kind=dp)), &
        cmplx(-0.59422063930623703_dp, -2.4912991845187826_dp, kind=dp), 4.0e-15_dp)
    call check_complex_rel('spence -10-5i', spence(cmplx(-10.0_dp, -5.0_dp, kind=dp)), &
        cmplx(-1.1404439840810280_dp, 6.8007592404935260_dp, kind=dp), 5.0e-15_dp)

    call check_complex_rel('spence cut upper', spence(cmplx(-2.0_dp, 0.0_dp, kind=dp)), &
        cmplx(2.3201804233130980_dp, -3.4513922952232030_dp, kind=dp), 5.0e-15_dp)
    call check_complex_rel('spence cut lower signed zero', &
        spence(cmplx(-2.0_dp, -0.0_dp, kind=dp)), &
        cmplx(2.3201804233130980_dp, 3.4513922952232030_dp, kind=dp), 5.0e-15_dp)
    call check_complex_rel('spence near cut upper', &
        spence(cmplx(-2.0_dp, 1.0e-8_dp, kind=dp)), &
        cmplx(2.3201804128411228_dp, -3.4513922929127125_dp, kind=dp), 5.0e-15_dp)
    call check_complex_rel('spence near cut lower', &
        spence(cmplx(-2.0_dp, -1.0e-8_dp, kind=dp)), &
        cmplx(2.3201804128411228_dp, 3.4513922929127125_dp, kind=dp), 5.0e-15_dp)

    call check_complex_rel('spence real complex zero', spence(cmplx(0.0_dp, 0.0_dp, kind=dp)), &
        cmplx(scifort_pi * scifort_pi / 6.0_dp, 0.0_dp, kind=dp), 3.0e-15_dp)
    call check_complex_close('spence real complex one', spence(cmplx(1.0_dp, 0.0_dp, kind=dp)), &
        cmplx(0.0_dp, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex_rel('spence real complex two', spence(cmplx(2.0_dp, 0.0_dp, kind=dp)), &
        cmplx(-0.82246703342411322_dp, 0.0_dp, kind=dp), 4.0e-15_dp)

    z = cmplx(3.0_dp, 4.0_dp, kind=dp)
    lhs = spence(z) + spence(1.0_dp - z)
    rhs = cmplx(scifort_pi * scifort_pi / 6.0_dp, 0.0_dp, kind=dp) - &
        log(z) * log(1.0_dp - z)
    call check_complex_rel('spence reflection identity', lhs, rhs, 8.0e-15_dp)

    values = spence([cmplx(1.0_dp, -1.0_dp, kind=dp), &
        cmplx(-2.0_dp, 0.0_dp, kind=dp), cmplx(1.5_dp, 2.0_dp, kind=dp), &
        cmplx(0.0_dp, 3.0_dp, kind=dp)])
    call check_complex_rel('spence elemental first', values(1), &
        cmplx(-0.20561675835602825_dp, 0.91596559417721902_dp, kind=dp), 4.0e-15_dp)
    call check_complex_rel('spence elemental last', values(4), &
        cmplx(-0.59422063930623703_dp, -2.4912991845187826_dp, kind=dp), 4.0e-15_dp)

    z = spence(cmplx(1.0e300_dp, 1.0e300_dp, kind=dp))
    call check_true('spence huge finite real', ieee_is_finite(real(z, dp)))
    call check_true('spence huge finite imag', ieee_is_finite(aimag(z)))

    z = spence(cmplx(pos_inf, 0.0_dp, kind=dp))
    call check_true('spence complex infinity real nan', ieee_is_nan(real(z, dp)))
    call check_true('spence complex infinity imag nan', ieee_is_nan(aimag(z)))
    z = spence(cmplx(nan_value, 1.0_dp, kind=dp))
    call check_true('spence complex nan real', ieee_is_nan(real(z, dp)))
    call check_true('spence complex nan imag', ieee_is_nan(aimag(z)))

    print *, 'complex Spence special-function tests passed'

contains

    subroutine check_complex_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum absolute complex error

        if (actual == expected) return
        if (ieee_is_nan(real(actual, dp)) .or. ieee_is_nan(aimag(actual)) .or. &
                abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_close

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

end program test_special_complex_spence
