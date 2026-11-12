! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_faddeeva
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : dawsn, dp, erfcx, wofz
    implicit none

    complex(dp) :: actual
    complex(dp) :: array_actual(3)
    complex(dp) :: array_input(3)
    complex(dp) :: reflected
    complex(dp) :: z
    real(dp) :: nan_value
    real(dp) :: pinf
    real(dp) :: x

    pinf = ieee_value(1.0_dp, ieee_positive_inf)
    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)

    call check_complex('wofz zero', wofz(cmplx(0.0_dp, 0.0_dp, kind=dp)), &
        cmplx(1.0_dp, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex('wofz series', wofz(cmplx(0.25_dp, 0.5_dp, kind=dp)), &
        cmplx(0.5937088515052158_dp, 0.12476492296402969_dp, kind=dp), 2.0e-15_dp)
    call check_complex('wofz Taylor', wofz(cmplx(2.0_dp, 0.1_dp, kind=dp)), &
        cmplx(0.040201398161451296_dp, 0.3315826873345632_dp, kind=dp), 2.0e-15_dp)
    call check_complex('wofz quadrature', wofz(cmplx(3.0_dp, 2.0_dp, kind=dp)), &
        cmplx(0.09271076642644344_dp, 0.12831696222826169_dp, kind=dp), 3.0e-15_dp)
    call check_complex('wofz continued fraction', wofz(cmplx(8.0_dp, 3.0_dp, kind=dp)), &
        cmplx(0.02359322812772812_dp, 0.06204131036433501_dp, kind=dp), 2.0e-15_dp)
    call check_complex('wofz lower reflection', wofz(cmplx(-3.0_dp, -1.25_dp, kind=dp)), &
        cmplx(-0.07506171466819871_dp, -0.1636550559585121_dp, kind=dp), 3.0e-14_dp)

    call check_complex('complex erfcx', erfcx(cmplx(0.25_dp, 0.5_dp, kind=dp)), &
        cmplx(0.6383373967914268_dp, -0.32581483982634513_dp, kind=dp), 3.0e-15_dp)
    call check_close('real erfcx generic preserved', erfcx(0.5_dp), &
        0.6156903441929258_dp, 4.0e-16_dp)

    x = 1.25_dp
    actual = wofz(cmplx(x, 0.0_dp, kind=dp))
    call check_close('wofz real-axis real', real(actual, dp), exp(-x * x), 3.0e-16_dp)
    call check_close('wofz real-axis imaginary', aimag(actual), &
        2.0_dp * dawsn(x) / sqrt(acos(-1.0_dp)), 5.0e-16_dp)

    z = cmplx(2.0_dp, 0.75_dp, kind=dp)
    reflected = 2.0_dp * exp(-z * z) - wofz(z)
    call check_complex('wofz reflection identity', wofz(-z), reflected, 5.0e-15_dp)

    array_input = [cmplx(0.1_dp, 0.2_dp, kind=dp), &
        cmplx(-1.0_dp, 0.4_dp, kind=dp), cmplx(4.0_dp, 1.0_dp, kind=dp)]
    array_actual = wofz(array_input)
    call check_complex('wofz elemental 1', array_actual(1), wofz(array_input(1)), 0.0_dp)
    call check_complex('wofz elemental 2', array_actual(2), wofz(array_input(2)), 0.0_dp)
    call check_complex('wofz elemental 3', array_actual(3), wofz(array_input(3)), 0.0_dp)

    actual = wofz(cmplx(nan_value, 0.0_dp, kind=dp))
    call check_true('wofz nan real', ieee_is_nan(real(actual, dp)))
    call check_true('wofz nan imaginary', ieee_is_nan(aimag(actual)))
    call check_complex('wofz +imag infinity', wofz(cmplx(0.0_dp, pinf, kind=dp)), &
        cmplx(0.0_dp, 0.0_dp, kind=dp), 0.0_dp)

    print *, 'special Faddeeva tests passed'

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
        complex(dp), intent(in) :: expected !! complex reference value
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

end program test_special_faddeeva
