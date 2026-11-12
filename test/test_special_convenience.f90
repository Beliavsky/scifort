! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_convenience
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : cbrt, cosdg, cotdg, diric, dp, exp10, exp2, &
        expm1, log1p, powm1, radian, round, sinc, sindg, softplus, tandg
    implicit none

    real(dp) :: array_actual(3)
    real(dp) :: array_input(3)
    real(dp) :: nan_value
    real(dp) :: pinf

    pinf = ieee_value(1.0_dp, ieee_positive_inf)
    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)

    call check_close('cbrt negative', cbrt(-8.0_dp), -2.0_dp, 2.0e-15_dp)
    call check_relative('cbrt large', cbrt(1.0e300_dp), 1.0e100_dp, 3.0e-16_dp)
    call check_close('cbrt signed zero', cbrt(-0.0_dp), -0.0_dp, 0.0_dp)
    call check_true('cbrt infinity', cbrt(pinf) > huge(1.0_dp))

    call check_close('exp10 quarter', exp10(0.25_dp), &
        1.7782794100389228_dp, 3.0e-16_dp)
    call check_close('exp2 half', exp2(0.5_dp), &
        1.4142135623730949_dp, 3.0e-16_dp)
    call check_true('exp10 overflow', exp10(400.0_dp) > huge(1.0_dp))
    call check_true('exp2 overflow', exp2(2000.0_dp) > huge(1.0_dp))
    call check_close('exp2 negative infinity', exp2(-pinf), 0.0_dp, 0.0_dp)

    call check_close('expm1 tiny', expm1(1.0e-12_dp), &
        1.0000000000005e-12_dp, 2.0e-27_dp)
    call check_close('log1p tiny', log1p(1.0e-12_dp), &
        9.999999999995e-13_dp, 2.0e-27_dp)
    call check_close('expm1 negative infinity', expm1(-pinf), -1.0_dp, 0.0_dp)
    call check_true('log1p lower pole', log1p(-1.0_dp) < -huge(1.0_dp))
    call check_true('log1p invalid', ieee_is_nan(log1p(-2.0_dp)))

    call check_close('powm1 near unity', powm1(1.0_dp + 1.0e-12_dp, 0.5_dp), &
        5.000444502910455e-13_dp, 2.0e-28_dp)
    call check_close('powm1 negative even', powm1(-2.0_dp, 4.0_dp), 15.0_dp, 0.0_dp)
    call check_close('powm1 negative odd', powm1(-2.0_dp, 3.0_dp), -9.0_dp, 0.0_dp)
    call check_close('powm1 zero exponent', powm1(nan_value, 0.0_dp), 0.0_dp, 0.0_dp)
    call check_true('powm1 invalid real branch', ieee_is_nan(powm1(-2.0_dp, 0.5_dp)))

    call check_close('radian dms', radian(12.0_dp, 34.0_dp, 56.0_dp), &
        0.21960120499537542_dp, 3.0e-17_dp)
    call check_close('radian signed components', radian(-1.0_dp, 30.0_dp, 0.0_dp), &
        -0.008726646259971648_dp, 2.0e-18_dp)
    call check_true('radian infinity', radian(pinf, 0.0_dp, 0.0_dp) > huge(1.0_dp))

    call check_close('round tie even 0', round(0.5_dp), 0.0_dp, 0.0_dp)
    call check_close('round tie even 2', round(1.5_dp), 2.0_dp, 0.0_dp)
    call check_close('round negative tie', round(-2.5_dp), -2.0_dp, 0.0_dp)
    call check_close('round ordinary', round(-1.9_dp), -2.0_dp, 0.0_dp)
    call check_true('round infinity', round(pinf) > huge(1.0_dp))

    call check_close('sinc zero', sinc(0.0_dp), 1.0_dp, 0.0_dp)
    call check_close('sinc half', sinc(0.5_dp), &
        0.63661977236758138_dp, 2.0e-16_dp)
    call check_close('sinc integer zero', sinc(3.0_dp), 0.0_dp, 0.0_dp)
    call check_true('sinc infinity', ieee_is_nan(sinc(pinf)))

    call check_close('sindg 30', sindg(30.0_dp), 0.5_dp, 2.0e-16_dp)
    call check_close('sindg 90 exact', sindg(90.0_dp), 1.0_dp, 0.0_dp)
    call check_close('cosdg 60', cosdg(60.0_dp), 0.5_dp, 2.0e-16_dp)
    call check_close('cosdg 180 exact', cosdg(180.0_dp), -1.0_dp, 0.0_dp)
    call check_close('tandg 60', tandg(60.0_dp), &
        1.7320508075688767_dp, 2.0e-15_dp)
    call check_true('tandg 90 pole', tandg(90.0_dp) > huge(1.0_dp))
    call check_close('cotdg 60', cotdg(60.0_dp), &
        0.57735026918962573_dp, 7.0e-16_dp)
    call check_true('cotdg zero pole', cotdg(0.0_dp) > huge(1.0_dp))
    call check_true('cotdg tiny overflow guard', cotdg(tiny(1.0_dp)) > huge(1.0_dp))
    call check_true('degree trig infinity nan', ieee_is_nan(sindg(pinf)))

    call check_close('softplus zero', softplus(0.0_dp), &
        0.69314718055994531_dp, 3.0e-16_dp)
    call check_close('softplus negative tail', softplus(-20.0_dp), &
        2.0611536203143808e-9_dp, 2.0e-24_dp)
    call check_close('softplus negative infinity', softplus(-pinf), 0.0_dp, 0.0_dp)
    call check_true('softplus positive infinity', softplus(pinf) > huge(1.0_dp))

    call check_close('diric ordinary', diric(1.0_dp, 5), &
        0.24966218772839896_dp, 4.0e-16_dp)
    call check_close('diric origin', diric(0.0_dp, 4), 1.0_dp, 0.0_dp)
    call check_close('diric two pi even', diric(2.0_dp * acos(-1.0_dp), 2), &
        -1.0_dp, 2.0e-15_dp)
    call check_close('diric pi odd', diric(acos(-1.0_dp), 3), &
        -1.0_dp / 3.0_dp, 3.0e-16_dp)
    call check_true('diric invalid order', ieee_is_nan(diric(1.0_dp, 0)))

    array_input = [-1.0_dp, 0.0_dp, 1.0_dp]
    array_actual = softplus(array_input)
    call check_close('softplus elemental 1', array_actual(1), &
        0.31326168751822286_dp, 3.0e-16_dp)
    call check_close('softplus elemental 2', array_actual(2), &
        0.69314718055994531_dp, 3.0e-16_dp)
    call check_close('softplus elemental 3', array_actual(3), &
        1.3132616875182228_dp, 4.0e-16_dp)

    print *, 'special convenience-function tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
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

    subroutine check_relative(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! nonzero reference value
        real(dp), intent(in) :: tolerance !! maximum relative error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_relative

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_convenience
