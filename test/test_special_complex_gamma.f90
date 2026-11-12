! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_complex_gamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_special, only : gamma, loggamma, rgamma
    implicit none

    complex(dp) :: array_actual(3)
    complex(dp) :: array_input(3)
    complex(dp) :: lower_cut
    complex(dp) :: upper_cut
    complex(dp) :: value
    complex(dp) :: z
    real(dp) :: nan_value

    call check_complex('loggamma ordinary', loggamma(cmplx(1.5_dp, 2.0_dp, kind=dp)), &
        cmplx(-1.4991963725850939_dp, 0.7332806816909994_dp, kind=dp), 8.0e-15_dp)
    call check_complex('gamma ordinary', gamma(cmplx(2.5_dp, 1.0_dp, kind=dp)), &
        cmplx(0.7747621045510842_dp, 0.7076312043795936_dp, kind=dp), 8.0e-15_dp)
    call check_complex('rgamma ordinary', rgamma(cmplx(2.5_dp, 1.0_dp, kind=dp)), &
        cmplx(0.7036905931010070_dp, -0.6427178342636960_dp, kind=dp), 8.0e-15_dp)

    upper_cut = loggamma(cmplx(-3.2_dp, 0.0_dp, kind=dp))
    lower_cut = loggamma(cmplx(-3.2_dp, -0.0_dp, kind=dp))
    call check_close('loggamma upper cut real', real(upper_cut, dp), &
        -0.37243213612996806_dp, 8.0e-15_dp)
    call check_close('loggamma upper cut imaginary', aimag(upper_cut), &
        -4.0_dp * scifort_pi, 4.0e-15_dp)
    call check_close('loggamma lower cut real', real(lower_cut, dp), &
        -0.37243213612996806_dp, 8.0e-15_dp)
    call check_close('loggamma lower cut imaginary', aimag(lower_cut), &
        4.0_dp * scifort_pi, 4.0e-15_dp)
    call check_complex('loggamma negative half upper', &
        loggamma(cmplx(-0.5_dp, 0.0_dp, kind=dp)), &
        cmplx(1.265512123484647_dp, -scifort_pi, kind=dp), 5.0e-15_dp)
    call check_complex('loggamma negative half lower', &
        loggamma(cmplx(-0.5_dp, -0.0_dp, kind=dp)), &
        cmplx(1.265512123484647_dp, scifort_pi, kind=dp), 5.0e-15_dp)

    call check_complex('loggamma negative complex', &
        loggamma(cmplx(-3.2_dp, 0.7_dp, kind=dp)), &
        cmplx(-2.3406078939632646_dp, -10.713635915626586_dp, kind=dp), 2.0e-14_dp)
    call check_complex('loggamma large complex', &
        loggamma(cmplx(-100.0_dp, 50.0_dp, kind=dp)), &
        cmplx(-507.0115475078980_dp, -83.29671223031121_dp, kind=dp), 2.0e-14_dp)
    call check_complex('gamma large complex', gamma(cmplx(100.0_dp, 200.0_dp, kind=dp)), &
        cmplx(3.1120557973638144e94_dp, -2.1496300650288736e94_dp, kind=dp), 3.0e-14_dp)
    call check_complex('rgamma large complex', rgamma(cmplx(100.0_dp, 200.0_dp, kind=dp)), &
        cmplx(2.175380366403182e-95_dp, 1.5026282762844462e-95_dp, kind=dp), 3.0e-14_dp)
    value = gamma(cmplx(200.0_dp, 1.0_dp, kind=dp))
    call check_true('gamma overflow real sign', real(value, dp) > huge(1.0_dp))
    call check_true('gamma overflow imaginary sign', aimag(value) < -huge(1.0_dp))
    value = rgamma(cmplx(-200.0_dp, 1.0_dp, kind=dp))
    call check_true('rgamma overflow real sign', real(value, dp) < -huge(1.0_dp))
    call check_true('rgamma overflow imaginary sign', aimag(value) > huge(1.0_dp))
    call check_complex('gamma near pole', gamma(cmplx(-3.0_dp, 1.0e-12_dp, kind=dp)), &
        cmplx(-0.20938548530594636_dp, 1.6666666666666660e11_dp, kind=dp), 3.0e-14_dp)

    z = cmplx(1.25_dp, -0.75_dp, kind=dp)
    call check_complex('loggamma recurrence', loggamma(z + 1.0_dp), &
        log(z) + loggamma(z), 2.0e-14_dp)
    call check_complex('gamma recurrence', gamma(z + 1.0_dp), z * gamma(z), 2.0e-14_dp)
    call check_complex('gamma reciprocal identity', gamma(z) * rgamma(z), &
        cmplx(1.0_dp, 0.0_dp, kind=dp), 2.0e-14_dp)
    call check_complex('loggamma conjugation', loggamma(conjg(z)), conjg(loggamma(z)), 2.0e-14_dp)
    call check_complex('gamma conjugation', gamma(conjg(z)), conjg(gamma(z)), 2.0e-14_dp)

    call check_close('real gamma generic preserved', gamma(-0.5_dp), &
        -3.5449077018110318_dp, 3.0e-15_dp)
    call check_true('real loggamma branch preserved', ieee_is_nan(loggamma(-0.5_dp)))
    call check_close('real rgamma generic preserved', rgamma(20.0_dp), &
        8.22063524662433e-18_dp, 5.0e-15_dp)

    value = gamma(cmplx(-2.0_dp, 0.0_dp, kind=dp))
    call check_true('complex gamma pole real nan', ieee_is_nan(real(value, dp)))
    call check_true('complex gamma pole imaginary nan', ieee_is_nan(aimag(value)))
    value = loggamma(cmplx(0.0_dp, 0.0_dp, kind=dp))
    call check_true('complex loggamma zero real nan', ieee_is_nan(real(value, dp)))
    call check_true('complex loggamma zero imaginary nan', ieee_is_nan(aimag(value)))
    call check_complex('complex rgamma pole zero', rgamma(cmplx(-3.0_dp, 0.0_dp, kind=dp)), &
        cmplx(0.0_dp, 0.0_dp, kind=dp), 0.0_dp)

    array_input = [cmplx(0.25_dp, 0.5_dp, kind=dp), &
        cmplx(-1.25_dp, 0.4_dp, kind=dp), cmplx(4.0_dp, -2.0_dp, kind=dp)]
    array_actual = loggamma(array_input)
    call check_complex('loggamma elemental 1', array_actual(1), loggamma(array_input(1)), 0.0_dp)
    call check_complex('loggamma elemental 2', array_actual(2), loggamma(array_input(2)), 0.0_dp)
    call check_complex('loggamma elemental 3', array_actual(3), loggamma(array_input(3)), 0.0_dp)

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    value = gamma(cmplx(nan_value, 0.0_dp, kind=dp))
    call check_true('complex gamma nan real', ieee_is_nan(real(value, dp)))
    call check_true('complex gamma nan imaginary', ieee_is_nan(aimag(value)))

    print *, 'special complex gamma tests passed'

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

end program test_special_complex_gamma
