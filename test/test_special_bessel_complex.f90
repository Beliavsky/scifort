! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_bessel_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, iv, ive, ivp, jn, jv, jve, jvp, kv, kve, kvp, &
        yv, yve, yvp
    implicit none

    real(dp), parameter :: orders(3) = [0.0_dp, -0.5_dp, 2.3_dp]
    complex(dp), parameter :: arguments(3) = [ &
        cmplx(2.0_dp, 1.0_dp, kind=dp), &
        cmplx(0.5_dp, -3.0_dp, kind=dp), &
        cmplx(-2.0_dp, 1.0_dp, kind=dp)]
    complex(dp), parameter :: j_refs(3) = [ &
        cmplx(0.18785372808246176_dp, -0.6461694351539806_dp, kind=dp), &
        cmplx(1.6639848742985486_dp, 4.2894276110328440_dp, kind=dp), &
        cmplx(0.3788230216137206_dp, 0.07504467420501187_dp, kind=dp)]
    complex(dp) :: values(3)
    real(dp) :: nan_value
    real(dp) :: pos_inf

    call check_complex_rel('jv upper right', jv(2.3_dp, cmplx(1.2_dp, 0.7_dp, kind=dp)), &
        cmplx(0.07010731829731956_dp, 0.13283504662956458_dp, kind=dp), 6.0e-10_dp)
    call check_complex_rel('jv negative order upper left', &
        jv(-0.5_dp, cmplx(-2.0_dp, 1.0_dp, kind=dp)), &
        cmplx(0.47621033866924845_dp, 0.46447163413219034_dp, kind=dp), 6.0e-10_dp)
    call check_complex_rel('jn complex alias', jn(2.3_dp, cmplx(1.2_dp, 0.7_dp, kind=dp)), &
        jv(2.3_dp, cmplx(1.2_dp, 0.7_dp, kind=dp)), 0.0_dp)
    call check_complex_rel('jve larger complex', jve(10.2_dp, cmplx(20.0_dp, 10.0_dp, kind=dp)), &
        cmplx(0.022037046705856562_dp, 0.019758059971445453_dp, kind=dp), 8.0e-10_dp)

    call check_complex_rel('jv large-order small complex series', &
        jv(20.25_dp, cmplx(0.2_dp, 0.3_dp, kind=dp)), &
        cmplx(8.217467251114453e-35_dp, 1.4345746781724917e-34_dp, kind=dp), 1.0e-44_dp)
    call check_complex_rel('jve large-order small complex series', &
        jve(20.25_dp, cmplx(0.2_dp, 0.3_dp, kind=dp)), &
        cmplx(6.087649467480913e-35_dp, 1.0627590605187962e-34_dp, kind=dp), 1.0e-44_dp)

    call check_complex_rel('yv upper left', yv(2.3_dp, cmplx(-2.0_dp, 1.0_dp, kind=dp)), &
        cmplx(-0.3578546039355087_dp, 0.6052005957571243_dp, kind=dp), 7.0e-10_dp)
    call check_complex_rel('yve lower right negative order', &
        yve(-0.5_dp, cmplx(0.5_dp, -3.0_dp, kind=dp)), &
        cmplx(0.21332958424204593_dp, -0.08173411068365298_dp, kind=dp), 6.0e-10_dp)

    call check_complex_rel('iv upper left', iv(2.3_dp, cmplx(-2.0_dp, 1.0_dp, kind=dp)), &
        cmplx(0.5559535894373280_dp, -0.24266929258099262_dp, kind=dp), 7.0e-10_dp)
    call check_complex_rel('ive lower large negative order', &
        ive(-0.5_dp, cmplx(20.0_dp, -10.0_dp, kind=dp)), &
        cmplx(-0.07944039122408875_dp, 0.028405099477655304_dp, kind=dp), 8.0e-10_dp)

    call check_complex_rel('iv high-order turning regime', &
        iv(100.5_dp, cmplx(100.0_dp, 20.0_dp, kind=dp)), &
        cmplx(-1.1691049462104622e22_dp, 2.371818855679421e21_dp, kind=dp), 2.0e-10_dp)
    call check_complex_abs('ive high-order turning regime', &
        ive(100.5_dp, cmplx(100.0_dp, 20.0_dp, kind=dp)), &
        cmplx(-4.349159223844675e-22_dp, 8.823346344486091e-23_dp, kind=dp), 2.0e-34_dp)

    call check_complex_rel('iv large-order small complex series', &
        iv(20.25_dp, cmplx(0.2_dp, 0.3_dp, kind=dp)), &
        cmplx(8.167314651162235e-35_dp, 1.4351997188422873e-34_dp, kind=dp), 1.0e-44_dp)
    call check_complex_rel('ive large-order small complex series', &
        ive(20.25_dp, cmplx(0.2_dp, 0.3_dp, kind=dp)), &
        cmplx(6.686831674970872e-35_dp, 1.1750421466250504e-34_dp, kind=dp), 1.0e-44_dp)

    call check_complex_rel('kv upper left', kv(2.3_dp, cmplx(-2.0_dp, 1.0_dp, kind=dp)), &
        cmplx(-1.6440432817314512_dp, -0.25168948674181313_dp, kind=dp), 8.0e-10_dp)
    call check_complex_rel('kve lower large negative order', &
        kve(-0.5_dp, cmplx(20.0_dp, -10.0_dp, kind=dp)), &
        cmplx(0.25795348454925277_dp, 0.060894557386565344_dp, kind=dp), 8.0e-10_dp)

    call check_complex_rel('jvp complex n2', jvp(0.5_dp, cmplx(1.2_dp, 0.7_dp, kind=dp), 2), &
        cmplx(-0.5770925139144379_dp, 0.14657839934284167_dp, kind=dp), 9.0e-10_dp)
    call check_complex_rel('yvp upper left n3', &
        yvp(2.3_dp, cmplx(-1.5_dp, 0.8_dp, kind=dp), 3), &
        cmplx(0.15060340227315494_dp, -3.2935908386357293_dp, kind=dp), 1.2e-9_dp)
    call check_complex_rel('ivp lower left n2', &
        ivp(-0.7_dp, cmplx(-1.5_dp, -0.8_dp, kind=dp), 2), &
        cmplx(-0.7219737077108694_dp, 0.28069692759693404_dp, kind=dp), 1.0e-9_dp)
    call check_complex_rel('kvp lower left n3', &
        kvp(2.3_dp, cmplx(-1.5_dp, -0.8_dp, kind=dp), 3), &
        cmplx(0.5693235561139090_dp, -6.6083824197948795_dp, kind=dp), 1.5e-9_dp)

    call check_complex_rel('jvp large-order small complex n2', &
        jvp(20.25_dp, cmplx(0.2_dp, 0.3_dp, kind=dp), 2), &
        cmplx(3.022242541110661e-31_dp, -3.930395761696471e-31_dp, kind=dp), 1.0e-40_dp)
    call check_complex_rel('ivp large-order small complex n2', &
        ivp(20.25_dp, cmplx(0.2_dp, 0.3_dp, kind=dp), 2), &
        cmplx(3.0313565728971083e-31_dp, -3.914432630278655e-31_dp, kind=dp), 1.0e-40_dp)

    call check_complex_rel('negative-axis jv branch', &
        jv(2.3_dp, cmplx(-2.0_dp, 0.0_dp, kind=dp)), &
        cmplx(0.1599170783911904_dp, 0.2201069754722054_dp, kind=dp), 7.0e-10_dp)
    call check_complex_rel('negative-axis yv branch', &
        yv(2.3_dp, cmplx(-2.0_dp, 0.0_dp, kind=dp)), &
        cmplx(-0.43481956090839935_dp, 0.9183119390839153_dp, kind=dp), 8.0e-10_dp)
    call check_complex_rel('negative-axis iv branch', &
        iv(2.3_dp, cmplx(-2.0_dp, 0.0_dp, kind=dp)), &
        cmplx(0.2936416476491177_dp, 0.4041630549216118_dp, kind=dp), 8.0e-10_dp)
    call check_complex_rel('negative-axis kv branch', &
        kv(2.3_dp, cmplx(-2.0_dp, 0.0_dp, kind=dp)), &
        cmplx(0.1910940681243293_dp, -1.8324733182372310_dp, kind=dp), 9.0e-10_dp)

    call check_complex_rel('jvp n0', jvp(0.7_dp, cmplx(1.1_dp, 0.4_dp, kind=dp), 0), &
        jv(0.7_dp, cmplx(1.1_dp, 0.4_dp, kind=dp)), 0.0_dp)
    call check_complex_rel('kvp n0', kvp(0.7_dp, cmplx(1.1_dp, 0.4_dp, kind=dp), 0), &
        kv(0.7_dp, cmplx(1.1_dp, 0.4_dp, kind=dp)), 0.0_dp)

    values = jv(orders, arguments)
    call check_complex_array_rel('jv complex elemental arrays', values, j_refs, 8.0e-10_dp)

    call check_complex_rel('jv complex zero order zero', jv(0.0_dp, (0.0_dp, 0.0_dp)), &
        (1.0_dp, 0.0_dp), 0.0_dp)
    call check_complex_rel('iv complex positive order zero', iv(2.0_dp, (0.0_dp, 0.0_dp)), &
        (0.0_dp, 0.0_dp), 0.0_dp)
    call check_true('yve complex zero is nan', complex_is_nan(yve(0.0_dp, (0.0_dp, 0.0_dp))))
    call check_true('kve complex zero is nan', complex_is_nan(kve(0.0_dp, (0.0_dp, 0.0_dp))))

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)
    call check_true('complex nan argument', complex_is_nan( &
        jv(0.0_dp, cmplx(nan_value, 1.0_dp, kind=dp))))
    call check_true('complex infinite argument', complex_is_nan( &
        kv(0.0_dp, cmplx(1.0_dp, pos_inf, kind=dp))))
    call check_true('complex negative derivative order', complex_is_nan( &
        ivp(0.0_dp, cmplx(1.0_dp, 1.0_dp, kind=dp), -1)))

    print *, 'complex Bessel J/Y/I/K special-function tests passed'

contains

    subroutine check_complex_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! finite complex reference value
        real(dp), intent(in) :: tolerance !! maximum scaled absolute error

        if (complex_is_nan(actual) .or. &
                abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' scaled tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_rel

    subroutine check_complex_abs(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! finite complex reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (complex_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' absolute tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_abs

    subroutine check_complex_array_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual(:) !! computed complex values
        complex(dp), intent(in) :: expected(:) !! finite complex reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise scaled absolute error

        integer :: i

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', name, ': size mismatch'
            error stop 1
        end if
        do i = 1, size(actual)
            call check_complex_rel(name, actual(i), expected(i), tolerance)
        end do
    end subroutine check_complex_array_rel

    pure logical function complex_is_nan(z) result(is_nan)
        use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
        complex(dp), intent(in) :: z !! complex value to classify

        is_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
    end function complex_is_nan

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_bessel_complex
