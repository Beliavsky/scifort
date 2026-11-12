! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_bessel_integrals
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, it2i0k0, it2j0y0, iti0k0, itj0y0
    implicit none

    real(dp) :: a
    real(dp) :: b
    real(dp) :: nan_value
    real(dp) :: pos_inf
    real(dp) :: values_a(3)
    real(dp) :: values_b(3)
    real(dp), parameter :: x_values(3) = [0.1_dp, 1.0_dp, 20.0_dp]

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call itj0y0(0.1_dp, a, b)
    call check_rel('itj0y0 J series', a, &
        0.09991669791046702330006932168851845_dp, 3.0e-14_dp)
    call check_rel('itj0y0 Y series', b, &
        -0.2174305666201166107559790741138575_dp, 3.0e-14_dp)

    call itj0y0(20.0_dp, a, b)
    call check_rel('itj0y0 J quadrature', a, &
        1.058378821421127758487940219255186_dp, 3.0e-14_dp)
    call check_close('itj0y0 Y quadrature', b, &
        -0.1682159767721502961135076318084084_dp, 2.0e-14_dp)

    call itj0y0(40.0_dp, a, b)
    call check_rel('itj0y0 J asymptotic', a, &
        1.125776150359991460304483362480519_dp, 3.0e-14_dp)
    call check_rel('itj0y0 Y asymptotic', b, &
        -0.008932466273627416184069357915951210_dp, 3.0e-14_dp)

    call it2j0y0(1.0_dp, a, b)
    call check_rel('it2j0y0 J series', a, &
        0.1211652469950687059698058645623951_dp, 3.0e-14_dp)
    call check_rel('it2j0y0 Y series', b, &
        0.3952729016992932814252300382729309_dp, 3.0e-14_dp)

    call it2j0y0(20.0_dp, a, b)
    call check_rel('it2j0y0 J quadrature', a, &
        2.877310607840809640475111866125831_dp, 3.0e-14_dp)
    call check_rel('it2j0y0 Y quadrature', b, &
        0.008503152668301641945861580230087621_dp, 2.0e-12_dp)

    call it2j0y0(40.0_dp, a, b)
    call check_rel('it2j0y0 J asymptotic', a, &
        3.569813944072160865613961110204289_dp, 3.0e-14_dp)
    call check_rel('it2j0y0 Y asymptotic', b, &
        0.0003011239815299519001769589302971811_dp, 3.0e-13_dp)

    call iti0k0(1.0_dp, a, b)
    call check_rel('iti0k0 I series', a, &
        1.086521097023589815837941923492506_dp, 3.0e-14_dp)
    call check_rel('iti0k0 K series', b, &
        1.242509848623778266223432356721264_dp, 3.0e-14_dp)

    call iti0k0(10.0_dp, a, b)
    call check_rel('iti0k0 I large series', a, &
        2993.044542811447950490029594281177_dp, 3.0e-14_dp)
    call check_rel('iti0k0 K tail', b, &
        1.570779311615978859831228128596967_dp, 3.0e-13_dp)

    call it2i0k0(1.0_dp, a, b)
    call check_rel('it2i0k0 I series', a, &
        0.1289794424945685369518314009417339_dp, 3.0e-14_dp)
    call check_rel('it2i0k0 K series', b, &
        0.2085182909001294927382536398229784_dp, 3.0e-14_dp)

    call it2i0k0(10.0_dp, a, b)
    call check_rel('it2i0k0 I large series', a, &
        340.8153668040809672287908756626073_dp, 3.0e-14_dp)
    call check_rel('it2i0k0 K tail', b, &
        0.000001562928143076405436669322822113825_dp, 3.0e-12_dp)

    call itj0y0(-1.0_dp, a, b)
    call check_rel('itj0y0 negative J parity', a, &
        -0.9197304100897602393144211940806200_dp, 3.0e-14_dp)
    call check_true('itj0y0 negative Y domain', ieee_is_nan(b))

    call it2j0y0(-1.0_dp, a, b)
    call check_rel('it2j0y0 negative J parity', a, &
        0.1211652469950687059698058645623951_dp, 3.0e-14_dp)
    call check_true('it2j0y0 negative Y domain', ieee_is_nan(b))

    call iti0k0(-1.0_dp, a, b)
    call check_rel('iti0k0 negative I parity', a, &
        -1.086521097023589815837941923492506_dp, 3.0e-14_dp)
    call check_true('iti0k0 negative K domain', ieee_is_nan(b))

    call it2i0k0(-1.0_dp, a, b)
    call check_rel('it2i0k0 negative I parity', a, &
        0.1289794424945685369518314009417339_dp, 3.0e-14_dp)
    call check_true('it2i0k0 negative K domain', ieee_is_nan(b))

    call it2j0y0(0.0_dp, a, b)
    call check_close('it2j0y0 zero J', a, 0.0_dp, 0.0_dp)
    call check_true('it2j0y0 zero Y divergence', .not. ieee_is_finite(b) .and. b < 0.0_dp)

    call it2i0k0(0.0_dp, a, b)
    call check_close('it2i0k0 zero I', a, 0.0_dp, 0.0_dp)
    call check_true('it2i0k0 zero K divergence', .not. ieee_is_finite(b) .and. b > 0.0_dp)

    call itj0y0(pos_inf, a, b)
    call check_close('itj0y0 positive infinity J', a, 1.0_dp, 0.0_dp)
    call check_close('itj0y0 positive infinity Y', b, 0.0_dp, 0.0_dp)

    call iti0k0(pos_inf, a, b)
    call check_true('iti0k0 positive infinity I', .not. ieee_is_finite(a) .and. a > 0.0_dp)
    call check_rel('iti0k0 positive infinity K', b, &
        1.570796326794896619231321691639751_dp, 3.0e-16_dp)

    call itj0y0(nan_value, a, b)
    call check_true('itj0y0 NaN propagation', ieee_is_nan(a) .and. ieee_is_nan(b))

    call itj0y0(x_values, values_a, values_b)
    call check_rel('itj0y0 elemental J 1', values_a(1), &
        0.09991669791046702330006932168851845_dp, 3.0e-14_dp)
    call check_rel('itj0y0 elemental J 2', values_a(2), &
        0.9197304100897602393144211940806200_dp, 3.0e-14_dp)
    call check_rel('itj0y0 elemental J 3', values_a(3), &
        1.058378821421127758487940219255186_dp, 3.0e-14_dp)

    print *, 'Bessel integral special-function tests passed'

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

    subroutine check_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! nonzero reference value
        real(dp), intent(in) :: tolerance !! maximum relative error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' relative tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_rel

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_bessel_integrals
