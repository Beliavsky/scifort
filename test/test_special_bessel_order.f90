! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_bessel_order
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, iv, ive, kn, kv, kve, yn
    implicit none

    integer, parameter :: kn_orders(3) = [-2, 0, 3]
    real(dp), parameter :: iv_orders(3) = [-1.5_dp, 0.0_dp, 2.5_dp]
    real(dp), parameter :: iv_array_ref(3) = [ &
        -0.29352532634747990_dp, 1.2660658777520084_dp, 0.05709890920304825_dp]
    real(dp), parameter :: kn_array_ref(3) = [ &
        1.6248388986351774_dp, 0.42102443824070834_dp, 7.101262824737944_dp]
    real(dp), parameter :: x_array(3) = [0.1_dp, 1.0_dp, 10.0_dp]
    real(dp), parameter :: yn_array_ref(3) = [ &
        -127.64478324269018_dp, -1.6506826068162546_dp, -0.005868082442208836_dp]
    real(dp) :: nan_value
    real(dp) :: pos_inf
    real(dp) :: values(3)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call check_rel('yn -3 0.1', yn(-3, 0.1_dp), 5099.332378612905_dp, 1.0e-13_dp)
    call check_rel('yn -2 10', yn(-2, 10.0_dp), -0.005868082442208836_dp, 1.0e-13_dp)
    call check_rel('yn 4 1', yn(4, 1.0_dp), -33.27842302897212_dp, 1.0e-13_dp)
    call check_rel('yn 20 100', yn(20, 100.0_dp), 0.051247973076187975_dp, 1.0e-13_dp)
    call check_true('yn -1 zero positive infinity', yn(-1, 0.0_dp) > huge(1.0_dp))
    call check_true('yn -2 zero negative infinity', yn(-2, 0.0_dp) < -huge(1.0_dp))
    call check_true('yn negative x', ieee_is_nan(yn(2, -1.0_dp)))
    call check_true('yn infinity', ieee_is_nan(yn(2, pos_inf)))
    call check_true('yn nan', ieee_is_nan(yn(2, nan_value)))

    call check_rel('kn -5 1', kn(-5, 1.0_dp), 360.96058960124066_dp, 4.0e-12_dp)
    call check_rel('kn 2 0.1', kn(2, 0.1_dp), 199.5039646421141_dp, 4.0e-12_dp)
    call check_rel('kn 4 10', kn(4, 10.0_dp), 3.786143716089198e-5_dp, 4.0e-12_dp)
    call check_rel('kn 20 100', kn(20, 100.0_dp), 3.3852054148901705e-44_dp, 4.0e-12_dp)
    call check_true('kn zero', kn(-2, 0.0_dp) > huge(1.0_dp))
    call check_true('kn negative x', ieee_is_nan(kn(2, -1.0_dp)))
    call check_close('kn infinity', kn(2, pos_inf), 0.0_dp, 0.0_dp)

    call check_rel('iv 2.5 1.2', iv(2.5_dp, 1.2_dp), 0.09289108546788319_dp, 2.0e-13_dp)
    call check_rel('iv -0.5 1', iv(-0.5_dp, 1.0_dp), 1.2312002145929675_dp, 2.0e-13_dp)
    call check_rel('iv -1.5 1', iv(-1.5_dp, 1.0_dp), -0.2935253263474799_dp, 2.0e-13_dp)
    call check_rel('iv -1.5 10', iv(-1.5_dp, 10.0_dp), 2500.9061423416138_dp, 2.0e-13_dp)
    call check_rel('iv 2 -1.2', iv(2.0_dp, -1.2_dp), 0.2025956815463259_dp, 2.0e-13_dp)
    call check_rel('iv -3 -1.2', iv(-3.0_dp, -1.2_dp), -0.039359003064890014_dp, 2.0e-13_dp)
    call check_rel('iv 20 100', iv(20.0_dp, 100.0_dp), 1.4483461256427172e41_dp, 3.0e-13_dp)
    call check_close('iv negative-order zero crossing', iv(-1.5_dp, 1.1996786402577342_dp), &
        0.0_dp, 8.0e-15_dp)
    call check_true('iv negative noninteger zero', iv(-1.5_dp, 0.0_dp) > huge(1.0_dp))
    call check_true('iv negative x noninteger', ieee_is_nan(iv(2.5_dp, -1.0_dp)))
    call check_true('iv infinity', ieee_is_nan(iv(2.0_dp, pos_inf)))
    call check_true('iv opposite infinities', ieee_is_nan(iv(pos_inf, -pos_inf)))
    call check_true('iv nan', ieee_is_nan(iv(nan_value, 1.0_dp)))

    call check_rel('ive 2.5 1.2', ive(2.5_dp, 1.2_dp), 0.027978257281168106_dp, 2.0e-13_dp)
    call check_rel('ive -0.5 1', ive(-0.5_dp, 1.0_dp), 0.4529332469146212_dp, 2.0e-13_dp)
    call check_rel('ive -1.5 1', ive(-1.5_dp, 1.0_dp), -0.10798193302637601_dp, 2.0e-13_dp)
    call check_rel('ive -1.5 10', ive(-1.5_dp, 10.0_dp), 0.11354096320487644_dp, 2.0e-13_dp)
    call check_rel('ive 2 -1.2', ive(2.0_dp, -1.2_dp), 0.06102064664016112_dp, 2.0e-13_dp)
    call check_rel('ive -3 -1.2', ive(-3.0_dp, -1.2_dp), -0.011854703909779501_dp, 2.0e-13_dp)
    call check_rel('ive 20 100', ive(20.0_dp, 100.0_dp), 0.005387957626966431_dp, 3.0e-13_dp)
    call check_true('ive negative noninteger zero', ieee_is_nan(ive(-1.5_dp, 0.0_dp)))
    call check_true('ive infinity', ieee_is_nan(ive(2.0_dp, pos_inf)))
    call check_true('ive opposite infinities', ieee_is_nan(ive(-pos_inf, pos_inf)))

    call check_rel('kv 2.5 1.2', kv(2.5_dp, 1.2_dp), 1.924019608752178_dp, 4.0e-12_dp)
    call check_rel('kv negative order', kv(-2.5_dp, 1.2_dp), 1.924019608752178_dp, 4.0e-12_dp)
    call check_rel('kv 5.5 0.1', kv(5.5_dp, 0.1_dp), 374326429.2282701_dp, 4.0e-12_dp)
    call check_rel('kv 20 100', kv(20.0_dp, 100.0_dp), 3.3852054148901705e-44_dp, 4.0e-12_dp)
    call check_true('kv zero', kv(2.0_dp, 0.0_dp) > huge(1.0_dp))
    call check_true('kv negative x', ieee_is_nan(kv(2.0_dp, -1.0_dp)))
    call check_close('kv infinity', kv(2.0_dp, pos_inf), 0.0_dp, 0.0_dp)
    call check_true('kv infinite order', ieee_is_nan(kv(pos_inf, 1.0_dp)))
    call check_true('kv opposite infinities', ieee_is_nan(kv(pos_inf, -pos_inf)))

    call check_rel('kve 2.5 1.2', kve(2.5_dp, 1.2_dp), 6.3879700626950555_dp, 4.0e-12_dp)
    call check_rel('kve negative order', kve(-2.5_dp, 1.2_dp), 6.3879700626950555_dp, 4.0e-12_dp)
    call check_rel('kve 5.5 0.1', kve(5.5_dp, 0.1_dp), 413694683.45018625_dp, 4.0e-12_dp)
    call check_rel('kve 20 100', kve(20.0_dp, 100.0_dp), 0.909982870433507_dp, 4.0e-12_dp)
    call check_rel('kve 2.5 700', kve(2.5_dp, 700.0_dp), 0.04757412957545403_dp, 4.0e-12_dp)
    call check_true('kve zero', kve(2.0_dp, 0.0_dp) > huge(1.0_dp))
    call check_true('kve infinity', ieee_is_nan(kve(2.0_dp, pos_inf)))
    call check_true('kve opposite infinities', ieee_is_nan(kve(-pos_inf, pos_inf)))

    values = yn(2, x_array)
    call check_array('yn elemental array', values, yn_array_ref, 2.0e-13_dp)
    values = kn(kn_orders, 1.0_dp)
    call check_array_rel('kn elemental order array', values, kn_array_ref, 4.0e-12_dp)
    values = iv(iv_orders, 1.0_dp)
    call check_array_rel('iv elemental order array', values, iv_array_ref, 2.0e-13_dp)

    print *, 'ordered Bessel special-function tests passed'

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
        if (ieee_is_nan(actual) .or. &
                abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' relative tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_rel

    subroutine check_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (size(actual) /= size(expected) .or. &
                any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array

    subroutine check_array_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! nonzero reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise relative error

        if (size(actual) /= size(expected) .or. &
                any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance * abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array_rel

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_bessel_order
