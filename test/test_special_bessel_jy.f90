! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_bessel_jy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, jn, jv, jve, yv, yve
    implicit none

    real(dp), parameter :: order_array(3) = [0.25_dp, 2.5_dp, -1.5_dp]
    real(dp), parameter :: x_array(3) = [0.5_dp, 3.0_dp, 2.0_dp]
    real(dp), parameter :: j_array_ref(3) = [ &
        0.7416565701571465_dp, 0.4127100322097165_dp, -0.39562328135870356_dp]
    real(dp), parameter :: y_array_ref(3) = [ &
        -0.7568435456944942_dp, -0.369040730073798_dp, -0.49129377868716273_dp]
    real(dp) :: nan_value
    real(dp) :: pos_inf
    real(dp) :: values(3)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call check_rel('jv half order', jv(0.5_dp, 1.0_dp), &
        0.67139670714180388_dp, 3.0e-13_dp)
    call check_rel('jv middle recurrence', jv(20.75_dp, 5.0_dp), &
        5.6989999348203777e-12_dp, 5.0e-13_dp)
    call check_rel('jv large argument', jv(100.5_dp, 1000.0_dp), &
        -0.0063901639529748908_dp, 4.0e-12_dp)
    call check_rel('jv negative half order', jv(-1.5_dp, 1.0_dp), &
        -1.1024955751601793_dp, 5.0e-13_dp)
    call check_rel('jv near negative integer below', jv(-0.999999_dp, 1.0e-6_dp), &
        1.9999716373727692_dp, 2.0e-12_dp)
    call check_rel('jv near negative integer above', jv(-1.000001_dp, 1.0e-6_dp), &
        -2.0000283629190165_dp, 2.0e-12_dp)
    call check_rel('jv negative x positive integer order', jv(3.0_dp, -2.0_dp), &
        -0.12894324947440208_dp, 3.0e-13_dp)
    call check_rel('jv negative x negative integer order', jv(-3.0_dp, -2.0_dp), &
        0.12894324947440208_dp, 3.0e-13_dp)
    call check_true('jv negative x noninteger', ieee_is_nan(jv(2.5_dp, -2.0_dp)))

    call check_rel('yv small fractional', yv(0.25_dp, 0.1_dp), &
        -1.911768321207173_dp, 8.0e-13_dp)
    call check_rel('yv middle fractional', yv(5.5_dp, 10.0_dp), &
        0.23675446066584144_dp, 8.0e-13_dp)
    call check_rel('yv large argument', yv(100.5_dp, 1000.0_dp), &
        -0.024474988493919471_dp, 4.0e-12_dp)
    call check_rel('yv negative half order', yv(-1.5_dp, 1.0_dp), &
        -0.24029783912342725_dp, 8.0e-13_dp)
    call check_rel('yv near negative integer', yv(-0.999999_dp, 1.0e-6_dp), &
        636610.90339985467_dp, 2.0e-12_dp)
    call check_rel('yv half-order cancellation', yv(-100.5_dp, 0.1_dp), &
        1.8829894563603392e-290_dp, 2.0e-12_dp)
    call check_true('jv reflected overflow', jv(-100.5_dp, 1.0e-12_dp) > huge(1.0_dp))
    call check_close('yv reflected underflow', yv(-100.5_dp, 1.0e-12_dp), 0.0_dp, 0.0_dp)
    call check_true('yv negative x', ieee_is_nan(yv(2.0_dp, -1.0_dp)))

    call check_rel('jve real scaling', jve(2.5_dp, 10.0_dp), &
        0.19665848358181753_dp, 5.0e-13_dp)
    call check_rel('yve real scaling', yve(0.75_dp, 1.0_dp), &
        -0.62186941744297486_dp, 8.0e-13_dp)
    call check_close('jn aliases jv', jn(2.5_dp, 10.0_dp), jv(2.5_dp, 10.0_dp), 0.0_dp)

    call check_close('jv zero order at zero', jv(0.0_dp, 0.0_dp), 1.0_dp, 0.0_dp)
    call check_close('jv positive order at zero', jv(0.5_dp, 0.0_dp), 0.0_dp, 0.0_dp)
    call check_true('jv negative half order at zero', jv(-0.5_dp, 0.0_dp) > huge(1.0_dp))
    call check_true('jve negative noninteger at zero', ieee_is_nan(jve(-0.5_dp, 0.0_dp)))
    call check_true('yv zero singular', yv(2.0_dp, 0.0_dp) < -huge(1.0_dp))
    call check_true('yve zero is nan', ieee_is_nan(yve(2.0_dp, 0.0_dp)))
    call check_true('jv infinite argument', ieee_is_nan(jv(0.5_dp, pos_inf)))
    call check_close('jv infinite order', jv(pos_inf, 1.0_dp), 0.0_dp, 0.0_dp)
    call check_true('jve infinite order', ieee_is_nan(jve(pos_inf, 1.0_dp)))
    call check_true('yv infinite order', ieee_is_nan(yv(pos_inf, 1.0_dp)))
    call check_true('jv nan order', ieee_is_nan(jv(nan_value, 1.0_dp)))
    call check_true('yv nan argument', ieee_is_nan(yv(1.0_dp, nan_value)))

    values = jv(order_array, x_array)
    call check_array_rel('jv elemental arrays', values, j_array_ref, 8.0e-13_dp)
    values = yv(order_array, x_array)
    call check_array_rel('yv elemental arrays', values, y_array_ref, 8.0e-13_dp)

    print *, 'real-order J/Y Bessel special-function tests passed'

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

end program test_special_bessel_jy
