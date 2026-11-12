! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_bessel_derivatives
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, iv, ivp, jv, jvp, kv, kvp, yv, yvp
    implicit none

    real(dp), parameter :: order_array(3) = [0.25_dp, 2.5_dp, -1.5_dp]
    real(dp), parameter :: x_array(3) = [0.5_dp, 3.0_dp, 2.0_dp]
    real(dp), parameter :: j_array_ref(3) = [ &
        -0.99443430764124674_dp, -0.17070357259041180_dp, -0.092666400118205261_dp]
    real(dp), parameter :: y_array_ref(3) = [ &
        -2.4187042551713418_dp, -0.018751565516006208_dp, 0.28721392944886170_dp]
    real(dp), parameter :: i_array_ref(3) = [ &
        -0.11588261581372922_dp, 1.9554249787142197_dp, 0.84702748009914020_dp]
    real(dp), parameter :: k_array_ref(3) = [ &
        4.7442113635325640_dp, 0.18179779534084911_dp, 0.40853803576620934_dp]
    real(dp) :: nan_value
    real(dp) :: pos_inf
    real(dp) :: values(3)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call check_close('jvp n0', jvp(0.0_dp, 1.0_dp, 0), jv(0.0_dp, 1.0_dp), 0.0_dp)
    call check_rel('jvp default n1', jvp(0.0_dp, 1.0_dp), &
        -0.44005058574493355_dp, 4.0e-13_dp)
    call check_rel('jvp n2', jvp(0.0_dp, 1.0_dp, 2), &
        -0.32514710081303311_dp, 4.0e-13_dp)
    call check_rel('jvp fractional n4', jvp(2.5_dp, 3.0_dp, 4), &
        0.12489610571047538_dp, 2.0e-12_dp)
    call check_rel('jvp negative order', jvp(-1.5_dp, 2.0_dp, 2), &
        -0.092666400118205261_dp, 2.0e-12_dp)
    call check_rel('jvp negative x integer order', jvp(3.0_dp, -2.0_dp, 2), &
        -0.081469484640985287_dp, 1.0e-12_dp)
    call check_rel('jvp n20', jvp(0.0_dp, 1.0_dp, 20), &
        0.098607264095423350_dp, 2.0e-12_dp)
    call check_rel('jvp n1000 normalized weights', jvp(0.0_dp, 1.0_dp, 1000), &
        0.013639730196363853_dp, 5.0e-11_dp)

    call check_close('yvp n0', yvp(0.0_dp, 1.0_dp, 0), yv(0.0_dp, 1.0_dp), 0.0_dp)
    call check_rel('yvp default n1', yvp(0.0_dp, 1.0_dp), &
        0.78121282130028891_dp, 8.0e-13_dp)
    call check_rel('yvp n2', yvp(0.0_dp, 1.0_dp, 2), &
        -0.86946978551596588_dp, 8.0e-13_dp)
    call check_rel('yvp fractional n4', yvp(2.5_dp, 3.0_dp, 4), &
        -0.59729366029768727_dp, 3.0e-12_dp)
    call check_rel('yvp negative order', yvp(-1.5_dp, 2.0_dp, 2), &
        0.28721392944886170_dp, 3.0e-12_dp)
    call check_rel('yvp large argument n3', yvp(5.5_dp, 10.0_dp, 3), &
        0.093248131373937782_dp, 3.0e-12_dp)

    call check_close('ivp n0', ivp(0.0_dp, 1.0_dp, 0), iv(0.0_dp, 1.0_dp), 0.0_dp)
    call check_rel('ivp default n1', ivp(0.0_dp, 1.0_dp), &
        0.56515910399248503_dp, 4.0e-13_dp)
    call check_rel('ivp n2', ivp(0.0_dp, 1.0_dp, 2), &
        0.70090677375952326_dp, 4.0e-13_dp)
    call check_rel('ivp fractional n4', ivp(2.5_dp, 3.0_dp, 4), &
        1.9584398055585808_dp, 2.0e-12_dp)
    call check_rel('ivp negative order', ivp(-1.5_dp, 2.0_dp, 2), &
        0.84702748009914020_dp, 2.0e-12_dp)
    call check_rel('ivp negative x integer order', ivp(3.0_dp, -2.0_dp, 2), &
        -0.50648561311004148_dp, 2.0e-12_dp)
    call check_rel('ivp n20', ivp(0.0_dp, 1.0_dp, 20), &
        0.26722594163768459_dp, 3.0e-12_dp)

    call check_close('kvp n0', kvp(0.0_dp, 1.0_dp, 0), kv(0.0_dp, 1.0_dp), 0.0_dp)
    call check_rel('kvp default n1', kvp(0.0_dp, 1.0_dp), &
        -0.60190723019723458_dp, 4.0e-12_dp)
    call check_rel('kvp n2', kvp(0.0_dp, 1.0_dp, 2), &
        1.0229316684379428_dp, 4.0e-12_dp)
    call check_rel('kvp fractional n4', kvp(2.5_dp, 3.0_dp, 4), &
        0.60968976226465488_dp, 5.0e-12_dp)
    call check_rel('kvp negative order', kvp(-1.5_dp, 2.0_dp, 2), &
        0.40853803576620934_dp, 5.0e-12_dp)
    call check_rel('kvp large argument n3', kvp(5.5_dp, 10.0_dp, 3), &
        -0.00012835232117489479_dp, 5.0e-12_dp)
    call check_rel('kvp n20', kvp(0.0_dp, 1.0_dp, 20), &
        1.2182343645771726e17_dp, 8.0e-12_dp)

    values = jvp(order_array, x_array, 2)
    call check_array_rel('jvp elemental arrays', values, j_array_ref, 4.0e-12_dp)
    values = yvp(order_array, x_array, 2)
    call check_array_rel('yvp elemental arrays', values, y_array_ref, 5.0e-12_dp)
    values = ivp(order_array, x_array, 2)
    call check_array_rel('ivp elemental arrays', values, i_array_ref, 5.0e-12_dp)
    values = kvp(order_array, x_array, 2)
    call check_array_rel('kvp elemental arrays', values, k_array_ref, 6.0e-12_dp)

    call check_close('jvp zero first', jvp(0.0_dp, 0.0_dp, 1), 0.0_dp, 0.0_dp)
    call check_close('jvp zero second', jvp(0.0_dp, 0.0_dp, 2), -0.5_dp, 0.0_dp)
    call check_true('yvp zero first positive infinity', yvp(0.0_dp, 0.0_dp, 1) > huge(1.0_dp))
    call check_true('yvp zero second nan', ieee_is_nan(yvp(0.0_dp, 0.0_dp, 2)))
    call check_close('ivp zero first', ivp(0.0_dp, 0.0_dp, 1), 0.0_dp, 0.0_dp)
    call check_close('ivp zero second', ivp(0.0_dp, 0.0_dp, 2), 0.5_dp, 0.0_dp)
    call check_true('kvp zero first negative infinity', kvp(0.0_dp, 0.0_dp, 1) < -huge(1.0_dp))
    call check_true('kvp zero second positive infinity', kvp(0.0_dp, 0.0_dp, 2) > huge(1.0_dp))

    call check_true('jvp negative derivative order', ieee_is_nan(jvp(0.0_dp, 1.0_dp, -1)))
    call check_true('yvp negative derivative order', ieee_is_nan(yvp(0.0_dp, 1.0_dp, -1)))
    call check_true('ivp negative derivative order', ieee_is_nan(ivp(0.0_dp, 1.0_dp, -1)))
    call check_true('kvp negative derivative order', ieee_is_nan(kvp(0.0_dp, 1.0_dp, -1)))
    call check_true('jvp negative x noninteger', ieee_is_nan(jvp(0.5_dp, -1.0_dp, 1)))
    call check_true('yvp negative x', ieee_is_nan(yvp(1.0_dp, -1.0_dp, 1)))
    call check_true('ivp negative x noninteger', ieee_is_nan(ivp(0.5_dp, -1.0_dp, 1)))
    call check_true('kvp negative x', ieee_is_nan(kvp(1.0_dp, -1.0_dp, 1)))
    call check_true('jvp infinity', ieee_is_nan(jvp(0.0_dp, pos_inf, 1)))
    call check_true('yvp infinity', ieee_is_nan(yvp(0.0_dp, pos_inf, 1)))
    call check_true('ivp infinity', ieee_is_nan(ivp(0.0_dp, pos_inf, 1)))
    call check_close('kvp infinity odd', kvp(0.0_dp, pos_inf, 1), -0.0_dp, 0.0_dp)
    call check_true('jvp nan order', ieee_is_nan(jvp(nan_value, 1.0_dp, 1)))
    call check_true('kvp nan argument', ieee_is_nan(kvp(0.0_dp, nan_value, 1)))

    print *, 'Bessel derivative special-function tests passed'

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

end program test_special_bessel_derivatives
