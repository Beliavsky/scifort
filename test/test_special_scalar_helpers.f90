! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_scalar_helpers
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : agm, dp, huber, log_ndtr, ndtri_exp, owens_t, &
        pseudo_huber, voigt_profile
    implicit none

    real(dp) :: array_actual(3)
    real(dp) :: array_input(3)
    real(dp) :: nan_value
    real(dp) :: pinf

    pinf = ieee_value(1.0_dp, ieee_positive_inf)
    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)

    call check_close('agm ordinary', agm(24.0_dp, 6.0_dp), &
        13.458171481725614_dp, 4.0e-15_dp)
    call check_close('agm negative pair', agm(-24.0_dp, -6.0_dp), &
        -13.458171481725614_dp, 4.0e-15_dp)
    call check_relative('agm disparate scales', agm(1.0e-300_dp, 1.0e300_dp), &
        1.1358405546107696e297_dp, 2.0e-15_dp)
    call check_close('agm zero', agm(0.0_dp, 2.0_dp), 0.0_dp, 0.0_dp)
    call check_true('agm mixed signs', ieee_is_nan(agm(-1.0_dp, 1.0_dp)))
    call check_true('agm infinity', agm(pinf, 2.0_dp) > huge(1.0_dp))

    call check_close('owens_t ordinary', owens_t(1.0_dp, 1.0_dp), &
        0.06674188216570097_dp, 3.0e-16_dp)
    call check_close('owens_t large slope', owens_t(0.25_dp, 1000.0_dp), &
        0.20064683715853815_dp, 5.0e-16_dp)
    call check_relative('owens_t tail', owens_t(20.0_dp, 10.0_dp), &
        1.3768120593031561e-89_dp, 7.0e-14_dp)
    call check_close('owens_t infinite slope', owens_t(0.0_dp, pinf), &
        0.25_dp, 3.0e-16_dp)
    call check_close('owens_t signs', owens_t(-1.0_dp, -2.0_dp), &
        -0.07846818699308411_dp, 4.0e-16_dp)
    call check_true('owens_t nan', ieee_is_nan(owens_t(nan_value, 1.0_dp)))

    call check_close('voigt ordinary', voigt_profile(0.0_dp, 1.0_dp, 1.0_dp), &
        0.20870928052036772_dp, 6.0e-15_dp)
    call check_close('voigt Gaussian limit', voigt_profile(1.0_dp, 1.0_dp, 0.0_dp), &
        0.24197072451914337_dp, 6.0e-15_dp)
    call check_close('voigt Cauchy limit', voigt_profile(1.0_dp, 0.0_dp, 1.0_dp), &
        0.15915494309189535_dp, 4.0e-16_dp)
    call check_close('voigt signed sigma', voigt_profile(0.0_dp, -1.0_dp, 1.0_dp), &
        -1.1067799664385465_dp, 2.0e-13_dp)
    call check_close('voigt signed gamma', voigt_profile(0.0_dp, 1.0_dp, -1.0_dp), &
        1.1067799664385465_dp, 2.0e-13_dp)
    call check_true('voigt degenerate peak', &
        voigt_profile(0.0_dp, 0.0_dp, 0.0_dp) > huge(1.0_dp))
    call check_close('voigt degenerate away', &
        voigt_profile(1.0_dp, 0.0_dp, 0.0_dp), 0.0_dp, 0.0_dp)

    call check_close('ndtri_exp ordinary', ndtri_exp(-10.0_dp), &
        -3.9139462405318932_dp, 4.0e-14_dp)
    call check_close('ndtri_exp deep tail', ndtri_exp(-1000.0_dp), &
        -44.615747731969407_dp, 4.0e-13_dp)
    call check_close('ndtri_exp extreme tail', ndtri_exp(-1.0e6_dp), &
        -1414.2077829910173_dp, 4.0e-13_dp)
    call check_close('ndtri_exp median', ndtri_exp(log(0.5_dp)), 0.0_dp, 2.0e-15_dp)
    call check_close('ndtri_exp upper', ndtri_exp(-1.0e-12_dp), &
        7.0344838253012023_dp, 3.0e-13_dp)
    call check_close('ndtri_exp round trip', log_ndtr(ndtri_exp(-100.0_dp)), &
        -100.0_dp, 3.0e-13_dp)
    call check_true('ndtri_exp zero', ndtri_exp(0.0_dp) > huge(1.0_dp))
    call check_true('ndtri_exp invalid', ieee_is_nan(ndtri_exp(1.0_dp)))

    call check_close('huber quadratic', huber(1.0_dp, 0.5_dp), 0.125_dp, 0.0_dp)
    call check_close('huber linear', huber(1.0_dp, 2.0_dp), 1.5_dp, 0.0_dp)
    call check_close('huber wide', huber(2.0_dp, 10.0_dp), 18.0_dp, 0.0_dp)
    call check_true('huber negative delta', huber(-1.0_dp, 2.0_dp) > huge(1.0_dp))

    call check_close('pseudo_huber small', pseudo_huber(1.0_dp, 0.5_dp), &
        0.11803398874989485_dp, 3.0e-16_dp)
    call check_close('pseudo_huber ordinary', pseudo_huber(1.0_dp, 2.0_dp), &
        1.2360679774997896_dp, 5.0e-16_dp)
    call check_close('pseudo_huber wide', pseudo_huber(2.0_dp, 10.0_dp), &
        16.396078054371138_dp, 4.0e-15_dp)
    call check_close('pseudo_huber tiny residual', pseudo_huber(2.0_dp, 1.0e-12_dp), &
        5.0e-25_dp, 2.0e-40_dp)
    call check_true('pseudo_huber negative delta', &
        pseudo_huber(-1.0_dp, 2.0_dp) > huge(1.0_dp))

    array_input = [-1.0_dp, 0.0_dp, 1.0_dp]
    array_actual = huber(1.0_dp, array_input)
    call check_close('huber elemental 1', array_actual(1), 0.5_dp, 0.0_dp)
    call check_close('huber elemental 2', array_actual(2), 0.0_dp, 0.0_dp)
    call check_close('huber elemental 3', array_actual(3), 0.5_dp, 0.0_dp)

    print *, 'special scalar-helper tests passed'

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

end program test_special_scalar_helpers
