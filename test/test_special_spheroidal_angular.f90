! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_spheroidal_angular
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_math, only : positive_infinity
    use scifort_special, only : dp, obl_ang1, obl_ang1_cv, obl_cv, pro_ang1, &
        pro_ang1_cv, pro_cv
    implicit none

    real(dp) :: characteristic
    real(dp) :: derivative
    real(dp) :: derivative_cv
    real(dp) :: value
    real(dp) :: value_cv

    call pro_ang1(0, 0, 1.0_dp, 0.5_dp, value, derivative)
    call check_close(value, 0.9606110837921175_dp, 2.0e-12_dp, 'pro 0 0 value')
    call check_close(derivative, -0.15562274047847688_dp, 2.0e-12_dp, &
        'pro 0 0 derivative')

    call obl_ang1(0, 0, 1.0_dp, 0.5_dp, value, derivative)
    call check_close(value, 1.0441336969542285_dp, 2.0e-12_dp, 'obl 0 0 value')
    call check_close(derivative, 0.17878182305745155_dp, 2.0e-12_dp, &
        'obl 0 0 derivative')

    call pro_ang1(1, 2, 1.0_dp, 0.2_dp, value, derivative)
    call check_close(value, 0.5862149019789092_dp, 2.0e-12_dp, 'pro 1 2 value')
    call check_close(derivative, 2.7923384075869953_dp, 2.0e-12_dp, &
        'pro 1 2 derivative')

    call obl_ang1(2, 4, 3.0_dp, -0.6_dp, value, derivative)
    call check_close(value, 5.710137628039639_dp, 3.0e-12_dp, 'obl 2 4 value')
    call check_close(derivative, -29.784852753167264_dp, 3.0e-12_dp, &
        'obl 2 4 derivative')

    call pro_ang1(3, 5, 10.0_dp, 0.9_dp, value, derivative)
    call check_close(value, 4.482912834543368_dp, 5.0e-12_dp, 'pro c10 value')
    call check_close(derivative, -84.72150464959435_dp, 5.0e-12_dp, &
        'pro c10 derivative')

    call obl_ang1(3, 5, 10.0_dp, 0.9_dp, value, derivative)
    call check_close(value, 91.54576572629422_dp, 5.0e-12_dp, 'obl c10 value')
    call check_close(derivative, -69.68039095960876_dp, 5.0e-12_dp, &
        'obl c10 derivative')

    characteristic = pro_cv(0, 3, 5.0_dp)
    call pro_ang1(0, 3, 5.0_dp, -0.75_dp, value, derivative)
    call pro_ang1_cv(0, 3, 5.0_dp, characteristic, -0.75_dp, value_cv, derivative_cv)
    call check_close(value_cv, value, 2.0e-13_dp, 'pro cv value agreement')
    call check_close(derivative_cv, derivative, 2.0e-13_dp, 'pro cv derivative agreement')

    characteristic = obl_cv(4, 7, 0.1_dp)
    call obl_ang1(4, 7, 0.1_dp, 0.4_dp, value, derivative)
    call obl_ang1_cv(4, 7, 0.1_dp, characteristic, 0.4_dp, value_cv, derivative_cv)
    call check_close(value_cv, value, 2.0e-13_dp, 'obl cv value agreement')
    call check_close(derivative_cv, derivative, 2.0e-13_dp, 'obl cv derivative agreement')

    call pro_ang1(1, 2, 0.0_dp, 0.2_dp, value, derivative)
    call check_close(value, 0.5878775382679627_dp, 2.0e-14_dp, 'pro c0 value')
    call check_close(derivative, 2.8169132042006546_dp, 2.0e-14_dp, &
        'pro c0 derivative')

    call pro_ang1(2, 4, 3.0_dp, 0.6_dp, value, derivative)
    call pro_ang1(2, 4, 3.0_dp, -0.6_dp, value_cv, derivative_cv)
    call check_close(value_cv, value, 3.0e-12_dp, 'even-family parity value')
    call check_close(derivative_cv, -derivative, 3.0e-12_dp, &
        'even-family parity derivative')

    call pro_ang1(1, 2, 3.0_dp, 0.6_dp, value, derivative)
    call pro_ang1(1, 2, 3.0_dp, -0.6_dp, value_cv, derivative_cv)
    call check_close(value_cv, -value, 3.0e-12_dp, 'odd-family parity value')
    call check_close(derivative_cv, derivative, 3.0e-12_dp, &
        'odd-family parity derivative')

    call pro_ang1(0, 0, 20.0_dp, 0.0_dp, value, derivative)
    call check_close(value, 1.0_dp, 2.0e-12_dp, 'pro central normalization')
    call check_close(derivative, 0.0_dp, 0.0_dp, 'pro central derivative')

    call obl_ang1(0, 0, 20.0_dp, 0.0_dp, value, derivative)
    call check_close(value, 1.0_dp, 2.0e-12_dp, 'obl central normalization')
    call check_close(derivative, 0.0_dp, 0.0_dp, 'obl central derivative')

    call pro_ang1(-1, 2, 1.0_dp, 0.2_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'negative m')
    call pro_ang1(2, 1, 1.0_dp, 0.2_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'n below m')
    call pro_ang1(0, 0, 1.0_dp, 1.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'endpoint x')
    call obl_ang1(0, 0, positive_infinity(1.0_dp), 0.2_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'infinite c')

    print *, 'test_special_spheroidal_angular: PASS'

contains

    subroutine check_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! independent reference value
        real(dp), intent(in) :: tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must hold
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

end program test_special_spheroidal_angular
