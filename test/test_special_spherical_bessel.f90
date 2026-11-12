! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_spherical_bessel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, spherical_in, spherical_jn, spherical_kn, &
        spherical_yn
    implicit none

    integer, parameter :: orders(3) = [0, 2, 5]
    real(dp), parameter :: x_values(3) = [0.1_dp, 1.0_dp, 10.0_dp]
    real(dp), parameter :: jn_array_ref(3) = [ &
        0.9983341664682815_dp, 0.062035052011373916_dp, &
        -0.05553451162145216_dp]
    real(dp), parameter :: in_array_ref(3) = [ &
        1.0016675001984405_dp, 0.07156287012947449_dp, &
        236.83958270651144_dp]
    real(dp) :: nan_value
    real(dp) :: pos_inf
    real(dp) :: values(3)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call check_rel('spherical_jn 0 0.1', spherical_jn(0, 0.1_dp), &
        0.9983341664682815_dp, 2.0e-14_dp)
    call check_rel('spherical_jn 5 1', spherical_jn(5, 1.0_dp), &
        9.256115861125825e-5_dp, 2.0e-13_dp)
    call check_rel('spherical_jn 20 5', spherical_jn(20, 5.0_dp), &
        5.4277267607932615e-12_dp, 3.0e-13_dp)
    call check_rel('spherical_jn 20 100', spherical_jn(20, 100.0_dp), &
        0.010107671283873054_dp, 2.0e-13_dp)
    call check_rel('spherical_jn 100 100', spherical_jn(100, 100.0_dp), &
        0.010880477011438346_dp, 3.0e-13_dp)
    call check_rel('spherical_jn derivative', spherical_jn(100, 100.0_dp, .true.), &
        0.002287300435009257_dp, 4.0e-13_dp)
    call check_rel('spherical_jn negative argument', spherical_jn(5, -10.0_dp), &
        0.05553451162145216_dp, 2.0e-13_dp)
    call check_rel('spherical_jn negative derivative', spherical_jn(5, -10.0_dp, .true.), &
        -0.07226857814482038_dp, 3.0e-13_dp)
    call check_rel('spherical_jn tiny derivative n1', &
        spherical_jn(1, 1.0e-300_dp, .true.), 1.0_dp / 3.0_dp, 8.0e-14_dp)
    call check_rel('spherical_jn tiny derivative n2', &
        spherical_jn(2, 1.0e-200_dp, .true.), 1.3333333333333333e-201_dp, 1.5e-13_dp)
    call check_close('spherical_jn zero n0', spherical_jn(0, 0.0_dp), 1.0_dp, 0.0_dp)
    call check_close('spherical_jn zero n1 derivative', spherical_jn(1, 0.0_dp, .true.), &
        1.0_dp / 3.0_dp, 0.0_dp)
    call check_true('spherical_jn invalid order', ieee_is_nan(spherical_jn(-1, 1.0_dp)))
    call check_close('spherical_jn infinity', spherical_jn(4, pos_inf), 0.0_dp, 0.0_dp)

    call check_rel('spherical_yn 0 0.1', spherical_yn(0, 0.1_dp), &
        -9.950041652780257_dp, 2.0e-14_dp)
    call check_rel('spherical_yn 5 1', spherical_yn(5, 1.0_dp), &
        -999.4403433922365_dp, 2.0e-13_dp)
    call check_rel('spherical_yn 20 5', spherical_yn(20, 5.0_dp), &
        -926795140.3057541_dp, 3.0e-13_dp)
    call check_rel('spherical_yn 100 100', spherical_yn(100, 100.0_dp), &
        -0.02298385049156228_dp, 4.0e-13_dp)
    call check_rel('spherical_yn derivative', spherical_yn(20, 100.0_dp, .true.), &
        0.009892593037261125_dp, 3.0e-13_dp)
    call check_rel('spherical_yn negative argument', spherical_yn(5, -1.0_dp), &
        -999.4403433922365_dp, 2.0e-13_dp)
    call check_rel('spherical_yn negative derivative', spherical_yn(5, -1.0_dp, .true.), &
        -5883.743876138712_dp, 3.0e-13_dp)
    call check_true('spherical_yn zero', spherical_yn(2, 0.0_dp) < -huge(1.0_dp))
    call check_true('spherical_yn zero derivative', &
        spherical_yn(2, 0.0_dp, .true.) > huge(1.0_dp))
    call check_true('spherical_yn tiny overflow', &
        spherical_yn(2, 1.0e-300_dp) < -huge(1.0_dp))

    call check_rel('spherical_in 0 0.1', spherical_in(0, 0.1_dp), &
        1.0016675001984405_dp, 2.0e-14_dp)
    call check_rel('spherical_in 5 1', spherical_in(5, 1.0_dp), &
        9.996237520068258e-5_dp, 3.0e-13_dp)
    call check_rel('spherical_in 20 5', spherical_in(20, 5.0_dp), &
        9.708266441147003e-12_dp, 4.0e-13_dp)
    call check_rel('spherical_in 20 100', spherical_in(20, 100.0_dp), &
        1.6407475505991204e40_dp, 4.0e-13_dp)
    call check_rel('spherical_in 100 100', spherical_in(100, 100.0_dp), &
        3.735988741596959e20_dp, 6.0e-13_dp)
    call check_rel('spherical_in derivative', spherical_in(100, 100.0_dp, .true.), &
        5.268762728123326e20_dp, 8.0e-13_dp)
    call check_rel('spherical_in negative argument', spherical_in(5, -10.0_dp), &
        -236.83958270651144_dp, 4.0e-13_dp)
    call check_rel('spherical_in negative derivative', spherical_in(5, -10.0_dp, .true.), &
        249.4166740551066_dp, 5.0e-13_dp)
    call check_rel('spherical_in tiny derivative n1', &
        spherical_in(1, 1.0e-200_dp, .true.), 1.0_dp / 3.0_dp, 8.0e-14_dp)
    call check_rel('spherical_in tiny derivative n2', &
        spherical_in(2, 1.0e-100_dp, .true.), 1.3333333333333333e-101_dp, 8.0e-14_dp)
    call check_rel('spherical_in tiny derivative n22', &
        spherical_in(22, 1.0e-12_dp, .true.), 8.6703637264923585e-280_dp, 5.0e-14_dp)
    call check_close('spherical_in tiny n0', spherical_in(0, 1.0e-300_dp), &
        1.0_dp, 2.0e-14_dp)
    call check_true('spherical_in overflow', spherical_in(0, 1000.0_dp) > huge(1.0_dp))

    call check_rel('spherical_kn 0 0.1', spherical_kn(0, 0.1_dp), &
        14.213152925974637_dp, 3.0e-14_dp)
    call check_rel('spherical_kn 5 1', spherical_kn(5, 1.0_dp), &
        1404.7865936708652_dp, 3.0e-13_dp)
    call check_rel('spherical_kn 20 5', spherical_kn(20, 5.0_dp), &
        766744623.5383387_dp, 5.0e-13_dp)
    call check_rel('spherical_kn 20 100', spherical_kn(20, 100.0_dp), &
        4.68935421826122e-45_dp, 5.0e-13_dp)
    call check_rel('spherical_kn 100 100', spherical_kn(100, 100.0_dp), &
        1.4827965292343432e-25_dp, 1.0e-12_dp)
    call check_rel('spherical_kn derivative', spherical_kn(100, 100.0_dp, .true.), &
        -2.1133522415149233e-25_dp, 1.5e-12_dp)
    call check_rel('spherical_kn high-order finite 1', spherical_kn(150, 1.0_dp), &
        5.8857785702867412e306_dp, 1.0e-13_dp)
    call check_rel('spherical_kn high-order finite 0.99', spherical_kn(150, 0.99_dp), &
        2.6847508422762158e307_dp, 1.0e-13_dp)
    call check_rel('spherical_kn negative argument', spherical_kn(5, -10.0_dp), &
        -744.0534640570654_dp, 8.0e-13_dp)
    call check_rel('spherical_kn negative derivative', spherical_kn(5, -10.0_dp, .true.), &
        783.5656266217329_dp, 1.0e-12_dp)
    call check_close('spherical_kn cancellation zero', spherical_kn(1, -1.0_dp), &
        0.0_dp, 2.0e-14_dp)
    call check_true('spherical_kn zero', spherical_kn(2, 0.0_dp) > huge(1.0_dp))
    call check_true('spherical_kn zero derivative', &
        spherical_kn(2, 0.0_dp, .true.) < -huge(1.0_dp))
    call check_close('spherical_kn positive infinity', spherical_kn(2, pos_inf), &
        0.0_dp, 0.0_dp)

    values = spherical_jn(orders, x_values)
    call check_array_rel('spherical_jn elemental', values, jn_array_ref, 3.0e-13_dp)
    values = spherical_in(orders, x_values)
    call check_array_rel('spherical_in elemental', values, in_array_ref, 5.0e-13_dp)
    call check_true('spherical NaN', ieee_is_nan(spherical_kn(2, nan_value)))

    print *, 'spherical Bessel special-function tests passed'

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

        if (size(actual) /= size(expected) .or. any(ieee_is_nan(actual)) .or. &
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

end program test_special_spherical_bessel
