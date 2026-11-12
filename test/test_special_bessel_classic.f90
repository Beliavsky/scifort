! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_bessel_classic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, i0, i1, j0, j1, k0, k0e, k1, k1e, y0, y1
    implicit none

    real(dp), parameter :: x_array(3) = [0.0_dp, 1.0_dp, 10.0_dp]
    real(dp), parameter :: j0_array_ref(3) = [ &
        1.0_dp, 0.76519768655796649_dp, -0.24593576445134832_dp]
    real(dp) :: nan_value, pos_inf
    real(dp) :: values(3)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call check_close('j0 0.1', j0(0.1_dp), 0.99750156206604002_dp, 5.0e-15_dp)
    call check_close('j0 100', j0(100.0_dp), 0.019985850304223330_dp, 5.0e-15_dp)
    call check_close('j1 -1', j1(-1.0_dp), -0.44005058574493355_dp, 5.0e-15_dp)
    call check_close('j1 700', j1(700.0_dp), 0.029489824084030461_dp, 5.0e-14_dp)
    call check_true('j0 infinity', ieee_is_nan(j0(pos_inf)))
    call check_true('j1 negative infinity', ieee_is_nan(j1(-pos_inf)))
    call check_true('j0 nan', ieee_is_nan(j0(nan_value)))

    call check_close('y0 0.1', y0(0.1_dp), -1.5342386513503667_dp, 8.0e-15_dp)
    call check_close('y0 100', y0(100.0_dp), -0.077244313365083098_dp, 8.0e-15_dp)
    call check_close('y1 1', y1(1.0_dp), -0.78121282130028880_dp, 8.0e-15_dp)
    call check_close('y1 700', y1(700.0_dp), 0.0063093414214519597_dp, 5.0e-14_dp)
    call check_true('y0 zero', y0(0.0_dp) < -huge(1.0_dp))
    call check_true('y1 negative', ieee_is_nan(y1(-1.0_dp)))
    call check_true('y0 infinity', ieee_is_nan(y0(pos_inf)))

    call check_rel('i0 1', i0(1.0_dp), 1.2660658777520082_dp, 4.0e-15_dp)
    call check_rel('i0 100', i0(100.0_dp), 1.0737517071310738e42_dp, 8.0e-14_dp)
    call check_rel('i1 -10', i1(-10.0_dp), -2670.9883037012551_dp, 8.0e-14_dp)
    call check_rel('i1 700', i1(700.0_dp), 1.5285003902339006e302_dp, 2.0e-13_dp)
    call check_rel('i0 713 high precision', i0(713.0_dp), &
        6.7051282636709967e307_dp, 3.0e-13_dp)
    call check_rel('i1 713 high precision', i1(713.0_dp), &
        6.7004245591864025e307_dp, 3.0e-13_dp)
    call check_true('i0 713 finite', ieee_is_finite(i0(713.0_dp)))
    call check_true('i0 infinity', ieee_is_nan(i0(pos_inf)))
    call check_true('i1 negative infinity', ieee_is_nan(i1(-pos_inf)))

    call check_rel('k0 0.1', k0(0.1_dp), 2.4270690247020164_dp, 3.0e-12_dp)
    call check_rel('k0 10', k0(10.0_dp), 1.7780062316167650e-5_dp, 3.0e-12_dp)
    call check_rel('k1 1', k1(1.0_dp), 0.60190723019723458_dp, 3.0e-12_dp)
    call check_rel('k1 100', k1(100.0_dp), 4.6798537356369101e-45_dp, 3.0e-12_dp)
    call check_rel('k0e 1', k0e(1.0_dp), 1.1444630798068947_dp, 3.0e-12_dp)
    call check_rel('k0e 700', k0e(700.0_dp), 0.047362369454613570_dp, 3.0e-12_dp)
    call check_rel('k1e 10', k1e(10.0_dp), 0.41076657059578869_dp, 3.0e-12_dp)
    call check_rel('k1e 700', k1e(700.0_dp), 0.047396187653494549_dp, 3.0e-12_dp)
    call check_true('k0 zero', k0(0.0_dp) > huge(1.0_dp))
    call check_true('k1 negative', ieee_is_nan(k1(-1.0_dp)))
    call check_close('k0 infinity', k0(pos_inf), 0.0_dp, 0.0_dp)
    call check_close('k1e infinity', k1e(pos_inf), 0.0_dp, 0.0_dp)

    values = j0(x_array)
    call check_array('j0 elemental array', values, j0_array_ref, 8.0e-15_dp)

    print *, 'classic Bessel special-function tests passed'

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

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_bessel_classic
