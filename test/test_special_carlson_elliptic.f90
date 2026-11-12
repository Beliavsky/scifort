! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_carlson_elliptic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, elliprc, elliprd, elliprf, elliprg, elliprj
    implicit none

    real(dp), parameter :: x_array(3) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: y_array(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: z_array(3) = [2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: rf_array_ref(3) = [ &
        1.3110287771460598_dp, 0.7269459354689082_dp, 0.5840828416771517_dp]
    real(dp), parameter :: scale_large = 1.0e200_dp
    real(dp), parameter :: scale_small = 1.0e-200_dp

    real(dp) :: base
    real(dp) :: nan_value
    real(dp) :: pos_inf
    real(dp) :: values(3)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call check_rel('elliprc(0,1)', elliprc(0.0_dp, 1.0_dp), &
        1.5707963267948963_dp, 5.0e-15_dp)
    call check_rel('elliprc positive', elliprc(1.0_dp, 2.0_dp), &
        0.7853981633974482_dp, 5.0e-15_dp)
    call check_rel('elliprc principal value', elliprc(2.0_dp, -1.0_dp), &
        0.6617680207599845_dp, 8.0e-15_dp)
    call check_rel('elliprc equal args', elliprc(7.0_dp, 7.0_dp), &
        1.0_dp / sqrt(7.0_dp), 2.0e-15_dp)

    call check_rel('elliprf(0,1,2)', elliprf(0.0_dp, 1.0_dp, 2.0_dp), &
        1.3110287771460598_dp, 5.0e-15_dp)
    call check_rel('elliprf(1,2,3)', elliprf(1.0_dp, 2.0_dp, 3.0_dp), &
        0.7269459354689082_dp, 5.0e-15_dp)
    call check_close('elliprf equal args', elliprf(4.0_dp, 4.0_dp, 4.0_dp), &
        0.5_dp, 3.0e-15_dp)
    call check_rel('elliprf symmetry', elliprf(1.0_dp, 2.0_dp, 3.0_dp), &
        elliprf(3.0_dp, 1.0_dp, 2.0_dp), 3.0e-15_dp)

    call check_rel('elliprd(0,1,2)', elliprd(0.0_dp, 1.0_dp, 2.0_dp), &
        1.0679379896673960_dp, 8.0e-15_dp)
    call check_rel('elliprd(0,2,1)', elliprd(0.0_dp, 2.0_dp, 1.0_dp), &
        1.7972103521033884_dp, 8.0e-15_dp)
    call check_rel('elliprd(1,2,3)', elliprd(1.0_dp, 2.0_dp, 3.0_dp), &
        0.2904602810289907_dp, 8.0e-15_dp)
    call check_close('elliprd equal args', elliprd(4.0_dp, 4.0_dp, 4.0_dp), &
        0.125_dp, 5.0e-15_dp)

    call check_rel('elliprj(0,1,2,3)', elliprj(0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp), &
        0.7768862377858233_dp, 2.0e-14_dp)
    call check_rel('elliprj positive', elliprj(1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp), &
        0.23984809974956783_dp, 2.0e-14_dp)
    call check_rel('elliprj principal value', elliprj(1.0_dp, 2.0_dp, 3.0_dp, -1.0_dp), &
        -0.09324045243867635_dp, 5.0e-13_dp)
    call check_close('elliprj equal args', elliprj(4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp), &
        0.125_dp, 5.0e-15_dp)
    call check_rel('elliprd equals elliprj z', elliprd(1.0_dp, 2.0_dp, 3.0_dp), &
        elliprj(1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp), 2.0e-14_dp)

    call check_rel('elliprg(0,1,2)', elliprg(0.0_dp, 1.0_dp, 2.0_dp), &
        0.9550494472569279_dp, 2.0e-14_dp)
    call check_close('elliprg two zeros', elliprg(0.0_dp, 0.0_dp, 1.0_dp), &
        0.5_dp, 0.0_dp)
    call check_rel('elliprg(1,2,3)', elliprg(1.0_dp, 2.0_dp, 3.0_dp), &
        1.4018470999908950_dp, 2.0e-14_dp)
    call check_close('elliprg equal args', elliprg(4.0_dp, 4.0_dp, 4.0_dp), &
        2.0_dp, 5.0e-15_dp)
    call check_rel('elliprg symmetry', elliprg(1.0_dp, 2.0_dp, 3.0_dp), &
        elliprg(3.0_dp, 1.0_dp, 2.0_dp), 3.0e-14_dp)

    values = elliprf(x_array, y_array, z_array)
    call check_array_rel('elliprf elemental arrays', values, rf_array_ref, 2.0e-14_dp)

    base = elliprf(1.0_dp, 2.0_dp, 3.0_dp)
    call check_rel('elliprf large homogeneity', &
        elliprf(scale_large, 2.0_dp * scale_large, 3.0_dp * scale_large), &
        base / sqrt(scale_large), 3.0e-14_dp)
    call check_rel('elliprf small homogeneity', &
        elliprf(scale_small, 2.0_dp * scale_small, 3.0_dp * scale_small), &
        base / sqrt(scale_small), 3.0e-14_dp)

    base = elliprd(1.0_dp, 2.0_dp, 3.0_dp)
    call check_rel('elliprd large homogeneity', &
        elliprd(scale_large, 2.0_dp * scale_large, 3.0_dp * scale_large), &
        base / (scale_large * sqrt(scale_large)), 5.0e-14_dp)
    call check_rel('elliprd small homogeneity', &
        elliprd(scale_small, 2.0_dp * scale_small, 3.0_dp * scale_small), &
        base / (scale_small * sqrt(scale_small)), 5.0e-14_dp)

    base = elliprj(1.0_dp, 2.0_dp, 3.0_dp, -0.5_dp)
    call check_rel('elliprj principal large homogeneity', &
        elliprj(scale_large, 2.0_dp * scale_large, 3.0_dp * scale_large, &
            -0.5_dp * scale_large), &
        base / (scale_large * sqrt(scale_large)), 2.0e-12_dp)
    call check_rel('elliprj principal small homogeneity', &
        elliprj(scale_small, 2.0_dp * scale_small, 3.0_dp * scale_small, &
            -0.5_dp * scale_small), &
        base / (scale_small * sqrt(scale_small)), 2.0e-12_dp)

    call check_true('elliprd tiny scale overflow', &
        elliprd(1.0e-250_dp, 2.0e-250_dp, 3.0e-250_dp) > huge(1.0_dp))
    call check_true('elliprj principal tiny negative overflow', &
        elliprj(1.0e-250_dp, 2.0e-250_dp, 3.0e-250_dp, -1.0e-250_dp) < -huge(1.0_dp))

    base = elliprg(1.0_dp, 2.0_dp, 3.0_dp)
    call check_rel('elliprg large homogeneity', &
        elliprg(scale_large, 2.0_dp * scale_large, 3.0_dp * scale_large), &
        base * sqrt(scale_large), 5.0e-14_dp)
    call check_rel('elliprg small homogeneity', &
        elliprg(scale_small, 2.0_dp * scale_small, 3.0_dp * scale_small), &
        base * sqrt(scale_small), 5.0e-14_dp)

    call check_true('elliprc negative x nan', ieee_is_nan(elliprc(-1.0_dp, 1.0_dp)))
    call check_true('elliprc zero y nan', ieee_is_nan(elliprc(1.0_dp, 0.0_dp)))
    call check_close('elliprc principal zero x', elliprc(0.0_dp, -1.0_dp), 0.0_dp, 0.0_dp)
    call check_true('elliprf negative nan', ieee_is_nan(elliprf(-1.0_dp, 1.0_dp, 2.0_dp)))
    call check_true('elliprf two zeros inf', elliprf(0.0_dp, 0.0_dp, 1.0_dp) > huge(1.0_dp))
    call check_true('elliprd zero z inf', elliprd(1.0_dp, 2.0_dp, 0.0_dp) > huge(1.0_dp))
    call check_true('elliprd two zeros inf', elliprd(0.0_dp, 0.0_dp, 1.0_dp) > huge(1.0_dp))
    call check_true('elliprj zero p nan', ieee_is_nan(elliprj(1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp)))
    call check_true('elliprj two zeros inf', elliprj(0.0_dp, 0.0_dp, 1.0_dp, -2.0_dp) > huge(1.0_dp))
    call check_close('elliprg all zero', elliprg(0.0_dp, 0.0_dp, 0.0_dp), 0.0_dp, 0.0_dp)
    call check_true('elliprg negative nan', ieee_is_nan(elliprg(-1.0_dp, 1.0_dp, 2.0_dp)))

    call check_close('elliprc infinity', elliprc(pos_inf, 1.0_dp), 0.0_dp, 0.0_dp)
    call check_close('elliprc negative infinity y', elliprc(1.0_dp, -pos_inf), 0.0_dp, 0.0_dp)
    call check_close('elliprf infinity', elliprf(pos_inf, 1.0_dp, 2.0_dp), 0.0_dp, 0.0_dp)
    call check_close('elliprd infinity', elliprd(1.0_dp, pos_inf, 2.0_dp), 0.0_dp, 0.0_dp)
    call check_close('elliprj infinity', elliprj(1.0_dp, 2.0_dp, 3.0_dp, pos_inf), 0.0_dp, 0.0_dp)
    call check_true('elliprg infinity', elliprg(pos_inf, 1.0_dp, 2.0_dp) > huge(1.0_dp))

    call check_true('elliprc nan', ieee_is_nan(elliprc(nan_value, 1.0_dp)))
    call check_true('elliprf nan', ieee_is_nan(elliprf(1.0_dp, nan_value, 2.0_dp)))
    call check_true('elliprd nan', ieee_is_nan(elliprd(1.0_dp, 2.0_dp, nan_value)))
    call check_true('elliprj nan', ieee_is_nan(elliprj(1.0_dp, 2.0_dp, 3.0_dp, nan_value)))
    call check_true('elliprg nan', ieee_is_nan(elliprg(1.0_dp, nan_value, 3.0_dp)))

    print *, 'Carlson elliptic special-function tests passed'

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

end program test_special_carlson_elliptic
