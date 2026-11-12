! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_exponential_integrals
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_special, only : exp1, expi, expn, shichi, sici
    implicit none

    real(dp) :: chi
    real(dp) :: ci
    real(dp) :: nan_value
    real(dp) :: pos_inf
    real(dp) :: shi
    real(dp) :: si
    real(dp) :: x_values(4)
    real(dp) :: y_values(4)
    real(dp) :: z_values(4)

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_dp, ieee_positive_inf)

    call check_rel('exp1 tiny', exp1(1.0e-300_dp), &
        690.1983122333121_dp, 3.0e-15_dp)
    call check_rel('exp1 ordinary', exp1(0.5_dp), &
        0.5597735947761608_dp, 3.0e-15_dp)
    call check_rel('exp1 large', exp1(100.0_dp), &
        3.683597761682032e-46_dp, 5.0e-15_dp)

    call check_rel('expi negative', expi(-10.0_dp), &
        -4.156968929685325e-6_dp, 5.0e-15_dp)
    call check_rel('expi ordinary', expi(10.0_dp), &
        2492.2289762418773_dp, 5.0e-15_dp)
    call check_rel('expi asymptotic', expi(500.0_dp), &
        2.8128213978862943375e214_dp, 8.0e-14_dp)
    call check_rel('expi extended finite range', expi(715.0_dp), &
        4.6436256703705711004e307_dp, 8.0e-14_dp)

    call check_rel('expn n0', expn(0, 0.5_dp), &
        1.2130613194252668_dp, 3.0e-15_dp)
    call check_rel('expn n2', expn(2, 0.5_dp), &
        0.3266438623245532_dp, 3.0e-15_dp)
    call check_rel('expn high order', expn(5000, 0.5_dp), &
        1.2131826134156755e-4_dp, 5.0e-15_dp)
    call check_rel('expn large x', expn(10, 20.0_dp), &
        6.944390555412316e-11_dp, 5.0e-15_dp)

    call sici(1.0_dp, si, ci)
    call check_rel('sici si ordinary', si, 0.9460830703671831_dp, 3.0e-15_dp)
    call check_rel('sici ci ordinary', ci, 0.33740392290096816_dp, 5.0e-15_dp)
    call sici(-10.0_dp, si, ci)
    call check_rel('sici si negative', si, -1.658347594218874_dp, 3.0e-15_dp)
    call check_rel('sici ci negative', ci, -0.04545643300445537_dp, 8.0e-14_dp)
    call sici(16.0_dp, si, ci)
    call check_rel('sici si breakpoint', si, 1.6313022682700327_dp, 3.0e-15_dp)
    call check_rel('sici ci breakpoint', ci, -0.01420019012019002_dp, 8.0e-14_dp)
    call sici(1.0e6_dp, si, ci)
    call check_rel('sici si large', si, 1.570795390043119_dp, 3.0e-15_dp)
    call check_rel('sici ci large', ci, -3.499944389227205e-7_dp, 8.0e-14_dp)

    call shichi(1.0_dp, shi, chi)
    call check_rel('shichi shi ordinary', shi, 1.0572508753757286_dp, 3.0e-15_dp)
    call check_rel('shichi chi ordinary', chi, 0.8378669409802082_dp, 3.0e-15_dp)
    call shichi(-10.0_dp, shi, chi)
    call check_rel('shichi shi negative', shi, -1246.1144901994235_dp, 5.0e-15_dp)
    call check_rel('shichi chi negative', chi, 1246.1144860424545_dp, 5.0e-15_dp)
    call shichi(700.0_dp, shi, chi)
    call check_rel('shichi shi large', shi, 7.254893680262804e300_dp, 1.0e-13_dp)
    call check_rel('shichi chi large', chi, 7.254893680262804e300_dp, 1.0e-13_dp)

    call check_true('exp1 negative invalid', ieee_is_nan(exp1(-1.0_dp)))
    call check_true('exp1 zero infinity', exp1(0.0_dp) > huge(1.0_dp))
    call check_close('exp1 infinity limit', exp1(pos_inf), 0.0_dp, 0.0_dp)
    call check_true('exp1 nan', ieee_is_nan(exp1(nan_value)))

    call check_true('expi zero negative infinity', expi(0.0_dp) < -huge(1.0_dp))
    call check_close('expi negative infinity limit', expi(-pos_inf), -0.0_dp, 0.0_dp)
    call check_true('expi positive infinity', expi(pos_inf) > huge(1.0_dp))
    call check_true('expi nan', ieee_is_nan(expi(nan_value)))

    call check_true('expn negative order invalid', ieee_is_nan(expn(-1, 1.0_dp)))
    call check_true('expn negative x invalid', ieee_is_nan(expn(2, -1.0_dp)))
    call check_true('expn n1 zero infinity', expn(1, 0.0_dp) > huge(1.0_dp))
    call check_close('expn n2 zero', expn(2, 0.0_dp), 1.0_dp, 0.0_dp)
    call check_close('expn infinity limit', expn(20, pos_inf), 0.0_dp, 0.0_dp)

    call sici(0.0_dp, si, ci)
    call check_close('sici zero si', si, 0.0_dp, 0.0_dp)
    call check_true('sici zero ci', ci < -huge(1.0_dp))
    call sici(pos_inf, si, ci)
    call check_close('sici positive infinity si', si, 0.5_dp * scifort_pi, 0.0_dp)
    call check_close('sici positive infinity ci', ci, 0.0_dp, 0.0_dp)
    call sici(-pos_inf, si, ci)
    call check_close('sici negative infinity si', si, -0.5_dp * scifort_pi, 0.0_dp)
    call check_true('sici negative infinity ci scipy behavior', ieee_is_nan(ci))
    call sici(nan_value, si, ci)
    call check_true('sici nan si', ieee_is_nan(si))
    call check_true('sici nan ci', ieee_is_nan(ci))

    call shichi(0.0_dp, shi, chi)
    call check_close('shichi zero shi', shi, 0.0_dp, 0.0_dp)
    call check_true('shichi zero chi', chi < -huge(1.0_dp))
    call shichi(pos_inf, shi, chi)
    call check_true('shichi infinity shi', shi > huge(1.0_dp))
    call check_true('shichi infinity chi', chi > huge(1.0_dp))

    x_values = [0.1_dp, 0.5_dp, 1.0_dp, 10.0_dp]
    y_values = exp1(x_values)
    call check_rel('exp1 elemental first', y_values(1), 1.8229239584193906_dp, 3.0e-15_dp)
    call check_rel('exp1 elemental last', y_values(4), 4.156968929685325e-6_dp, 5.0e-15_dp)

    y_values = expn([1, 2, 3, 4], x_values)
    call check_rel('expn elemental first', y_values(1), 1.8229239584193906_dp, 3.0e-15_dp)
    call check_rel('expn elemental second', y_values(2), 0.3266438623245532_dp, 3.0e-15_dp)

    call sici(x_values, y_values, z_values)
    call check_rel('sici elemental si', y_values(3), 0.9460830703671831_dp, 3.0e-15_dp)
    call check_rel('sici elemental ci', z_values(3), 0.33740392290096816_dp, 5.0e-15_dp)

    print *, 'Exponential/trigonometric integral special-function tests passed'

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

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_exponential_integrals
