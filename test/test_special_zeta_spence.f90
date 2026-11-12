! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_zeta_spence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_negative_inf, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_special, only : spence, zeta, zetac
    implicit none

    real(dp) :: nan_value
    real(dp) :: neg_inf
    real(dp) :: pos_inf
    real(dp) :: values(4)

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_dp, ieee_positive_inf)
    neg_inf = ieee_value(1.0_dp, ieee_negative_inf)

    call check_rel('spence zero', spence(0.0_dp), &
        1.64493406684822643647_dp, 3.0e-15_dp)
    call check_rel('spence tiny positive', spence(1.0e-12_dp), &
        1.64493406681959541536_dp, 4.0e-15_dp)
    call check_rel('spence tenth', spence(0.1_dp), &
        1.29971472300495872517_dp, 4.0e-15_dp)
    call check_rel('spence half', spence(0.5_dp), &
        0.58224052646501250590_dp, 4.0e-15_dp)
    call check_close('spence one', spence(1.0_dp), 0.0_dp, 0.0_dp)
    call check_rel('spence two', spence(2.0_dp), &
        -0.82246703342411321824_dp, 4.0e-15_dp)
    call check_rel('spence ten', spence(10.0_dp), &
        -3.95066377824415773045_dp, 4.0e-15_dp)
    call check_rel('spence million', spence(1.0e6_dp), &
        -97.0790852399416749023_dp, 4.0e-15_dp)
    call check_rel('spence huge finite', spence(1.0e100_dp), &
        -26511.135486458838279_dp, 4.0e-15_dp)
    call check_true('spence negative branch cut', ieee_is_nan(spence(-1.0_dp)))
    call check_true('spence positive infinity', ieee_is_nan(spence(pos_inf)))
    call check_true('spence negative infinity', ieee_is_nan(spence(neg_inf)))
    call check_true('spence nan', ieee_is_nan(spence(nan_value)))

    call check_rel('zeta two', zeta(2.0_dp), &
        1.64493406684822643647_dp, 3.0e-15_dp)
    call check_rel('zeta ten', zeta(10.0_dp), &
        1.00099457512781808534_dp, 3.0e-15_dp)
    call check_close('zeta positive infinity', zeta(pos_inf), 1.0_dp, 0.0_dp)
    call check_true('zeta pole', zeta(1.0_dp) > huge(1.0_dp))
    call check_rel('zeta half', zeta(0.5_dp), &
        -1.46035450880958681289_dp, 4.0e-15_dp)
    call check_rel('zeta tenth', zeta(0.1_dp), &
        -0.60303751985624171525_dp, 5.0e-15_dp)
    call check_close('zeta zero', zeta(0.0_dp), -0.5_dp, 0.0_dp)
    call check_rel('zeta minus half', zeta(-0.5_dp), &
        -0.20788622497735456602_dp, 5.0e-15_dp)
    call check_rel('zeta minus one', zeta(-1.0_dp), &
        -1.0_dp / 12.0_dp, 5.0e-15_dp)
    call check_close('zeta trivial zero', zeta(-20.0_dp), 0.0_dp, 0.0_dp)
    call check_rel('zeta negative noninteger', zeta(-19.5_dp), &
        33.168325785694607879_dp, 8.0e-15_dp)
    call check_rel('zeta far negative', zeta(-100.5_dp), &
        -1.2790431911215158384e78_dp, 4.0e-14_dp)
    call check_true('zeta negative infinity', ieee_is_nan(zeta(neg_inf)))
    call check_true('zeta nan', ieee_is_nan(zeta(nan_value)))

    call check_rel('zetac two', zetac(2.0_dp), &
        0.64493406684822643647_dp, 3.0e-15_dp)
    call check_rel('zetac ten', zetac(10.0_dp), &
        9.9457512781808533715e-4_dp, 4.0e-15_dp)
    call check_rel('zetac hundred', zetac(100.0_dp), &
        7.8886090522101180735e-31_dp, 4.0e-15_dp)
    call check_rel('zetac thousand range extension', zetac(1000.0_dp), &
        9.3326361850321887899e-302_dp, 5.0e-15_dp)
    call check_close('zetac trivial zero', zetac(-20.0_dp), -1.0_dp, 0.0_dp)
    call check_true('zetac pole', zetac(1.0_dp) > huge(1.0_dp))
    call check_close('zetac positive infinity', zetac(pos_inf), 0.0_dp, 0.0_dp)

    call check_rel('hurwitz zeta half', zeta(2.0_dp, 0.5_dp), &
        4.9348022005446793094_dp, 4.0e-15_dp)
    call check_rel('hurwitz zeta shifted', zeta(10.0_dp, 2.0_dp), &
        9.9457512781808533715e-4_dp, 4.0e-15_dp)
    call check_rel('hurwitz negative q even', zeta(2.0_dp, -0.5_dp), &
        8.9348022005446793094_dp, 4.0e-15_dp)
    call check_rel('hurwitz negative q odd', zeta(3.0_dp, -0.5_dp), &
        0.41439832211715999780_dp, 8.0e-15_dp)
    call check_rel('hurwitz negative q distant', zeta(4.0_dp, -3.5_dp), &
        32.464643259910412_dp, 2.0e-14_dp)
    call check_true('hurwitz noninteger negative q', ieee_is_nan(zeta(1.5_dp, -0.5_dp)))
    call check_true('hurwitz nonpositive integer q', zeta(2.0_dp, -1.0_dp) > huge(1.0_dp))
    call check_true('hurwitz s below one', ieee_is_nan(zeta(0.5_dp, 1.0_dp)))
    call check_close('hurwitz positive infinity q', zeta(2.0_dp, pos_inf), 0.0_dp, 0.0_dp)

    values = spence([0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp])
    call check_rel('spence elemental first', values(1), &
        1.64493406684822643647_dp, 3.0e-15_dp)
    call check_rel('spence elemental last', values(4), &
        -0.82246703342411321824_dp, 4.0e-15_dp)

    values = zeta([2.0_dp, 0.5_dp, -1.0_dp, -2.0_dp])
    call check_rel('zeta elemental first', values(1), &
        1.64493406684822643647_dp, 3.0e-15_dp)
    call check_close('zeta elemental last', values(4), 0.0_dp, 0.0_dp)

    print *, 'Zeta/Spence special-function tests passed'

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

end program test_special_zeta_spence
