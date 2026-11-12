! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_gamma_combinatorics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_negative_inf, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_special, only : beta, betaln, binom, poch, polygamma
    implicit none

    real(dp) :: nan_value
    real(dp) :: neg_inf
    real(dp) :: pos_inf
    real(dp) :: values(4)

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_dp, ieee_positive_inf)
    neg_inf = ieee_value(1.0_dp, ieee_negative_inf)

    call check_rel('polygamma trigamma one', polygamma(1, 1.0_dp), &
        1.6449340668482266_dp, 3.0e-15_dp)
    call check_rel('polygamma second one', polygamma(2, 1.0_dp), &
        -2.404113806319188_dp, 4.0e-15_dp)
    call check_rel('polygamma third half', polygamma(3, 0.5_dp), &
        97.40909103400242_dp, 4.0e-15_dp)
    call check_rel('polygamma negative trigamma', polygamma(1, -0.5_dp), &
        8.934802200544679_dp, 8.0e-15_dp)
    call check_close('polygamma negative second', polygamma(2, -0.5_dp), &
        -0.82879664423432_dp, 2.0e-13_dp)
    call check_close('polygamma reflected fourth', polygamma(4, -1.5_dp), &
        -0.3137559995067314_dp, 2.0e-12_dp)
    call check_rel('polygamma high positive order', polygamma(10, 3.25_dp), &
        -8.98832431208088_dp, 3.0e-13_dp)
    call check_rel('polygamma far positive', polygamma(20, 100.0_dp), &
        -1.3423500975739656e-23_dp, 2.0e-12_dp)
    call check_close('polygamma far negative reflection', polygamma(8, -20.5_dp), &
        -1.3235280212547434e-7_dp, 2.0e-20_dp)
    call check_true('polygamma invalid order', ieee_is_nan(polygamma(-1, 1.0_dp)))
    call check_true('polygamma odd pole positive inf', polygamma(1, -2.0_dp) > huge(1.0_dp))
    call check_true('polygamma even pole negative inf', polygamma(2, -2.0_dp) < -huge(1.0_dp))
    call check_close('polygamma infinity signed zero', polygamma(2, pos_inf), -0.0_dp, 0.0_dp)
    call check_true('polygamma nan', ieee_is_nan(polygamma(1, nan_value)))

    call check_rel('poch half half', poch(0.5_dp, 0.5_dp), &
        0.5641895835477564_dp, 4.0e-15_dp)
    call check_rel('poch negative noninteger', poch(-0.5_dp, 1.5_dp), &
        -0.28209479177387814_dp, 4.0e-15_dp)
    call check_close('poch negative integer one', poch(-2.5_dp, 1.0_dp), -2.5_dp, 0.0_dp)
    call check_close('poch negative integer two', poch(-2.5_dp, 2.0_dp), 3.75_dp, 0.0_dp)
    call check_close('poch integer cancellation', poch(-2.0_dp, 3.0_dp), 0.0_dp, 0.0_dp)
    call check_close('poch negative shift', poch(-3.0_dp, -2.0_dp), 0.05_dp, 2.0e-16_dp)
    call check_close('poch positive negative shift', poch(3.0_dp, -2.0_dp), 0.5_dp, 0.0_dp)
    call check_close('poch pole denominator', poch(-3.0_dp, 0.5_dp), 0.0_dp, 0.0_dp)
    call check_rel('poch general negative order', poch(10.25_dp, -3.5_dp), &
        7.086790738887486e-4_dp, 2.0e-14_dp)
    call check_rel('poch large', poch(100.5_dp, 20.25_dp), &
        2.1666995375026756e41_dp, 3.0e-13_dp)
    call check_true('poch numerator pole', poch(-0.5_dp, 0.5_dp) > huge(1.0_dp))
    call check_true('poch negative-order pole', poch(3.0_dp, -3.0_dp) > huge(1.0_dp))
    call check_close('poch zero order nan convention', poch(nan_value, 0.0_dp), 1.0_dp, 0.0_dp)

    call check_rel('beta ordinary', beta(2.0_dp, 3.0_dp), &
        1.0_dp / 12.0_dp, 3.0e-15_dp)
    call check_rel('beta half half', beta(0.5_dp, 0.5_dp), &
        3.1415926535897927_dp, 3.0e-15_dp)
    call check_rel('beta negative continuation', beta(-0.5_dp, 1.5_dp), &
        -3.1415926535897927_dp, 5.0e-15_dp)
    call check_true('beta signed zero continuation', sign(1.0_dp, beta(-0.5_dp, 0.5_dp)) < 0.0_dp)
    call check_close('beta integer cancellation one', beta(-2.0_dp, 1.0_dp), -0.5_dp, 0.0_dp)
    call check_close('beta integer cancellation two', beta(-2.0_dp, 2.0_dp), 0.5_dp, 0.0_dp)
    call check_rel('beta integer cancellation three', beta(-4.0_dp, 3.0_dp), &
        -1.0_dp / 12.0_dp, 3.0e-15_dp)
    call check_rel('beta negative general', beta(-2.25_dp, 3.5_dp), &
        -6.390081251923158_dp, 2.0e-14_dp)
    call check_true('beta uncancelled pole', beta(-2.0_dp, 3.0_dp) > huge(1.0_dp))
    call check_close('beta positive infinity limit', beta(pos_inf, 2.0_dp), 0.0_dp, 0.0_dp)
    call check_true('beta nan', ieee_is_nan(beta(nan_value, 1.0_dp)))

    call check_rel('betaln negative continuation', betaln(-0.5_dp, 1.5_dp), &
        1.1447298858494_dp, 4.0e-15_dp)
    call check_true('betaln zero continuation', betaln(-0.5_dp, 0.5_dp) < -huge(1.0_dp))
    call check_rel('betaln integer cancellation', betaln(-4.0_dp, 3.0_dp), &
        -2.4849066497880004_dp, 4.0e-15_dp)
    call check_true('betaln uncancelled pole', betaln(-2.0_dp, 3.0_dp) > huge(1.0_dp))
    call check_true('betaln positive infinity limit', betaln(pos_inf, 2.0_dp) < -huge(1.0_dp))
    call check_true('betaln negative infinity', betaln(neg_inf, 2.0_dp) > huge(1.0_dp))

    call check_close('binom ordinary integer', binom(5.0_dp, 2.0_dp), 10.0_dp, 0.0_dp)
    call check_rel('binom generalized ordinary', binom(5.5_dp, 6.0_dp), &
        0.2255859375_dp, 3.0e-15_dp)
    call check_rel('binom negative upper noninteger', binom(-2.5_dp, 3.0_dp), &
        -6.5625_dp, 4.0e-15_dp)
    call check_rel('binom general real', binom(10.25_dp, 3.5_dp), &
        184.21552615039906_dp, 3.0e-14_dp)
    call check_rel('binom large real', binom(100.5_dp, 20.25_dp), &
        8.424574429798476e20_dp, 4.0e-13_dp)
    call check_close('binom negative lower integer', binom(3.0_dp, -1.0_dp), 0.0_dp, 0.0_dp)
    call check_close('binom integer beyond support', binom(3.0_dp, 4.0_dp), 0.0_dp, 0.0_dp)
    call check_true('binom negative integer upper', ieee_is_nan(binom(-3.0_dp, 2.0_dp)))
    call check_true('binom positive infinity', binom(pos_inf, 2.0_dp) > huge(1.0_dp))

    values = beta([-0.5_dp, 0.5_dp, 2.0_dp, 5.0_dp], &
        [1.5_dp, 0.5_dp, 3.0_dp, 1.0_dp])
    call check_rel('beta elemental first', values(1), -3.1415926535897927_dp, 5.0e-15_dp)
    call check_close('beta elemental last', values(4), 0.2_dp, 2.0e-16_dp)

    values = poch([0.5_dp, -2.5_dp, 3.0_dp, 10.25_dp], &
        [0.5_dp, 2.0_dp, -2.0_dp, -3.5_dp])
    call check_rel('poch elemental first', values(1), 0.5641895835477564_dp, 4.0e-15_dp)
    call check_rel('poch elemental last', values(4), 7.086790738887486e-4_dp, 2.0e-14_dp)

    values = polygamma(1, [0.5_dp, 1.0_dp, 2.0_dp, 10.0_dp])
    call check_rel('polygamma elemental first', values(1), 4.93480220054468_dp, 4.0e-15_dp)
    call check_rel('polygamma elemental last', values(4), 0.10516633568168575_dp, 4.0e-15_dp)

    print *, 'Gamma combinatorics special-function tests passed'

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

end program test_special_gamma_combinatorics
