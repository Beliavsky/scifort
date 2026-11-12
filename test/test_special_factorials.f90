! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_factorials
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_special, only : comb, factorial, factorial2, factorialk, perm
    implicit none

    real(dp) :: nan_value
    real(dp) :: pos_inf
    real(dp) :: values(4)

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_dp, ieee_positive_inf)

    call check_close('factorial integer', factorial(5), 120.0_dp, 0.0_dp)
    call check_rel('factorial noninteger', factorial(3.5_dp), &
        11.63172839656745_dp, 3.0e-15_dp)
    call check_close('factorial negative', factorial(-0.5_dp), 0.0_dp, 0.0_dp)
    call check_close('factorial negative integer', factorial(-3), 0.0_dp, 0.0_dp)
    call check_close('factorial zero', factorial(0), 1.0_dp, 0.0_dp)
    call check_true('factorial infinity', factorial(pos_inf) > huge(1.0_dp))
    call check_true('factorial nan', ieee_is_nan(factorial(nan_value)))
    call check_rel('factorial near overflow', factorial(170.5_dp), &
        9.483367566824803e307_dp, 8.0e-14_dp)
    call check_true('factorial overflow', factorial(171) > huge(1.0_dp))

    call check_close('factorial2 negative', factorial2(-5), 0.0_dp, 0.0_dp)
    call check_close('factorial2 zero', factorial2(0), 1.0_dp, 0.0_dp)
    call check_close('factorial2 odd', factorial2(7), 105.0_dp, 0.0_dp)
    call check_close('factorial2 even', factorial2(8), 384.0_dp, 0.0_dp)
    call check_rel('factorial2 large finite', factorial2(300), &
        8.154414069380595e307_dp, 5.0e-15_dp)
    call check_true('factorial2 overflow', factorial2(301) > huge(1.0_dp))

    call check_close('factorialk ordinary', factorialk(17, 4), 9945.0_dp, 0.0_dp)
    call check_close('factorialk step three', factorialk(5, 3), 10.0_dp, 0.0_dp)
    call check_close('factorialk zero', factorialk(0, 3), 1.0_dp, 0.0_dp)
    call check_close('factorialk negative', factorialk(-3, 3), 0.0_dp, 0.0_dp)
    call check_close('factorialk step beyond n', factorialk(1000, 1001), 1000.0_dp, 0.0_dp)
    call check_close('factorialk near n', factorialk(1000, 999), 1000.0_dp, 0.0_dp)
    call check_true('factorialk invalid step', ieee_is_nan(factorialk(5, 0)))
    call check_true('factorialk overflow', factorialk(171, 1) > huge(1.0_dp))

    call check_close('comb integer', comb(10, 3), 120.0_dp, 0.0_dp)
    call check_rel('comb real', comb(10.0_dp, 3.5_dp), &
        166.71933427325396_dp, 2.0e-14_dp)
    call check_close('comb k beyond n', comb(3, 4), 0.0_dp, 0.0_dp)
    call check_close('comb negative n', comb(-3, 2), 0.0_dp, 0.0_dp)
    call check_close('comb negative k', comb(3, -1), 0.0_dp, 0.0_dp)
    call check_close('comb nan nonrepetition', comb(nan_value, 2.0_dp), 0.0_dp, 0.0_dp)
    call check_true('comb infinity', comb(pos_inf, 2.0_dp) > huge(1.0_dp))
    call check_true('comb inf inf', ieee_is_nan(comb(pos_inf, pos_inf)))
    call check_close('comb repetition', comb(10, 3, repetition=.true.), 220.0_dp, 0.0_dp)
    call check_close('comb repetition zero zero', comb(0, 0, repetition=.true.), &
        1.0_dp, 0.0_dp)
    call check_close('comb repetition one two', comb(1, 2, repetition=.true.), &
        1.0_dp, 0.0_dp)
    call check_rel('comb repetition real', comb(2.2_dp, 3.1_dp, repetition=.true.), &
        5.072893716317679_dp, 2.0e-14_dp)
    call check_true('comb repetition negative integer upper', &
        ieee_is_nan(comb(-3, 2, repetition=.true.)))
    call check_true('comb repetition nan', &
        ieee_is_nan(comb(nan_value, 2.0_dp, repetition=.true.)))

    call check_close('perm integer', perm(10, 3), 720.0_dp, 0.0_dp)
    call check_rel('perm real', perm(10.0_dp, 3.5_dp), &
        1939.2340147230286_dp, 2.0e-14_dp)
    call check_close('perm k beyond n', perm(3, 4), 0.0_dp, 0.0_dp)
    call check_close('perm negative n', perm(-3, 2), 0.0_dp, 0.0_dp)
    call check_close('perm negative k', perm(3, -1), 0.0_dp, 0.0_dp)
    call check_close('perm nan', perm(nan_value, 2.0_dp), 0.0_dp, 0.0_dp)
    call check_true('perm infinity', perm(pos_inf, 2.0_dp) > huge(1.0_dp))
    call check_true('perm inf inf', ieee_is_nan(perm(pos_inf, pos_inf)))

    values = factorial([3, 4, 5, 10])
    call check_close('factorial elemental first', values(1), 6.0_dp, 0.0_dp)
    call check_close('factorial elemental last', values(4), 3628800.0_dp, 0.0_dp)

    values = factorial2([3, 4, 7, 8])
    call check_close('factorial2 elemental first', values(1), 3.0_dp, 0.0_dp)
    call check_close('factorial2 elemental last', values(4), 384.0_dp, 0.0_dp)

    values = comb([10, 10, 3, -1], [3, 4, 4, 0])
    call check_close('comb elemental first', values(1), 120.0_dp, 0.0_dp)
    call check_close('comb elemental second', values(2), 210.0_dp, 0.0_dp)
    call check_close('comb elemental last', values(4), 0.0_dp, 0.0_dp)

    values = perm([10, 10, 3, -1], [3, 4, 4, 0])
    call check_close('perm elemental first', values(1), 720.0_dp, 0.0_dp)
    call check_close('perm elemental second', values(2), 5040.0_dp, 0.0_dp)
    call check_close('perm elemental last', values(4), 0.0_dp, 0.0_dp)

    print *, 'Factorial/combinatorial special-function tests passed'

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

end program test_special_factorials
