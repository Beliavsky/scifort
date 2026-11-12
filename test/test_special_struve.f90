! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_struve
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : dp, it2struve0, itmodstruve0, itstruve0, modstruve, struve
    implicit none

    real(dp) :: values(3)

    call check_rel('struve H0(1)', struve(0.0_dp, 1.0_dp), &
        0.5686566270482879509864228863223533_dp, 3.0e-14_dp)
    call check_rel('struve H1(2)', struve(1.0_dp, 2.0_dp), &
        0.6467637282835621171228332858530974_dp, 3.0e-14_dp)
    call check_rel('struve H-4(15)', struve(-4.0_dp, 15.0_dp), &
        0.1726934661655210428735302335615137_dp, 3.0e-12_dp)
    call check_rel('struve half order', struve(-0.5_dp, 3.0_dp), &
        0.06500818287737577811400469640462895_dp, 5.0e-14_dp)
    call check_rel('struve noninteger order', struve(2.3_dp, 5.0_dp), &
        1.631478056297301711743234057470048_dp, 8.0e-14_dp)
    call check_rel('struve large argument', struve(10.0_dp, 30.0_dp), &
        19567.71527019647542541333392445306_dp, 2.0e-13_dp)
    call check_rel('struve negative-order transition', struve(-7.0_dp, 20.0_dp), &
        0.004427956994919247042212481160820535_dp, 2.0e-11_dp)
    call check_rel('struve negative x parity', struve(1.0_dp, -2.0_dp), &
        0.6467637282835621171228332858530974_dp, 3.0e-14_dp)
    call check_true('struve noninteger negative x', ieee_is_nan(struve(0.5_dp, -1.0_dp)))
    call check_close('struve corrected half-integer zero', struve(-1.5_dp, 0.0_dp), 0.0_dp, 0.0_dp)

    call check_rel('modstruve L0(1)', modstruve(0.0_dp, 1.0_dp), &
        0.7102431859378908887385266778116507_dp, 3.0e-14_dp)
    call check_rel('modstruve L1(2)', modstruve(1.0_dp, 2.0_dp), &
        1.102759787367715817589721650181909_dp, 3.0e-14_dp)
    call check_rel('modstruve noninteger order', modstruve(2.3_dp, 5.0_dp), &
        14.26495178849319367651486762379355_dp, 8.0e-14_dp)
    call check_rel('modstruve cancellation-sensitive series', &
        modstruve(1.5_dp, 13.358100117952025_dp), &
        63912.25001765464985927163968398475_dp, 3.0e-13_dp)
    call check_rel('modstruve large argument', modstruve(10.0_dp, 30.0_dp), &
        145831791220.2592481964515700172662_dp, 3.0e-13_dp)

    call check_rel('itstruve0 1', itstruve0(1.0_dp), &
        0.3010904267080554841008306098194878_dp, 3.0e-14_dp)
    call check_rel('itstruve0 20', itstruve0(20.0_dp), &
        2.548451692293957187034415718094268_dp, 2.0e-13_dp)
    call check_rel('itstruve0 corrected large x', itstruve0(40.0_dp), &
        3.148417627447185790283349638852635_dp, 3.0e-13_dp)
    call check_rel('itstruve0 symmetry', itstruve0(-20.0_dp), &
        2.548451692293957187034415718094268_dp, 2.0e-13_dp)

    call check_rel('it2struve0 1', it2struve0(1.0_dp), &
        0.9571973506383524735836015484408061_dp, 3.0e-14_dp)
    call check_rel('it2struve0 20', it2struve0(20.0_dp), &
        0.04030795889403845267409310978008914_dp, 2.0e-12_dp)
    call check_rel('it2struve0 40', it2struve0(40.0_dp), &
        0.01621331363128439754017503786745591_dp, 3.0e-12_dp)
    call check_rel('it2struve0 negative x', it2struve0(-1.0_dp), &
        2.184395302951440764879041834838697_dp, 3.0e-14_dp)

    call check_rel('itmodstruve0 1', itmodstruve0(1.0_dp), &
        0.3364726286440383596450803764566356_dp, 3.0e-14_dp)
    call check_rel('itmodstruve0 10', itmodstruve0(10.0_dp), &
        2990.773297132525600100548338903940_dp, 5.0e-14_dp)
    call check_rel('itmodstruve0 symmetry', itmodstruve0(-2.0_dp), &
        1.588284995334598312102565340583887_dp, 3.0e-14_dp)

    values = struve([0.0_dp, 1.0_dp, 2.0_dp], [1.0_dp, 2.0_dp, 5.0_dp])
    call check_rel('struve elemental 1', values(1), 0.5686566270482879509864228863223533_dp, 3.0e-14_dp)
    call check_rel('struve elemental 2', values(2), 0.6467637282835621171228332858530974_dp, 3.0e-14_dp)
    call check_rel('struve elemental 3', values(3), 1.569374548040279571833068323894561_dp, 8.0e-14_dp)

    print *, 'Struve special-function tests passed'

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
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance * abs(expected)) then
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

end program test_special_struve
