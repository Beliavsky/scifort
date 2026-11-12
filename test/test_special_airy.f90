! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_airy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : airy, airye, dp
    implicit none

    real(dp), parameter :: x_array(4) = [-2.0_dp, -0.1_dp, 0.1_dp, 2.0_dp]
    real(dp), parameter :: ai_array_ref(4) = [ &
        0.22740742820168564_dp, 0.38084866812012153_dp, &
        0.32920312994353801_dp, 0.034924130423274358_dp]
    real(dp), parameter :: aip_array_ref(4) = [ &
        0.61825902074169103_dp, -0.25695811232364618_dp, &
        -0.25713042190758623_dp, -0.053090384433653881_dp]
    real(dp), parameter :: bi_array_ref(4) = [ &
        -0.41230258795639851_dp, 0.56999904300295490_dp, &
        0.65986169019418905_dp, 3.2980949999782143_dp]
    real(dp), parameter :: bip_array_ref(4) = [ &
        0.27879516692116974_dp, 0.45121336229346126_dp, &
        0.45151263114964657_dp, 4.1006820499328898_dp]
    real(dp) :: ai
    real(dp) :: ai_values(4)
    real(dp) :: aip
    real(dp) :: aip_values(4)
    real(dp) :: bi
    real(dp) :: bi_values(4)
    real(dp) :: bip
    real(dp) :: bip_values(4)
    real(dp) :: eai
    real(dp) :: eaip
    real(dp) :: ebi
    real(dp) :: ebip
    real(dp) :: nan_value
    real(dp) :: pos_inf

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call airy(0.0_dp, ai, aip, bi, bip)
    call check_close('airy Ai(0)', ai, 0.35502805388781722_dp, 2.0e-16_dp)
    call check_close('airy Aip(0)', aip, -0.25881940379280682_dp, 2.0e-16_dp)
    call check_close('airy Bi(0)', bi, 0.61492662744600068_dp, 2.0e-16_dp)
    call check_close('airy Bip(0)', bip, 0.44828835735382638_dp, 2.0e-16_dp)

    call airy(-100.0_dp, ai, aip, bi, bip)
    call check_rel('airy negative large Ai', ai, 0.17675339323955203_dp, 3.0e-12_dp)
    call check_rel('airy negative large Aip', aip, -0.24229703166065122_dp, 3.0e-12_dp)
    call check_rel('airy negative large Bi', bi, 0.024273887680166775_dp, 3.0e-12_dp)
    call check_rel('airy negative large Bip', bip, 1.7675948932340515_dp, 3.0e-12_dp)

    call airy(1.00000001_dp, ai, aip, bi, bip)
    call check_rel('airy switch Ai', ai, 0.13529241472140713_dp, 3.0e-12_dp)
    call check_rel('airy switch Aip', aip, -0.15914743994386918_dp, 3.0e-12_dp)
    call check_rel('airy switch Bi', bi, 1.2074236042772308_dp, 3.0e-12_dp)
    call check_rel('airy switch Bip', bip, 0.93243594546701158_dp, 3.0e-12_dp)

    call airy(100.0_dp, ai, aip, bi, bip)
    call check_rel('airy positive tail Ai', ai, 2.6344821520883423e-291_dp, 5.0e-12_dp)
    call check_rel('airy positive tail Aip', aip, -2.6351403616043710e-290_dp, 5.0e-12_dp)
    call check_rel('airy positive tail Bi', bi, 6.0412239966699719e288_dp, 5.0e-12_dp)
    call check_rel('airy positive tail Bip', bip, 6.0397127453103741e289_dp, 5.0e-12_dp)

    call airye(2.0_dp, eai, eaip, ebi, ebip)
    call check_rel('airye positive eAi', eai, 0.23016491865250383_dp, 5.0e-12_dp)
    call check_rel('airye positive eAip', eaip, -0.34988828258010268_dp, 5.0e-12_dp)
    call check_rel('airye positive eBi', ebi, 0.50043725430409480_dp, 5.0e-12_dp)
    call check_rel('airye positive eBip', ebip, 0.62221799731543759_dp, 5.0e-12_dp)

    call airye(1000.0_dp, eai, eaip, ebi, ebip)
    call check_rel('airye large eAi', eai, 0.050164170749970864_dp, 5.0e-12_dp)
    call check_rel('airye large eAip', eaip, -1.5863429058298848_dp, 5.0e-12_dp)
    call check_rel('airye large eBi', ebi, 0.10032900247310518_dp, 5.0e-12_dp)
    call check_rel('airye large eBip', ebip, 3.1726565491304126_dp, 5.0e-12_dp)

    call airye(-2.0_dp, eai, eaip, ebi, ebip)
    call check_true('airye negative eAi is nan', ieee_is_nan(eai))
    call check_true('airye negative eAip is nan', ieee_is_nan(eaip))
    call check_rel('airye negative eBi', ebi, -0.41230258795639885_dp, 5.0e-12_dp)
    call check_rel('airye negative eBip', ebip, 0.27879516692116863_dp, 5.0e-12_dp)

    call airy(x_array, ai_values, aip_values, bi_values, bip_values)
    call check_array_rel('airy elemental Ai', ai_values, ai_array_ref, 5.0e-12_dp)
    call check_array_rel('airy elemental Aip', aip_values, aip_array_ref, 5.0e-12_dp)
    call check_array_rel('airy elemental Bi', bi_values, bi_array_ref, 5.0e-12_dp)
    call check_array_rel('airy elemental Bip', bip_values, bip_array_ref, 5.0e-12_dp)

    call airy(1000.0_dp, ai, aip, bi, bip)
    call check_close('airy underflow Ai', ai, 0.0_dp, 0.0_dp)
    call check_close('airy underflow Aip', aip, -0.0_dp, 0.0_dp)
    call check_true('airy overflow Bi', bi > huge(1.0_dp))
    call check_true('airy overflow Bip', bip > huge(1.0_dp))

    call airye(1.0e300_dp, eai, eaip, ebi, ebip)
    call check_true('airye huge eAi finite', ieee_is_finite(eai) .and. eai > 0.0_dp)
    call check_true('airye huge eAip finite', ieee_is_finite(eaip) .and. eaip < 0.0_dp)
    call check_true('airye huge eBi finite', ieee_is_finite(ebi) .and. ebi > 0.0_dp)
    call check_true('airye huge eBip finite', ieee_is_finite(ebip) .and. ebip > 0.0_dp)
    call check_rel('airye huge eBi/eAi ratio', ebi / eai, 2.0_dp, 2.0e-14_dp)
    call check_rel('airye huge eBip/-eAip ratio', ebip / (-eaip), 2.0_dp, 2.0e-14_dp)

    call airy(pos_inf, ai, aip, bi, bip)
    call check_true('airy positive infinity Ai nan', ieee_is_nan(ai))
    call check_true('airy positive infinity Aip nan', ieee_is_nan(aip))
    call check_true('airy positive infinity Bi nan', ieee_is_nan(bi))
    call check_true('airy positive infinity Bip nan', ieee_is_nan(bip))
    call airy(-pos_inf, ai, aip, bi, bip)
    call check_true('airy negative infinity Ai nan', ieee_is_nan(ai))
    call check_true('airy negative infinity Aip nan', ieee_is_nan(aip))
    call check_true('airy negative infinity Bi nan', ieee_is_nan(bi))
    call check_true('airy negative infinity Bip nan', ieee_is_nan(bip))
    call airye(nan_value, eai, eaip, ebi, ebip)
    call check_true('airye nan eAi', ieee_is_nan(eai))
    call check_true('airye nan eAip', ieee_is_nan(eaip))
    call check_true('airye nan eBi', ieee_is_nan(ebi))
    call check_true('airye nan eBip', ieee_is_nan(ebip))

    print *, 'Airy special-function tests passed'

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

        if (size(actual) /= size(expected) .or. &
                any(ieee_is_nan(actual)) .or. &
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

end program test_special_airy
