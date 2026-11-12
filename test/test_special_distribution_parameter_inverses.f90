! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_distribution_parameter_inverses
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_special, only : bdtrik, bdtrin, betainc, betaincc, btdtria, btdtrib, &
        chdtr, chdtriv, dp, fdtr, fdtridfd, gammainc, gammaincc, gdtr, gdtrib, &
        nbdtrik, nbdtrin, ndtr, nrdtrimn, nrdtrisd, pdtrik, stdtr, stdtridf
    implicit none

    call check_close('bdtrik ordinary', bdtrik(0.2_dp, 10.0_dp, 0.4_dp), &
        2.1810661490739176_dp, 2.0e-12_dp)
    call check_close('bdtrik upper', bdtrik(0.9_dp, 10.0_dp, 0.4_dp), &
        5.4966425430249455_dp, 3.0e-12_dp)
    call check_close('bdtrik round trip', &
        betaincc(bdtrik(0.2_dp, 10.0_dp, 0.4_dp) + 1.0_dp, &
        10.0_dp - bdtrik(0.2_dp, 10.0_dp, 0.4_dp), 0.4_dp), &
        0.2_dp, 2.0e-13_dp)

    call check_close('bdtrin ordinary', bdtrin(3.0_dp, 0.2_dp, 0.4_dp), &
        12.417517722119964_dp, 6.0e-12_dp)
    call check_close('bdtrin round trip', &
        betaincc(4.0_dp, bdtrin(3.0_dp, 0.2_dp, 0.4_dp) - 3.0_dp, 0.4_dp), &
        0.2_dp, 2.0e-13_dp)

    call check_close('btdtria ordinary', btdtria(0.3_dp, 2.5_dp, 0.4_dp), &
        2.6800228127146104_dp, 2.0e-12_dp)
    call check_close('btdtrib ordinary', btdtrib(2.5_dp, 0.3_dp, 0.4_dp), &
        2.291216994664174_dp, 2.0e-12_dp)
    call check_close('btdtria round trip', &
        betainc(btdtria(0.3_dp, 2.5_dp, 0.4_dp), 2.5_dp, 0.4_dp), &
        0.3_dp, 2.0e-13_dp)
    call check_close('btdtrib round trip', &
        betainc(2.5_dp, btdtrib(2.5_dp, 0.3_dp, 0.4_dp), 0.4_dp), &
        0.3_dp, 2.0e-13_dp)

    call check_close('chdtriv ordinary', chdtriv(0.2_dp, 5.0_dp), &
        8.518942437722675_dp, 5.0e-12_dp)
    call check_close('chdtriv round trip', chdtr(chdtriv(0.2_dp, 5.0_dp), 5.0_dp), &
        0.2_dp, 2.0e-13_dp)

    call check_close('fdtridfd ordinary', fdtridfd(3.0_dp, 0.7_dp, 1.5_dp), &
        6.600824923384406_dp, 8.0e-11_dp)
    call check_close('fdtridfd round trip', &
        fdtr(3.0_dp, fdtridfd(3.0_dp, 0.7_dp, 1.5_dp), 1.5_dp), &
        0.7_dp, 3.0e-13_dp)

    call check_close('gdtrib ordinary', gdtrib(2.0_dp, 0.7_dp, 1.5_dp), &
        2.472602993054606_dp, 2.0e-12_dp)
    call check_close('gdtrib round trip', &
        gdtr(2.0_dp, gdtrib(2.0_dp, 0.7_dp, 1.5_dp), 1.5_dp), &
        0.7_dp, 2.0e-13_dp)

    call check_close('nbdtrik ordinary', nbdtrik(0.2_dp, 5.0_dp, 0.4_dp), &
        3.296971876451575_dp, 2.0e-11_dp)
    call check_close('nbdtrik round trip', &
        betainc(5.0_dp, nbdtrik(0.2_dp, 5.0_dp, 0.4_dp) + 1.0_dp, 0.4_dp), &
        0.2_dp, 2.0e-13_dp)
    call check_close('nbdtrin ordinary', nbdtrin(3.0_dp, 0.2_dp, 0.4_dp), &
        4.734052621648643_dp, 2.0e-11_dp)
    call check_close('nbdtrin round trip', &
        betainc(nbdtrin(3.0_dp, 0.2_dp, 0.4_dp), 4.0_dp, 0.4_dp), &
        0.2_dp, 2.0e-13_dp)

    call check_close('pdtrik ordinary', pdtrik(0.2_dp, 5.0_dp), &
        2.579728414025716_dp, 3.0e-12_dp)
    call check_close('pdtrik round trip', &
        gammaincc(pdtrik(0.2_dp, 5.0_dp) + 1.0_dp, 5.0_dp), &
        0.2_dp, 2.0e-13_dp)

    call check_close('stdtridf ordinary', stdtridf(0.9_dp, 2.0_dp), &
        1.7673464458150185_dp, 3.0e-12_dp)
    call check_close('stdtridf symmetry', stdtridf(0.1_dp, -2.0_dp), &
        1.7673464458150172_dp, 3.0e-12_dp)
    call check_close('stdtridf round trip', &
        stdtr(stdtridf(0.9_dp, 2.0_dp), 2.0_dp), 0.9_dp, 2.0e-13_dp)

    call check_close('nrdtrimn ordinary', nrdtrimn(0.1_dp, 2.0_dp, 3.0_dp), &
        5.563103131089201_dp, 3.0e-14_dp)
    call check_close('nrdtrisd ordinary', nrdtrisd(0.0_dp, 0.9_dp, 1.0_dp), &
        0.780304146072379_dp, 3.0e-14_dp)
    call check_close('nrdtrimn round trip', &
        ndtr((3.0_dp - nrdtrimn(0.1_dp, 2.0_dp, 3.0_dp)) / 2.0_dp), &
        0.1_dp, 2.0e-15_dp)
    call check_close('nrdtrisd round trip', &
        ndtr(1.0_dp / nrdtrisd(0.0_dp, 0.9_dp, 1.0_dp)), &
        0.9_dp, 2.0e-15_dp)

    call check_true('btdtria lower endpoint infinite', .not. ieee_is_finite(btdtria(0.0_dp, 2.0_dp, 0.5_dp)))
    call check_close('btdtrib lower endpoint', btdtrib(2.0_dp, 0.0_dp, 0.5_dp), 0.0_dp, 0.0_dp)
    call check_true('nrdtrisd invalid sign', ieee_is_nan(nrdtrisd(0.0_dp, 0.1_dp, 1.0_dp)))
    call check_true('stdtridf invalid range', ieee_is_nan(stdtridf(0.99_dp, 2.0_dp)))

    print *, 'special distribution-parameter inverse tests passed'

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

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_distribution_parameter_inverses
