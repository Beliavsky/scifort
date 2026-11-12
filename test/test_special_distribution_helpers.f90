! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_distribution_helpers
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_value
    use scifort_special, only : bdtr, bdtrc, bdtri, chdtr, chdtrc, chdtri, dp, &
        fdtr, fdtrc, fdtri, gdtr, gdtrc, gdtria, gdtrix, nbdtr, nbdtrc, nbdtri, &
        pdtr, pdtrc, pdtri, stdtr, stdtrit
    implicit none

    real(dp) :: pinf

    pinf = ieee_value(1.0_dp, ieee_positive_inf)

    call check_close('bdtr ordinary', bdtr(2.2_dp, 5, 0.3_dp), &
        0.8369199999999999_dp, 2.0e-15_dp)
    call check_close('bdtrc ordinary', bdtrc(2.2_dp, 5, 0.3_dp), &
        0.16307999999999995_dp, 2.0e-15_dp)
    call check_close('bdtri ordinary', bdtri(2.2_dp, 5, 0.7_dp), &
        0.38981835053095826_dp, 3.0e-14_dp)
    call check_close('bdtri round trip', bdtr(2.2_dp, 5, bdtri(2.2_dp, 5, 0.7_dp)), &
        0.7_dp, 4.0e-14_dp)
    call check_close('bdtr upper support', bdtr(5.9_dp, 5, 0.3_dp), 1.0_dp, 0.0_dp)
    call check_true('bdtr invalid support', ieee_is_nan(bdtr(6.0_dp, 5, 0.3_dp)))

    call check_close('nbdtr ordinary', nbdtr(2, 3, 0.4_dp), &
        0.3174400000000003_dp, 2.0e-15_dp)
    call check_close('nbdtrc ordinary', nbdtrc(2, 3, 0.4_dp), &
        0.6825599999999997_dp, 2.0e-15_dp)
    call check_close('nbdtri ordinary', nbdtri(2, 3, 0.7_dp), &
        0.6101816494690417_dp, 3.0e-14_dp)
    call check_close('nbdtri round trip', nbdtr(2, 3, nbdtri(2, 3, 0.7_dp)), &
        0.7_dp, 4.0e-14_dp)

    call check_close('pdtr ordinary', pdtr(2.2_dp, 3.0_dp), &
        0.42319008112684364_dp, 2.0e-15_dp)
    call check_close('pdtrc ordinary', pdtrc(2.2_dp, 3.0_dp), &
        0.5768099188731566_dp, 2.0e-15_dp)
    call check_close('pdtri ordinary', pdtri(2.2_dp, 0.7_dp), &
        1.9137757941270623_dp, 3.0e-14_dp)
    call check_close('pdtri round trip', pdtr(2.2_dp, pdtri(2.2_dp, 0.7_dp)), &
        0.7_dp, 4.0e-14_dp)
    call check_close('pdtri exact upper endpoint', pdtri(2.0_dp, 1.0_dp), 0.0_dp, 0.0_dp)
    call check_true('pdtri lower endpoint', .not. ieee_is_finite(pdtri(2.0_dp, 0.0_dp)))
    call check_true('pdtr both infinities', ieee_is_nan(pdtr(pinf, pinf)))

    call check_close('chdtr ordinary', chdtr(3.0_dp, 2.0_dp), &
        0.42759329552912023_dp, 2.0e-15_dp)
    call check_close('chdtrc ordinary', chdtrc(3.0_dp, 2.0_dp), &
        0.5724067044708798_dp, 2.0e-15_dp)
    call check_close('chdtri ordinary', chdtri(3.0_dp, 0.7_dp), &
        1.4236522430352798_dp, 3.0e-14_dp)
    call check_close('chdtri round trip', chdtrc(3.0_dp, chdtri(3.0_dp, 0.7_dp)), &
        0.7_dp, 4.0e-14_dp)
    call check_true('chdtr invalid df', ieee_is_nan(chdtr(0.0_dp, 2.0_dp)))

    call check_close('fdtr ordinary', fdtr(3.0_dp, 4.0_dp, 2.0_dp), &
        0.7436128024718243_dp, 4.0e-15_dp)
    call check_close('fdtrc ordinary', fdtrc(3.0_dp, 4.0_dp, 2.0_dp), &
        0.256387197528176_dp, 4.0e-15_dp)
    call check_close('fdtri ordinary', fdtri(3.0_dp, 4.0_dp, 0.7_dp), &
        1.7214950496159525_dp, 5.0e-14_dp)
    call check_close('fdtri round trip', fdtr(3.0_dp, 4.0_dp, fdtri(3.0_dp, 4.0_dp, 0.7_dp)), &
        0.7_dp, 5.0e-14_dp)
    call check_close('fdtr asymmetric tail', fdtr(0.2_dp, 1000.0_dp, 1.0e-8_dp), &
        0.13231839644830626_dp, 2.0e-13_dp)

    call check_close('gdtr ordinary', gdtr(2.0_dp, 3.0_dp, 4.0_dp), &
        0.986246032255997_dp, 3.0e-15_dp)
    call check_close('gdtrc ordinary', gdtrc(2.0_dp, 3.0_dp, 4.0_dp), &
        0.013753967744002971_dp, 3.0e-15_dp)
    call check_close('gdtria ordinary', gdtria(0.7_dp, 3.0_dp, 4.0_dp), &
        0.9038919164664978_dp, 4.0e-14_dp)
    call check_close('gdtrix ordinary', gdtrix(2.0_dp, 3.0_dp, 0.7_dp), &
        1.8077838329329956_dp, 4.0e-14_dp)
    call check_close('gdtr scaled product', gdtr(1.0e200_dp, 2.5_dp, 1.0e-200_dp), &
        0.15085496391539038_dp, 4.0e-15_dp)
    call check_close('gdtria round trip', gdtr(gdtria(0.7_dp, 3.0_dp, 4.0_dp), 3.0_dp, 4.0_dp), &
        0.7_dp, 5.0e-14_dp)
    call check_close('gdtrix round trip', gdtr(2.0_dp, 3.0_dp, gdtrix(2.0_dp, 3.0_dp, 0.7_dp)), &
        0.7_dp, 5.0e-14_dp)

    call check_close('stdtr ordinary', stdtr(3.0_dp, 2.0_dp), &
        0.9303370157205785_dp, 4.0e-15_dp)
    call check_close('stdtrit ordinary', stdtrit(3.0_dp, 0.7_dp), &
        0.5843897274398182_dp, 4.0e-14_dp)
    call check_relative('stdtr tail', stdtr(2.5_dp, -100.0_dp), &
        7.191649211666189e-6_dp, 3.0e-13_dp)
    call check_relative('stdtrit deep tail', stdtrit(2.5_dp, 1.0e-12_dp), &
        -55306.17407651582_dp, 3.0e-12_dp)
    call check_close('stdtrit round trip', stdtr(3.0_dp, stdtrit(3.0_dp, 0.7_dp)), &
        0.7_dp, 5.0e-14_dp)
    call check_true('stdtrit lower endpoint sign', stdtrit(3.0_dp, 0.0_dp) < 0.0_dp)
    call check_true('stdtrit lower endpoint infinite', &
        .not. ieee_is_finite(stdtrit(3.0_dp, 0.0_dp)))
    call check_true('stdtrit upper endpoint infinite', &
        .not. ieee_is_finite(stdtrit(3.0_dp, 1.0_dp)))

    print *, 'special distribution-helper tests passed'

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

    subroutine check_relative(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! nonzero reference value
        real(dp), intent(in) :: tolerance !! maximum relative error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_relative

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_distribution_helpers
