! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_statistical_legacy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, ieee_value
    use scifort_special, only : dp, kolmogi, kolmogorov, smirnov, smirnovi, tklmbda
    implicit none

    real(dp) :: p
    real(dp) :: x

    call check_close('smirnov n5 d0.2', smirnov(5, 0.2_dp), &
        0.58528_dp, 2.0e-15_dp)
    call check_close('smirnov n24 d0.2', smirnov(24, 0.2_dp), &
        0.12826831885365378_dp, 3.0e-14_dp)
    call check_close('smirnov n100 d0.1', smirnov(100, 0.1_dp), &
        0.1265906584562817_dp, 3.0e-14_dp)
    call check_close('smirnov endpoint zero', smirnov(5, 0.0_dp), 1.0_dp, 0.0_dp)
    call check_close('smirnov endpoint one', smirnov(5, 1.0_dp), 0.0_dp, 0.0_dp)
    call check_true('smirnov invalid d', ieee_is_nan(smirnov(5, -0.1_dp)))
    call check_true('smirnov invalid n', ieee_is_nan(smirnov(0, 0.2_dp)))

    call check_close('smirnovi n5 p0.1', smirnovi(5, 0.1_dp), &
        0.44698006121219586_dp, 3.0e-15_dp)
    call check_close('smirnovi n24 p0.5', smirnovi(24, 0.5_dp), &
        0.11376168733165032_dp, 3.0e-15_dp)
    call check_close('smirnovi n100 p0.01', smirnovi(100, 0.01_dp), &
        0.14986798620592912_dp, 3.0e-15_dp)
    call check_close('smirnovi endpoint zero', smirnovi(5, 0.0_dp), 1.0_dp, 0.0_dp)
    call check_close('smirnovi endpoint one', smirnovi(5, 1.0_dp), 0.0_dp, 0.0_dp)
    call check_true('smirnovi invalid p', ieee_is_nan(smirnovi(5, -0.1_dp)))

    p = kolmogorov(1.0_dp)
    call check_close('kolmogorov 1', p, 0.26999967167735456_dp, 3.0e-15_dp)
    call check_close('kolmogorov 2', kolmogorov(2.0_dp), &
        0.0006709252557796953_dp, 3.0e-17_dp)
    call check_close('kolmogorov negative', kolmogorov(-0.1_dp), 1.0_dp, 0.0_dp)
    call check_close('kolmogi 0.5', kolmogi(0.5_dp), &
        0.8275735551899059_dp, 3.0e-15_dp)
    call check_close('kolmogi 0.1', kolmogi(0.1_dp), &
        1.2238478702170823_dp, 3.0e-15_dp)
    call check_true('kolmogi invalid', ieee_is_nan(kolmogi(1.1_dp)))
    x = kolmogi(p)
    call check_close('kolmogorov inverse round trip', x, 1.0_dp, 3.0e-15_dp)

    call check_close('tklmbda logistic', tklmbda(1.0_dp, 0.0_dp), &
        0.7310585786300049_dp, 3.0e-15_dp)
    call check_close('tklmbda negative lambda', tklmbda(-1.0_dp, -1.0_dp), &
        0.3819660112501069_dp, 3.0e-15_dp)
    call check_close('tklmbda positive lambda', tklmbda(0.5_dp, 0.25_dp), &
        0.6458122937336483_dp, 8.0e-15_dp)
    call check_close('tklmbda support low', tklmbda(-3.0_dp, 0.5_dp), 0.0_dp, 0.0_dp)
    call check_close('tklmbda support high', tklmbda(3.0_dp, 0.5_dp), 1.0_dp, 0.0_dp)
    call check_true('tklmbda invalid lambda', &
        ieee_is_nan(tklmbda(0.0_dp, ieee_value(0.0_dp, ieee_positive_inf))))

    print *, 'special statistical-legacy tests passed'

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

end program test_special_statistical_legacy
