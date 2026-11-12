! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_error_fresnel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dawsn, dp, erfcx, erfi, fresnel
    implicit none

    real(dp) :: c
    real(dp) :: nan_value
    real(dp) :: pinf
    real(dp) :: s

    pinf = ieee_value(1.0_dp, ieee_positive_inf)
    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)

    call check_close('dawsn 0.1', dawsn(0.1_dp), &
        0.09933599239785290_dp, 3.0e-16_dp)
    call check_close('dawsn 1', dawsn(1.0_dp), &
        0.5380795069127684_dp, 5.0e-16_dp)
    call check_close('dawsn 6', dawsn(6.0_dp), &
        0.08454268897454388_dp, 3.0e-16_dp)
    call check_close('dawsn 10', dawsn(10.0_dp), &
        0.05025384718759854_dp, 2.0e-16_dp)
    call check_close('dawsn odd', dawsn(-2.0_dp), &
        -0.3013403889237920_dp, 5.0e-16_dp)
    call check_true('dawsn +inf', dawsn(pinf) == 0.0_dp)
    call check_true('dawsn nan', ieee_is_nan(dawsn(nan_value)))

    call check_close('erfi 0.5', erfi(0.5_dp), &
        0.6149520946965109_dp, 3.0e-15_dp)
    call check_close('erfi 2', erfi(2.0_dp), &
        18.564802414575553_dp, 8.0e-14_dp)
    call check_relative('erfi 8', erfi(8.0_dp), &
        4.432449746002334e26_dp, 2.0e-14_dp)
    call check_true('erfi overflow', .not. ieee_is_finite(erfi(30.0_dp)))
    call check_true('erfi +inf', .not. ieee_is_finite(erfi(pinf)))
    call check_true('erfi nan', ieee_is_nan(erfi(nan_value)))

    call check_close('erfcx 0.5', erfcx(0.5_dp), &
        0.6156903441929258_dp, 4.0e-16_dp)
    call check_close('erfcx 6', erfcx(6.0_dp), &
        0.09277656780053836_dp, 4.0e-16_dp)
    call check_close('erfcx 20', erfcx(20.0_dp), &
        0.028174348741051323_dp, 2.0e-16_dp)
    call check_relative('erfcx negative', erfcx(-8.0_dp), &
        1.2470298161623233e28_dp, 2.0e-14_dp)
    call check_true('erfcx +inf', erfcx(pinf) == 0.0_dp)
    call check_true('erfcx -inf', .not. ieee_is_finite(erfcx(-pinf)))
    call check_true('erfcx nan', ieee_is_nan(erfcx(nan_value)))

    call fresnel(0.5_dp, s, c)
    call check_close('fresnel S 0.5', s, 0.06473243285999929_dp, 5.0e-16_dp)
    call check_close('fresnel C 0.5', c, 0.4923442258714464_dp, 5.0e-16_dp)
    call fresnel(2.0_dp, s, c)
    call check_close('fresnel S 2', s, 0.34341567836369824_dp, 4.0e-15_dp)
    call check_close('fresnel C 2', c, 0.48825340607534073_dp, 4.0e-15_dp)
    call fresnel(6.0_dp, s, c)
    call check_close('fresnel S 6', s, 0.4469607612369303_dp, 5.0e-16_dp)
    call check_close('fresnel C 6', c, 0.4995314678555011_dp, 6.0e-16_dp)
    call fresnel(-4.0_dp, s, c)
    call check_close('fresnel S odd', s, -0.42051575424692844_dp, 4.0e-15_dp)
    call check_close('fresnel C odd', c, -0.4984260330381776_dp, 4.0e-15_dp)

    ! Exact binary64 x = 100000000.25. High-precision evaluation of the DLMF
    ! asymptotic form keeps the quadratic phase that SciPy 1.17.0 loses here.
    call fresnel(100000000.25_dp, s, c)
    call check_close('fresnel large-phase S', s, &
        0.49999999683222864_dp, 2.0e-16_dp)
    call check_close('fresnel large-phase C', c, &
        0.50000000031199825_dp, 2.0e-16_dp)

    ! At 2**52 + 1 the binary64 input is an exact odd integer, so x**2 mod 4
    ! is exactly one and the C correction has a known positive sign.
    call fresnel(4503599627370497.0_dp, s, c)
    call check_close('fresnel integer-phase S', s, 0.5_dp, 0.0_dp)
    call check_close('fresnel integer-phase C', c, &
        0.5000000000000001_dp, 0.0_dp)

    call fresnel(pinf, s, c)
    call check_close('fresnel +inf S', s, 0.5_dp, 0.0_dp)
    call check_close('fresnel +inf C', c, 0.5_dp, 0.0_dp)
    call fresnel(-pinf, s, c)
    call check_close('fresnel -inf S', s, -0.5_dp, 0.0_dp)
    call check_close('fresnel -inf C', c, -0.5_dp, 0.0_dp)
    call fresnel(nan_value, s, c)
    call check_true('fresnel nan S', ieee_is_nan(s))
    call check_true('fresnel nan C', ieee_is_nan(c))

    print *, 'special error/fresnel tests passed'

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

        if (ieee_is_nan(actual) .or. &
                abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' relative tolerance=', tolerance
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

end program test_special_error_fresnel
