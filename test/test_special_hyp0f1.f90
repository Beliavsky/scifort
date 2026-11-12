program test_special_hyp0f1
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, hyp0f1
    implicit none

    complex(dp) :: complex_values(2)
    real(dp) :: pinf
    real(dp) :: qnan
    real(dp) :: real_values(2)

    pinf = ieee_value(0.0_dp, ieee_positive_inf)
    qnan = ieee_value(0.0_dp, ieee_quiet_nan)

    call check_real_rel('hyp0f1 ordinary real', hyp0f1(2.5_dp, 3.0_dp), &
        2.84569534228708099_dp, 2.0e-14_dp)
    call check_real_rel('hyp0f1 oscillatory real', hyp0f1(3.0_dp, -100.0_dp), &
        -3.20682703845996415e-3_dp, 2.0e-14_dp)
    call check_real_rel('hyp0f1 large positive real', hyp0f1(3.0_dp, 10000.0_dp), &
        4.03868271583280647e81_dp, 3.0e-14_dp)
    call check_real_rel('hyp0f1 negative parameter real', hyp0f1(-0.5_dp, -100.0_dp), &
        1.86669870763659524e1_dp, 3.0e-14_dp)
    call check_real_rel('hyp0f1 negative parameter positive argument', &
        hyp0f1(-0.5_dp, 100.0_dp), -4.60906935639301014e9_dp, 4.0e-14_dp)

    call check_complex_rel('hyp0f1 ordinary complex', &
        hyp0f1(2.5_dp, cmplx(3.0_dp, 4.0_dp, kind=dp)), &
        cmplx(1.20368363795208455_dp, 3.07327074391132182_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('hyp0f1 continued complex', &
        hyp0f1(3.0_dp, cmplx(-100.0_dp, 20.0_dp, kind=dp)), &
        cmplx(-9.01249619191309831e-3_dp, -9.03137917082272833e-3_dp, kind=dp), &
        3.0e-14_dp)
    call check_complex_rel('hyp0f1 negative parameter complex', &
        hyp0f1(-0.5_dp, cmplx(-100.0_dp, 20.0_dp, kind=dp)), &
        cmplx(6.99827037066951618e1_dp, -2.64363462672023957e1_dp, kind=dp), 4.0e-13_dp)
    call check_complex_rel('hyp0f1 small oscillatory complex result', &
        hyp0f1(10.0_dp, cmplx(-250.0_dp, 15.0_dp, kind=dp)), &
        cmplx(9.55556503794538380e-7_dp, -4.07587536951161428e-7_dp, kind=dp), 2.0e-14_dp)

    real_values = hyp0f1(2.5_dp, [3.0_dp, -2.0_dp])
    call check_real_rel('hyp0f1 elemental real first', real_values(1), &
        2.84569534228708099_dp, 2.0e-14_dp)
    call check_real_rel('hyp0f1 elemental real second', real_values(2), &
        3.97606101444104800e-1_dp, 2.0e-14_dp)

    complex_values = hyp0f1(2.5_dp, [cmplx(3.0_dp, 4.0_dp, kind=dp), &
        cmplx(-5.0_dp, 2.0_dp, kind=dp)])
    call check_complex_rel('hyp0f1 elemental complex first', complex_values(1), &
        cmplx(1.20368363795208455_dp, 3.07327074391132182_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('hyp0f1 elemental complex second', complex_values(2), &
        cmplx(-6.03702814395699060e-2_dp, 1.19637138955561220e-1_dp, kind=dp), &
        3.0e-14_dp)

    call check_real_nan('hyp0f1 real pole', hyp0f1(-2.0_dp, 1.0_dp))
    call check_complex_nan_realpart('hyp0f1 complex pole', &
        hyp0f1(-2.0_dp, cmplx(1.0_dp, 2.0_dp, kind=dp)))
    call check_real_nan('hyp0f1 infinite real argument', hyp0f1(2.5_dp, pinf))
    call check_real_nan('hyp0f1 NaN real argument', hyp0f1(2.5_dp, qnan))
    call check_real_rel('hyp0f1 positive infinite parameter', hyp0f1(pinf, 3.0_dp), &
        1.0_dp, 0.0_dp)
    call check_real_rel('hyp0f1 NaN parameter at zero', hyp0f1(qnan, 0.0_dp), &
        1.0_dp, 0.0_dp)

    print *, 'hyp0f1 special-function tests passed'

contains

    subroutine check_real_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed real value
        real(dp), intent(in) :: expected !! finite real reference value
        real(dp), intent(in) :: tolerance !! maximum error scaled by max(1,abs(reference))

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_real_rel

    subroutine check_complex_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! finite complex reference value
        real(dp), intent(in) :: tolerance !! maximum error scaled by max(1,abs(reference))

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_rel

    subroutine check_real_nan(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value expected to be NaN

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_real_nan

    subroutine check_complex_nan_realpart(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed value expected to have a NaN real part

        if (.not. ieee_is_nan(real(actual, dp))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_complex_nan_realpart

end program test_special_hyp0f1
