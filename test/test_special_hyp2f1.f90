program test_special_hyp2f1
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, hyp2f1
    implicit none

    complex(dp) :: complex_values(2)
    real(dp) :: pinf
    real(dp) :: qnan
    real(dp) :: real_values(2)

    pinf = ieee_value(0.0_dp, ieee_positive_inf)
    qnan = ieee_value(0.0_dp, ieee_quiet_nan)

    call check_real_rel('hyp2f1 ordinary real', hyp2f1(2.3_dp, 1.1_dp, 3.7_dp, 0.5_dp), &
        1.5515729233917759_dp, 3.0e-14_dp)
    call check_real_rel('hyp2f1 near one', hyp2f1(0.25_dp, 1.75_dp, 2.5_dp, 0.99_dp), &
        1.5690469526417128_dp, 3.0e-13_dp)
    call check_real_rel('hyp2f1 negative continuation', hyp2f1(3.2_dp, -0.4_dp, 1.3_dp, -10.0_dp), &
        3.8457177923179713_dp, 8.0e-14_dp)
    call check_real_rel('hyp2f1 finite beyond cut', hyp2f1(1.0_dp, 2.0_dp, 1.0_dp, 1.2_dp), &
        25.0_dp, 2.0e-13_dp)
    call check_real_rel('hyp2f1 Gauss value at one', hyp2f1(1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp), &
        2.0_dp, 3.0e-14_dp)
    call check_real_rel('hyp2f1 terminating polynomial', hyp2f1(-3.0_dp, 2.5_dp, 4.0_dp, 2.0_dp), &
        -0.125_dp, 3.0e-14_dp)
    call check_real_rel('hyp2f1 terminating before denominator pole', &
        hyp2f1(-1.0_dp, 2.0_dp, -2.0_dp, 0.5_dp), 1.5_dp, 3.0e-14_dp)
    call check_real_rel('hyp2f1 cancelled denominator pole', &
        hyp2f1(-2.0_dp, 2.0_dp, -2.0_dp, 0.5_dp), 4.0_dp, 3.0e-14_dp)

    call check_complex_rel('hyp2f1 ordinary complex', &
        hyp2f1(2.2_dp, 3.1_dp, 4.5_dp, cmplx(0.9_dp, 0.8_dp, kind=dp)), &
        cmplx(-0.9012861841487853_dp, 1.4460875039582420_dp, kind=dp), 8.0e-14_dp)
    call check_complex_rel('hyp2f1 upper branch cut', &
        hyp2f1(1.0_dp, 1.0_dp, 2.0_dp, cmplx(1.2_dp, 0.0_dp, kind=dp)), &
        cmplx(1.3411982603617505_dp, 2.6179938779914945_dp, kind=dp), 8.0e-14_dp)
    call check_complex_rel('hyp2f1 lower branch cut', &
        hyp2f1(1.0_dp, 1.0_dp, 2.0_dp, cmplx(1.2_dp, -0.0_dp, kind=dp)), &
        cmplx(1.3411982603617505_dp, -2.6179938779914945_dp, kind=dp), 8.0e-14_dp)
    call check_complex_rel('hyp2f1 far complex', &
        hyp2f1(0.5_dp, 1.0_dp, 1.5_dp, cmplx(4.0_dp, 2.0_dp, kind=dp)), &
        cmplx(0.3809013397134929_dp, 0.6069882990957904_dp, kind=dp), 1.0e-13_dp)
    call check_complex_rel('hyp2f1 difficult complex continuation', &
        hyp2f1(4.012622322850383_dp, 4.329347914335271_dp, -3.2550282172290967_dp, &
            cmplx(-1.3824143923222665_dp, 0.5547167482664364_dp, kind=dp)), &
        cmplx(0.1228145191442493_dp, -2.9898622217910607_dp, kind=dp), 3.0e-12_dp)

    real_values = hyp2f1(1.0_dp, 1.0_dp, 2.0_dp, [0.5_dp, -2.0_dp])
    call check_real_rel('hyp2f1 elemental real first', real_values(1), &
        1.3862943611198901_dp, 5.0e-14_dp)
    call check_real_rel('hyp2f1 elemental real second', real_values(2), &
        0.5493061443340546_dp, 5.0e-14_dp)

    complex_values = hyp2f1(1.2_dp, 2.3_dp, 3.4_dp, &
        [cmplx(-3.0_dp, 4.0_dp, kind=dp), cmplx(0.2_dp, -0.3_dp, kind=dp)])
    call check_complex_rel('hyp2f1 elemental complex first', complex_values(1), &
        cmplx(0.14657416768786263_dp, 0.15026255508862367_dp, kind=dp), 8.0e-14_dp)
    call check_complex_rel('hyp2f1 elemental complex second', complex_values(2), &
        cmplx(1.0983480472193027_dp, -0.3203993819367633_dp, kind=dp), 8.0e-14_dp)

    call check_real_inf('hyp2f1 singular denominator pole', hyp2f1(-3.0_dp, 2.0_dp, -2.0_dp, 0.5_dp))
    call check_real_inf('hyp2f1 real branch cut', hyp2f1(1.0_dp, 1.0_dp, 2.0_dp, 1.2_dp))
    call check_real_inf('hyp2f1 divergent at one', hyp2f1(1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp))
    call check_real_negative_inf('hyp2f1 negative divergence at one', &
        hyp2f1(-0.5_dp, 3.1_dp, 0.4_dp, 1.0_dp))
    call check_real_rel('hyp2f1 zero numerator with pole', hyp2f1(0.0_dp, 1.0_dp, -2.0_dp, 0.5_dp), &
        1.0_dp, 0.0_dp)
    call check_real_rel('hyp2f1 positive infinite denominator', hyp2f1(1.0_dp, 1.0_dp, pinf, 0.5_dp), &
        1.0_dp, 0.0_dp)
    call check_real_nan('hyp2f1 NaN numerator', hyp2f1(qnan, 1.0_dp, 2.0_dp, 0.5_dp))
    call check_real_nan('hyp2f1 NaN argument', hyp2f1(1.0_dp, 1.0_dp, 2.0_dp, qnan))
    call check_real_nan('hyp2f1 infinite argument limitation', hyp2f1(1.0_dp, 1.0_dp, 2.0_dp, pinf))

    print *, 'hyp2f1 special-function tests passed'

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

    subroutine check_real_inf(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value expected to be positive infinity

        if (ieee_is_finite(actual) .or. actual <= 0.0_dp) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_real_inf

    subroutine check_real_negative_inf(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value expected to be negative infinity

        if (ieee_is_finite(actual) .or. actual >= 0.0_dp) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_real_negative_inf

    subroutine check_real_nan(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value expected to be NaN

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_real_nan

end program test_special_hyp2f1
