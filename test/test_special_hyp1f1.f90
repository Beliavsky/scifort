program test_special_hyp1f1
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, hyp1f1
    implicit none

    complex(dp) :: complex_values(2)
    real(dp) :: pinf
    real(dp) :: qnan
    real(dp) :: real_values(2)

    pinf = ieee_value(0.0_dp, ieee_positive_inf)
    qnan = ieee_value(0.0_dp, ieee_quiet_nan)

    call check_real_rel('hyp1f1 ordinary real', hyp1f1(2.5_dp, 3.5_dp, 4.0_dp), &
        2.32554446729673430e1_dp, 3.0e-14_dp)
    call check_real_rel('hyp1f1 large negative real', hyp1f1(2.3_dp, 4.7_dp, -100.0_dp), &
        3.0206934965050350e-4_dp, 2.0e-15_dp)
    call check_real_rel('hyp1f1 large positive real', hyp1f1(1.0_dp, 2.0_dp, 700.0_dp), &
        1.4489029353357207e301_dp, 6.0e-14_dp)
    call check_real_rel('hyp1f1 negative numerator real', hyp1f1(-0.5_dp, 2.5_dp, -100.0_dp), &
        6.713002792756121_dp, 3.0e-14_dp)
    call check_real_rel('hyp1f1 terminating polynomial', hyp1f1(-3.0_dp, 2.5_dp, 4.0_dp), &
        6.031746031746035e-2_dp, 3.0e-14_dp)
    call check_real_rel('hyp1f1 removable real denominator pole', &
        hyp1f1(-2.0_dp, -3.0_dp, 2.0_dp), 3.0_dp, 3.0e-14_dp)

    call check_complex_rel('hyp1f1 ordinary complex', &
        hyp1f1(2.5_dp, 3.5_dp, cmplx(3.0_dp, 4.0_dp, kind=dp)), &
        cmplx(-8.200797833110100_dp, -1.595858400403759_dp, kind=dp), 5.0e-14_dp)
    call check_complex_rel('hyp1f1 continued left half-plane', &
        hyp1f1(2.3_dp, 4.7_dp, cmplx(-50.0_dp, 30.0_dp, kind=dp)), &
        cmplx(3.595309665207361e-4_dp, 9.637923287187143e-4_dp, kind=dp), 2.0e-12_dp)
    call check_complex_rel('hyp1f1 large complex asymptotic', &
        hyp1f1(2.3_dp, 4.7_dp, cmplx(100.0_dp, 50.0_dp, kind=dp)), &
        cmplx(8.568858375957164e38_dp, -4.1160128231354073e39_dp, kind=dp), 3.0e-13_dp)
    call check_complex_rel('hyp1f1 large complex left half-plane', &
        hyp1f1(2.3_dp, 4.7_dp, cmplx(-100.0_dp, 50.0_dp, kind=dp)), &
        cmplx(1.1637638940827755e-4_dp, 2.0444663762627922e-4_dp, kind=dp), 3.0e-14_dp)
    call check_complex_rel('hyp1f1 negative numerator complex', &
        hyp1f1(-1.5_dp, 3.2_dp, cmplx(20.0_dp, 10.0_dp, kind=dp)), &
        cmplx(1.3498804231111174e2_dp, 3.6035423666644800e2_dp, kind=dp), 5.0e-14_dp)

    real_values = hyp1f1(2.5_dp, 3.5_dp, [4.0_dp, -2.0_dp])
    call check_real_rel('hyp1f1 elemental real first', real_values(1), &
        2.32554446729673430e1_dp, 3.0e-14_dp)
    call check_real_rel('hyp1f1 elemental real second', real_values(2), &
        2.6471407416488240e-1_dp, 5.0e-14_dp)

    complex_values = hyp1f1(2.5_dp, 3.5_dp, [cmplx(3.0_dp, 4.0_dp, kind=dp), &
        cmplx(-2.0_dp, 1.0_dp, kind=dp)])
    call check_complex_rel('hyp1f1 elemental complex first', complex_values(1), &
        cmplx(-8.200797833110100_dp, -1.595858400403759_dp, kind=dp), 5.0e-14_dp)
    call check_complex_rel('hyp1f1 elemental complex second', complex_values(2), &
        cmplx(2.1064043912150807e-1_dp, 1.4770440647027590e-1_dp, kind=dp), 5.0e-14_dp)

    call check_real_inf('hyp1f1 singular real denominator pole', hyp1f1(0.5_dp, -1.0_dp, 0.0_dp))
    call check_complex_inf('hyp1f1 singular complex denominator pole', &
        hyp1f1(-1.0_dp, -1.0_dp, cmplx(1.0_dp, 0.0_dp, kind=dp)))
    call check_real_rel('hyp1f1 a equals b', hyp1f1(2.0_dp, 2.0_dp, 4.0_dp), &
        exp(4.0_dp), 3.0e-15_dp)
    call check_real_rel('hyp1f1 zero numerator', hyp1f1(0.0_dp, 2.0_dp, -20.0_dp), &
        1.0_dp, 0.0_dp)
    call check_real_rel('hyp1f1 positive infinite denominator', hyp1f1(1.0_dp, pinf, 3.0_dp), &
        1.0_dp, 0.0_dp)
    call check_real_nan('hyp1f1 NaN numerator', hyp1f1(qnan, 2.0_dp, 1.0_dp))
    call check_real_nan('hyp1f1 NaN argument', hyp1f1(1.0_dp, 2.0_dp, qnan))

    print *, 'hyp1f1 special-function tests passed'

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

    subroutine check_complex_inf(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed value expected to have +Inf real part

        if (ieee_is_finite(real(actual, dp)) .or. real(actual, dp) <= 0.0_dp .or. &
                aimag(actual) /= 0.0_dp) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_complex_inf

    subroutine check_real_nan(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value expected to be NaN

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_real_nan

end program test_special_hyp1f1
