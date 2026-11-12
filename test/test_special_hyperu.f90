program test_special_hyperu
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, hyperu
    implicit none

    real(dp) :: pinf
    real(dp) :: qnan
    real(dp) :: transformed
    real(dp) :: values(2)

    pinf = ieee_value(0.0_dp, ieee_positive_inf)
    qnan = ieee_value(0.0_dp, ieee_quiet_nan)

    call check_rel('hyperu ordinary', hyperu(1.0_dp, 1.0_dp, 1.0_dp), &
        0.59634736232319407434_dp, 3.0e-12_dp)
    call check_rel('hyperu second ordinary', hyperu(2.0_dp, 1.0_dp, 1.0_dp), &
        0.19269472464638814868_dp, 3.0e-12_dp)
    call check_rel('hyperu noninteger parameters', hyperu(2.3_dp, 1.7_dp, 4.2_dp), &
        0.020134246423386408244_dp, 3.0e-12_dp)
    call check_rel('hyperu negative a', hyperu(-0.5_dp, 2.0_dp, 1.0_dp), &
        0.17023014757496953553_dp, 3.0e-12_dp)
    call check_rel('hyperu small x noninteger b', hyperu(1.2_dp, 2.5_dp, 1.0e-6_dp), &
        965211965.69109443075_dp, 5.0e-10_dp)
    call check_rel('hyperu small x integer b', hyperu(1.2_dp, 2.0_dp, 1.0e-6_dp), &
        1089121.3823733429142_dp, 2.0e-8_dp)
    call check_rel('hyperu near integer b', hyperu(1.2_dp, 2.000001_dp, 0.01_dp), &
        107.87334873774195118_dp, 2.0e-8_dp)
    call check_rel('hyperu small x negative b', hyperu(5.0_dp, 0.0_dp, 0.01_dp), &
        0.0071380662380992986486_dp, 2.0e-8_dp)
    call check_rel('hyperu large x', hyperu(3.5_dp, -2.0_dp, 100.0_dp), &
        8.0570645987601123776e-8_dp, 2.0e-12_dp)
    call check_rel('hyperu difficult continuation', hyperu(20.0_dp, 30.0_dp, 2.0_dp), &
        8739.3477568626403809_dp, 3.0e-10_dp)
    call check_rel('hyperu terminating polynomial', hyperu(-3.0_dp, 2.5_dp, 0.7_dp), &
        -12.572000000000001324_dp, 3.0e-14_dp)
    call check_rel('hyperu negative noninteger small x', hyperu(-2.3_dp, 2.5_dp, 0.001_dp), &
        -19505.080643650935504_dp, 3.0e-10_dp)

    transformed = 2.3_dp**(1.0_dp - 1.7_dp) * hyperu(2.3_dp - 1.7_dp + 1.0_dp, &
        2.0_dp - 1.7_dp, 2.3_dp)
    call check_rel('hyperu Kummer transformation', hyperu(2.3_dp, 1.7_dp, 2.3_dp), &
        transformed, 5.0e-12_dp)

    values = hyperu(1.0_dp, 1.0_dp, [1.0_dp, 10.0_dp])
    call check_rel('hyperu elemental first', values(1), 0.59634736232319407434_dp, 3.0e-12_dp)
    call check_rel('hyperu elemental second', values(2), 0.091563333939788081876_dp, 3.0e-12_dp)

    call check_rel('hyperu zero finite limit', hyperu(1.0_dp, 0.5_dp, 0.0_dp), 2.0_dp, 3.0e-14_dp)
    call check_rel('hyperu zero from reciprocal gamma pole', hyperu(-0.5_dp, 0.5_dp, 0.0_dp), &
        0.0_dp, 0.0_dp)
    call check_rel('hyperu zero polynomial', hyperu(-2.0_dp, 0.5_dp, 0.0_dp), 0.75_dp, 3.0e-14_dp)
    call check_negative_inf('hyperu signed zero divergence', hyperu(-0.5_dp, 2.0_dp, 0.0_dp))
    call check_rel('hyperu a zero identity', hyperu(0.0_dp, 3.0_dp, 0.0_dp), 1.0_dp, 0.0_dp)

    call check_rel('hyperu positive infinity decay', hyperu(2.0_dp, 1.0_dp, pinf), 0.0_dp, 0.0_dp)
    call check_positive_inf('hyperu negative a at infinity', hyperu(-0.5_dp, 1.0_dp, pinf))
    call check_nan('hyperu negative branch cut', hyperu(1.0_dp, 1.0_dp, -1.0_dp))
    call check_nan('hyperu NaN parameter', hyperu(qnan, 1.0_dp, 1.0_dp))
    call check_nan('hyperu NaN argument', hyperu(1.0_dp, 1.0_dp, qnan))
    call check_nan('hyperu infinite parameter', hyperu(pinf, 1.0_dp, 1.0_dp))

    print *, 'hyperu special-function tests passed'

contains

    subroutine check_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed finite value
        real(dp), intent(in) :: expected !! finite reference value
        real(dp), intent(in) :: tolerance !! scaled absolute tolerance

        if (.not. ieee_is_finite(actual) .or. &
                abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_rel

    subroutine check_positive_inf(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! value expected to be positive infinity

        if (ieee_is_finite(actual) .or. actual <= 0.0_dp) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_positive_inf

    subroutine check_negative_inf(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! value expected to be negative infinity

        if (ieee_is_finite(actual) .or. actual >= 0.0_dp) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_negative_inf

    subroutine check_nan(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! value expected to be NaN

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_nan

end program test_special_hyperu
