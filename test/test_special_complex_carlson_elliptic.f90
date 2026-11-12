program test_special_complex_carlson_elliptic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, elliprc, elliprd, elliprf, elliprg, elliprj
    implicit none

    complex(dp), parameter :: scale = cmplx(0.3_dp, 0.4_dp, kind=dp)
    complex(dp) :: base
    complex(dp) :: value
    complex(dp) :: x_values(2)
    complex(dp) :: y_values(2)
    complex(dp) :: z_values(2)
    complex(dp) :: rf_values(2)
    real(dp) :: pinf
    real(dp) :: qnan

    pinf = ieee_value(0.0_dp, ieee_positive_inf)
    qnan = ieee_value(0.0_dp, ieee_quiet_nan)

    call check_complex_rel('elliprc ordinary complex', &
        elliprc(cmplx(1.2_dp, 3.4_dp, kind=dp), cmplx(5.0_dp, 0.0_dp, kind=dp)), &
        cmplx(4.78722846608336916e-1_dp, -9.02710542226324991e-2_dp, kind=dp), 3.0e-13_dp)
    call check_complex_rel('elliprc zero first argument', &
        elliprc(cmplx(0.0_dp, 0.0_dp, kind=dp), cmplx(1.2_dp, 3.4_dp, kind=dp)), &
        cmplx(6.75312534611681392e-1_dp, -4.77793802638808607e-1_dp, kind=dp), 3.0e-13_dp)
    call check_complex_rel('elliprc negative real second principal value', &
        elliprc(cmplx(1.0_dp, 2.0_dp, kind=dp), cmplx(-2.0_dp, 0.0_dp, kind=dp)), &
        cmplx(4.99394805608676329e-1_dp, 7.34481761309273562e-2_dp, kind=dp), 3.0e-13_dp)
    call check_complex_rel('elliprc near cut', &
        elliprc(cmplx(-1.0_dp, 1.0e-8_dp, kind=dp), cmplx(2.0_dp, 0.0_dp, kind=dp)), &
        cmplx(9.06899681084063825e-1_dp, -3.80172996638973637e-1_dp, kind=dp), 4.0e-13_dp)

    call check_complex_rel('elliprf ordinary complex', &
        elliprf(cmplx(1.2_dp, 3.4_dp, kind=dp), cmplx(5.0_dp, 0.0_dp, kind=dp), &
            cmplx(6.0_dp, 0.0_dp, kind=dp)), &
        cmplx(4.63739354934839865e-1_dp, -8.50398414569596556e-2_dp, kind=dp), 3.0e-13_dp)
    call check_complex_rel('elliprf conjugate pair', &
        elliprf(cmplx(-1.0_dp, 2.0_dp, kind=dp), cmplx(-1.0_dp, -2.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp)), &
        cmplx(8.67085315389173328e-1_dp, 0.0_dp, kind=dp), 3.0e-13_dp)
    call check_complex_rel('elliprf one zero', &
        elliprf(cmplx(0.0_dp, 0.0_dp, kind=dp), cmplx(1.0_dp, 2.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp)), &
        cmplx(9.58593760618404067e-1_dp, -2.61653747787980229e-1_dp, kind=dp), 3.0e-13_dp)

    call check_complex_rel('elliprd ordinary complex', &
        elliprd(cmplx(1.2_dp, 3.4_dp, kind=dp), cmplx(5.0_dp, 0.0_dp, kind=dp), &
            cmplx(6.0_dp, 0.0_dp, kind=dp)), &
        cmplx(8.28846236219513560e-2_dp, -2.83768097451232992e-2_dp, kind=dp), 4.0e-13_dp)
    call check_complex_rel('elliprd conjugate pair', &
        elliprd(cmplx(-1.0_dp, 2.0_dp, kind=dp), cmplx(-1.0_dp, -2.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp)), &
        cmplx(3.44314242713237517e-1_dp, 0.0_dp, kind=dp), 4.0e-13_dp)

    call check_complex_rel('elliprg ordinary complex', &
        elliprg(cmplx(1.2_dp, 3.4_dp, kind=dp), cmplx(5.0_dp, 0.0_dp, kind=dp), &
            cmplx(6.0_dp, 0.0_dp, kind=dp)), &
        cmplx(2.04850611118344750_dp, 3.15127800675830849e-1_dp, kind=dp), 5.0e-13_dp)
    call check_complex_rel('elliprg two zeros', &
        elliprg(cmplx(0.0_dp, 0.0_dp, kind=dp), cmplx(0.0_dp, 0.0_dp, kind=dp), &
            cmplx(1.0_dp, 2.0_dp, kind=dp)), &
        cmplx(6.36009824757034492e-1_dp, 3.93075688878711638e-1_dp, kind=dp), 5.0e-13_dp)

    call check_complex_rel('elliprj ordinary complex', &
        elliprj(cmplx(1.2_dp, 3.4_dp, kind=dp), cmplx(5.0_dp, 0.0_dp, kind=dp), &
            cmplx(6.0_dp, 0.0_dp, kind=dp), cmplx(7.0_dp, 0.0_dp, kind=dp)), &
        cmplx(7.53615890063274729e-2_dp, -2.52330972492828057e-2_dp, kind=dp), 3.0e-12_dp)
    call check_complex_rel('elliprj unbalanced complex', &
        elliprj(cmplx(0.013156576809215587_dp, 2.542745457113952_dp, kind=dp), &
            cmplx(6.913260494980437_dp, 0.033462634677567764_dp, kind=dp), &
            cmplx(0.005273402991439643_dp, 0.035473201706190303_dp, kind=dp), &
            cmplx(0.001075933916054898_dp, 0.04410193429902219_dp, kind=dp)), &
        cmplx(7.56167663812429801e-3_dp, -3.07090337225544507_dp, kind=dp), 5.0e-12_dp)
    call check_complex_rel('elliprj conjugate-domain complex', &
        elliprj(cmplx(-1.0_dp, 2.0_dp, kind=dp), cmplx(-1.0_dp, -2.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp), cmplx(4.0_dp, 1.0_dp, kind=dp)), &
        cmplx(2.79308381070909273e-1_dp, -4.48136963785417189e-2_dp, kind=dp), 5.0e-12_dp)
    call check_complex_rel('elliprj real xyz complex pole', &
        elliprj(cmplx(1.0_dp, 0.0_dp, kind=dp), cmplx(2.0_dp, 0.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp), cmplx(-0.5_dp, 0.25_dp, kind=dp)), &
        cmplx(3.3658765910387920e-1_dp, -1.0321040234370540_dp, kind=dp), 5.0e-12_dp)
    call check_complex_rel('elliprj real-axis principal value dispatch', &
        elliprj(cmplx(1.0_dp, 0.0_dp, kind=dp), cmplx(2.0_dp, 0.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp), cmplx(-0.5_dp, 0.0_dp, kind=dp)), &
        cmplx(2.0722001115871863e-1_dp, 0.0_dp, kind=dp), 5.0e-12_dp)

    base = elliprf(cmplx(1.2_dp, 3.4_dp, kind=dp), cmplx(5.0_dp, 0.0_dp, kind=dp), &
        cmplx(6.0_dp, 0.0_dp, kind=dp))
    call check_complex_rel('elliprf complex homogeneity', &
        elliprf(scale * cmplx(1.2_dp, 3.4_dp, kind=dp), scale * 5.0_dp, scale * 6.0_dp), &
        base / sqrt(scale), 8.0e-13_dp)

    x_values = [cmplx(1.2_dp, 3.4_dp, kind=dp), cmplx(0.0_dp, 0.0_dp, kind=dp)]
    y_values = [cmplx(5.0_dp, 0.0_dp, kind=dp), cmplx(1.0_dp, 2.0_dp, kind=dp)]
    z_values = [cmplx(6.0_dp, 0.0_dp, kind=dp), cmplx(3.0_dp, 0.0_dp, kind=dp)]
    rf_values = elliprf(x_values, y_values, z_values)
    call check_complex_rel('elliprf elemental array first', rf_values(1), &
        cmplx(4.63739354934839865e-1_dp, -8.50398414569596556e-2_dp, kind=dp), 3.0e-13_dp)
    call check_complex_rel('elliprf elemental array second', rf_values(2), &
        cmplx(9.58593760618404067e-1_dp, -2.61653747787980229e-1_dp, kind=dp), 3.0e-13_dp)

    value = elliprf(cmplx(-1.0_dp, 0.0_dp, kind=dp), cmplx(1.0_dp, 0.0_dp, kind=dp), &
        cmplx(2.0_dp, 0.0_dp, kind=dp))
    call check_complex_nan('elliprf negative-real cut', value)
    value = elliprj(cmplx(1.0_dp, 1.0_dp, kind=dp), cmplx(2.0_dp, 0.5_dp, kind=dp), &
        cmplx(3.0_dp, -0.2_dp, kind=dp), cmplx(-0.1_dp, 1.0_dp, kind=dp))
    call check_complex_nan('elliprj rejected insufficient domain', value)

    call check_complex_rel('elliprf positive-infinity limit', &
        elliprf(cmplx(pinf, 0.0_dp, kind=dp), cmplx(2.0_dp, 0.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp)), cmplx(0.0_dp, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex_positive_infinity('elliprg imaginary-infinity limit', &
        elliprg(cmplx(0.0_dp, pinf, kind=dp), cmplx(2.0_dp, 0.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp)))
    call check_complex_rel('elliprj infinite-pole limit', &
        elliprj(cmplx(1.0_dp, 1.0_dp, kind=dp), cmplx(2.0_dp, 0.5_dp, kind=dp), &
            cmplx(3.0_dp, -0.2_dp, kind=dp), cmplx(pinf, pinf, kind=dp)), &
        cmplx(0.0_dp, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex_nan('elliprf ambiguous infinity', &
        elliprf(cmplx(pinf, pinf, kind=dp), cmplx(2.0_dp, 0.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp)))
    call check_complex_nan('elliprf NaN input', &
        elliprf(cmplx(qnan, 0.0_dp, kind=dp), cmplx(2.0_dp, 0.0_dp, kind=dp), &
            cmplx(3.0_dp, 0.0_dp, kind=dp)))

    print *, 'complex Carlson elliptic special-function tests passed'

contains

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

    subroutine check_complex_positive_infinity(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! value expected to be positive infinity plus zero imaginary part

        if (ieee_is_finite(real(actual, dp)) .or. real(actual, dp) <= 0.0_dp .or. &
                aimag(actual) /= 0.0_dp) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_complex_positive_infinity

    subroutine check_complex_nan(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed value expected to contain NaNs

        if (.not. ieee_is_nan(real(actual, dp)) .or. .not. ieee_is_nan(aimag(actual))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_complex_nan

end program test_special_complex_carlson_elliptic
