! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_hankel_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, h1vp, h2vp, hankel1, hankel1e, hankel2, hankel2e
    implicit none

    real(dp), parameter :: orders(3) = [0.0_dp, -0.5_dp, 2.3_dp]
    complex(dp), parameter :: arguments(3) = [ &
        cmplx(2.0_dp, 1.0_dp, kind=dp), &
        cmplx(0.5_dp, -3.0_dp, kind=dp), &
        cmplx(-2.0_dp, 1.0_dp, kind=dp)]
    complex(dp), parameter :: h1_refs(3) = [ &
        cmplx(0.11221517779606792_dp, 0.15428168525601327_dp, kind=dp), &
        cmplx(3.3056583723189688_dp, 8.5742668521347320_dp, kind=dp), &
        cmplx(-0.22637757414340362_dp, -0.28280992973049685_dp, kind=dp)]
    complex(dp) :: values(3)
    complex(dp) :: z
    real(dp) :: nan_value
    real(dp) :: pos_inf

    call check_complex_rel('hankel1 upper right', &
        hankel1(0.0_dp, cmplx(2.0_dp, 1.0_dp, kind=dp)), &
        cmplx(0.11221517779606792_dp, 0.15428168525601327_dp, kind=dp), 2.0e-10_dp)
    call check_complex_rel('hankel2 lower right', &
        hankel2(2.3_dp, cmplx(2.0_dp, -1.0_dp, kind=dp)), &
        cmplx(-0.09573663979871791_dp, 0.34937481052484493_dp, kind=dp), 2.0e-10_dp)
    call check_complex_rel('hankel1e upper', &
        hankel1e(2.3_dp, cmplx(0.5_dp, 3.0_dp, kind=dp)), &
        cmplx(0.26323821838171790_dp, 0.88362438757353401_dp, kind=dp), 2.0e-10_dp)
    call check_complex_rel('hankel2e lower negative order', &
        hankel2e(-0.5_dp, cmplx(0.5_dp, -3.0_dp, kind=dp)), &
        cmplx(0.34909251734965030_dp, 0.29572572779871398_dp, kind=dp), 2.0e-10_dp)

    call check_complex_rel('hankel1 upper left', &
        hankel1(2.3_dp, cmplx(-2.0_dp, 1.0_dp, kind=dp)), &
        cmplx(-0.22637757414340362_dp, -0.28280992973049685_dp, kind=dp), 2.0e-10_dp)
    call check_complex_rel('hankel1 lower left', &
        hankel1(2.3_dp, cmplx(-2.0_dp, -1.0_dp, kind=dp)), &
        cmplx(0.98402361737084498_dp, -0.43289927814052065_dp, kind=dp), 2.0e-10_dp)
    call check_complex_rel('hankel2 upper left', &
        hankel2(2.3_dp, cmplx(-2.0_dp, 1.0_dp, kind=dp)), &
        cmplx(0.98402361737084498_dp, 0.43289927814052065_dp, kind=dp), 2.0e-10_dp)
    call check_complex_rel('hankel2 lower left', &
        hankel2(2.3_dp, cmplx(-2.0_dp, -1.0_dp, kind=dp)), &
        cmplx(-0.22637757414340362_dp, 0.28280992973049685_dp, kind=dp), 2.0e-10_dp)

    call check_complex_rel('h1vp complex n1', &
        h1vp(0.5_dp, cmplx(1.2_dp, 0.7_dp, kind=dp), 1), &
        cmplx(0.15136502000683200_dp, 0.38172577383574835_dp, kind=dp), 3.0e-10_dp)
    call check_complex_rel('h1vp complex n2', &
        h1vp(0.5_dp, cmplx(1.2_dp, 0.7_dp, kind=dp), 2), &
        cmplx(-0.50844870224235894_dp, -0.02631616256298830_dp, kind=dp), 3.0e-10_dp)
    call check_complex_rel('h2vp upper left n2', &
        h2vp(2.3_dp, cmplx(-1.5_dp, 0.8_dp, kind=dp), 2), &
        cmplx(-1.4919619869238390_dp, 0.92705071750606505_dp, kind=dp), 4.0e-10_dp)
    call check_complex_rel('h2vp lower left n3', &
        h2vp(2.3_dp, cmplx(-1.5_dp, -0.8_dp, kind=dp), 3), &
        cmplx(3.5278798647613883_dp, -0.28783758986931895_dp, kind=dp), 5.0e-10_dp)

    call check_complex_rel('hankel1e larger complex', &
        hankel1e(10.2_dp, cmplx(20.0_dp, 10.0_dp, kind=dp)), &
        cmplx(-0.38878345237351725_dp, -0.32956417483372769_dp, kind=dp), 6.0e-10_dp)
    call check_complex_rel('hankel2e larger complex', &
        hankel2e(10.2_dp, cmplx(20.0_dp, -10.0_dp, kind=dp)), &
        cmplx(-0.38878345237351725_dp, 0.32956417483372769_dp, kind=dp), 6.0e-10_dp)

    z = cmplx(1.7_dp, 0.6_dp, kind=dp)
    call check_complex_rel('Schwarz reflection H1/H2', hankel1(1.3_dp, conjg(z)), &
        conjg(hankel2(1.3_dp, z)), 2.0e-10_dp)
    call check_complex_rel('h1vp n0', h1vp(0.7_dp, z, 0), hankel1(0.7_dp, z), 0.0_dp)
    call check_complex_rel('h2vp n0', h2vp(0.7_dp, z, 0), hankel2(0.7_dp, z), 0.0_dp)
    call check_complex_rel('complex negative axis H1 matches real', &
        hankel1(2.3_dp, cmplx(-2.0_dp, 0.0_dp, kind=dp)), hankel1(2.3_dp, -2.0_dp), 0.0_dp)
    call check_complex_rel('complex negative axis H2 matches real', &
        hankel2(2.3_dp, cmplx(-2.0_dp, 0.0_dp, kind=dp)), hankel2(2.3_dp, -2.0_dp), 0.0_dp)

    values = hankel1(orders, arguments)
    call check_complex_array_rel('hankel1 complex elemental arrays', values, h1_refs, 3.0e-10_dp)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)
    call check_true('complex nan argument', complex_is_nan( &
        hankel1(0.0_dp, cmplx(nan_value, 1.0_dp, kind=dp))))
    call check_true('complex infinite argument', complex_is_nan( &
        hankel2(0.0_dp, cmplx(1.0_dp, pos_inf, kind=dp))))
    call check_true('complex negative derivative order', complex_is_nan( &
        h1vp(0.0_dp, cmplx(1.0_dp, 1.0_dp, kind=dp), -1)))
    call check_true('complex zero hankel1', complex_is_nan( &
        hankel1(1.0_dp, cmplx(0.0_dp, 0.0_dp, kind=dp))))

    print *, 'Hankel complex-argument special-function tests passed'

contains

    subroutine check_complex_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! finite complex reference value
        real(dp), intent(in) :: tolerance !! maximum scaled absolute error

        if (complex_is_nan(actual) .or. &
                abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' scaled tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_rel

    subroutine check_complex_array_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual(:) !! computed complex values
        complex(dp), intent(in) :: expected(:) !! finite complex reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise scaled absolute error

        integer :: i

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', name, ': size mismatch'
            error stop 1
        end if
        do i = 1, size(actual)
            call check_complex_rel(name, actual(i), expected(i), tolerance)
        end do
    end subroutine check_complex_array_rel

    pure logical function complex_is_nan(z) result(is_nan)
        use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
        complex(dp), intent(in) :: z !! complex value to classify

        is_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
    end function complex_is_nan

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_hankel_complex
