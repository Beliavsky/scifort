! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_spherical_bessel_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, &
        ieee_value
    use scifort_special, only : dp, spherical_in, spherical_jn, spherical_kn, &
        spherical_yn
    implicit none

    integer, parameter :: orders(3) = [0, 2, 5]
    complex(dp), parameter :: arguments(3) = [ &
        cmplx(1.0_dp, 2.0_dp, kind=dp), &
        cmplx(-1.0_dp, 2.0_dp, kind=dp), &
        cmplx(-1.0_dp, -2.0_dp, kind=dp)]
    complex(dp), parameter :: j_refs(3) = [ &
        cmplx(1.416996119211877_dp, -0.8743911970021471_dp, kind=dp), &
        cmplx(-0.15013127056923803_dp, -0.38603772637600986_dp, kind=dp), &
        cmplx(-0.0037590213911662413_dp, 0.004726747071524348_dp, kind=dp)]
    complex(dp) :: values(3)
    real(dp) :: nan_value
    real(dp) :: pos_inf

    call check_complex_rel('spherical_jn upper right', &
        spherical_jn(2, cmplx(1.0_dp, 2.0_dp, kind=dp)), &
        cmplx(-0.15013127056923795_dp, 0.38603772637600986_dp, kind=dp), 8.0e-10_dp)
    call check_complex_rel('spherical_jn derivative upper left', &
        spherical_jn(5, cmplx(-1.0_dp, 2.0_dp, kind=dp), .true.), &
        cmplx(-0.006690738198206036_dp, 0.012475882681904919_dp, kind=dp), 8.0e-10_dp)
    call check_complex_abs('spherical_jn high-order small complex', &
        spherical_jn(20, cmplx(0.2_dp, 0.3_dp, kind=dp)), &
        cmplx(7.291539797292776e-35_dp, 7.581983243979208e-35_dp, kind=dp), 2.0e-43_dp)

    call check_complex_rel('spherical_yn upper right', &
        spherical_yn(2, cmplx(1.0_dp, 2.0_dp, kind=dp)), &
        cmplx(-0.22499497388743211_dp, -0.10171042469570665_dp, kind=dp), 9.0e-10_dp)
    call check_complex_rel('spherical_yn derivative lower left', &
        spherical_yn(5, cmplx(-1.0_dp, -2.0_dp, kind=dp), .true.), &
        cmplx(-4.541291214618294_dp, 17.577003001532034_dp, kind=dp), 9.0e-10_dp)
    call check_complex_rel('spherical_yn high order', &
        spherical_yn(20, cmplx(10.0_dp, 5.0_dp, kind=dp)), &
        cmplx(39.61730833641062_dp, 67.28294580051607_dp, kind=dp), 1.0e-9_dp)

    call check_complex_rel('spherical_in upper right', &
        spherical_in(2, cmplx(1.0_dp, 2.0_dp, kind=dp)), &
        cmplx(-0.21890731036371977_dp, 0.15881574297707554_dp, kind=dp), 8.0e-10_dp)
    call check_complex_rel('spherical_in derivative upper left', &
        spherical_in(5, cmplx(-1.0_dp, 2.0_dp, kind=dp), .true.), &
        cmplx(-0.0006270233083748914_dp, 0.0102381161060989_dp, kind=dp), 9.0e-10_dp)
    call check_complex_abs('spherical_in high-order small complex', &
        spherical_in(20, cmplx(0.2_dp, 0.3_dp, kind=dp)), &
        cmplx(7.261903423623431e-35_dp, 7.593467427786096e-35_dp, kind=dp), 2.0e-43_dp)
    call check_complex_rel('spherical_in derivative high order', &
        spherical_in(20, cmplx(10.0_dp, 5.0_dp, kind=dp), .true.), &
        cmplx(-0.00027670488239090027_dp, -0.00019028822793018592_dp, kind=dp), 2.0e-9_dp)

    call check_complex_rel('spherical_kn upper right', &
        spherical_kn(2, cmplx(1.0_dp, 2.0_dp, kind=dp)), &
        cmplx(-0.335212477123592_dp, 0.4228669286177394_dp, kind=dp), 9.0e-10_dp)
    call check_complex_rel('spherical_kn derivative lower left', &
        spherical_kn(5, cmplx(-1.0_dp, -2.0_dp, kind=dp), .true.), &
        cmplx(-1.5588657711281817_dp, -35.651284714124074_dp, kind=dp), 1.0e-9_dp)
    call check_complex_rel('spherical_kn near negative axis high order', &
        spherical_kn(20, cmplx(-5.0_dp, 0.1_dp, kind=dp)), &
        cmplx(-693294691.8354329_dp, -320101996.0297102_dp, kind=dp), 2.0e-9_dp)

    call check_complex_rel('spherical_jn order100 turning regime', &
        spherical_jn(100, cmplx(100.0_dp, 20.0_dp, kind=dp)), &
        cmplx(2.8446301802590535_dp, -1.2791328632664498_dp, kind=dp), 2.0e-9_dp)
    call check_complex_rel('spherical_yn order100 turning regime', &
        spherical_yn(100, cmplx(100.0_dp, 20.0_dp, kind=dp)), &
        cmplx(1.2791184077558073_dp, 2.8446100775767937_dp, kind=dp), 2.0e-9_dp)
    call check_complex_rel('spherical_in order100 turning regime', &
        spherical_in(100, cmplx(100.0_dp, 20.0_dp, kind=dp)), &
        cmplx(-1.4148917044277168e21_dp, 4.3590439517947114e20_dp, kind=dp), 2.0e-9_dp)
    call check_complex_rel('spherical_in derivative order100 turning regime', &
        spherical_in(100, cmplx(100.0_dp, 20.0_dp, kind=dp), .true.), &
        cmplx(-1.8890218324487352e21_dp, 7.895773709998732e20_dp, kind=dp), 2.0e-9_dp)

    call check_complex_rel('complex negative real spherical_jn', &
        spherical_jn(2, cmplx(-1.0_dp, 0.0_dp, kind=dp)), &
        cmplx(spherical_jn(2, -1.0_dp), 0.0_dp, kind=dp), 0.0_dp)
    call check_complex_rel('complex negative real spherical_kn', &
        spherical_kn(2, cmplx(-1.0_dp, 0.0_dp, kind=dp)), &
        cmplx(spherical_kn(2, -1.0_dp), 0.0_dp, kind=dp), 0.0_dp)

    values = spherical_jn(orders, arguments)
    call check_complex_array_rel('spherical_jn elemental complex arrays', values, j_refs, 1.0e-9_dp)

    call check_complex_rel('spherical_jn complex zero n0', &
        spherical_jn(0, cmplx(0.0_dp, 0.0_dp, kind=dp)), &
        cmplx(1.0_dp, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex_rel('spherical_jn complex zero n1 derivative', &
        spherical_jn(1, cmplx(0.0_dp, 0.0_dp, kind=dp), .true.), &
        cmplx(1.0_dp / 3.0_dp, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex_rel('spherical_in complex zero n0', &
        spherical_in(0, cmplx(0.0_dp, 0.0_dp, kind=dp)), &
        cmplx(1.0_dp, 0.0_dp, kind=dp), 0.0_dp)
    call check_true('spherical_yn complex zero is nan', &
        complex_is_nan(spherical_yn(0, cmplx(0.0_dp, 0.0_dp, kind=dp))))
    call check_true('spherical_kn complex zero is nan', &
        complex_is_nan(spherical_kn(0, cmplx(0.0_dp, 0.0_dp, kind=dp))))

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)
    call check_true('spherical complex invalid order', complex_is_nan( &
        spherical_jn(-1, cmplx(1.0_dp, 1.0_dp, kind=dp))))
    call check_true('spherical complex nan argument', complex_is_nan( &
        spherical_in(0, cmplx(nan_value, 1.0_dp, kind=dp))))
    call check_true('spherical complex infinite argument', complex_is_nan( &
        spherical_kn(0, cmplx(1.0_dp, pos_inf, kind=dp))))

    print *, 'complex spherical Bessel special-function tests passed'

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

    subroutine check_complex_abs(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! finite complex reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (complex_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' absolute tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_abs

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

end program test_special_spherical_bessel_complex
