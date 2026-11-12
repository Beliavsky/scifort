! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_complex_elementary
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, expm1, log1p, log_ndtr, ndtr, xlog1py, xlogy
    implicit none

    complex(dp) :: actual(3)
    complex(dp) :: input(3)
    complex(dp) :: nan_result
    real(dp) :: nan_value

    call check_complex('ndtr ordinary', ndtr(cmplx(0.25_dp, 0.5_dp, kind=dp)), &
        cmplx(0.6115602918247813_dp, 0.20115839591423415_dp, kind=dp), 2.0e-13_dp)
    call check_complex('ndtr large complex', ndtr(cmplx(-10.0_dp, -30.0_dp, kind=dp)), &
        cmplx(6.2052682526467586e171_dp, 2.22615980800883e171_dp, kind=dp), 3.0e-13_dp)
    call check_complex('log_ndtr far left', log_ndtr(cmplx(-20.0_dp, 20.0_dp, kind=dp)), &
        cmplx(-4.261248303073217_dp, -1.3397114720099623_dp, kind=dp), 3.0e-13_dp)
    call check_complex('log_ndtr right half-plane', log_ndtr(cmplx(30.0_dp, -30.0_dp, kind=dp)), &
        cmplx(0.00619568293739492_dp, -0.007034681681139299_dp, kind=dp), 3.0e-13_dp)

    call check_complex('expm1 ordinary', expm1(cmplx(0.25_dp, 0.5_dp, kind=dp)), &
        cmplx(0.1268383147091815_dp, 0.6155945769770066_dp, kind=dp), 3.0e-15_dp)
    call check_complex('expm1 tiny', expm1(cmplx(1.0e-20_dp, 2.0e-20_dp, kind=dp)), &
        cmplx(1.0e-20_dp, 2.0e-20_dp, kind=dp), 2.0e-35_dp)
    call check_complex('log1p ordinary', log1p(cmplx(0.25_dp, 0.5_dp, kind=dp)), &
        cmplx(0.2973535538733464_dp, 0.3805063771123649_dp, kind=dp), 3.0e-15_dp)
    call check_complex('log1p tiny', log1p(cmplx(1.0e-20_dp, 2.0e-20_dp, kind=dp)), &
        cmplx(1.0e-20_dp, 2.0e-20_dp, kind=dp), 2.0e-35_dp)
    call check_complex('log1p upper cut', log1p(cmplx(-2.0_dp, 0.0_dp, kind=dp)), &
        cmplx(0.0_dp, acos(-1.0_dp), kind=dp), 2.0e-15_dp)
    call check_complex('log1p lower cut', log1p(cmplx(-2.0_dp, -0.0_dp, kind=dp)), &
        cmplx(0.0_dp, -acos(-1.0_dp), kind=dp), 2.0e-15_dp)

    call check_complex('xlogy complex', xlogy(cmplx(2.0_dp, 3.0_dp, kind=dp), &
        cmplx(-1.0_dp, 2.0_dp, kind=dp)), &
        cmplx(-4.493893894953008_dp, 6.483044740242557_dp, kind=dp), 4.0e-15_dp)
    call check_complex('xlogy real complex', xlogy(2.0_dp, cmplx(-1.0_dp, 2.0_dp, kind=dp)), &
        cmplx(1.6094379124341005_dp, 4.068887871591405_dp, kind=dp), 4.0e-15_dp)
    call check_complex('xlogy complex real', xlogy(cmplx(2.0_dp, 3.0_dp, kind=dp), 4.0_dp), &
        cmplx(2.772588722239781_dp, 4.1588830833596715_dp, kind=dp), 4.0e-15_dp)
    call check_complex('xlogy zero convention', xlogy(cmplx(0.0_dp, 0.0_dp, kind=dp), &
        cmplx(0.0_dp, 0.0_dp, kind=dp)), cmplx(0.0_dp, 0.0_dp, kind=dp), 0.0_dp)

    call check_complex('xlog1py complex', xlog1py(cmplx(2.0_dp, 3.0_dp, kind=dp), &
        cmplx(-2.0_dp, 2.0_dp, kind=dp)), &
        cmplx(-4.493893894953008_dp, 6.483044740242557_dp, kind=dp), 4.0e-15_dp)
    call check_complex('xlog1py real complex', xlog1py(2.0_dp, cmplx(-2.0_dp, 2.0_dp, kind=dp)), &
        cmplx(1.6094379124341005_dp, 4.068887871591405_dp, kind=dp), 4.0e-15_dp)
    call check_complex('xlog1py complex real', xlog1py(cmplx(2.0_dp, 3.0_dp, kind=dp), 4.0_dp), &
        cmplx(3.2188758248682006_dp, 4.828313737302301_dp, kind=dp), 4.0e-15_dp)

    call check_close('real ndtr preserved', ndtr(1.0_dp), 0.8413447460685429_dp, 3.0e-16_dp)
    call check_close('real log_ndtr preserved', log_ndtr(-10.0_dp), -53.23128515051247_dp, 2.0e-14_dp)
    call check_close('real expm1 preserved', expm1(1.0e-12_dp), 1.0000000000005e-12_dp, 2.0e-28_dp)
    call check_close('real log1p preserved', log1p(1.0e-12_dp), 9.999999999995e-13_dp, 2.0e-28_dp)

    input = [cmplx(0.1_dp, 0.2_dp, kind=dp), cmplx(-2.0_dp, 0.5_dp, kind=dp), &
        cmplx(3.0_dp, -1.0_dp, kind=dp)]
    actual = ndtr(input)
    call check_complex('ndtr elemental 1', actual(1), ndtr(input(1)), 0.0_dp)
    call check_complex('ndtr elemental 2', actual(2), ndtr(input(2)), 0.0_dp)
    call check_complex('ndtr elemental 3', actual(3), ndtr(input(3)), 0.0_dp)

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    nan_result = log1p(cmplx(nan_value, 0.0_dp, kind=dp))
    call check_true('complex nan real', ieee_is_nan(real(nan_result, dp)))
    call check_true('complex nan imaginary', ieee_is_nan(aimag(nan_result)))

    print *, 'special complex elementary tests passed'

contains

    subroutine check_close(name, actual_value, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual_value !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual_value == expected) return
        if (ieee_is_nan(actual_value) .or. abs(actual_value - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual_value
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_complex(name, actual_value, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual_value !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum scaled absolute error

        if (actual_value == expected) return
        if (ieee_is_nan(real(actual_value, dp)) .or. ieee_is_nan(aimag(actual_value)) .or. &
                abs(actual_value - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual_value
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_complex_elementary
