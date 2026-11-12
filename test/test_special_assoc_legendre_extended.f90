! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_assoc_legendre_extended
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : assoc_legendre_p, assoc_legendre_p_derivatives, dp, &
        sph_harm_y_derivatives, sph_legendre_p_derivatives
    implicit none

    complex(dp) :: assoc_derivs(0:2)
    complex(dp) :: assoc_high_derivs(0:4)
    complex(dp) :: harmonic_derivs(0:2, 0:2)
    real(dp) :: sph_derivs(0:2)

    call check_close('assoc real branch2', assoc_legendre_p(4, 2, 0.3_dp), &
        -2.5252500000000007_dp, 3.0e-14_dp)
    call check_close('assoc real branch3', assoc_legendre_p(4, 1, 2.0_dp, 3), &
        216.50635094610968_dp, 3.0e-13_dp)
    call check_close('assoc normalized negative order', assoc_legendre_p(4, -2, 0.3_dp, 2, .true.), &
        -0.28233153300906716_dp, 4.0e-15_dp)
    call check_true('assoc real branch2 cut', ieee_is_nan(assoc_legendre_p(4, 1, 2.0_dp)))
    call check_true('assoc real branch3 cut', ieee_is_nan(assoc_legendre_p(4, 1, 0.3_dp, 3)))

    call check_complex_close('assoc complex branch2', &
        assoc_legendre_p(4, 2, cmplx(0.3_dp, 0.2_dp, dp)), &
        cmplx(-3.8752499999999994_dp, 6.57_dp, dp), 2.0e-14_dp)
    call check_complex_close('assoc complex branch3', &
        assoc_legendre_p(3, 1, cmplx(2.0_dp, 0.4_dp, dp), 3), &
        cmplx(42.17911923464781_dp, 33.462434152503576_dp, dp), 8.0e-14_dp)
    call check_complex_close('assoc branch3 negative real', &
        assoc_legendre_p(2, 1, cmplx(-2.0_dp, 0.0_dp, dp), 3), &
        cmplx(10.392304845413264_dp, 0.0_dp, dp), 3.0e-14_dp)

    call assoc_legendre_p_derivatives(4, 2, cmplx(0.3_dp, 0.2_dp, dp), 2, assoc_derivs)
    call check_complex_close('assoc derivative 0', assoc_derivs(0), &
        cmplx(-3.8752499999999994_dp, 6.57_dp, dp), 2.0e-14_dp)
    call check_complex_close('assoc derivative 1', assoc_derivs(1), &
        cmplx(37.89_dp, 14.339999999999996_dp, dp), 8.0e-14_dp)
    call check_complex_close('assoc derivative 2', assoc_derivs(2), &
        cmplx(88.49999999999999_dp, -75.60000000000001_dp, dp), 2.0e-13_dp)

    call assoc_legendre_p_derivatives(4, 2, cmplx(0.3_dp, 0.2_dp, dp), 4, assoc_high_derivs)
    call check_complex_close('assoc derivative 3', assoc_high_derivs(3), &
        cmplx(-378.0_dp, -252.0_dp, dp), 8.0e-13_dp)
    call check_complex_close('assoc derivative 4', assoc_high_derivs(4), &
        cmplx(-1260.0_dp, 0.0_dp, dp), 2.0e-12_dp)

    call assoc_legendre_p_derivatives(4, -2, cmplx(0.3_dp, 0.2_dp, dp), 2, assoc_derivs, 2, .true.)
    call check_complex_close('assoc normalized derivative 0', assoc_derivs(0), &
        cmplx(-0.4332661214903031_dp, 0.7345483306086809_dp, dp), 5.0e-15_dp)
    call check_complex_close('assoc normalized derivative 1', assoc_derivs(1), &
        cmplx(4.236230783373352_dp, 1.6032607398673495_dp, dp), 2.0e-14_dp)
    call check_complex_close('assoc normalized derivative 2', assoc_derivs(2), &
        cmplx(9.89460080043657_dp, -8.452336954949203_dp, dp), 5.0e-14_dp)

    call sph_legendre_p_derivatives(3, 1, 0.7_dp, 2, sph_derivs)
    call check_close('sph derivative 0', sph_derivs(0), -0.4007648002460704_dp, 4.0e-15_dp)
    call check_close('sph derivative 1', sph_derivs(1), 0.5500405041508836_dp, 8.0e-15_dp)
    call check_close('sph derivative 2', sph_derivs(2), 3.1904864206724945_dp, 2.0e-14_dp)

    call sph_legendre_p_derivatives(4, -2, 1.2_dp, 2, sph_derivs)
    call check_close('sph negative order derivative 0', sph_derivs(0), &
        -0.02350309267922801_dp, 4.0e-15_dp)
    call check_close('sph negative order derivative 1', sph_derivs(1), &
        -1.392299007746487_dp, 1.0e-14_dp)
    call check_close('sph negative order derivative 2', sph_derivs(2), &
        0.9031369636683038_dp, 2.0e-14_dp)

    call sph_harm_y_derivatives(3, 1, 0.7_dp, 1.1_dp, 2, harmonic_derivs)
    call check_complex_close('harm value', harmonic_derivs(0, 0), &
        cmplx(-0.1817853589955138_dp, -0.35716453963284894_dp, dp), 5.0e-15_dp)
    call check_complex_close('harm theta', harmonic_derivs(1, 0), &
        cmplx(0.24949623930980994_dp, 0.49020014563117_dp, dp), 8.0e-15_dp)
    call check_complex_close('harm phi', harmonic_derivs(0, 1), &
        cmplx(0.35716453963284894_dp, -0.1817853589955138_dp, dp), 8.0e-15_dp)
    call check_complex_close('harm theta theta', harmonic_derivs(2, 0), &
        cmplx(1.4471922658780163_dp, 2.8433849802793922_dp, dp), 3.0e-14_dp)
    call check_complex_close('harm theta phi', harmonic_derivs(1, 1), &
        cmplx(-0.49020014563117_dp, 0.24949623930980994_dp, dp), 1.0e-14_dp)
    call check_complex_close('harm phi phi', harmonic_derivs(0, 2), &
        cmplx(0.1817853589955138_dp, 0.35716453963284894_dp, dp), 8.0e-15_dp)

    print *, 'special extended associated-Legendre tests passed'

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

    subroutine check_complex_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum absolute complex error

        if (actual == expected) return
        if (complex_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

    pure function complex_is_nan(z) result(is_nan)
        complex(dp), intent(in) :: z !! complex value to inspect
        logical :: is_nan

        is_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
    end function complex_is_nan

end program test_special_assoc_legendre_extended
