! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_lambert_w
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    use scifort_special, only : dp, lambertw, wrightomega
    implicit none

    complex(dp) :: value
    complex(dp) :: z
    real(dp), parameter :: pi = acos(-1.0_dp)

    call check_complex_close('W0 1', lambertw(1.0_dp), &
        cmplx(0.56714329040978384_dp, 0.0_dp, dp), 3.0e-15_dp)
    call check_complex_close('W0 -0.1', lambertw(-0.1_dp), &
        cmplx(-0.11183255915896297_dp, 0.0_dp, dp), 3.0e-15_dp)
    call check_complex_close('W-1 -0.1', lambertw(-0.1_dp, -1), &
        cmplx(-3.5771520639572971_dp, 0.0_dp, dp), 5.0e-15_dp)
    call check_complex_close('W0 -2', lambertw(-2.0_dp), &
        cmplx(0.17281600283999998_dp, 1.6736864137408427_dp, dp), 5.0e-15_dp)
    call check_complex_close('W-1 -2', lambertw(-2.0_dp, -1), &
        cmplx(0.17281600283999998_dp, -1.6736864137408427_dp, dp), 5.0e-15_dp)
    call check_complex_close('W3 complex', lambertw(cmplx(1.0_dp, 2.0_dp, dp), 3), &
        cmplx(-2.1072075066264553_dp, 18.27108543113687_dp, dp), 3.0e-14_dp)
    call check_complex_close('W-4 complex', lambertw(cmplx(1.0_dp, 2.0_dp, dp), -4), &
        cmplx(-2.3074940367904766_dp, -22.351925909251833_dp, dp), 4.0e-14_dp)

    value = lambertw(-exp(-1.0_dp))
    call check_complex_close('W branch point', value, cmplx(-1.0_dp, 0.0_dp, dp), 5.0e-15_dp)
    call check_complex_close('W near branch point', &
        lambertw(cmplx(-exp(-1.0_dp), 1.0e-8_dp, dp)), &
        cmplx(-0.9998351278743025_dp, 0.00016485400656392362_dp, dp), 8.0e-13_dp)

    z = cmplx(-0.7_dp, 2.3_dp, dp)
    value = lambertw(z, 2)
    call check_complex_close('W defining equation', value * exp(value), z, 3.0e-14_dp)
    call check_true('W invalid tol', complex_is_nan(lambertw(1.0_dp, 0, -1.0_dp)))
    call check_true('W nonprincipal zero real is -inf', &
        .not. ieee_is_finite(real(lambertw(0.0_dp, 1), dp)) .and. &
        real(lambertw(0.0_dp, 1), dp) < 0.0_dp)

    call check_close('omega 0', wrightomega(0.0_dp), &
        0.56714329040978384_dp, 3.0e-15_dp)
    call check_close('omega 1', wrightomega(1.0_dp), 1.0_dp, 3.0e-15_dp)
    call check_close('omega -1', wrightomega(-1.0_dp), &
        0.27846454276107374_dp, 3.0e-15_dp)
    call check_complex_close('omega upper strip', wrightomega(cmplx(1.0_dp, 4.0_dp, dp)), &
        cmplx(0.09583285722030857_dp, 2.4680141222274354_dp, dp), 6.0e-15_dp)
    call check_complex_close('omega lower strip', wrightomega(cmplx(1.0_dp, -4.0_dp, dp)), &
        cmplx(0.09583285722030857_dp, -2.4680141222274354_dp, dp), 6.0e-15_dp)
    call check_complex_close('omega large positive', &
        wrightomega(cmplx(1000.0_dp, 100.0_dp, dp)), &
        cmplx(993.0941403447679_dp, 99.89974283382134_dp, dp), 3.0e-12_dp)
    call check_complex_close('omega large negative noncentral', &
        wrightomega(cmplx(-1000.0_dp, 4.0_dp, dp)), &
        cmplx(-1006.914646493052_dp, 0.8592607062263946_dp, dp), 3.0e-12_dp)
    call check_complex_close('omega branch point upper', &
        wrightomega(cmplx(-1.0_dp, pi, dp)), cmplx(-1.0_dp, 0.0_dp, dp), 5.0e-15_dp)
    call check_complex_close('omega branch point lower', &
        wrightomega(cmplx(-1.0_dp, -pi, dp)), cmplx(-1.0_dp, 0.0_dp, dp), 5.0e-15_dp)

    z = cmplx(3.0_dp, 8.0_dp, dp)
    value = wrightomega(z)
    call check_complex_close('omega defining equation', value + log(value), z, 4.0e-14_dp)

    print *, 'special Lambert-W tests passed'

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
        complex(dp), intent(in) :: expected !! complex reference value
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

end program test_special_lambert_w
