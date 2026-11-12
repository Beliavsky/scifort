program test_special_complex_digamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_special, only : digamma, psi
    implicit none

    complex(dp), parameter :: z = cmplx(3.0_dp, 4.0_dp, kind=dp)
    complex(dp), parameter :: zn = cmplx(-2.75_dp, 0.625_dp, kind=dp)
    complex(dp) :: values(3)
    complex(dp) :: expected(3)
    complex(dp) :: near_pole
    complex(dp) :: lhs
    complex(dp) :: rhs

    call check_complex_rel('digamma 3+4i', digamma(z), &
        cmplx(1.55035981733341_dp, 1.0105022091860445_dp, kind=dp), 8.0e-15_dp)
    call check_complex_rel('psi alias 3+4i', psi(z), digamma(z), 2.0e-16_dp)
    call check_complex_rel('digamma negative complex', digamma(zn), &
        cmplx(1.0765535364330212_dp, 2.9505376308648430_dp, kind=dp), 2.0e-14_dp)

    lhs = digamma(z + 1.0_dp)
    rhs = digamma(z) + 1.0_dp / z
    call check_complex_rel('digamma recurrence', lhs, rhs, 8.0e-15_dp)

    lhs = digamma(zn)
    rhs = digamma(1.0_dp - zn) - acos(-1.0_dp) / tan(acos(-1.0_dp) * zn)
    call check_complex_rel('digamma reflection', lhs, rhs, 2.0e-14_dp)

    call check_complex_rel('digamma conjugation', digamma(conjg(z)), conjg(digamma(z)), 4.0e-15_dp)
    call check_complex_rel('complex real-axis dispatch positive', &
        digamma(cmplx(0.5_dp, 0.0_dp, kind=dp)), &
        cmplx(digamma(0.5_dp), 0.0_dp, kind=dp), 2.0e-16_dp)
    call check_complex_rel('complex real-axis dispatch negative', &
        digamma(cmplx(-0.5_dp, -0.0_dp, kind=dp)), &
        cmplx(digamma(-0.5_dp), -0.0_dp, kind=dp), 2.0e-16_dp)

    near_pole = digamma(cmplx(-3.0_dp, 1.0e-8_dp, kind=dp))
    call check_complex_rel('digamma near pole', near_pole, &
        cmplx(1.2561176684318005_dp, 1.0000000000000001e8_dp, kind=dp), 4.0e-14_dp)

    call check_complex_rel('digamma ultra-near pole', &
        digamma(cmplx(-3.0_dp, 1.0e-300_dp, kind=dp)), &
        cmplx(1.2561176684318005_dp, 1.0e300_dp, kind=dp), 4.0e-14_dp)

    call check_true('complex digamma pole zero', complex_is_nan(digamma(cmplx(0.0_dp, 0.0_dp, kind=dp))))
    call check_true('complex digamma pole negative integer', &
        complex_is_nan(digamma(cmplx(-4.0_dp, -0.0_dp, kind=dp))))
    call check_true('complex digamma +inf real', &
        real(digamma(cmplx(ieee_value(0.0_dp, ieee_positive_inf), 0.0_dp, kind=dp)), dp) > huge(1.0_dp))

    values = [z, zn, cmplx(0.5_dp, 0.25_dp, kind=dp)]
    expected = [ &
        cmplx(1.55035981733341_dp, 1.0105022091860445_dp, kind=dp), &
        cmplx(1.0765535364330212_dp, 2.9505376308648430_dp, kind=dp), &
        cmplx(-1.5381612557092352_dp, 1.0301191246287900_dp, kind=dp)]
    call check_complex_array('digamma elemental array', digamma(values), expected, 2.0e-14_dp)

    call check_complex_rel('digamma large complex', digamma(cmplx(100.0_dp, 200.0_dp, kind=dp)), &
        cmplx(5.408890142204208_dp, 1.1091500511306238_dp, kind=dp), 2.0e-14_dp)

    print *, 'special complex digamma tests passed'

contains

    pure logical function complex_is_nan(value) result(is_nan)
        complex(dp), intent(in) :: value !! complex value to classify

        is_nan = ieee_is_nan(real(value, dp)) .and. ieee_is_nan(aimag(value))
    end function complex_is_nan

    subroutine check_complex_rel(name, actual, reference, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: reference !! expected complex value
        real(dp), intent(in) :: tolerance !! maximum scaled absolute error

        real(dp) :: scale

        scale = max(1.0_dp, abs(reference))
        if (abs(actual - reference) <= tolerance * scale) return
        print *, 'FAIL ', trim(name)
        print *, ' actual  = ', actual
        print *, ' expected= ', reference
        print *, ' scaled tolerance= ', tolerance
        error stop 1
    end subroutine check_complex_rel

    subroutine check_complex_array(name, actual, reference, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual(:) !! computed complex values
        complex(dp), intent(in) :: reference(:) !! expected complex values
        real(dp), intent(in) :: tolerance !! maximum scaled absolute error

        integer :: i

        do i = 1, size(actual)
            call check_complex_rel(name, actual(i), reference(i), tolerance)
        end do
    end subroutine check_complex_array

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (condition) return
        print *, 'FAIL ', trim(name)
        error stop 1
    end subroutine check_true

end program test_special_complex_digamma
