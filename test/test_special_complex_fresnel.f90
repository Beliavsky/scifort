! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_complex_fresnel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_constants, only : scifort_pi
    use scifort_special, only : dp, erf, fresnel
    implicit none

    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, kind=dp)
    complex(dp) :: array_c(3)
    complex(dp) :: array_s(3)
    complex(dp) :: array_z(3)
    complex(dp) :: c
    complex(dp) :: rhs
    complex(dp) :: s
    complex(dp) :: z
    real(dp) :: c_real
    real(dp) :: nan_value
    real(dp) :: pinf
    real(dp) :: s_real
    integer :: k

    z = cmplx(0.5_dp, 0.25_dp, kind=dp)
    call fresnel(z, s, c)
    call check_complex('complex fresnel S ordinary', s, &
        cmplx(0.017932781442058978_dp, 0.09014567263112916_dp, kind=dp), 3.0e-14_dp)
    call check_complex('complex fresnel C ordinary', c, &
        cmplx(0.5090783640851632_dp, 0.23999205598010134_dp, kind=dp), 3.0e-14_dp)

    z = cmplx(1.0_dp, 2.0_dp, kind=dp)
    call fresnel(z, s, c)
    call check_complex('complex fresnel S growing', s, &
        cmplx(36.725464883991435_dp, 15.587751104404592_dp, kind=dp), 5.0e-14_dp)
    call check_complex('complex fresnel C growing', c, &
        cmplx(16.08787137412547_dp, -36.22568799288166_dp, kind=dp), 5.0e-14_dp)

    z = cmplx(3.0_dp, -2.0_dp, kind=dp)
    call fresnel(z, s, c)
    call check_complex('complex fresnel S large', s, &
        cmplx(3788100.7001828905_dp, -5815898.602940464_dp, kind=dp), 6.0e-14_dp)
    call check_complex('complex fresnel C large', c, &
        cmplx(5815899.102940465_dp, 3788100.2001828905_dp, kind=dp), 6.0e-14_dp)

    z = cmplx(5.0_dp, 5.0_dp, kind=dp)
    call fresnel(z, s, c)
    call check_complex('complex fresnel S diagonal', s, &
        cmplx(-2.0608160801207947e32_dp, 2.0608160801207947e32_dp, kind=dp), 6.0e-14_dp)
    call check_complex('complex fresnel C diagonal', c, &
        cmplx(2.0608160801207947e32_dp, 2.0608160801207947e32_dp, kind=dp), 6.0e-14_dp)

    call fresnel(cmplx(0.0_dp, 2.0_dp, kind=dp), s, c)
    call check_complex('complex fresnel S imaginary axis', s, &
        cmplx(0.0_dp, -0.34341567836369824_dp, kind=dp), 5.0e-15_dp)
    call check_complex('complex fresnel C imaginary axis', c, &
        cmplx(0.0_dp, 0.48825340607534073_dp, kind=dp), 5.0e-15_dp)

    call fresnel(cmplx(0.0_dp, 4.0_dp, kind=dp), s, c)
    call check_complex('complex fresnel S imaginary-axis continuation', s, &
        cmplx(0.0_dp, -0.42051575424692844_dp, kind=dp), 5.0e-15_dp)
    call check_complex('complex fresnel C imaginary-axis continuation', c, &
        cmplx(0.0_dp, 0.4984260330381776_dp, kind=dp), 5.0e-15_dp)

    call fresnel(2.0_dp, s_real, c_real)
    call check_close('real fresnel S preserved', s_real, 0.34341567836369824_dp, 4.0e-15_dp)
    call check_close('real fresnel C preserved', c_real, 0.48825340607534073_dp, 4.0e-15_dp)

    z = cmplx(1.25_dp, -0.75_dp, kind=dp)
    call fresnel(z, s, c)
    rhs = 0.5_dp * cmplx(1.0_dp, 1.0_dp, kind=dp) * &
        erf(0.5_dp * sqrt(scifort_pi) * cmplx(1.0_dp, -1.0_dp, kind=dp) * z)
    call check_complex('complex fresnel erf identity', c + imaginary_unit * s, rhs, 5.0e-14_dp)

    call fresnel(-z, array_s(1), array_c(1))
    call check_complex('complex fresnel S odd', array_s(1), -s, 5.0e-14_dp)
    call check_complex('complex fresnel C odd', array_c(1), -c, 5.0e-14_dp)
    call fresnel(conjg(z), array_s(1), array_c(1))
    call check_complex('complex fresnel S conjugate', array_s(1), conjg(s), 5.0e-14_dp)
    call check_complex('complex fresnel C conjugate', array_c(1), conjg(c), 5.0e-14_dp)

    array_z = [cmplx(0.1_dp, 0.2_dp, kind=dp), &
        cmplx(-1.0_dp, 0.4_dp, kind=dp), cmplx(2.0_dp, -0.5_dp, kind=dp)]
    call fresnel(array_z, array_s, array_c)
    do k = 1, 3
        call check_elemental(k, array_z(k), array_s(k), array_c(k))
    end do

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    call fresnel(cmplx(nan_value, 0.0_dp, kind=dp), s, c)
    call check_nan_complex('complex fresnel nan S', s)
    call check_nan_complex('complex fresnel nan C', c)

    pinf = ieee_value(1.0_dp, ieee_positive_inf)
    call fresnel(cmplx(pinf, 0.0_dp, kind=dp), s, c)
    call check_nan_complex('complex fresnel infinite S', s)
    call check_nan_complex('complex fresnel infinite C', c)

    print *, 'special complex Fresnel tests passed'

contains

    subroutine check_elemental(k, input, s_value, c_value)
        integer, intent(in) :: k !! array position used in the check label
        complex(dp), intent(in) :: input !! complex Fresnel argument at this position
        complex(dp), intent(in) :: s_value !! elemental-array Fresnel S result
        complex(dp), intent(in) :: c_value !! elemental-array Fresnel C result

        complex(dp) :: c_scalar
        complex(dp) :: s_scalar
        character(len=40) :: label

        call fresnel(input, s_scalar, c_scalar)
        write(label, '(a,i0)') 'complex fresnel elemental S ', k
        call check_complex(trim(label), s_value, s_scalar, 0.0_dp)
        write(label, '(a,i0)') 'complex fresnel elemental C ', k
        call check_complex(trim(label), c_value, c_scalar, 0.0_dp)
    end subroutine check_elemental

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

    subroutine check_nan_complex(name, actual_value)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual_value !! result expected to have NaN components

        if (.not. ieee_is_nan(real(actual_value, dp)) .or. .not. ieee_is_nan(aimag(actual_value))) then
            print *, 'FAIL ', name
            print *, ' actual = ', actual_value
            error stop 1
        end if
    end subroutine check_nan_complex

end program test_special_complex_fresnel
