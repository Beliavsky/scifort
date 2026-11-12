! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_associated_legendre
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_special, only : dp, lpmv, sph_harm_y, sph_legendre_p
    implicit none

    complex(dp) :: harmonics(3)
    real(dp) :: values(3)

    call check_close('lpmv noninteger m1', lpmv(1, 2.5_dp, 0.3_dp), &
        -0.12922999435771895_dp, 2.0e-13_dp)
    call check_close('lpmv noninteger m2', lpmv(2, 5.3_dp, -0.2_dp), &
        11.168937127392732_dp, 3.0e-12_dp)
    call check_close('lpmv negative order', lpmv(-2, 2.5_dp, -0.2_dp), &
        0.02132223825959121_dp, 2.0e-13_dp)
    call check_rel('lpmv difficult continuation', lpmv(5, 10.25_dp, -0.95_dp), &
        -113163.20899160246_dp, 2.0e-12_dp)

    call check_close('lpmv integer degree', lpmv(2, 3.0_dp, 0.3_dp), &
        4.095_dp, 2.0e-14_dp)
    call check_close('lpmv integer negative order', lpmv(-2, 3.0_dp, 0.3_dp), &
        0.034125_dp, 2.0e-15_dp)
    call check_close('lpmv negative integer degree symmetry', lpmv(4, -5.0_dp, 0.2_dp), &
        96.768_dp, 3.0e-13_dp)
    call check_close('lpmv positive order above degree', lpmv(5, 3.0_dp, 0.3_dp), &
        0.0_dp, 0.0_dp)
    call check_true('lpmv negative order below degree domain', ieee_is_nan(lpmv(-5, 3.0_dp, 0.3_dp)))
    call check_true('lpmv x outside Ferrers interval', ieee_is_nan(lpmv(1, 2.5_dp, 1.1_dp)))
    call check_true('lpmv noninteger x=-1 negative infinity', &
        .not. ieee_is_finite(lpmv(0, 2.5_dp, -1.0_dp)) .and. lpmv(0, 2.5_dp, -1.0_dp) < 0.0_dp)
    call check_true('lpmv noninteger m1 x=-1 positive infinity', &
        .not. ieee_is_finite(lpmv(1, 2.5_dp, -1.0_dp)) .and. lpmv(1, 2.5_dp, -1.0_dp) > 0.0_dp)

    call check_close('sph legendre 00', sph_legendre_p(0, 0, 0.7_dp), &
        0.28209479177387814_dp, 2.0e-15_dp)
    call check_close('sph legendre 11', sph_legendre_p(1, 1, 0.7_dp), &
        -0.22257344192657688_dp, 2.0e-15_dp)
    call check_close('sph legendre negative order', sph_legendre_p(1, -1, 0.7_dp), &
        0.22257344192657688_dp, 2.0e-15_dp)
    call check_close('sph legendre 21', sph_legendre_p(2, 1, 0.7_dp), &
        -0.3806538080852601_dp, 3.0e-15_dp)
    call check_close('sph legendre 3-2', sph_legendre_p(3, -2, 1.9_dp), &
        -0.29586539627892683_dp, 4.0e-15_dp)
    call check_close('sph legendre high degree', sph_legendre_p(12, 5, 0.7_dp), &
        -0.09795458270632562_dp, 5.0e-15_dp)
    call check_close('sph legendre order above degree', sph_legendre_p(2, 3, 0.7_dp), 0.0_dp, 0.0_dp)
    call check_close('sph legendre negative degree', sph_legendre_p(-1, 0, 0.7_dp), 0.0_dp, 0.0_dp)

    call check_complex_close('sph harmonic 00', sph_harm_y(0, 0, 0.7_dp, 1.2_dp), &
        cmplx(0.28209479177387814_dp, 0.0_dp, dp), 2.0e-15_dp)
    call check_complex_close('sph harmonic 11', sph_harm_y(1, 1, 0.7_dp, 1.2_dp), &
        cmplx(-0.08065121262265872_dp, -0.20744714737382625_dp, dp), 3.0e-15_dp)
    call check_complex_close('sph harmonic 1-1', sph_harm_y(1, -1, 0.7_dp, 1.2_dp), &
        cmplx(0.08065121262265872_dp, -0.20744714737382625_dp, dp), 3.0e-15_dp)
    call check_complex_close('sph harmonic 3-2', sph_harm_y(3, -2, 1.9_dp, 5.0_dp), &
        cmplx(0.2482522304565697_dp, -0.1609570215573854_dp, dp), 4.0e-15_dp)

    values = lpmv(1, 2.5_dp, [-0.2_dp, 0.3_dp, 0.8_dp])
    call check_close('lpmv elemental middle', values(2), -0.12922999435771895_dp, 2.0e-13_dp)
    values = sph_legendre_p(2, 1, [-0.7_dp, 0.7_dp, 1.9_dp])
    call check_close('sph legendre elemental symmetry', values(1), values(2), 3.0e-15_dp)
    harmonics = sph_harm_y(1, 1, 0.7_dp, [-1.2_dp, 0.0_dp, 1.2_dp])
    call check_complex_close('sph harmonic conjugate phi', harmonics(1), conjg(harmonics(3)), 3.0e-15_dp)

    print *, 'special associated-Legendre tests passed'

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

    subroutine check_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! nonzero reference value
        real(dp), intent(in) :: tolerance !! maximum relative error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' relative tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_rel

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

end program test_special_associated_legendre
