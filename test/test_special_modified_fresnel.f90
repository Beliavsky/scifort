! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_modified_fresnel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_special, only : modfresnelm, modfresnelp
    implicit none

    complex(dp) :: fm
    complex(dp) :: fp
    complex(dp) :: km
    complex(dp) :: kp
    real(dp) :: nan_value
    real(dp) :: pinf

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    pinf = ieee_value(1.0_dp, ieee_positive_inf)

    call modfresnelp(0.0_dp, fp, kp)
    call check_complex('modfresnelp F 0', fp, &
        cmplx(0.6266570686577501_dp, 0.6266570686577501_dp, kind=dp), 4.0e-16_dp)
    call check_complex('modfresnelp K 0', kp, cmplx(0.5_dp, 0.0_dp, kind=dp), 0.0_dp)

    call modfresnelp(1.0_dp, fp, kp)
    call check_complex('modfresnelp F 1', fp, &
        cmplx(-0.27786716924252196_dp, 0.31638876693436902_dp, kind=dp), 7.0e-16_dp)
    call check_complex('modfresnelp K 1', kp, &
        cmplx(0.20779404795392433_dp, 0.11515989377745532_dp, kind=dp), 7.0e-16_dp)

    call modfresnelp(10.0_dp, fp, kp)
    call check_complex('modfresnelp F 10', fp, &
        cmplx(0.025531883844305777_dp, 0.042986168728126783_dp, kind=dp), 2.0e-15_dp)
    call check_complex('modfresnelp K 10', kp, &
        cmplx(0.020045317520283938_dp, 0.019845921064366676_dp, kind=dp), 2.0e-15_dp)

    ! The documented integral definition extends to negative lower limits.
    ! These high-precision values intentionally differ from SciPy 1.17.0's
    ! first returned component for negative x; its scaled K value is correct.
    call modfresnelp(-1.0_dp, fp, kp)
    call check_complex('modfresnelp F -1', fp, &
        cmplx(1.5311813065580222_dp, 0.93692537038113123_dp, kind=dp), 2.0e-15_dp)
    call check_complex('modfresnelp K -1', kp, &
        cmplx(0.33250825791421539_dp, -0.95663087858535183_dp, kind=dp), 2.0e-15_dp)

    ! Exact binary64 x = 100000000.25. The compensated square phase keeps
    ! the oscillatory F value accurate where SciPy 1.17.0 loses phase bits.
    call modfresnelp(100000000.25_dp, fp, kp)
    call check_complex('modfresnelp F large phase', fp, &
        cmplx(6.9613199638473871e-10_dp, -4.9513028708219212e-9_dp, kind=dp), 2.0e-23_dp)
    call check_complex('modfresnelp K large phase', kp, &
        cmplx(1.9947113970203850e-9_dp, 1.9947113970203848e-9_dp, kind=dp), 2.0e-23_dp)

    call modfresnelm(3.0_dp, fm, km)
    call check_complex('modfresnelm F 3', fm, &
        cmplx(-0.076206489072518605_dp, 0.14690545823601889_dp, kind=dp), 2.0e-15_dp)
    call check_complex('modfresnelm K 3', km, &
        cmplx(0.069474878530462927_dp, -0.062380421395377479_dp, kind=dp), 2.0e-15_dp)

    call modfresnelp(nan_value, fp, kp)
    call require(ieee_is_nan(real(fp, dp)) .and. ieee_is_nan(aimag(fp)), &
        'modfresnelp NaN F')
    call require(ieee_is_nan(real(kp, dp)) .and. ieee_is_nan(aimag(kp)), &
        'modfresnelp NaN K')

    call modfresnelm(pinf, fm, km)
    call require(ieee_is_nan(real(fm, dp)) .and. ieee_is_nan(aimag(fm)), &
        'modfresnelm infinity F')
    call require(ieee_is_nan(real(km, dp)) .and. ieee_is_nan(aimag(km)), &
        'modfresnelm infinity K')

    print *, 'special modified Fresnel tests passed'

contains

    subroutine check_complex(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum scaled absolute error

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex

    subroutine require(condition, name)
        logical, intent(in) :: condition !! assertion condition expected to be true
        character(len=*), intent(in) :: name !! check label printed on failure

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine require

end program test_special_modified_fresnel
