! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_hankel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, h1vp, h2vp, hankel1, hankel1e, hankel2, &
        hankel2e, jv, yv
    implicit none

    real(dp), parameter :: orders(3) = [0.0_dp, -0.5_dp, 2.3_dp]
    real(dp), parameter :: arguments(3) = [1.0_dp, 2.0_dp, -2.0_dp]
    complex(dp), parameter :: h1_refs(3) = [ &
        cmplx(0.76519768655796638_dp, 0.08825696421567700_dp, kind=dp), &
        cmplx(-0.23478571040624860_dp, 0.51301613656182787_dp, kind=dp), &
        cmplx(-0.75839486069272488_dp, -0.21471258543619393_dp, kind=dp)]
    complex(dp) :: nan_complex
    complex(dp) :: values(3)
    real(dp) :: nan_value
    real(dp) :: pos_inf

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)
    nan_complex = cmplx(nan_value, nan_value, kind=dp)

    call check_complex_rel('hankel1 positive', hankel1(2.3_dp, 3.5_dp), &
        cmplx(0.46694776765903079_dp, -0.08889279748125363_dp, kind=dp), 3.0e-12_dp)
    call check_complex_rel('hankel2 positive', hankel2(2.3_dp, 3.5_dp), &
        cmplx(0.46694776765903079_dp, 0.08889279748125363_dp, kind=dp), 3.0e-12_dp)
    call check_complex_rel('hankel1 negative axis', hankel1(2.3_dp, -2.0_dp), &
        cmplx(-0.75839486069272488_dp, -0.21471258543619393_dp, kind=dp), 5.0e-12_dp)
    call check_complex_rel('hankel2 negative axis', hankel2(2.3_dp, -2.0_dp), &
        cmplx(1.0782290174751057_dp, 0.65492653638060472_dp, kind=dp), 7.0e-12_dp)
    call check_complex_rel('hankel1 negative order', hankel1(-0.5_dp, -1.0_dp), &
        cmplx(-0.67139670714180311_dp, -0.43109886801837627_dp, kind=dp), 4.0e-12_dp)
    call check_complex_rel('hankel2 negative order', hankel2(-0.5_dp, -1.0_dp), &
        cmplx(0.67139670714180455_dp, -0.43109886801837616_dp, kind=dp), 4.0e-12_dp)

    call check_complex_rel('hankel1 J/Y identity', hankel1(1.0_dp, 0.1_dp), &
        cmplx(jv(1.0_dp, 0.1_dp), yv(1.0_dp, 0.1_dp), kind=dp), 3.0e-12_dp)
    call check_complex_rel('hankel2 J/Y identity', hankel2(1.0_dp, 0.1_dp), &
        cmplx(jv(1.0_dp, 0.1_dp), -yv(1.0_dp, 0.1_dp), kind=dp), 3.0e-12_dp)

    call check_complex_rel('hankel1e positive', hankel1e(2.3_dp, 3.5_dp), &
        cmplx(-0.40609425722097447_dp, 0.24704169975520393_dp, kind=dp), 4.0e-12_dp)
    call check_complex_rel('hankel2e positive', hankel2e(2.3_dp, 3.5_dp), &
        cmplx(-0.40609425722097447_dp, -0.24704169975520393_dp, kind=dp), 4.0e-12_dp)
    call check_complex_rel('hankel1e negative axis', hankel1e(2.3_dp, -2.0_dp), &
        cmplx(0.51084122357510875_dp, -0.60025453214958602_dp, kind=dp), 7.0e-12_dp)
    call check_complex_rel('hankel2e negative axis', hankel2e(2.3_dp, -2.0_dp), &
        cmplx(0.14682141959514100_dp, -1.2529764774044623_dp, kind=dp), 7.0e-12_dp)

    call check_complex_rel('h1vp default', h1vp(0.0_dp, 1.0_dp), &
        cmplx(-0.44005058574493355_dp, 0.78121282130028891_dp, kind=dp), 5.0e-12_dp)
    call check_complex_rel('h2vp n2', h2vp(2.5_dp, 3.0_dp, 2), &
        cmplx(-0.17070357259041172_dp, 0.01875156551600621_dp, kind=dp), 5.0e-12_dp)
    call check_complex_rel('h1vp negative axis n1', h1vp(2.3_dp, -2.0_dp, 1), &
        cmplx(-0.32250194855090547_dp, -0.51102022633878308_dp, kind=dp), 8.0e-12_dp)
    call check_complex_rel('h1vp negative axis n2', h1vp(2.3_dp, -2.0_dp, 2), &
        cmplx(-0.40583331684885654_dp, -0.32475492197256378_dp, kind=dp), 8.0e-12_dp)
    call check_complex_rel('h2vp negative axis n2', h2vp(2.3_dp, -2.0_dp, 2), &
        cmplx(0.37739695716693289_dp, 0.28561563062234896_dp, kind=dp), 8.0e-12_dp)
    call check_complex_rel('h1vp n20', h1vp(0.0_dp, 1.0_dp, 20), &
        cmplx(5.2723352338201819_dp, -7.7328768259473408e16_dp, kind=dp), 7.0e-11_dp)
    call check_complex_rel('h2vp n0', h2vp(-0.5_dp, -1.0_dp, 0), &
        hankel2(-0.5_dp, -1.0_dp), 0.0_dp)

    values = hankel1(orders, arguments)
    call check_complex_array_rel('hankel1 elemental arrays', values, h1_refs, 7.0e-12_dp)

    call check_true('hankel1 zero nan', complex_is_nan(hankel1(1.0_dp, 0.0_dp)))
    call check_true('hankel2 zero order0 real nan', &
        ieee_is_nan(real(hankel2(0.0_dp, 0.0_dp), dp)))
    call check_true('hankel2 zero order0 imag +inf', &
        aimag(hankel2(0.0_dp, 0.0_dp)) > huge(1.0_dp))
    call check_true('hankel2 zero nonzero order nan', &
        complex_is_nan(hankel2(1.0_dp, 0.0_dp)))
    call check_true('scaled zero nan', complex_is_nan(hankel1e(0.0_dp, 0.0_dp)))
    call check_true('negative derivative order nan', complex_is_nan(h1vp(0.0_dp, 1.0_dp, -1)))
    call check_true('nan order', complex_is_nan(hankel1(nan_value, 1.0_dp)))
    call check_true('nan argument', complex_is_nan(hankel2(0.0_dp, nan_value)))
    call check_true('infinite argument', complex_is_nan(hankel1(0.0_dp, pos_inf)))
    call check_true('infinite order', complex_is_nan(hankel2(pos_inf, 1.0_dp)))
    call check_true('nan complex helper sanity', complex_is_nan(nan_complex))

    print *, 'Hankel real-axis special-function tests passed'

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

end program test_special_hankel
