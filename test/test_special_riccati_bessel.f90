! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_riccati_bessel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_special, only : dp, riccati_jn, riccati_yn, &
        spherical_bessel_status_invalid_argument, spherical_bessel_status_success
    implicit none

    real(dp), allocatable :: derivatives(:)
    real(dp), allocatable :: values(:)
    integer :: status

    call riccati_jn(3, 1.0_dp, values, derivatives, status)
    call check_true('riccati_jn status', status == spherical_bessel_status_success)
    call check_true('riccati_jn bounds', lbound(values, 1) == 0 .and. ubound(values, 1) == 3)
    call check_rel('riccati_jn n0', values(0), 0.8414709848078965_dp, 3.0e-14_dp)
    call check_rel('riccati_jn n1', values(1), 0.30116867893975674_dp, 4.0e-14_dp)
    call check_rel('riccati_jn n3', values(3), 0.009006581117112517_dp, 4.0e-13_dp)
    call check_rel('riccati_jn derivative n0', derivatives(0), 0.5403023058681398_dp, 3.0e-14_dp)
    call check_rel('riccati_jn derivative n2', derivatives(2), 0.17709857491700903_dp, 5.0e-13_dp)

    call riccati_jn(3, -1.0_dp, values, derivatives, status)
    call check_rel('riccati_jn negative n0', values(0), -0.8414709848078965_dp, 3.0e-14_dp)
    call check_rel('riccati_jn negative n1', values(1), 0.30116867893975674_dp, 4.0e-14_dp)
    call check_rel('riccati_jn negative derivative n1', derivatives(1), -0.5403023058681398_dp, 4.0e-14_dp)

    call riccati_jn(3, 0.0_dp, values, derivatives, status)
    call check_close('riccati_jn zero n0', values(0), 0.0_dp, 0.0_dp)
    call check_close('riccati_jn zero derivative n0', derivatives(0), 1.0_dp, 0.0_dp)
    call check_close('riccati_jn zero derivative n2', derivatives(2), 0.0_dp, 0.0_dp)

    call riccati_yn(3, 1.0_dp, values, derivatives, status)
    call check_true('riccati_yn status', status == spherical_bessel_status_success)
    call check_rel('riccati_yn n0', values(0), -0.5403023058681398_dp, 3.0e-14_dp)
    call check_rel('riccati_yn n1', values(1), -1.3817732906760363_dp, 4.0e-14_dp)
    call check_rel('riccati_yn n3', values(3), -16.64331454012381_dp, 4.0e-13_dp)
    call check_rel('riccati_yn derivative n0', derivatives(0), 0.8414709848078965_dp, 3.0e-14_dp)
    call check_rel('riccati_yn derivative n2', derivatives(2), 5.8282618416439025_dp, 5.0e-13_dp)

    call riccati_yn(3, 0.0_dp, values, derivatives, status)
    call check_close('riccati_yn zero n0', values(0), -1.0_dp, 0.0_dp)
    call check_close('riccati_yn zero derivative n0', derivatives(0), 0.0_dp, 0.0_dp)
    call check_true('riccati_yn zero n1', .not. ieee_is_finite(values(1)) .and. values(1) < 0.0_dp)
    call check_true('riccati_yn zero derivative n1', &
        .not. ieee_is_finite(derivatives(1)) .and. derivatives(1) > 0.0_dp)

    call riccati_jn(-1, 1.0_dp, values, derivatives, status)
    call check_true('riccati invalid order', status == spherical_bessel_status_invalid_argument)
    call check_true('riccati invalid empty', size(values) == 0 .and. size(derivatives) == 0)

    print *, 'Riccati-Bessel special-function tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual == expected) return
        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! nonzero reference value
        real(dp), intent(in) :: tolerance !! maximum relative error

        if (actual == expected) return
        if (abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_rel

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_riccati_bessel
