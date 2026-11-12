program test_special_spheroidal_radial2
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_special, only : obl_cv, obl_rad1, obl_rad2, obl_rad2_cv, &
        pro_cv, pro_rad1, pro_rad2, pro_rad2_cv
    implicit none

    real(dp) :: cv
    real(dp) :: derivative
    real(dp) :: derivative_cv
    real(dp) :: value
    real(dp) :: value_cv

    ! SciPy 1.17.0 references in the stable large-argument expansion region.
    call pro_rad2(0, 0, 2.0_dp, 1.5_dp, value, derivative)
    call check_close(value, 0.3068766159404849_dp, 2.0e-10_dp, 'pro 0 0 value')
    call check_close(derivative, 0.33944957398859427_dp, 2.0e-10_dp, &
        'pro 0 0 derivative')

    call pro_rad2(1, 2, 2.0_dp, 2.0_dp, value, derivative)
    call check_close(value, -0.06895436445035619_dp, 2.0e-11_dp, 'pro 1 2 value')
    call check_close(derivative, 0.5751025563707028_dp, 2.0e-11_dp, &
        'pro 1 2 derivative')

    call pro_rad2(2, 4, 3.0_dp, 2.0_dp, value, derivative)
    call check_close(value, -0.10388988124556839_dp, 2.0e-11_dp, 'pro 2 4 value')
    call check_close(derivative, 0.5303964082662372_dp, 2.0e-11_dp, &
        'pro 2 4 derivative')

    call pro_rad2(0, 3, 10.0_dp, 1.2_dp, value, derivative)
    call check_close(value, -0.05229478087130487_dp, 3.0e-10_dp, &
        'pro high-c value')
    call check_close(derivative, -1.5011595561304383_dp, 3.0e-10_dp, &
        'pro high-c derivative')

    call obl_rad2(0, 0, 2.0_dp, 1.5_dp, value, derivative)
    call check_close(value, 0.2851395143556178_dp, 2.0e-10_dp, 'obl 0 0 value')
    call check_close(derivative, -0.26212372708626147_dp, 2.0e-10_dp, &
        'obl 0 0 derivative')

    call obl_rad2(1, 2, 2.0_dp, 2.0_dp, value, derivative)
    call check_close(value, 0.061557357045144964_dp, 2.0e-11_dp, 'obl 1 2 value')
    call check_close(derivative, 0.35454761268210117_dp, 2.0e-11_dp, &
        'obl 1 2 derivative')

    call obl_rad2(2, 4, 3.0_dp, 2.0_dp, value, derivative)
    call check_close(value, 0.02908347659751293_dp, 2.0e-11_dp, 'obl 2 4 value')
    call check_close(derivative, 0.35103597959983474_dp, 2.0e-11_dp, &
        'obl 2 4 derivative')

    call obl_rad2(3, 5, 5.0_dp, 3.0_dp, value, derivative)
    call check_close(value, 0.04717545263511711_dp, 3.0e-11_dp, 'obl 3 5 value')
    call check_close(derivative, 0.196301190885361_dp, 3.0e-11_dp, &
        'obl 3 5 derivative')

    ! Near the singular/small-argument region, verify the independent radial
    ! Wronskian rather than legacy SciPy values that can lose normalization.
    call check_pro_wronskian(4, 4, 2.0_dp, 1.05_dp, 3.0e-10_dp)
    call check_pro_wronskian(5, 10, 1.0_dp, 1.001_dp, 3.0e-10_dp)
    call check_obl_wronskian(1, 2, 2.0_dp, 0.5_dp, 3.0e-10_dp)
    call check_obl_wronskian(5, 10, 0.2_dp, 0.0_dp, 3.0e-10_dp)

    cv = pro_cv(2, 4, 3.0_dp)
    call pro_rad2(2, 4, 3.0_dp, 1.4_dp, value, derivative)
    call pro_rad2_cv(2, 4, 3.0_dp, cv, 1.4_dp, value_cv, derivative_cv)
    call check_close(value_cv, value, 5.0e-13_dp, 'pro cv value agreement')
    call check_close(derivative_cv, derivative, 5.0e-13_dp, &
        'pro cv derivative agreement')

    cv = obl_cv(2, 4, 3.0_dp)
    call obl_rad2(2, 4, 3.0_dp, 0.7_dp, value, derivative)
    call obl_rad2_cv(2, 4, 3.0_dp, cv, 0.7_dp, value_cv, derivative_cv)
    call check_close(value_cv, value, 5.0e-13_dp, 'obl cv value agreement')
    call check_close(derivative_cv, derivative, 5.0e-13_dp, &
        'obl cv derivative agreement')

    call pro_rad2(-1, 0, 1.0_dp, 2.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'negative m')
    call pro_rad2(1, 0, 1.0_dp, 2.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'n below m')
    call pro_rad2(0, 0, 0.0_dp, 2.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'zero c')
    call pro_rad2(0, 0, -1.0_dp, 2.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'negative c')
    call pro_rad2(0, 0, 1.0_dp, 1.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'pro endpoint x')
    call obl_rad2(0, 0, 1.0_dp, -0.1_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'negative obl x')

    print *, 'test_special_spheroidal_radial2: PASS'

contains

    subroutine check_pro_wronskian(m, n, c, x, tolerance)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree
        real(dp), intent(in) :: c !! positive spheroidal parameter
        real(dp), intent(in) :: x !! prolate radial argument
        real(dp), intent(in) :: tolerance !! absolute tolerance for scaled Wronskian

        real(dp) :: d1
        real(dp) :: d2
        real(dp) :: r1
        real(dp) :: r2
        real(dp) :: scaled_wronskian

        call pro_rad1(m, n, c, x, r1, d1)
        call pro_rad2(m, n, c, x, r2, d2)
        scaled_wronskian = c * (x * x - 1.0_dp) * (r1 * d2 - d1 * r2)
        call check_close(scaled_wronskian, 1.0_dp, tolerance, 'pro Wronskian')
    end subroutine check_pro_wronskian

    subroutine check_obl_wronskian(m, n, c, x, tolerance)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree
        real(dp), intent(in) :: c !! positive spheroidal parameter
        real(dp), intent(in) :: x !! oblate radial argument
        real(dp), intent(in) :: tolerance !! absolute tolerance for scaled Wronskian

        real(dp) :: d1
        real(dp) :: d2
        real(dp) :: r1
        real(dp) :: r2
        real(dp) :: scaled_wronskian

        call obl_rad1(m, n, c, x, r1, d1)
        call obl_rad2(m, n, c, x, r2, d2)
        scaled_wronskian = c * (1.0_dp + x * x) * (r1 * d2 - d1 * r2)
        call check_close(scaled_wronskian, 1.0_dp, tolerance, 'obl Wronskian')
    end subroutine check_obl_wronskian

    subroutine check_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! independent reference value
        real(dp), intent(in) :: tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must hold
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

end program test_special_spheroidal_radial2
