program test_special_spheroidal_radial
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_special, only : obl_cv, obl_rad1, obl_rad1_cv, pro_cv, pro_rad1, &
        pro_rad1_cv
    implicit none

    real(dp) :: cv
    real(dp) :: derivative
    real(dp) :: derivative_cv
    real(dp) :: value
    real(dp) :: value_cv

    call pro_rad1(0, 0, 2.0_dp, 1.5_dp, value, derivative)
    call check_close(value, 0.2559589154341927_dp, 3.0e-13_dp, 'pro 0 0 value')
    call check_close(derivative, -1.020328167403485_dp, 3.0e-13_dp, &
        'pro 0 0 derivative')

    call pro_rad1(1, 2, 2.0_dp, 1.5_dp, value, derivative)
    call check_close(value, 0.2708016141903853_dp, 3.0e-13_dp, 'pro 1 2 value')
    call check_close(derivative, 0.2438429619275685_dp, 3.0e-13_dp, &
        'pro 1 2 derivative')

    call pro_rad1(2, 4, 3.0_dp, 1.2_dp, value, derivative)
    call check_close(value, 0.03509616565750883_dp, 5.0e-13_dp, 'pro 2 4 value')
    call check_close(derivative, 0.2187619759796699_dp, 5.0e-13_dp, &
        'pro 2 4 derivative')

    call pro_rad1(0, 3, 10.0_dp, 1.05_dp, value, derivative)
    call check_close(value, 0.09119813981004403_dp, 3.0e-11_dp, &
        'pro near-boundary value')
    call check_close(derivative, -4.721059439815246_dp, 3.0e-11_dp, &
        'pro near-boundary derivative')

    call obl_rad1(0, 0, 2.0_dp, 0.5_dp, value, derivative)
    call check_close(value, 0.5318688411864536_dp, 3.0e-13_dp, 'obl 0 0 value')
    call check_close(derivative, -0.4686740830333379_dp, 3.0e-13_dp, &
        'obl 0 0 derivative')

    call obl_rad1(1, 2, 2.0_dp, 0.5_dp, value, derivative)
    call check_close(value, 0.11790749181902452_dp, 3.0e-13_dp, 'obl 1 2 value')
    call check_close(derivative, 0.24743390845149757_dp, 3.0e-13_dp, &
        'obl 1 2 derivative')

    call obl_rad1(3, 5, 5.0_dp, 2.0_dp, value, derivative)
    call check_close(value, -0.09214297282228393_dp, 5.0e-12_dp, 'obl 3 5 value')
    call check_close(derivative, -0.10649605955466138_dp, 5.0e-12_dp, &
        'obl 3 5 derivative')

    call obl_rad1(0, 0, 10.0_dp, 0.0_dp, value, derivative)
    call check_close(value, 0.10577359301793428_dp, 5.0e-11_dp, &
        'obl origin even value')
    call check_close(derivative, 0.0_dp, 0.0_dp, 'obl origin even derivative')

    call obl_rad1(2, 3, 2.0_dp, 0.0_dp, value, derivative)
    call check_close(value, 0.0_dp, 0.0_dp, 'obl origin odd value')
    call check_close(derivative, 0.06575955499605403_dp, 3.0e-13_dp, &
        'obl origin odd derivative')

    cv = pro_cv(2, 4, 3.0_dp)
    call pro_rad1(2, 4, 3.0_dp, 1.4_dp, value, derivative)
    call pro_rad1_cv(2, 4, 3.0_dp, cv, 1.4_dp, value_cv, derivative_cv)
    call check_close(value_cv, value, 2.0e-13_dp, 'pro cv value agreement')
    call check_close(derivative_cv, derivative, 2.0e-13_dp, &
        'pro cv derivative agreement')

    cv = obl_cv(2, 4, 3.0_dp)
    call obl_rad1(2, 4, 3.0_dp, 0.7_dp, value, derivative)
    call obl_rad1_cv(2, 4, 3.0_dp, cv, 0.7_dp, value_cv, derivative_cv)
    call check_close(value_cv, value, 2.0e-13_dp, 'obl cv value agreement')
    call check_close(derivative_cv, derivative, 2.0e-13_dp, &
        'obl cv derivative agreement')

    call pro_rad1(0, 0, 0.0_dp, 2.0_dp, value, derivative)
    call check_close(value, 1.0_dp, 0.0_dp, 'pro c0 base value')
    call check_close(derivative, 0.0_dp, 0.0_dp, 'pro c0 base derivative')
    call pro_rad1(0, 1, 0.0_dp, 2.0_dp, value, derivative)
    call check_close(value, 0.0_dp, 0.0_dp, 'pro c0 higher value')
    call check_close(derivative, 0.0_dp, 0.0_dp, 'pro c0 higher derivative')

    call obl_rad1(0, 0, 0.0_dp, 0.0_dp, value, derivative)
    call check_close(value, 1.0_dp, 0.0_dp, 'obl c0 base value')
    call check_close(derivative, 0.0_dp, 0.0_dp, 'obl c0 base derivative')

    call pro_rad1(-1, 0, 1.0_dp, 2.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'negative m')
    call pro_rad1(1, 0, 1.0_dp, 2.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'n below m')
    call pro_rad1(0, 0, -1.0_dp, 2.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'negative c')
    call pro_rad1(0, 0, 1.0_dp, 1.0_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'pro endpoint x')
    call obl_rad1(0, 0, 1.0_dp, -0.1_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'negative obl x')

    print *, 'test_special_spheroidal_radial: PASS'

contains

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

end program test_special_spheroidal_radial
