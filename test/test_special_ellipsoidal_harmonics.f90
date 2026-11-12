program test_special_ellipsoidal_harmonics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_special, only : ellip_harm, ellip_harm_2, ellip_normal
    implicit none

    real(dp) :: inf_value
    real(dp) :: nan_value
    real(dp) :: values(3)
    real(dp) :: x(3)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    inf_value = ieee_value(0.0_dp, ieee_positive_inf)

    ! SciPy 1.17.0 regression references spanning all four Lame classes.
    call check_close(ellip_harm(5.0_dp, 8.0_dp, 1, 1, 2.5_dp), &
        2.5_dp, 2.0e-14_dp, 'canonical first kind')
    call check_close(ellip_harm(5.0_dp, 8.0_dp, 2, 1, 0.7_dp), &
        -6.176666666666667_dp, 3.0e-14_dp, 'K class')
    call check_close(ellip_harm(5.0_dp, 8.0_dp, 2, 3, 2.4_dp), &
        2.0922714928995227_dp, 4.0e-14_dp, 'L class')
    call check_close(ellip_harm(5.0_dp, 8.0_dp, 2, 4, 3.2_dp), &
        4.789321455070647_dp, 5.0e-14_dp, 'M class')
    call check_close(ellip_harm(5.0_dp, 8.0_dp, 2, 5, -3.2_dp), &
        3.4260180968582192_dp, 5.0e-14_dp, 'N class')
    call check_close(ellip_harm(1.0_dp, 4.0_dp, 3, 2, 3.0_dp), &
        24.794733192202056_dp, 8.0e-14_dp, 'alternate focal parameters')

    call check_close(ellip_harm(5.0_dp, 8.0_dp, 3, 3, 2.5_dp, -1, 1), &
        1.1793413146046718_dp, 8.0e-14_dp, 'signm control')
    call check_close(ellip_harm(5.0_dp, 8.0_dp, 3, 5, 2.5_dp, 1, -1), &
        -0.5594365146293352_dp, 8.0e-14_dp, 'signn control')
    call check_close(ellip_harm(5.0_dp, 8.0_dp, 3, 7, 2.5_dp, -1, -1), &
        3.69754986443726_dp, 8.0e-14_dp, 'combined sign control')

    call check_close(ellip_harm_2(5.0_dp, 8.0_dp, 2, 1, 10.0_dp), &
        0.0010805685338243045_dp, 3.0e-13_dp, 'canonical second kind')
    call check_close(ellip_harm_2(5.0_dp, 8.0_dp, 0, 1, 3.2_dp), &
        0.43338468774846045_dp, 3.0e-13_dp, 'second kind degree zero')
    call check_close(ellip_harm_2(5.0_dp, 8.0_dp, 1, 2, 3.2_dp), &
        0.19169610382753502_dp, 4.0e-13_dp, 'second kind L class')
    call check_close(ellip_harm_2(5.0_dp, 8.0_dp, 2, 5, 4.5_dp), &
        0.017060684283441874_dp, 5.0e-13_dp, 'second kind N class')
    call check_close(ellip_harm_2(1.0_dp, 4.0_dp, 3, 2, 4.0_dp), &
        0.004587767708335918_dp, 6.0e-13_dp, 'second kind alternate focal parameters')
    call check_close(ellip_harm_2(5.0_dp, 8.0_dp, 2, 1, -10.0_dp), &
        -0.0010805685338243045_dp, 3.0e-13_dp, 'negative external coordinate')

    call check_close(ellip_normal(5.0_dp, 8.0_dp, 0, 1), &
        12.56637061435918_dp, 3.0e-13_dp, 'normalization degree zero')
    call check_close(ellip_normal(5.0_dp, 8.0_dp, 2, 1), &
        347.5144466193153_dp, 5.0e-13_dp, 'normalization K class')
    call check_close(ellip_normal(5.0_dp, 8.0_dp, 2, 5), &
        301.5928947446201_dp, 5.0e-13_dp, 'normalization N class')
    call check_close(ellip_normal(1.0_dp, 4.0_dp, 3, 2), &
        7.699843214081672_dp, 5.0e-13_dp, 'normalization alternate focal parameters')
    call check_close(ellip_normal(5.0_dp, 8.0_dp, 3, 7), &
        1723.3879699692575_dp, 8.0e-13_dp, 'canonical normalization')

    x = [0.7_dp, 2.4_dp, 3.2_dp]
    values = ellip_harm(5.0_dp, 8.0_dp, 2, [1, 3, 4], x)
    call check_close(values(1), ellip_harm(5.0_dp, 8.0_dp, 2, 1, x(1)), &
        0.0_dp, 'elemental first value')
    call check_close(values(2), ellip_harm(5.0_dp, 8.0_dp, 2, 3, x(2)), &
        0.0_dp, 'elemental second value')
    call check_close(values(3), ellip_harm(5.0_dp, 8.0_dp, 2, 4, x(3)), &
        0.0_dp, 'elemental third value')

    call require(ieee_is_nan(ellip_harm(5.0_dp, 4.0_dp, 2, 1, 3.0_dp)), &
        'k2 must exceed h2')
    call require(ieee_is_nan(ellip_harm(0.0_dp, 8.0_dp, 2, 1, 3.0_dp)), &
        'h2 must be positive')
    call require(ieee_is_nan(ellip_harm(5.0_dp, 8.0_dp, -1, 1, 3.0_dp)), &
        'negative degree')
    call require(ieee_is_nan(ellip_harm(5.0_dp, 8.0_dp, 2, 0, 3.0_dp)), &
        'low order')
    call require(ieee_is_nan(ellip_harm(5.0_dp, 8.0_dp, 2, 6, 3.0_dp)), &
        'high order')
    call require(ieee_is_nan(ellip_harm(5.0_dp, 8.0_dp, 2, 1, 3.0_dp, 0, 1)), &
        'invalid signm')
    call require(ieee_is_nan(ellip_harm(5.0_dp, 8.0_dp, 2, 1, nan_value)), &
        'NaN coordinate')
    call require(.not. ieee_is_nan(ellip_harm(5.0_dp, 8.0_dp, 2, 1, inf_value)), &
        'infinite first-kind coordinate')
    call require(ieee_is_nan(ellip_harm_2(5.0_dp, 8.0_dp, 2, 1, 2.0_dp)), &
        'second-kind coordinate inside k')
    call require(ieee_is_nan(ellip_harm_2(5.0_dp, 8.0_dp, 2, 1, inf_value)), &
        'second-kind infinite coordinate')
    call require(ieee_is_nan(ellip_normal(5.0_dp, 4.0_dp, 2, 1)), &
        'normalization invalid focal parameters')

    print *, 'special ellipsoidal-harmonic tests passed'

contains

    subroutine check_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! SciPy 1.17.0 reference value
        real(dp), intent(in) :: tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. &
                abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
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

end program test_special_ellipsoidal_harmonics
