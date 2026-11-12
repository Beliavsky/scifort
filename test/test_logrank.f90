! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_logrank
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_stats, only : censored_data, ecdf, ecdf_result, logrank, logrank_result
    implicit none

    type(censored_data) :: x
    type(censored_data) :: y
    type(ecdf_result) :: er
    type(logrank_result) :: res

    x = censored_data( &
        uncensored=[6.0_dp, 13.0_dp, 21.0_dp, 30.0_dp, 37.0_dp, 38.0_dp, &
            49.0_dp, 50.0_dp, 63.0_dp, 79.0_dp, 86.0_dp, 98.0_dp, 202.0_dp, 219.0_dp], &
        right=[31.0_dp, 47.0_dp, 80.0_dp, 82.0_dp, 82.0_dp, 149.0_dp])
    y = censored_data( &
        uncensored=[10.0_dp, 10.0_dp, 12.0_dp, 13.0_dp, 14.0_dp, 15.0_dp, &
            16.0_dp, 17.0_dp, 18.0_dp, 20.0_dp, 24.0_dp, 24.0_dp, 25.0_dp, &
            28.0_dp, 30.0_dp, 33.0_dp, 35.0_dp, 37.0_dp, 40.0_dp, 40.0_dp, &
            46.0_dp, 48.0_dp, 76.0_dp, 81.0_dp, 82.0_dp, 91.0_dp, 112.0_dp, 181.0_dp], &
        right=[34.0_dp, 40.0_dp, 70.0_dp])

    if (x%n_observations() /= 20) error stop 'censored_data observation count'
    er = ecdf(x)
    call check_scalar(er%sf%probabilities(size(er%sf%probabilities)), 0.0_dp, &
        2.0e-14_dp, 'ecdf censored-data overload')

    res = logrank(x, y)
    call check_scalar(res%statistic, -2.7379909000095908_dp, 3.0e-13_dp, &
        'logrank censored statistic')
    call check_scalar(res%pvalue, 0.006181578637461755_dp, 3.0e-13_dp, &
        'logrank censored pvalue')

    res = logrank(x, y, alternative='less')
    call check_scalar(res%pvalue, 0.0030907893187308777_dp, 3.0e-13_dp, &
        'logrank less pvalue')
    res = logrank(x, y, alternative='greater')
    call check_scalar(res%pvalue, 0.9969092106812691_dp, 3.0e-13_dp, &
        'logrank greater pvalue')

    res = logrank([1.0_dp], [2.0_dp])
    call check_scalar(res%statistic, 1.0_dp, 2.0e-14_dp, 'logrank arrays statistic')
    call check_scalar(res%pvalue, 0.31731050786291415_dp, 3.0e-13_dp, &
        'logrank arrays pvalue')
    res = logrank([2.0_dp], [1.0_dp], alternative='less')
    call check_scalar(res%pvalue, 0.15865525393145707_dp, 3.0e-13_dp, &
        'logrank arrays less')

    res = logrank(censored_data(uncensored=[1.0_dp], right=[2.0_dp]), [3.0_dp])
    if (.not. (res%pvalue >= 0.0_dp .and. res%pvalue <= 1.0_dp)) then
        error stop 'mixed censored/array overload'
    end if
    res = logrank([1.0_dp], censored_data(uncensored=[2.0_dp], right=[3.0_dp]))
    if (.not. (res%pvalue >= 0.0_dp .and. res%pvalue <= 1.0_dp)) then
        error stop 'mixed array/censored overload'
    end if

    call check_random_case( &
        censored_data(uncensored=[0.55_dp, 1.94_dp, 14.07_dp, 4.36_dp], &
            right=[1.26_dp, 1.53_dp, 3.97_dp]), &
        censored_data(uncensored=[0.8_dp, 7.54_dp, 6.25_dp, 0.4_dp, 6.68_dp, 2.19_dp], &
            right=[0.53_dp, 11.79_dp, 4.26_dp]), &
        -0.1761504693855175_dp, 0.8601757319332395_dp, 'random case 0')
    call check_random_case( &
        censored_data(uncensored=[7.12_dp, 3.84_dp, 3.26_dp, 6.44_dp], &
            right=[9.58_dp, 5.02_dp, 2.11_dp, 0.27_dp]), &
        censored_data(uncensored=[1.74_dp, 9.6_dp, 4.35_dp, 13.04_dp, 2.32_dp, &
            2.01_dp, 2.03_dp, 0.04_dp], right=[1.61_dp, 3.15_dp]), &
        -0.9346646317483472_dp, 0.3499611594412728_dp, 'random case 1')

    ! The final tied group must terminate without reading beyond event_times.
    res = logrank([1.0_dp, 2.0_dp], [2.0_dp, 2.0_dp])
    call check_scalar(res%statistic, 1.0_dp, 2.0e-14_dp, 'logrank terminal ties')
    call check_scalar(res%pvalue, 0.31731050786291415_dp, 3.0e-13_dp, &
        'logrank terminal ties pvalue')

    res = logrank([real(dp) ::], [real(dp) ::])
    if (.not. ieee_is_nan(res%statistic) .or. .not. ieee_is_nan(res%pvalue)) then
        error stop 'empty logrank should be NaN'
    end if

    res = logrank([1.0_dp], [1.0_dp])
    if (.not. ieee_is_nan(res%statistic) .or. .not. ieee_is_nan(res%pvalue)) then
        error stop 'zero-variance logrank should be NaN'
    end if

    print *, 'logrank survival-analysis tests passed'

contains

    subroutine check_random_case(a, b, expected_statistic, expected_pvalue, label)
        type(censored_data), intent(in) :: a
        type(censored_data), intent(in) :: b
        real(dp), intent(in) :: expected_statistic
        real(dp), intent(in) :: expected_pvalue
        character(len=*), intent(in) :: label
        type(logrank_result) :: result

        result = logrank(a, b)
        call check_scalar(result%statistic, expected_statistic, 4.0e-13_dp, &
            trim(label) // ' statistic')
        call check_scalar(result%pvalue, expected_pvalue, 4.0e-13_dp, &
            trim(label) // ' pvalue')
    end subroutine check_random_case

    subroutine check_scalar(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_scalar

end program test_logrank
