! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py.

program test_descriptive_scipy_compat
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_stats, only : dp, gmean, hmean, iqr, kurtosis, &
        median_abs_deviation, moment, pmean, sem, skew, variation
    implicit none

    real(dp), parameter :: x(6) = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]
    real(dp), parameter :: weights(6) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    real(dp), parameter :: zero_data(3) = [0.0_dp, 2.0_dp, 4.0_dp]
    real(dp), parameter :: negative_data(3) = [-1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: constant_data(4) = [3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp]
    real(dp), parameter :: zero_mean_data(2) = [-1.0_dp, 1.0_dp]
    real(dp), parameter :: huge_data(2) = [1.0e300_dp, 1.0e300_dp]
    real(dp), parameter :: tiny_data(2) = [1.0e-300_dp, 1.0e-300_dp]
    real(dp), parameter :: short_weights(2) = [1.0_dp, 2.0_dp]
    real(dp), parameter :: zero_weights(6) = 0.0_dp
    real(dp) :: empty(0)

    call check_close('gmean', gmean(x), 3.5635948725613567_dp, 2.0e-15_dp)
    call check_close('weighted gmean', gmean(x, weights), &
        5.564262476983317_dp, 5.0e-15_dp)
    call check_close('hmean', hmean(x), 2.4615384615384617_dp, 2.0e-15_dp)
    call check_close('weighted hmean', hmean(x, weights), &
        3.818181818181818_dp, 3.0e-15_dp)

    call check_close('pmean -2', pmean(x, -2.0_dp), &
        1.9474579822405906_dp, 3.0e-15_dp)
    call check_close('pmean 0', pmean(x, 0.0_dp), &
        3.5635948725613567_dp, 2.0e-15_dp)
    call check_close('pmean 0.5', pmean(x, 0.5_dp), &
        4.449887763691481_dp, 5.0e-15_dp)
    call check_close('pmean 1', pmean(x, 1.0_dp), 5.5_dp, 2.0e-15_dp)
    call check_close('pmean 2', pmean(x, 2.0_dp), &
        7.582875444051551_dp, 8.0e-15_dp)
    call check_close('weighted pmean 2', pmean(x, 2.0_dp, weights), &
        9.613977919080709_dp, 1.0e-14_dp)

    call check_close('sem ddof 0', sem(x, 0), 2.1311186420907373_dp, 3.0e-15_dp)
    call check_close('sem default', sem(x), 2.3345235059857505_dp, 3.0e-15_dp)
    call check_close('sem ddof 2', sem(x, 2), 2.6100766272276377_dp, 4.0e-15_dp)

    call check_close('variation default', variation(x), &
        0.9491187735373228_dp, 2.0e-15_dp)
    call check_close('variation ddof 1', variation(x, 1), &
        1.039707524036058_dp, 3.0e-15_dp)
    call check_close('variation ddof 2', variation(x, 2), &
        1.162428350231311_dp, 3.0e-15_dp)

    call check_close('iqr default', iqr(x), 5.0_dp, 0.0_dp)
    call check_close('iqr 10-90', iqr(x, [10.0_dp, 90.0_dp]), 10.5_dp, 0.0_dp)
    call check_close('iqr reversed range', iqr(x, [75.0_dp, 25.0_dp]), 5.0_dp, 0.0_dp)
    call check_close('iqr full range', iqr(x, [0.0_dp, 100.0_dp]), 15.0_dp, 0.0_dp)
    call check_close('iqr scale', iqr(x, scale=1.3489795003921634_dp), &
        3.706505546264005_dp, 4.0e-15_dp)

    call check_close('mad default', median_abs_deviation(x), 1.5_dp, 0.0_dp)
    call check_close('mad scale 2', median_abs_deviation(x, 2.0_dp), 0.75_dp, 0.0_dp)
    call check_close('mad normal scale', median_abs_deviation(x, 1.482602218505602_dp), &
        1.0117346252941224_dp, 2.0e-15_dp)

    call check_close('moment 0', moment(x, 0), 1.0_dp, 0.0_dp)
    call check_close('moment 1', moment(x), 0.0_dp, 0.0_dp)
    call check_close('moment 2', moment(x, 2), 27.25_dp, 5.0e-14_dp)
    call check_close('moment 3', moment(x, 3), 165.5_dp, 3.0e-13_dp)
    call check_close('moment 5', moment(x, 5), 20803.75_dp, 5.0e-11_dp)
    call check_close('moment center 0', moment(x, 2, 0.0_dp), 57.5_dp, 8.0e-14_dp)
    call check_close('moment center 3', moment(x, 2, 3.0_dp), 33.5_dp, 6.0e-14_dp)

    call check_close('skew biased', skew(x), 1.1634513776447746_dp, 4.0e-15_dp)
    call check_close('skew unbiased', skew(x, .false.), &
        1.593121410241262_dp, 5.0e-15_dp)
    call check_close('kurtosis fisher biased', kurtosis(x), &
        -0.10251662317986687_dp, 6.0e-15_dp)
    call check_close('kurtosis fisher unbiased', kurtosis(x, .true., .false.), &
        2.2009931823920548_dp, 8.0e-15_dp)
    call check_close('kurtosis pearson biased', kurtosis(x, .false.), &
        2.897483376820133_dp, 6.0e-15_dp)
    call check_close('kurtosis pearson unbiased', kurtosis(x, .false., .false.), &
        5.200993182392055_dp, 8.0e-15_dp)

    call check_close('gmean zero', gmean(zero_data), 0.0_dp, 0.0_dp)
    call check_close('hmean zero', hmean(zero_data), 0.0_dp, 0.0_dp)
    call check_close('pmean zero p2', pmean(zero_data, 2.0_dp), &
        2.581988897471611_dp, 4.0e-15_dp)
    call check_true('pmean overflow limit', pmean(huge_data, 2.0_dp) > huge(1.0_dp))
    call check_close('pmean underflow limit', pmean(tiny_data, -2.0_dp), 0.0_dp, 0.0_dp)
    call check_true('variation zero mean infinite', variation(zero_mean_data) > huge(1.0_dp))

    call check_nan('gmean negative', gmean(negative_data))
    call check_nan('hmean negative', hmean(negative_data))
    call check_nan('pmean negative', pmean(negative_data, 2.0_dp))
    call check_nan('gmean short weights', gmean(x, short_weights))
    call check_nan('hmean zero weight total', hmean(x, zero_weights))
    call check_nan('pmean zero weight total', pmean(x, 2.0_dp, zero_weights))
    call check_nan('iqr invalid range', iqr(x, [-1.0_dp, 75.0_dp]))
    call check_nan('skew constant', skew(constant_data))
    call check_nan('kurtosis constant', kurtosis(constant_data))
    call check_nan('empty gmean', gmean(empty))
    call check_nan('empty sem', sem(empty))
    call check_nan('empty iqr', iqr(empty))
    call check_nan('empty mad', median_abs_deviation(empty))
    call check_nan('empty moment', moment(empty, 0))

    print *, 'SciPy-compatible descriptive statistics tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! SciPy reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_nan(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value expected to be NaN

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', name, actual
            error stop 1
        end if
    end subroutine check_nan

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_scipy_compat
