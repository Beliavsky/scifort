! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_sample_distance
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : energy_distance, wasserstein_distance
    implicit none

    real(dp) :: nan_value
    real(dp) :: positive_inf
    real(dp) :: value

    call check_close(wasserstein_distance([0.0_dp, 1.0_dp, 3.0_dp], &
        [5.0_dp, 6.0_dp, 8.0_dp]), 5.0_dp, 2.0e-14_dp, 'wasserstein basic')
    call check_close(energy_distance([0.0_dp, 1.0_dp, 3.0_dp], &
        [5.0_dp, 6.0_dp, 8.0_dp]), 2.7080128015453204_dp, 2.0e-14_dp, &
        'energy basic')

    call check_close(wasserstein_distance([0.0_dp, 1.0_dp], [0.0_dp, 1.0_dp], &
        [3.0_dp, 1.0_dp], [2.0_dp, 2.0_dp]), 0.25_dp, 2.0e-14_dp, &
        'wasserstein weighted')
    call check_close(energy_distance([0.0_dp, 1.0_dp], [0.0_dp, 1.0_dp], &
        [3.0_dp, 1.0_dp], [2.0_dp, 2.0_dp]), 0.3535533905932738_dp, &
        2.0e-14_dp, 'energy weighted')

    call check_close(wasserstein_distance([3.4_dp, 3.9_dp, 7.5_dp, 7.8_dp], &
        [4.5_dp, 1.4_dp], [1.4_dp, 0.9_dp, 3.1_dp, 7.2_dp], &
        [3.2_dp, 3.5_dp]), 4.078133143804785_dp, 3.0e-14_dp, &
        'wasserstein scipy example')
    call check_close(energy_distance([0.7_dp, 7.4_dp, 2.4_dp, 6.8_dp], &
        [1.4_dp, 8.0_dp], [2.1_dp, 4.2_dp, 7.4_dp, 8.0_dp], &
        [7.6_dp, 8.8_dp]), 0.8800334097615822_dp, 3.0e-14_dp, &
        'energy scipy example')

    call check_close(wasserstein_distance([0.0_dp, 0.0_dp, 1.0_dp], &
        [0.0_dp, 2.0_dp]), 2.0_dp / 3.0_dp, 2.0e-14_dp, &
        'wasserstein repeated support')
    call check_close(energy_distance([0.0_dp, 0.0_dp, 1.0_dp], &
        [0.0_dp, 2.0_dp]), 0.7453559924999299_dp, 2.0e-14_dp, &
        'energy repeated support')

    call check_close(wasserstein_distance([0.0_dp, 2.0_dp], [0.0_dp, 2.0_dp], &
        [3.0_dp, 1.0_dp], [6.0_dp, 2.0_dp]), 0.0_dp, 0.0_dp, &
        'weight normalization')

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    value = wasserstein_distance([0.0_dp, nan_value], [1.0_dp, 2.0_dp])
    if (.not. ieee_is_nan(value)) error stop 'wasserstein NaN propagation mismatch'
    value = energy_distance([0.0_dp, nan_value], [1.0_dp, 2.0_dp])
    if (.not. ieee_is_nan(value)) error stop 'energy NaN propagation mismatch'

    positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
    value = wasserstein_distance([0.0_dp, positive_inf], [1.0_dp, 2.0_dp])
    if (ieee_is_finite(value) .or. value < 0.0_dp) error stop 'wasserstein infinity mismatch'
    value = energy_distance([0.0_dp, positive_inf], [1.0_dp, 2.0_dp])
    if (ieee_is_finite(value) .or. value < 0.0_dp) error stop 'energy infinity mismatch'
    value = wasserstein_distance([0.0_dp, positive_inf], [0.0_dp, positive_inf])
    if (.not. ieee_is_nan(value)) error stop 'duplicate infinity mismatch'

    print *, 'sample-distance tests passed'

contains

    subroutine check_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: atol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > atol) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' abs err = ', abs(actual - expected)
            error stop 1
        end if
    end subroutine check_close

end program test_sample_distance
