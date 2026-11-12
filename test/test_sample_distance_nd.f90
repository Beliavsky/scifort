! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_sample_distance_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : wasserstein_distance_nd
    implicit none

    real(dp) :: positive_inf
    real(dp) :: u3(2, 3)
    real(dp) :: v3(2, 3)
    real(dp) :: uw(3)
    real(dp) :: vw(2)
    real(dp) :: u2(3, 2)
    real(dp) :: v2(2, 2)
    real(dp) :: same(2, 2)
    real(dp), allocatable :: zero_u(:, :)
    real(dp), allocatable :: zero_v(:, :)
    real(dp) :: value

    u3(1, :) = [0.0_dp, 2.0_dp, 3.0_dp]
    u3(2, :) = [1.0_dp, 2.0_dp, 5.0_dp]
    v3(1, :) = [3.0_dp, 2.0_dp, 3.0_dp]
    v3(2, :) = [4.0_dp, 2.0_dp, 5.0_dp]
    call check_close(wasserstein_distance_nd(u3, v3), 3.0_dp, 2.0e-14_dp, &
        'wasserstein nd scipy 3d example')

    u2(1, :) = [0.0_dp, 2.75_dp]
    u2(2, :) = [2.0_dp, 209.3_dp]
    u2(3, :) = [0.0_dp, 0.0_dp]
    v2(1, :) = [0.2_dp, 0.322_dp]
    v2(2, :) = [4.5_dp, 25.1808_dp]
    uw = [0.4_dp, 5.2_dp, 0.114_dp]
    vw = [0.8_dp, 1.5_dp]
    call check_close(wasserstein_distance_nd(u2, v2, uw, vw), &
        174.15840245217169_dp, 5.0e-13_dp, 'wasserstein nd scipy weighted example')

    call check_close(wasserstein_distance_nd([0.0_dp, 1.0_dp, 3.0_dp], &
        [5.0_dp, 6.0_dp, 8.0_dp]), 5.0_dp, 2.0e-14_dp, &
        'wasserstein nd vector dispatch')

    same(1, :) = [-1.0_dp, 2.0_dp]
    same(2, :) = [3.0_dp, 4.0_dp]
    call check_close(wasserstein_distance_nd(same, same), 0.0_dp, 1.0e-15_dp, &
        'wasserstein nd identical distributions')

    allocate(zero_u(2, 0), zero_v(3, 0))
    call check_close(wasserstein_distance_nd(zero_u, zero_v), 0.0_dp, 0.0_dp, &
        'wasserstein nd zero-dimensional observations')

    same(1, :) = [1.0e200_dp, 0.0_dp]
    same(2, :) = [0.0_dp, 0.0_dp]
    value = wasserstein_distance_nd(same(1:1, :), same(2:2, :))
    if (.not. ieee_is_finite(value)) error stop 'wasserstein nd stable norm is not finite'
    call check_close(value / 1.0e200_dp, 1.0_dp, 2.0e-14_dp, &
        'wasserstein nd stable euclidean norm')

    positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
    same(1, :) = [0.0_dp, 0.0_dp]
    same(2, :) = [positive_inf, 1.0_dp]
    value = wasserstein_distance_nd(same, reshape([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], [2, 2]))
    if (ieee_is_finite(value) .or. value < 0.0_dp) then
        error stop 'wasserstein nd one-sided infinity mismatch'
    end if
    value = wasserstein_distance_nd(same, same)
    if (.not. ieee_is_nan(value)) error stop 'wasserstein nd two-sided infinity mismatch'

    print *, 'multidimensional sample-distance tests passed'

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

end program test_sample_distance_nd
