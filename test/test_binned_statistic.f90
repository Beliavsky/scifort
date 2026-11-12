! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_binned_statistic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats
    use test_binned_statistic_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_one_dimensional(failures)
    call test_two_dimensional(failures)
    call test_multi_dimensional(failures)
    call test_invalid_inputs(failures)

    if (failures > 0) then
        print '(a,i0)', 'test_binned_statistic: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_binned_statistic: PASS'

contains

    subroutine test_one_dimensional(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: i
        real(dp) :: nan_value
        type(binned_statistic_result) :: res

        res = binned_statistic(binned_x1_ref, binned_v1_ref, 'sum', bins=2)
        call check_true('binned 1d sum status', res%status == binned_statistic_status_success, failures)
        do i = 1, 2
            call check_close('binned 1d sum', res%statistic(i), binned_sum1_ref(i), &
                2.0e-15_dp, 2.0e-15_dp, failures)
        end do
        call check_true('binned 1d edges', all(res%bin_edges == binned_edges1_ref), failures)
        call check_true('binned 1d numbers', all(res%binnumber == binned_numbers1_ref), failures)

        res = binned_statistic(binned_x1_ref, binned_v1_ref, 'mean', bins=2)
        do i = 1, 2
            call check_close('binned 1d mean', res%statistic(i), binned_mean1_ref(i), &
                2.0e-15_dp, 2.0e-15_dp, failures)
        end do
        res = binned_statistic(binned_x1_ref, binned_v1_ref, 'std', bins=2)
        do i = 1, 2
            call check_close('binned 1d std', res%statistic(i), binned_std1_ref(i), &
                2.0e-15_dp, 2.0e-15_dp, failures)
        end do
        res = binned_statistic(binned_x1_ref, binned_v1_ref, 'median', bins=2)
        call check_true('binned 1d median', all(res%statistic == binned_median1_ref), failures)
        res = binned_statistic(binned_x1_ref, statistic='count', bins=2)
        call check_true('binned 1d count', all(res%statistic == binned_count1_ref), failures)
        res = binned_statistic(binned_x1_ref, binned_v1_ref, 'min', bins=2)
        call check_true('binned 1d min', all(res%statistic == binned_min1_ref), failures)
        res = binned_statistic(binned_x1_ref, binned_v1_ref, 'max', bins=2)
        call check_true('binned 1d max', all(res%statistic == binned_max1_ref), failures)

        res = binned_statistic(binned_xout_ref, binned_vout_ref, 'mean', &
            bin_edges=binned_out_edges_ref)
        call check_true('binned 1d explicit edges', all(res%bin_edges == binned_out_edges_ref), failures)
        call check_true('binned 1d outlier numbers', &
            all(res%binnumber == binned_out_numbers_ref), failures)
        call check_true('binned 1d explicit mean', all(res%statistic == binned_out_mean_ref), failures)

        res = binned_statistic(binned_xout_ref, binned_vout_ref, &
            bin_edges=binned_custom_edges_ref, stat_func=sample_range)
        call check_true('binned 1d custom numbers', &
            all(res%binnumber == binned_custom_numbers_ref), failures)
        call check_true('binned 1d custom statistic', all(res%statistic == binned_custom_ref), failures)

        res = binned_statistic([3.0_dp, 3.0_dp], statistic='count', bins=2)
        call check_true('binned 1d constant edges', all(res%bin_edges == binned_const_edges_ref), failures)
        call check_true('binned 1d constant count', &
            all(res%statistic == binned_const_count_ref), failures)
        call check_true('binned 1d constant numbers', &
            all(res%binnumber == binned_const_numbers_ref), failures)

        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        res = binned_statistic([0.0_dp, 0.0_dp], [nan_value, 1.0_dp], 'min', &
            bins=1, range=[0.0_dp, 1.0_dp])
        call check_close('binned min ignores nan with finite value', res%statistic(1), &
            1.0_dp, 0.0_dp, 0.0_dp, failures)
        res = binned_statistic([0.0_dp, 0.0_dp], [1.0_dp, nan_value], 'max', &
            bins=1, range=[0.0_dp, 1.0_dp])
        call check_true('binned max propagates nan', ieee_is_nan(res%statistic(1)), failures)
    end subroutine test_one_dimensional

    subroutine test_two_dimensional(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: i
        integer :: j
        real(dp) :: limits(2, 2)
        type(binned_statistic_2d_result) :: res

        limits = reshape([0.0_dp, 2.0_dp, 0.0_dp, 2.0_dp], [2, 2])
        res = binned_statistic_2d(binned_x2_ref, binned_y2_ref, binned_v2_ref, &
            'mean', bins=[2, 2], range=limits)
        call check_true('binned 2d status', res%status == binned_statistic_status_success, failures)
        do j = 1, 2
            do i = 1, 2
                call check_close('binned 2d mean', res%statistic(i, j), &
                    binned_2d_mean_ref(i, j), 2.0e-15_dp, 2.0e-15_dp, failures)
            end do
        end do
        call check_true('binned 2d x edges', all(res%x_edge == binned_2d_xedge_ref), failures)
        call check_true('binned 2d y edges', all(res%y_edge == binned_2d_yedge_ref), failures)
        call check_true('binned 2d flat numbers', &
            all(res%binnumber == binned_2d_numbers_ref), failures)
        call check_true('binned 2d expanded numbers', &
            all(res%binnumber_expanded == binned_2d_expanded_ref), failures)

        res = binned_statistic_2d(binned_x2_ref, binned_y2_ref, statistic='count', &
            bins=[2, 2], range=limits)
        call check_true('binned 2d count', all(res%statistic == binned_2d_count_ref), failures)
    end subroutine test_two_dimensional

    subroutine test_multi_dimensional(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: i
        real(dp) :: explicit_edges(3, 3)
        real(dp) :: limits(2, 3)
        type(binned_statistic_dd_result) :: res

        limits = reshape([0.0_dp, 2.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 2.0_dp], [2, 3])
        res = binned_statistic_dd(binned_dd_sample_ref, binned_dd_values_ref, &
            'sum', bins=[2, 2, 2], range=limits)
        call check_true('binned dd status', res%status == binned_statistic_status_success, failures)
        call check_true('binned dd shape', all(res%statistic_shape == [2, 2, 2]), failures)
        do i = 1, size(binned_dd_sum_ref)
            call check_close('binned dd sum', res%statistic(i), binned_dd_sum_ref(i), &
                0.0_dp, 0.0_dp, failures)
        end do
        call check_true('binned dd flat numbers', &
            all(res%binnumber == binned_dd_numbers_ref), failures)
        call check_true('binned dd expanded numbers', &
            all(res%binnumber_expanded == binned_dd_expanded_ref), failures)

        res = binned_statistic_dd(binned_dd_sample_ref, statistic='count', &
            bins=[2, 2, 2], range=limits)
        call check_true('binned dd count', all(res%statistic == binned_dd_count_ref), failures)

        explicit_edges(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp]
        explicit_edges(:, 2) = [0.0_dp, 0.5_dp, 2.0_dp]
        explicit_edges(:, 3) = [0.0_dp, 1.5_dp, 2.0_dp]
        res = binned_statistic_dd(binned_dd_sample_ref, statistic='count', &
            bin_edges=explicit_edges, edge_counts=[3, 3, 3])
        call check_true('binned dd explicit count', &
            all(res%statistic == binned_dd_explicit_count_ref), failures)
        call check_true('binned dd explicit numbers', &
            all(res%binnumber_expanded == binned_dd_explicit_expanded_ref), failures)
    end subroutine test_multi_dimensional

    subroutine test_invalid_inputs(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        type(binned_statistic_result) :: res1
        type(binned_statistic_2d_result) :: res2
        type(binned_statistic_dd_result) :: resd

        res1 = binned_statistic([1.0_dp], [1.0_dp], bins=0)
        call check_true('binned invalid zero bins', &
            res1%status == binned_statistic_status_invalid, failures)
        res1 = binned_statistic([1.0_dp, 2.0_dp], [1.0_dp], bins=2)
        call check_true('binned invalid value length', &
            res1%status == binned_statistic_status_invalid, failures)
        res1 = binned_statistic([1.0_dp], [1.0_dp], 'variance', bins=2)
        call check_true('binned invalid statistic', &
            res1%status == binned_statistic_status_invalid, failures)
        res1 = binned_statistic([1.0_dp], [1.0_dp], &
            bin_edges=[0.0_dp, 1.0_dp, 1.0_dp])
        call check_true('binned invalid repeated edge', &
            res1%status == binned_statistic_status_invalid, failures)

        res2 = binned_statistic_2d([1.0_dp], [1.0_dp, 2.0_dp], statistic='count', bins=[2])
        call check_true('binned 2d invalid lengths', &
            res2%status == binned_statistic_status_invalid, failures)

        resd = binned_statistic_dd(reshape([1.0_dp, 2.0_dp], [1, 2]), &
            statistic='count', bins=[2, 2, 2])
        call check_true('binned dd invalid bin shape', &
            resd%status == binned_statistic_status_invalid, failures)
    end subroutine test_invalid_inputs

    pure function sample_range(values) result(value)
        real(dp), intent(in) :: values(:) !! values in one bin
        real(dp) :: value

        if (size(values) == 0) then
            value = -99.0_dp
        else
            value = maxval(values) - minval(values)
        end if
    end function sample_range

end program test_binned_statistic
