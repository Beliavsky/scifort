! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_bayes_mvs
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats
    use test_bayes_mvs_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures = 0
    call test_documented_example(failures)
    call test_small_sample_edges(failures)
    call test_large_sample_approximation(failures)
    call test_flattening_and_invalid_inputs(failures)

    if (failures > 0) then
        print '(a,i0)', 'test_bayes_mvs: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_bayes_mvs: PASS'

contains

    subroutine test_documented_example(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: i
        integer :: status
        real(dp) :: intervals(2)
        type(bayes_mvs_result) :: summary
        type(mvsdist_result) :: distributions

        distributions = mvsdist(bayes_example_data_ref, status)
        call check_true('mvsdist documented status', status == mvs_status_success, failures)
        call check_close('mvsdist mean center', distributions%mean%mean(), &
            bayes_example_centers_ref(1), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist variance center', distributions%variance%mean(), &
            bayes_example_centers_ref(2), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist std center', distributions%std_dev%mean(), &
            bayes_example_centers_ref(3), 2.0e-15_dp, 2.0e-15_dp, failures)

        call check_close('mvsdist mean variance', distributions%mean%var(), &
            bayes_example_variances_ref(1), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('mvsdist variance variance', distributions%variance%var(), &
            bayes_example_variances_ref(2), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('mvsdist std variance', distributions%std_dev%var(), &
            bayes_example_variances_ref(3), 3.0e-15_dp, 3.0e-15_dp, failures)

        call check_close('mvsdist mean std', distributions%mean%std(), &
            bayes_example_stds_ref(1), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('mvsdist variance std', distributions%variance%std(), &
            bayes_example_stds_ref(2), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('mvsdist std std', distributions%std_dev%std(), &
            bayes_example_stds_ref(3), 3.0e-15_dp, 3.0e-15_dp, failures)

        intervals = distributions%mean%interval(0.9_dp)
        do i = 1, 2
            call check_close('mvsdist mean interval90', intervals(i), &
                bayes_example_interval90_ref(i), 5.0e-13_dp, 5.0e-14_dp, failures)
        end do
        intervals = distributions%variance%interval(0.9_dp)
        do i = 1, 2
            call check_close('mvsdist variance interval90', intervals(i), &
                bayes_example_interval90_ref(i + 2), 5.0e-13_dp, 5.0e-14_dp, failures)
        end do
        intervals = distributions%std_dev%interval(0.9_dp)
        do i = 1, 2
            call check_close('mvsdist std interval90', intervals(i), &
                bayes_example_interval90_ref(i + 4), 5.0e-13_dp, 5.0e-14_dp, failures)
        end do

        call check_close('mvsdist mean ppf25', distributions%mean%ppf(0.25_dp), &
            bayes_example_ppf25_ref(1), 2.0e-13_dp, 2.0e-14_dp, failures)
        call check_close('mvsdist variance ppf25', distributions%variance%ppf(0.25_dp), &
            bayes_example_ppf25_ref(2), 2.0e-13_dp, 2.0e-14_dp, failures)
        call check_close('mvsdist std ppf25', distributions%std_dev%ppf(0.25_dp), &
            bayes_example_ppf25_ref(3), 2.0e-13_dp, 2.0e-14_dp, failures)
        call check_close('mvsdist mean cdf', distributions%mean%cdf(1.0_dp), &
            bayes_example_cdf1_ref(1), 2.0e-13_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist variance cdf', distributions%variance%cdf(1.0_dp), &
            bayes_example_cdf1_ref(2), 2.0e-13_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist std cdf', distributions%std_dev%cdf(1.0_dp), &
            bayes_example_cdf1_ref(3), 2.0e-13_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist isf/ppf mean', distributions%mean%isf(0.75_dp), &
            bayes_example_ppf25_ref(1), 2.0e-13_dp, 2.0e-14_dp, failures)
        call check_close('mvsdist sf/cdf mean', distributions%mean%sf(1.0_dp), &
            1.0_dp - bayes_example_cdf1_ref(1), 2.0e-13_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist logcdf mean', exp(distributions%mean%logcdf(1.0_dp)), &
            bayes_example_cdf1_ref(1), 2.0e-13_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist logpdf/pdf mean', exp(distributions%mean%logpdf(9.0_dp)), &
            distributions%mean%pdf(9.0_dp), 2.0e-14_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist logsf/sf mean', exp(distributions%mean%logsf(9.0_dp)), &
            distributions%mean%sf(9.0_dp), 2.0e-14_dp, 2.0e-15_dp, failures)

        summary = bayes_mvs(bayes_example_data_ref, status=status)
        call check_true('bayes_mvs documented status', status == mvs_status_success, failures)
        call check_close('bayes_mvs mean statistic', summary%mean%statistic, &
            bayes_example_centers_ref(1), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('bayes_mvs variance statistic', summary%variance%statistic, &
            bayes_example_centers_ref(2), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('bayes_mvs std statistic', summary%std_dev%statistic, &
            bayes_example_centers_ref(3), 2.0e-15_dp, 2.0e-15_dp, failures)
        do i = 1, 2
            call check_close('bayes_mvs default mean interval', summary%mean%minmax(i), &
                bayes_example_interval90_ref(i), 5.0e-13_dp, 5.0e-14_dp, failures)
            call check_close('bayes_mvs default variance interval', summary%variance%minmax(i), &
                bayes_example_interval90_ref(i + 2), 5.0e-13_dp, 5.0e-14_dp, failures)
            call check_close('bayes_mvs default std interval', summary%std_dev%minmax(i), &
                bayes_example_interval90_ref(i + 4), 5.0e-13_dp, 5.0e-14_dp, failures)
        end do

        summary = bayes_mvs(bayes_example_data_ref, 0.95_dp, status)
        do i = 1, 2
            call check_close('bayes_mvs alpha95 mean interval', summary%mean%minmax(i), &
                bayes_example_interval95_ref(i), 5.0e-13_dp, 5.0e-14_dp, failures)
            call check_close('bayes_mvs alpha95 variance interval', summary%variance%minmax(i), &
                bayes_example_interval95_ref(i + 2), 5.0e-13_dp, 5.0e-14_dp, failures)
            call check_close('bayes_mvs alpha95 std interval', summary%std_dev%minmax(i), &
                bayes_example_interval95_ref(i + 4), 5.0e-13_dp, 5.0e-14_dp, failures)
        end do
    end subroutine test_documented_example

    subroutine test_small_sample_edges(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: i
        real(dp) :: limits(2)
        type(mvsdist_result) :: distributions

        distributions = mvsdist([1.0_dp, 2.0_dp])
        call check_true('mvsdist n2 mean infinite', &
            .not. ieee_is_finite(distributions%mean%mean()), failures)
        call check_true('mvsdist n2 mean variance nan', &
            ieee_is_nan(distributions%mean%var()), failures)
        call check_true('mvsdist n2 variance mean infinite', &
            .not. ieee_is_finite(distributions%variance%mean()), failures)
        call check_true('mvsdist n2 std mean infinite', &
            .not. ieee_is_finite(distributions%std_dev%mean()), failures)
        limits = distributions%mean%interval(0.9_dp)
        do i = 1, 2
            call check_close('mvsdist n2 mean interval', limits(i), &
                bayes_n2_interval90_ref(i), 8.0e-13_dp, 8.0e-14_dp, failures)
        end do
        limits = distributions%variance%interval(0.9_dp)
        do i = 1, 2
            call check_close('mvsdist n2 variance interval', limits(i), &
                bayes_n2_interval90_ref(i + 2), 8.0e-12_dp, 8.0e-14_dp, failures)
        end do
        limits = distributions%std_dev%interval(0.9_dp)
        do i = 1, 2
            call check_close('mvsdist n2 std interval', limits(i), &
                bayes_n2_interval90_ref(i + 4), 8.0e-13_dp, 8.0e-14_dp, failures)
        end do

        distributions = mvsdist([1.0_dp, 2.0_dp, 3.0_dp])
        call check_close('mvsdist n3 mean center', distributions%mean%mean(), &
            bayes_n3_finite_ref(1), 0.0_dp, 0.0_dp, failures)
        call check_true('mvsdist n3 mean variance infinite', &
            .not. ieee_is_finite(distributions%mean%var()), failures)
        call check_true('mvsdist n3 variance mean infinite', &
            .not. ieee_is_finite(distributions%variance%mean()), failures)
        call check_close('mvsdist n3 std center', distributions%std_dev%mean(), &
            bayes_n3_finite_ref(2), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_true('mvsdist n3 std variance infinite', &
            .not. ieee_is_finite(distributions%std_dev%var()), failures)
    end subroutine test_small_sample_edges

    subroutine test_large_sample_approximation(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: i
        real(dp) :: data(1001)
        real(dp) :: intervals(2)
        type(mvsdist_result) :: distributions

        do i = 1, size(data)
            data(i) = -2.0_dp + 5.0_dp * real(i - 1, dp) / 1000.0_dp
        end do
        distributions = mvsdist(data)

        call check_close('mvsdist large mean center', distributions%mean%mean(), &
            bayes_large_centers_ref(1), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist large variance center', distributions%variance%mean(), &
            bayes_large_centers_ref(2), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist large std center', distributions%std_dev%mean(), &
            bayes_large_centers_ref(3), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist large mean variance', distributions%mean%var(), &
            bayes_large_variances_ref(1), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist large variance variance', distributions%variance%var(), &
            bayes_large_variances_ref(2), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('mvsdist large std variance', distributions%std_dev%var(), &
            bayes_large_variances_ref(3), 2.0e-15_dp, 2.0e-15_dp, failures)

        intervals = distributions%mean%interval(0.9_dp)
        do i = 1, 2
            call check_close('mvsdist large mean interval', intervals(i), &
                bayes_large_interval90_ref(i), 2.0e-13_dp, 2.0e-14_dp, failures)
        end do
        intervals = distributions%variance%interval(0.9_dp)
        do i = 1, 2
            call check_close('mvsdist large variance interval', intervals(i), &
                bayes_large_interval90_ref(i + 2), 2.0e-13_dp, 2.0e-14_dp, failures)
        end do
        intervals = distributions%std_dev%interval(0.9_dp)
        do i = 1, 2
            call check_close('mvsdist large std interval', intervals(i), &
                bayes_large_interval90_ref(i + 4), 2.0e-13_dp, 2.0e-14_dp, failures)
        end do
        call check_close('mvsdist large mean ppf25', distributions%mean%ppf(0.25_dp), &
            bayes_large_ppf25_ref(1), 2.0e-13_dp, 2.0e-14_dp, failures)
        call check_close('mvsdist large variance ppf25', distributions%variance%ppf(0.25_dp), &
            bayes_large_ppf25_ref(2), 2.0e-13_dp, 2.0e-14_dp, failures)
        call check_close('mvsdist large std ppf25', distributions%std_dev%ppf(0.25_dp), &
            bayes_large_ppf25_ref(3), 2.0e-13_dp, 2.0e-14_dp, failures)
    end subroutine test_large_sample_approximation

    subroutine test_flattening_and_invalid_inputs(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: i
        integer :: status
        real(dp) :: matrix_data(2, 3)
        real(dp) :: nan_value
        type(bayes_mvs_result) :: summary
        type(mvsdist_result) :: distributions

        matrix_data = reshape([6.0_dp, 7.0_dp, 9.0_dp, 8.0_dp, 12.0_dp, 8.0_dp], [2, 3])
        summary = bayes_mvs(matrix_data, 0.9_dp, status)
        call check_true('bayes_mvs rank2 status', status == mvs_status_success, failures)
        call check_close('bayes_mvs rank2 mean center', summary%mean%statistic, &
            bayes_matrix_centers_ref(1), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('bayes_mvs rank2 variance center', summary%variance%statistic, &
            bayes_matrix_centers_ref(2), 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_close('bayes_mvs rank2 std center', summary%std_dev%statistic, &
            bayes_matrix_centers_ref(3), 3.0e-15_dp, 3.0e-15_dp, failures)
        do i = 1, 2
            call check_close('bayes_mvs rank2 mean interval', summary%mean%minmax(i), &
                bayes_matrix_interval90_ref(i), 5.0e-13_dp, 5.0e-14_dp, failures)
            call check_close('bayes_mvs rank2 variance interval', summary%variance%minmax(i), &
                bayes_matrix_interval90_ref(i + 2), 5.0e-13_dp, 5.0e-14_dp, failures)
            call check_close('bayes_mvs rank2 std interval', summary%std_dev%minmax(i), &
                bayes_matrix_interval90_ref(i + 4), 5.0e-13_dp, 5.0e-14_dp, failures)
        end do

        distributions = mvsdist([1.0_dp], status)
        call check_true('mvsdist short status', status == mvs_status_invalid_data, failures)
        call check_true('mvsdist short result nan', ieee_is_nan(distributions%mean%mean()), failures)

        summary = bayes_mvs([1.0_dp, 2.0_dp], 1.0_dp, status)
        call check_true('bayes_mvs invalid alpha status', status == mvs_status_invalid_alpha, failures)
        call check_true('bayes_mvs invalid alpha result', ieee_is_nan(summary%mean%statistic), failures)

        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        summary = bayes_mvs([1.0_dp, 2.0_dp], nan_value, status)
        call check_true('bayes_mvs nan alpha status', status == mvs_status_invalid_alpha, failures)
        call check_true('bayes_mvs nan alpha result', ieee_is_nan(summary%mean%statistic), failures)
        distributions = mvsdist([1.0_dp, 2.0_dp], status)
        call check_true('mvsdist nan ppf', ieee_is_nan(distributions%mean%ppf(nan_value)), failures)
        call check_true('mvsdist nan interval', &
            all(ieee_is_nan(distributions%mean%interval(nan_value))), failures)

        distributions = mvsdist([1.0_dp, 1.0_dp], status)
        call check_true('mvsdist constant status', status == mvs_status_success, failures)
        call check_true('mvsdist constant mean nan', ieee_is_nan(distributions%mean%mean()), failures)
        call check_true('mvsdist constant interval nan', &
            all(ieee_is_nan(distributions%mean%interval(0.9_dp))), failures)

        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        distributions = mvsdist([nan_value, 1.0_dp, 2.0_dp], status)
        call check_true('mvsdist nan status', status == mvs_status_success, failures)
        call check_true('mvsdist nan result', ieee_is_nan(distributions%variance%mean()), failures)
    end subroutine test_flattening_and_invalid_inputs

end program test_bayes_mvs
