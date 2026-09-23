program test_nonparametric_extended
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use nonparametric_extended_reference
    use scifort_kinds, only : dp
    use scifort_stats, only : anderson_ksamp, anderson_ksamp_result, ansari, &
        epps_singleton_2samp, make_sample_group, median_test, median_test_named, &
        median_test_result, mood, sample_group, significance_result
    implicit none

    real(dp), parameter :: x(5) = [1.0_dp, 2.0_dp, 5.0_dp, 7.0_dp, 9.0_dp]
    real(dp), parameter :: y(6) = [0.0_dp, 3.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp]
    real(dp), parameter :: xt(5) = [1.0_dp, 2.0_dp, 2.0_dp, 5.0_dp, 7.0_dp]
    real(dp), parameter :: yt(6) = [0.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, 8.0_dp]
    real(dp), parameter :: xe(8) = [1.0_dp, 2.0_dp, 2.5_dp, 3.5_dp, &
        5.0_dp, 8.0_dp, 13.0_dp, 21.0_dp]
    real(dp), parameter :: ye(9) = [0.5_dp, 1.5_dp, 2.2_dp, 4.0_dp, 4.5_dp, &
        7.0_dp, 11.0_dp, 18.0_dp, 29.0_dp]
    real(dp), parameter :: ax(10) = [-2.0_dp, -1.2_dp, -0.8_dp, -0.2_dp, 0.1_dp, &
        0.4_dp, 0.9_dp, 1.4_dp, 2.0_dp, 2.5_dp]
    real(dp), parameter :: ay(10) = [-0.3_dp, 0.4_dp, 1.0_dp, 1.5_dp, 1.8_dp, &
        2.2_dp, 2.6_dp, 3.1_dp, 3.7_dp, 4.3_dp]

    type(anderson_ksamp_result) :: ar
    integer :: failures
    type(sample_group) :: groups2(2)
    type(sample_group) :: groups3(3)
    type(median_test_result) :: mr
    type(significance_result) :: sr

    failures = 0

    sr = ansari(x, y)
    call check_pair('ansari exact two-sided', sr, ansari_ref(1:2), 2.0e-15_dp, failures)
    sr = ansari(x, y, 'less')
    call check_pair('ansari exact less', sr, ansari_ref(3:4), 2.0e-15_dp, failures)
    sr = ansari(x, y, 'greater')
    call check_pair('ansari exact greater', sr, ansari_ref(5:6), 2.0e-15_dp, failures)
    sr = ansari(xt, yt)
    call check_pair('ansari ties two-sided', sr, ansari_ref(7:8), 3.0e-15_dp, failures)
    sr = ansari(xt, yt, 'less')
    call check_pair('ansari ties less', sr, ansari_ref(9:10), 3.0e-15_dp, failures)
    sr = ansari(xt, yt, 'greater')
    call check_pair('ansari ties greater', sr, ansari_ref(11:12), 3.0e-15_dp, failures)

    sr = mood(x, y)
    call check_pair('mood no ties two-sided', sr, mood_ref(1:2), 3.0e-15_dp, failures)
    sr = mood(x, y, 'less')
    call check_pair('mood no ties less', sr, mood_ref(3:4), 3.0e-15_dp, failures)
    sr = mood(x, y, 'greater')
    call check_pair('mood no ties greater', sr, mood_ref(5:6), 3.0e-15_dp, failures)
    sr = mood(xt, yt)
    call check_pair('mood ties two-sided', sr, mood_ref(7:8), 3.0e-15_dp, failures)
    sr = mood(xt, yt, 'less')
    call check_pair('mood ties less', sr, mood_ref(9:10), 3.0e-15_dp, failures)
    sr = mood(xt, yt, 'greater')
    call check_pair('mood ties greater', sr, mood_ref(11:12), 3.0e-15_dp, failures)

    sr = epps_singleton_2samp(xe, ye)
    call check_pair('epps default', sr, epps_ref(1:2), 3.0e-12_dp, failures)
    sr = epps_singleton_2samp(xe, ye, [0.3_dp, 0.7_dp, 1.1_dp])
    call check_pair('epps custom t', sr, epps_ref(3:4), 5.0e-9_dp, failures)

    groups2(1) = make_sample_group(ax)
    groups2(2) = make_sample_group(ay)
    ar = anderson_ksamp(groups2, 'midrank')
    call check_anderson('anderson midrank', ar, anderson_ref(1:2), failures)
    ar = anderson_ksamp(groups2, 'right')
    call check_anderson('anderson right', ar, anderson_ref(3:4), failures)
    ar = anderson_ksamp(groups2, 'continuous')
    call check_anderson('anderson continuous', ar, anderson_ref(5:6), failures)
    call check_array('anderson critical', ar%critical_values, anderson_critical_ref, &
        2.0e-15_dp, failures)

    groups3(1) = make_sample_group([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 6.0_dp])
    groups3(2) = make_sample_group([0.0_dp, 2.0_dp, 3.0_dp, 5.0_dp])
    groups3(3) = make_sample_group([1.0_dp, 3.0_dp, 3.0_dp, 7.0_dp, 8.0_dp, 9.0_dp])
    ar = anderson_ksamp(groups3, 'midrank')
    call check_anderson('anderson tied midrank', ar, anderson_tie_ref(1:2), failures)
    ar = anderson_ksamp(groups3, 'right')
    call check_anderson('anderson tied right', ar, anderson_tie_ref(3:4), failures)

    call set_median_groups(groups3)
    mr = median_test(groups3)
    call check_median('median below', mr, median_ref(1:3), median_table_1, failures)
    mr = median_test(groups3, ties='above')
    call check_median('median above', mr, median_ref(4:6), median_table_2, failures)
    mr = median_test(groups3, ties='ignore')
    call check_median('median ignore', mr, median_ref(7:9), median_table_3, failures)
    mr = median_test_named(groups3, lambda_name='log-likelihood')
    call check_close('median log statistic', mr%statistic, median_ref(10), 3.0e-15_dp, failures)
    call check_close('median log pvalue', mr%pvalue, median_ref(11), 3.0e-15_dp, failures)

    sr = epps_singleton_2samp([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], ye)
    call check_true('epps short invalid', ieee_is_nan(sr%statistic), failures)
    sr = mood([1.0_dp], [2.0_dp], 'bad')
    call check_true('mood bad alternative', ieee_is_nan(sr%statistic), failures)
    ar = anderson_ksamp(groups2, 'bad')
    call check_true('anderson bad variant', ieee_is_nan(ar%statistic), failures)

    if (failures /= 0) then
        print *, 'nonparametric extended failures:', failures
        error stop 1
    end if
    print *, 'nonparametric extended tests passed'

contains

    subroutine set_median_groups(groups)
        type(sample_group), intent(out) :: groups(3) !! groups from the SciPy median-test example

        groups(1) = make_sample_group([10.0_dp, 14.0_dp, 14.0_dp, 18.0_dp, 20.0_dp, &
            22.0_dp, 24.0_dp, 25.0_dp, 31.0_dp, 31.0_dp, 32.0_dp, 39.0_dp, &
            43.0_dp, 43.0_dp, 48.0_dp, 49.0_dp])
        groups(2) = make_sample_group([28.0_dp, 30.0_dp, 31.0_dp, 33.0_dp, 34.0_dp, &
            35.0_dp, 36.0_dp, 40.0_dp, 44.0_dp, 55.0_dp, 57.0_dp, 61.0_dp, &
            91.0_dp, 92.0_dp, 99.0_dp])
        groups(3) = make_sample_group([0.0_dp, 3.0_dp, 9.0_dp, 22.0_dp, 23.0_dp, &
            25.0_dp, 25.0_dp, 33.0_dp, 34.0_dp, 34.0_dp, 40.0_dp, 45.0_dp, &
            46.0_dp, 48.0_dp, 62.0_dp, 67.0_dp, 84.0_dp])
    end subroutine set_median_groups

    subroutine check_pair(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        type(significance_result), intent(in) :: actual !! computed statistic and p-value
        real(dp), intent(in) :: expected_value(2) !! SciPy statistic and p-value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failure count

        call check_close(trim(name)//' statistic', actual%statistic, expected_value(1), &
            tolerance, count)
        call check_close(trim(name)//' pvalue', actual%pvalue, expected_value(2), &
            tolerance, count)
    end subroutine check_pair

    subroutine check_anderson(name, actual, expected_value, count)
        character(len=*), intent(in) :: name !! check label
        type(anderson_ksamp_result), intent(in) :: actual !! computed Anderson result
        real(dp), intent(in) :: expected_value(2) !! SciPy statistic and p-value
        integer, intent(inout) :: count !! accumulated failure count

        call check_close(trim(name)//' statistic', actual%statistic, expected_value(1), &
            5.0e-14_dp, count)
        call check_close(trim(name)//' pvalue', actual%pvalue, expected_value(2), &
            3.0e-13_dp, count)
    end subroutine check_anderson

    subroutine check_median(name, actual, expected_value, expected_table, count)
        character(len=*), intent(in) :: name !! check label
        type(median_test_result), intent(in) :: actual !! computed median-test result
        real(dp), intent(in) :: expected_value(3) !! statistic, p-value, and median
        integer, intent(in) :: expected_table(2, 3) !! expected contingency table
        integer, intent(inout) :: count !! accumulated failure count

        call check_close(trim(name)//' statistic', actual%statistic, expected_value(1), &
            3.0e-14_dp, count)
        call check_close(trim(name)//' pvalue', actual%pvalue, expected_value(2), &
            3.0e-14_dp, count)
        call check_close(trim(name)//' median', actual%median, expected_value(3), 0.0_dp, count)
        if (any(actual%table /= expected_table)) then
            print *, 'FAIL ', trim(name), ' table'
            count = count + 1
        end if
    end subroutine check_median

    subroutine check_close(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected_value !! SciPy 1.17 reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failure count

        if (ieee_is_nan(actual) .or. abs(actual - expected_value) > tolerance) then
            print *, 'FAIL ', trim(name), actual, expected_value
            count = count + 1
        end if
    end subroutine check_close

    subroutine check_array(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! computed array
        real(dp), intent(in) :: expected_value(:) !! SciPy 1.17 reference array
        real(dp), intent(in) :: tolerance !! maximum absolute element error
        integer, intent(inout) :: count !! accumulated failure count

        if (size(actual) /= size(expected_value)) then
            print *, 'FAIL ', trim(name), ' size'
            count = count + 1
        else if (any(abs(actual - expected_value) > tolerance)) then
            print *, 'FAIL ', trim(name)
            count = count + 1
        end if
    end subroutine check_array

    subroutine check_true(name, condition, count)
        character(len=*), intent(in) :: name !! check label
        logical, intent(in) :: condition !! condition expected to hold
        integer, intent(inout) :: count !! accumulated failure count

        if (.not. condition) then
            print *, 'FAIL ', trim(name)
            count = count + 1
        end if
    end subroutine check_true

end program test_nonparametric_extended
