program test_association_extended
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use association_extended_reference
    use scifort_kinds, only : dp
    use scifort_stats, only : association_result, brunnermunzel, kendalltau, &
        linregress, linregress_result, page_trend_result, page_trend_test, &
        pointbiserialr, siegelslopes, siegelslopes_result, theilslopes, theilslopes_result
    implicit none

    real(dp), parameter :: kx_tie(5) = [12.0_dp, 2.0_dp, 1.0_dp, 12.0_dp, 2.0_dp]
    real(dp), parameter :: ky_tie(5) = [1.0_dp, 4.0_dp, 7.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: kx(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: ky(5) = [5.0_dp, 1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: kx6(6) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    real(dp), parameter :: ky6(6) = [1.0_dp, 3.0_dp, 2.0_dp, 6.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: lx(6) = [0.2_dp, 1.1_dp, 2.7_dp, 3.0_dp, 5.2_dp, 8.1_dp]
    real(dp), parameter :: ly(6) = [1.0_dp, 1.7_dp, 4.2_dp, 4.0_dp, 8.9_dp, 12.3_dp]
    real(dp), parameter :: rx(7) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    real(dp), parameter :: ry(7) = [1.0_dp, 2.2_dp, 2.7_dp, 10.0_dp, 5.1_dp, 6.2_dp, 7.1_dp]
    real(dp), parameter :: bx(14) = [1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: by(11) = [3.0_dp, 3.0_dp, 4.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, &
        3.0_dp, 1.0_dp, 1.0_dp, 5.0_dp, 4.0_dp]

    type(association_result) :: ar
    type(linregress_result) :: lr
    type(page_trend_result) :: pr
    type(siegelslopes_result) :: sr
    type(theilslopes_result) :: tr
    real(dp) :: page_data(10, 3)
    real(dp) :: page_ranked(10, 3)
    real(dp) :: page4(3, 4)
    integer :: failures

    failures = 0

    ar = kendalltau(kx_tie, ky_tie)
    call check_pair('kendall tied two-sided', ar, kendall_tie_ref(1:2), 3.0e-15_dp, failures)
    ar = kendalltau(kx_tie, ky_tie, alternative='less')
    call check_pair('kendall tied less', ar, kendall_tie_ref(3:4), 3.0e-15_dp, failures)
    ar = kendalltau(kx_tie, ky_tie, alternative='greater')
    call check_pair('kendall tied greater', ar, kendall_tie_ref(5:6), 3.0e-15_dp, failures)
    ar = kendalltau(kx_tie, ky_tie, variant='c')
    call check_pair('kendall tau-c', ar, kendall_tie_ref(7:8), 3.0e-15_dp, failures)

    ar = kendalltau(kx, ky, method='exact')
    call check_pair('kendall exact two-sided', ar, kendall_exact_ref(1:2), 3.0e-15_dp, failures)
    ar = kendalltau(kx, ky, method='exact', alternative='less')
    call check_pair('kendall exact less', ar, kendall_exact_ref(3:4), 3.0e-15_dp, failures)
    ar = kendalltau(kx, ky, method='exact', alternative='greater')
    call check_pair('kendall exact greater', ar, kendall_exact_ref(5:6), 3.0e-15_dp, failures)
    ar = kendalltau(kx, ky, method='asymptotic')
    call check_pair('kendall asymptotic two-sided', ar, kendall_exact_ref(7:8), 3.0e-15_dp, failures)
    ar = kendalltau(kx, ky, method='asymptotic', alternative='less')
    call check_pair('kendall asymptotic less', ar, kendall_exact_ref(9:10), 3.0e-15_dp, failures)
    ar = kendalltau(kx, ky, method='asymptotic', alternative='greater')
    call check_pair('kendall asymptotic greater', ar, kendall_exact_ref(11:12), 3.0e-15_dp, failures)

    ar = kendalltau(kx6, ky6, method='exact')
    call check_pair('kendall exact n6 two-sided', ar, kendall_6_ref(1:2), 3.0e-15_dp, failures)
    ar = kendalltau(kx6, ky6, method='exact', alternative='less')
    call check_pair('kendall exact n6 less', ar, kendall_6_ref(3:4), 3.0e-15_dp, failures)
    ar = kendalltau(kx6, ky6, method='exact', alternative='greater')
    call check_pair('kendall exact n6 greater', ar, kendall_6_ref(5:6), 3.0e-15_dp, failures)

    ar = pointbiserialr([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
        [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp])
    call check_pair('point biserial', ar, pointbiserial_ref, 2.0e-14_dp, failures)

    lr = linregress(lx, ly)
    call check_linregress('linregress two-sided', lr, linregress_ref(1:6), failures)
    lr = linregress(lx, ly, 'less')
    call check_linregress('linregress less', lr, linregress_ref(7:12), failures)
    lr = linregress(lx, ly, 'greater')
    call check_linregress('linregress greater', lr, linregress_ref(13:18), failures)

    tr = theilslopes(ry, rx, 0.9_dp, 'separate')
    call check_theil('theil separate', tr, theil_ref(1:4), failures)
    tr = theilslopes(ry, rx, 0.9_dp, 'joint')
    call check_theil('theil joint', tr, theil_ref(5:8), failures)
    tr = theilslopes([1.0_dp, 1.3_dp, 2.1_dp, 3.2_dp, 4.5_dp, 5.0_dp], &
        [0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 0.95_dp, 'separate')
    call check_theil('theil repeated x', tr, theil_repeat_ref, failures)

    sr = siegelslopes(ry, rx, 'hierarchical')
    call check_siegel('siegel hierarchical', sr, siegel_ref(1:2), failures)
    sr = siegelslopes(ry, rx, 'separate')
    call check_siegel('siegel separate', sr, siegel_ref(3:4), failures)
    sr = siegelslopes([1.0_dp, 1.3_dp, 2.1_dp, 3.2_dp, 4.5_dp, 5.0_dp], &
        [0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 'hierarchical')
    call check_siegel('siegel repeated x', sr, siegel_repeat_ref, failures)

    ar = brunnermunzel(bx, by)
    call check_pair('brunner t two-sided', ar, brunnermunzel_ref(1:2), 2.0e-13_dp, failures)
    ar = brunnermunzel(bx, by, 'less')
    call check_pair('brunner t less', ar, brunnermunzel_ref(3:4), 2.0e-13_dp, failures)
    ar = brunnermunzel(bx, by, 'greater')
    call check_pair('brunner t greater', ar, brunnermunzel_ref(5:6), 2.0e-13_dp, failures)
    ar = brunnermunzel(bx, by, distribution='normal')
    call check_pair('brunner normal two-sided', ar, brunnermunzel_ref(7:8), 2.0e-13_dp, failures)
    ar = brunnermunzel(bx, by, 'less', 'normal')
    call check_pair('brunner normal less', ar, brunnermunzel_ref(9:10), 2.0e-13_dp, failures)
    ar = brunnermunzel(bx, by, 'greater', 'normal')
    call check_pair('brunner normal greater', ar, brunnermunzel_ref(11:12), 2.0e-13_dp, failures)

    call set_page_data(page_data)
    pr = page_trend_test(page_data, method='exact')
    call check_page('page exact', pr, page_ref(1:2), 'exact', failures)
    pr = page_trend_test(page_data, method='asymptotic')
    call check_page('page asymptotic', pr, page_ref(3:4), 'asymptotic', failures)
    call set_page_ranked(page_ranked)
    pr = page_trend_test(page_ranked, ranked=.true., predicted_ranks=[3, 1, 2], method='exact')
    call check_page('page ranked reorder', pr, page_ref(5:6), 'exact', failures)
    call set_page4(page4)
    pr = page_trend_test(page4, method='exact')
    call check_page('page four exact', pr, page4_ref(1:2), 'exact', failures)
    pr = page_trend_test(page4, method='asymptotic')
    call check_page('page four asymptotic', pr, page4_ref(3:4), 'asymptotic', failures)

    ar = kendalltau(kx_tie, ky_tie, method='exact')
    call check_true('kendall exact tied invalid', ieee_is_nan(ar%pvalue), failures)
    lr = linregress([1.0_dp, 1.0_dp, 1.0_dp], [2.0_dp, 3.0_dp, 4.0_dp])
    call check_true('linregress identical x invalid', ieee_is_nan(lr%slope), failures)
    pr = page_trend_test(page4, predicted_ranks=[1, 1, 2, 4])
    call check_true('page invalid predicted ranks', ieee_is_nan(pr%statistic), failures)

    if (failures /= 0) then
        print *, 'association extended failures:', failures
        error stop 1
    end if
    print *, 'association extended tests passed'

contains

    subroutine set_page_data(data)
        real(dp), intent(out) :: data(10, 3) !! unranked SciPy Page test example

        data(:, 1) = [3.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
        data(:, 2) = [4.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 2.0_dp, 4.0_dp, 4.0_dp, 3.0_dp]
        data(:, 3) = [3.0_dp, 4.0_dp, 5.0_dp, 2.0_dp, 2.0_dp, 5.0_dp, 4.0_dp, 4.0_dp, 5.0_dp, 4.0_dp]
    end subroutine set_page_data

    subroutine set_page_ranked(data)
        real(dp), intent(out) :: data(10, 3) !! ranked/reordered SciPy Page test example

        data(:, 1) = [1.5_dp, 3.0_dp, 3.0_dp, 2.0_dp, 1.5_dp, 3.0_dp, 3.0_dp, 2.5_dp, 3.0_dp, 3.0_dp]
        data(:, 2) = [1.5_dp, 1.5_dp, 1.5_dp, 1.0_dp, 1.5_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        data(:, 3) = [3.0_dp, 1.5_dp, 1.5_dp, 3.0_dp, 3.0_dp, 2.0_dp, 2.0_dp, 2.5_dp, 2.0_dp, 2.0_dp]
    end subroutine set_page_ranked

    subroutine set_page4(data)
        real(dp), intent(out) :: data(3, 4) !! four-condition exact/asymptotic Page case

        data(1, :) = [4.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]
        data(2, :) = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp]
        data(3, :) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    end subroutine set_page4

    subroutine check_pair(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        type(association_result), intent(in) :: actual !! statistic and p-value
        real(dp), intent(in) :: expected_value(2) !! SciPy reference statistic and p-value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failure count

        call check_close(trim(name)//' statistic', actual%statistic, expected_value(1), tolerance, count)
        call check_close(trim(name)//' pvalue', actual%pvalue, expected_value(2), tolerance, count)
    end subroutine check_pair

    subroutine check_linregress(name, actual, expected_value, count)
        character(len=*), intent(in) :: name !! check label
        type(linregress_result), intent(in) :: actual !! computed regression result
        real(dp), intent(in) :: expected_value(6) !! SciPy slope through intercept standard error
        integer, intent(inout) :: count !! accumulated failure count

        call check_close(trim(name)//' slope', actual%slope, expected_value(1), 2.0e-15_dp, count)
        call check_close(trim(name)//' intercept', actual%intercept, expected_value(2), 3.0e-15_dp, count)
        call check_close(trim(name)//' rvalue', actual%rvalue, expected_value(3), 2.0e-15_dp, count)
        call check_close(trim(name)//' pvalue', actual%pvalue, expected_value(4), 3.0e-14_dp, count)
        call check_close(trim(name)//' stderr', actual%stderr, expected_value(5), 2.0e-14_dp, count)
        call check_close(trim(name)//' intercept stderr', actual%intercept_stderr, expected_value(6), 2.0e-14_dp, count)
    end subroutine check_linregress

    subroutine check_theil(name, actual, expected_value, count)
        character(len=*), intent(in) :: name !! check label
        type(theilslopes_result), intent(in) :: actual !! computed Theil-Sen result
        real(dp), intent(in) :: expected_value(4) !! SciPy slope, intercept, and confidence limits
        integer, intent(inout) :: count !! accumulated failure count

        call check_close(trim(name)//' slope', actual%slope, expected_value(1), 3.0e-15_dp, count)
        call check_close(trim(name)//' intercept', actual%intercept, expected_value(2), 3.0e-15_dp, count)
        call check_close(trim(name)//' low', actual%low_slope, expected_value(3), 3.0e-15_dp, count)
        call check_close(trim(name)//' high', actual%high_slope, expected_value(4), 3.0e-15_dp, count)
    end subroutine check_theil

    subroutine check_siegel(name, actual, expected_value, count)
        character(len=*), intent(in) :: name !! check label
        type(siegelslopes_result), intent(in) :: actual !! computed Siegel result
        real(dp), intent(in) :: expected_value(2) !! SciPy slope and intercept
        integer, intent(inout) :: count !! accumulated failure count

        call check_close(trim(name)//' slope', actual%slope, expected_value(1), 3.0e-15_dp, count)
        call check_close(trim(name)//' intercept', actual%intercept, expected_value(2), 3.0e-15_dp, count)
    end subroutine check_siegel

    subroutine check_page(name, actual, expected_value, expected_method, count)
        character(len=*), intent(in) :: name !! check label
        type(page_trend_result), intent(in) :: actual !! computed Page result
        real(dp), intent(in) :: expected_value(2) !! SciPy statistic and p-value
        character(len=*), intent(in) :: expected_method !! expected selected method
        integer, intent(inout) :: count !! accumulated failure count

        call check_close(trim(name)//' statistic', actual%statistic, expected_value(1), 2.0e-15_dp, count)
        call check_close(trim(name)//' pvalue', actual%pvalue, expected_value(2), 3.0e-13_dp, count)
        if (trim(actual%method) /= trim(expected_method)) then
            print *, 'FAIL ', trim(name), ' method ', trim(actual%method)
            count = count + 1
        end if
    end subroutine check_page

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

    subroutine check_true(name, condition, count)
        character(len=*), intent(in) :: name !! check label
        logical, intent(in) :: condition !! condition expected to hold
        integer, intent(inout) :: count !! accumulated failure count

        if (.not. condition) then
            print *, 'FAIL ', trim(name)
            count = count + 1
        end if
    end subroutine check_true

end program test_association_extended
