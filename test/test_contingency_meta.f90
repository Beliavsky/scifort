program test_contingency_meta
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use contingency_meta_reference
    use scifort_kinds, only : dp
    use scifort_stats, only : association, barnard_exact, binomtest, binomtest_result, &
        boschloo_exact, combine_pvalues, combined_pvalue_result, confidence_interval, &
        exact_2x2_result, false_discovery_control, margins, odds_ratio, odds_ratio_result, &
        relative_risk, relative_risk_result
    implicit none

    integer, parameter :: vaccine(2, 2) = reshape([7, 8, 12, 3], [2, 2])
    integer, parameter :: bosch_table(2, 2) = reshape([74, 43, 31, 32], [2, 2])
    integer, parameter :: small_table(2, 2) = reshape([1, 3, 2, 4], [2, 2])
    integer, parameter :: or_table(2, 2) = reshape([7, 58, 15, 472], [2, 2])
    real(dp), parameter :: ps(5) = [0.01_dp, 0.2_dp, 0.03_dp, 0.5_dp, 0.001_dp]
    real(dp), parameter :: weights(5) = [1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp]

    type(binomtest_result) :: br
    type(combined_pvalue_result) :: cp
    type(confidence_interval) :: ci
    type(exact_2x2_result) :: er
    type(odds_ratio_result) :: orr
    type(relative_risk_result) :: rr
    integer :: assoc_table(4, 2)
    integer, allocatable :: m0(:, :)
    integer, allocatable :: m1(:, :)
    integer :: failures
    integer :: index
    real(dp), allocatable :: adjusted(:)

    failures = 0

    index = 1
    call check_binom('binom two-sided', 'two-sided', binom_ref(index:index + 7), failures)
    index = index + 8
    call check_binom('binom less', 'less', binom_ref(index:index + 7), failures)
    index = index + 8
    call check_binom('binom greater', 'greater', binom_ref(index:index + 7), failures)

    index = 1
    call check_barnard_set(vaccine, .true., barnard_ref(index:index + 5), failures)
    index = index + 6
    call check_barnard_set(vaccine, .false., barnard_ref(index:index + 5), failures)
    index = index + 6
    call check_barnard_set(small_table, .true., barnard_ref(index:index + 5), failures)
    index = index + 6
    call check_barnard_set(small_table, .false., barnard_ref(index:index + 5), failures)

    call check_boschloo_set(bosch_table, boschloo_ref(1:6), failures)
    call check_boschloo_set(small_table, boschloo_ref(7:12), failures)

    orr = odds_ratio(or_table, 'conditional')
    call check_close('conditional odds statistic', orr%statistic, odds_ref(1), 2.0e-13_dp, failures)
    ci = orr%confidence_interval()
    call check_close('conditional odds ci low', ci%low, odds_ref(2), 2.0e-13_dp, failures)
    call check_close('conditional odds ci high', ci%high, odds_ref(3), 2.0e-12_dp, failures)
    ci = orr%confidence_interval(alternative='less')
    call check_close('conditional odds less low', ci%low, odds_ref(4), 0.0_dp, failures)
    call check_close('conditional odds less high', ci%high, odds_ref(5), 2.0e-12_dp, failures)
    ci = orr%confidence_interval(alternative='greater')
    call check_close('conditional odds greater low', ci%low, odds_ref(6), 2.0e-13_dp, failures)
    call check_true('conditional odds greater high inf', .not. ieee_is_finite(ci%high), failures)

    orr = odds_ratio(or_table, 'sample')
    call check_close('sample odds statistic', orr%statistic, odds_ref(8), 2.0e-15_dp, failures)
    ci = orr%confidence_interval()
    call check_close('sample odds ci low', ci%low, odds_ref(9), 3.0e-15_dp, failures)
    call check_close('sample odds ci high', ci%high, odds_ref(10), 2.0e-14_dp, failures)
    ci = orr%confidence_interval(alternative='less')
    call check_close('sample odds less high', ci%high, odds_ref(12), 2.0e-14_dp, failures)
    ci = orr%confidence_interval(alternative='greater')
    call check_close('sample odds greater low', ci%low, odds_ref(13), 2.0e-14_dp, failures)
    call check_true('sample odds greater high inf', .not. ieee_is_finite(ci%high), failures)

    rr = relative_risk(27, 122, 44, 487)
    call check_close('relative risk', rr%relative_risk, relative_risk_ref(1), 2.0e-15_dp, failures)
    ci = rr%confidence_interval()
    call check_close('relative risk ci low', ci%low, relative_risk_ref(2), 2.0e-15_dp, failures)
    call check_close('relative risk ci high', ci%high, relative_risk_ref(3), 2.0e-15_dp, failures)

    assoc_table = reshape([100, 203, 420, 320, 150, 322, 700, 210], [4, 2])
    call check_close('association pearson', association(assoc_table, 'pearson'), association_ref(1), 2.0e-15_dp, failures)
    call check_close('association cramer', association(assoc_table, 'cramer'), association_ref(2), 2.0e-15_dp, failures)
    call check_close('association tschuprow', association(assoc_table, 'tschuprow'), association_ref(3), 2.0e-15_dp, failures)

    call check_combine('fisher', combine_ref(1:2), failures)
    call check_combine('pearson', combine_ref(3:4), failures)
    call check_combine('mudholkar_george', combine_ref(5:6), failures)
    call check_combine('tippett', combine_ref(7:8), failures)
    cp = combine_pvalues(ps, 'stouffer', weights)
    call check_close('stouffer statistic', cp%statistic, combine_ref(9), 3.0e-15_dp, failures)
    call check_close('stouffer pvalue', cp%pvalue, combine_ref(10), 3.0e-18_dp, failures)

    adjusted = false_discovery_control(ps, 'bh')
    call check_vector('FDR BH', adjusted, fdr_bh_ref, 2.0e-16_dp, failures)
    adjusted = false_discovery_control(ps, 'by')
    call check_vector('FDR BY', adjusted, fdr_by_ref, 2.0e-16_dp, failures)
    adjusted = false_discovery_control([0.5_dp, 0.1_dp, 0.1_dp, 0.01_dp], 'bh')
    call check_vector('FDR ties', adjusted, fdr_tie_ref, 2.0e-16_dp, failures)

    call margins(reshape([0, 6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11], [2, 6]), m0, m1)
    call check_true('margins row 1', all(m0(:, 1) == [15, 51]), failures)
    call check_true('margins columns', all(m1(1, :) == [6, 8, 10, 12, 14, 16]), failures)

    br = binomtest(0, 10, 0.0_dp)
    call check_close('binom degenerate pvalue', br%pvalue, 1.0_dp, 0.0_dp, failures)
    er = barnard_exact(reshape([0, 0, 0, 2], [2, 2]))
    call check_true('barnard zero column statistic nan', ieee_is_nan(er%statistic), failures)
    call check_close('barnard zero column pvalue', er%pvalue, 1.0_dp, 0.0_dp, failures)
    er = boschloo_exact(reshape([0, 0, 0, 2], [2, 2]))
    call check_true('boschloo zero column statistic nan', ieee_is_nan(er%statistic), failures)
    call check_true('boschloo zero column pvalue nan', ieee_is_nan(er%pvalue), failures)

    orr = odds_ratio(reshape([0, 1, 1, 0], [2, 2]), 'conditional')
    call check_close('conditional odds support low', orr%statistic, 0.0_dp, 0.0_dp, failures)
    ci = orr%confidence_interval()
    call check_close('conditional odds support low ci', ci%high, 39.0_dp, 5.0e-12_dp, failures)
    orr = odds_ratio(reshape([1, 0, 0, 1], [2, 2]), 'conditional')
    call check_true('conditional odds support high inf', .not. ieee_is_finite(orr%statistic), failures)

    rr = relative_risk(0, 10, 0, 20)
    call check_true('relative risk zero zero nan', ieee_is_nan(rr%relative_risk), failures)
    rr = relative_risk(0, 10, 2, 20)
    ci = rr%confidence_interval()
    call check_close('relative risk zero', rr%relative_risk, 0.0_dp, 0.0_dp, failures)
    call check_true('relative risk zero ci high nan', ieee_is_nan(ci%high), failures)

    cp = combine_pvalues([1.0_dp, 0.2_dp], 'pearson')
    call check_true('pearson endpoint statistic -inf', .not. ieee_is_finite(cp%statistic) .and. cp%statistic < 0.0_dp, failures)
    call check_close('pearson endpoint pvalue', cp%pvalue, 1.0_dp, 0.0_dp, failures)
    cp = combine_pvalues([0.0_dp, 1.0_dp], 'mudholkar_george')
    call check_true('mudholkar endpoint nan', ieee_is_nan(cp%statistic) .and. ieee_is_nan(cp%pvalue), failures)

    if (failures /= 0) then
        print *, 'contingency/meta failures:', failures
        error stop 1
    end if
    print *, 'contingency/meta tests passed'

contains

    subroutine check_binom(name, alternative, expected_value, count)
        character(len=*), intent(in) :: name !! check label
        character(len=*), intent(in) :: alternative !! binomial alternative
        real(dp), intent(in) :: expected_value(8) !! statistic, p-value, and three confidence intervals
        integer, intent(inout) :: count !! accumulated failure count
        type(binomtest_result) :: result_value
        type(confidence_interval) :: interval

        result_value = binomtest(7, 50, 0.1_dp, alternative)
        call check_close(trim(name)//' statistic', result_value%statistic, expected_value(1), 2.0e-15_dp, count)
        call check_close(trim(name)//' pvalue', result_value%pvalue, expected_value(2), 8.0e-16_dp, count)
        interval = result_value%proportion_ci(method='exact')
        call check_close(trim(name)//' exact low', interval%low, expected_value(3), 5.0e-13_dp, count)
        call check_close(trim(name)//' exact high', interval%high, expected_value(4), 5.0e-13_dp, count)
        interval = result_value%proportion_ci(method='wilson')
        call check_close(trim(name)//' wilson low', interval%low, expected_value(5), 3.0e-15_dp, count)
        call check_close(trim(name)//' wilson high', interval%high, expected_value(6), 3.0e-15_dp, count)
        interval = result_value%proportion_ci(method='wilsoncc')
        call check_close(trim(name)//' wilsoncc low', interval%low, expected_value(7), 3.0e-15_dp, count)
        call check_close(trim(name)//' wilsoncc high', interval%high, expected_value(8), 3.0e-15_dp, count)
    end subroutine check_binom

    subroutine check_barnard_set(table, pooled, expected_value, count)
        integer, intent(in) :: table(2, 2) !! 2x2 contingency table
        logical, intent(in) :: pooled !! pooled-Wald flag
        real(dp), intent(in) :: expected_value(6) !! statistic/p-value for three alternatives
        integer, intent(inout) :: count !! accumulated failure count
        type(exact_2x2_result) :: result_value
        character(len=9), parameter :: alternatives(3) = ['two-sided', 'less     ', 'greater  ']
        integer :: i

        do i = 1, 3
            result_value = barnard_exact(table, trim(alternatives(i)), pooled, 32)
            call check_close('barnard statistic', result_value%statistic, expected_value(2 * i - 1), 4.0e-15_dp, count)
            call check_close('barnard pvalue', result_value%pvalue, expected_value(2 * i), 5.0e-13_dp, count)
        end do
    end subroutine check_barnard_set

    subroutine check_boschloo_set(table, expected_value, count)
        integer, intent(in) :: table(2, 2) !! 2x2 contingency table
        real(dp), intent(in) :: expected_value(6) !! statistic/p-value for three alternatives
        integer, intent(inout) :: count !! accumulated failure count
        type(exact_2x2_result) :: result_value
        character(len=9), parameter :: alternatives(3) = ['two-sided', 'less     ', 'greater  ']
        integer :: i

        do i = 1, 3
            result_value = boschloo_exact(table, trim(alternatives(i)), 32)
            call check_close('boschloo statistic', result_value%statistic, expected_value(2 * i - 1), 2.0e-13_dp, count)
            call check_close('boschloo pvalue', result_value%pvalue, expected_value(2 * i), 8.0e-13_dp, count)
        end do
    end subroutine check_boschloo_set

    subroutine check_combine(method, expected_value, count)
        character(len=*), intent(in) :: method !! combination method
        real(dp), intent(in) :: expected_value(2) !! reference statistic and p-value
        integer, intent(inout) :: count !! accumulated failure count
        type(combined_pvalue_result) :: result_value

        result_value = combine_pvalues(ps, method)
        call check_close(trim(method)//' statistic', result_value%statistic, expected_value(1), 3.0e-14_dp, count)
        call check_close(trim(method)//' pvalue', result_value%pvalue, expected_value(2), 3.0e-15_dp, count)
    end subroutine check_combine

    subroutine check_vector(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! computed vector
        real(dp), intent(in) :: expected_value(:) !! SciPy reference vector
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failure count
        integer :: i

        if (size(actual) /= size(expected_value)) then
            print *, 'FAIL ', trim(name), ' size'
            count = count + 1
            return
        end if
        do i = 1, size(actual)
            call check_close(trim(name), actual(i), expected_value(i), tolerance, count)
        end do
    end subroutine check_vector

    subroutine check_close(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected_value !! SciPy reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failure count

        if (ieee_is_nan(actual) .or. abs(actual - expected_value) > tolerance) then
            print *, 'FAIL ', trim(name), actual, expected_value
            count = count + 1
        end if
    end subroutine check_close

    subroutine check_true(name, condition, count)
        character(len=*), intent(in) :: name !! check label
        logical, intent(in) :: condition !! condition expected to be true
        integer, intent(inout) :: count !! accumulated failure count

        if (.not. condition) then
            print *, 'FAIL ', trim(name)
            count = count + 1
        end if
    end subroutine check_true

end program test_contingency_meta
