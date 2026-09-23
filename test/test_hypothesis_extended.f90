! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_hypothesis_extended
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_stats, only : bartlett, bartlett_result, chi2_contingency, &
        chi2_contingency_result, chisquare, dp, fisher_exact, fisher_exact_result, &
        fligner, fligner_result, friedmanchisquare, friedmanchisquare_result, &
        f_oneway, f_oneway_result, kruskal, kruskal_result, levene, levene_result, &
        make_sample_group, power_divergence_named, power_divergence_result, ranksums, &
        ranksums_result, sample_group, wilcoxon, wilcoxon_result
    implicit none

    real(dp), parameter :: x(5) = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 9.0_dp]
    real(dp), parameter :: y(5) = [2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp]
    real(dp), parameter :: observed(6) = [16.0_dp, 18.0_dp, 16.0_dp, 14.0_dp, 12.0_dp, 12.0_dp]
    real(dp), parameter :: expected(6) = [16.0_dp, 16.0_dp, 16.0_dp, 16.0_dp, 16.0_dp, 8.0_dp]
    real(dp), parameter :: paired_x(10) = [125.0_dp, 115.0_dp, 130.0_dp, 140.0_dp, &
        140.0_dp, 115.0_dp, 140.0_dp, 125.0_dp, 140.0_dp, 135.0_dp]
    real(dp), parameter :: paired_y(10) = [110.0_dp, 122.0_dp, 125.0_dp, 120.0_dp, &
        140.0_dp, 124.0_dp, 123.0_dp, 137.0_dp, 135.0_dp, 145.0_dp]

    type(bartlett_result) :: bart
    type(chi2_contingency_result) :: contingency
    type(fisher_exact_result) :: fisher
    type(fligner_result) :: flig
    type(friedmanchisquare_result) :: friedman
    type(f_oneway_result) :: anova
    type(kruskal_result) :: kruskal_value
    type(levene_result) :: lev
    type(power_divergence_result) :: divergence
    type(ranksums_result) :: rank_sum
    type(sample_group) :: groups(3)
    type(wilcoxon_result) :: signed_rank
    integer :: table(2, 2)
    real(dp) :: friedman_data(4, 8)
    real(dp) :: contingency_table(2, 2)

    rank_sum = ranksums(x, y)
    call check_close('ranksums statistic', rank_sum%statistic, -0.5222329678670935_dp, 2.0e-15_dp)
    call check_close('ranksums pvalue', rank_sum%pvalue, 0.6015081344405899_dp, 3.0e-15_dp)
    rank_sum = ranksums(x, y, 'less')
    call check_close('ranksums less', rank_sum%pvalue, 0.30075406722029496_dp, 2.0e-15_dp)

    groups(1) = make_sample_group(x)
    groups(2) = make_sample_group(y)
    kruskal_value = kruskal(groups(1:2))
    call check_close('kruskal statistic', kruskal_value%statistic, 0.2727272727272734_dp, 2.0e-15_dp)
    call check_close('kruskal pvalue', kruskal_value%pvalue, 0.6015081344405895_dp, 3.0e-15_dp)
    groups(1) = make_sample_group([1.0_dp, 1.0_dp, 1.0_dp])
    groups(2) = make_sample_group([2.0_dp, 2.0_dp, 2.0_dp])
    groups(3) = make_sample_group([2.0_dp, 2.0_dp])
    kruskal_value = kruskal(groups)
    call check_close('kruskal ties statistic', kruskal_value%statistic, 7.0_dp, 0.0_dp)
    call check_close('kruskal ties pvalue', kruskal_value%pvalue, 0.0301973834223185_dp, 2.0e-15_dp)

    friedman_data(1, :) = [9.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 7.0_dp, 9.0_dp, 8.0_dp, 6.0_dp]
    friedman_data(2, :) = [7.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 6.0_dp, 7.0_dp, 7.0_dp, 5.0_dp]
    friedman_data(3, :) = [8.0_dp, 7.0_dp, 9.0_dp, 6.0_dp, 8.0_dp, 8.0_dp, 9.0_dp, 7.0_dp]
    friedman_data(4, :) = [6.0_dp, 8.0_dp, 5.0_dp, 9.0_dp, 5.0_dp, 6.0_dp, 6.0_dp, 8.0_dp]
    friedman = friedmanchisquare(friedman_data)
    call check_close('friedman statistic', friedman%statistic, 6.449999999999989_dp, 2.0e-14_dp)
    call check_close('friedman pvalue', friedman%pvalue, 0.09165537466946727_dp, 2.0e-15_dp)

    signed_rank = wilcoxon([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], method='exact')
    call check_close('wilcoxon exact statistic', signed_rank%statistic, 0.0_dp, 0.0_dp)
    call check_close('wilcoxon exact pvalue', signed_rank%pvalue, 0.0625_dp, 0.0_dp)
    signed_rank = wilcoxon(paired_x, paired_y, zero_method='wilcox', &
        correction=.true., method='auto')
    call check_close('wilcoxon auto statistic', signed_rank%statistic, 18.0_dp, 0.0_dp)
    call check_close('wilcoxon auto pvalue', signed_rank%pvalue, 0.6328125_dp, 0.0_dp)
    signed_rank = wilcoxon(paired_x, paired_y, zero_method='pratt', &
        correction=.true., method='asymptotic')
    call check_close('wilcoxon pratt pvalue', signed_rank%pvalue, 0.645818701677725_dp, 3.0e-15_dp)
    call check_close('wilcoxon pratt z', signed_rank%zstatistic, -0.4595786290693487_dp, 2.0e-15_dp)
    signed_rank = wilcoxon(paired_x, paired_y, zero_method='zsplit', &
        correction=.true., method='asymptotic')
    call check_close('wilcoxon zsplit pvalue', signed_rank%pvalue, 0.6462480913943793_dp, 3.0e-15_dp)

    divergence = power_divergence_named(observed, expected, lambda_name='log-likelihood')
    call check_close('g-test statistic', divergence%statistic, 3.3281031458963746_dp, 3.0e-15_dp)
    call check_close('g-test pvalue', divergence%pvalue, 0.6495419288047497_dp, 3.0e-15_dp)
    divergence = power_divergence_named(observed, expected, lambda_name='freeman-tukey')
    call check_close('freeman-tukey statistic', divergence%statistic, 3.2675381819469074_dp, 4.0e-15_dp)
    divergence = power_divergence_named(observed, expected, lambda_name='cressie-read')
    call check_close('cressie-read statistic', divergence%statistic, 3.434704733101811_dp, 4.0e-15_dp)
    divergence = chisquare(observed)
    call check_close('chisquare statistic', divergence%statistic, 2.0_dp, 0.0_dp)
    call check_close('chisquare pvalue', divergence%pvalue, 0.8491450360846096_dp, 3.0e-15_dp)
    divergence = chisquare(observed, ddof=1)
    call check_close('chisquare ddof pvalue', divergence%pvalue, 0.7357588823428847_dp, 3.0e-15_dp)

    table = reshape([6, 1, 2, 4], [2, 2])
    fisher = fisher_exact(table)
    call check_close('fisher odds ratio', fisher%statistic, 12.0_dp, 0.0_dp)
    call check_close('fisher two-sided pvalue', fisher%pvalue, 0.10256410256410256_dp, 1.0e-15_dp)
    fisher = fisher_exact(table, 'less')
    ! Exact lower tail: (7 + 70 + 175 + 140 + 35) / 429.
    ! Allow accumulated rounding in log-gamma terms and the five-term log sum.
    call check_close('fisher less pvalue', fisher%pvalue, 427.0_dp / 429.0_dp, &
        32.0_dp * epsilon(1.0_dp))
    fisher = fisher_exact(table, 'greater')
    call check_close('fisher greater pvalue', fisher%pvalue, 0.08624708624708625_dp, 2.0e-15_dp)

    contingency_table = reshape([12.0_dp, 29.0_dp, 5.0_dp, 2.0_dp], [2, 2])
    contingency = chi2_contingency(contingency_table)
    call check_close('chi2 contingency yates statistic', contingency%statistic, &
        2.9860164364723074_dp, 3.0e-15_dp)
    call check_close('chi2 contingency yates pvalue', contingency%pvalue, &
        0.08398654171499235_dp, 3.0e-15_dp)
    call check_true('chi2 contingency dof', contingency%dof == 1)
    call check_close('chi2 contingency expected 11', contingency%expected_freq(1, 1), &
        14.520833333333334_dp, 3.0e-15_dp)
    contingency = chi2_contingency(contingency_table, correction=.false.)
    call check_close('chi2 contingency statistic', contingency%statistic, &
        4.646430720203109_dp, 5.0e-15_dp)

    groups(1) = make_sample_group([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    groups(2) = make_sample_group([2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp])
    groups(3) = make_sample_group([7.0_dp, 8.0_dp, 9.0_dp])
    anova = f_oneway(groups)
    call check_close('anova statistic', anova%statistic, 18.848484848484837_dp, 2.0e-14_dp)
    call check_close('anova pvalue', anova%pvalue, 0.0009393130181095339_dp, 3.0e-17_dp)
    anova = f_oneway(groups, .false.)
    call check_close('welch anova statistic', anova%statistic, 21.004724409448816_dp, 2.0e-14_dp)
    call check_close('welch anova pvalue', anova%pvalue, 0.003145142232303294_dp, 8.0e-17_dp)

    groups(1) = make_sample_group([8.88_dp, 9.12_dp, 9.04_dp, 8.98_dp, 9.00_dp, &
        9.08_dp, 9.01_dp, 8.85_dp, 9.06_dp, 8.99_dp])
    groups(2) = make_sample_group([8.88_dp, 8.95_dp, 9.29_dp, 9.44_dp, 9.15_dp, &
        9.58_dp, 8.36_dp, 9.18_dp, 8.67_dp, 9.05_dp])
    groups(3) = make_sample_group([8.95_dp, 9.12_dp, 8.95_dp, 8.85_dp, 9.03_dp, &
        8.84_dp, 9.07_dp, 8.98_dp, 8.86_dp, 8.98_dp])
    bart = bartlett(groups)
    call check_close('bartlett statistic', bart%statistic, 22.789434813726768_dp, 2.0e-14_dp)
    call check_close('bartlett pvalue', bart%pvalue, 1.1254782518834628e-5_dp, 1.0e-18_dp)
    lev = levene(groups)
    call check_close('levene median statistic', lev%statistic, 7.584952754501659_dp, 3.0e-15_dp)
    call check_close('levene median pvalue', lev%pvalue, 0.002431505967249677_dp, 8.0e-17_dp)
    flig = fligner(groups)
    call check_close('fligner median statistic', flig%statistic, 10.803687663522238_dp, 3.0e-15_dp)
    call check_close('fligner median pvalue', flig%pvalue, 0.00450826080004775_dp, 8.0e-17_dp)
    lev = levene(groups, 'mean')
    call check_close('levene mean statistic', lev%statistic, 7.905194483442054_dp, 3.0e-15_dp)
    flig = fligner(groups, 'trimmed')
    call check_close('fligner trimmed statistic', flig%statistic, 10.141132565185508_dp, 1.0e-14_dp)

    groups(1) = make_sample_group([1.0_dp, 1.0_dp, 1.0_dp])
    groups(2) = make_sample_group([1.0_dp, 1.0_dp])
    kruskal_value = kruskal(groups(1:2))
    call check_true('kruskal all-identical invalid', ieee_is_nan(kruskal_value%statistic))
    divergence = power_divergence_named(observed, expected, lambda_name='bad')
    call check_true('invalid divergence name', ieee_is_nan(divergence%pvalue))

    print *, 'extended hypothesis tests passed'

contains

    subroutine check_close(name, actual, expected_value, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected_value !! SciPy 1.17.0 reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (ieee_is_nan(actual) .or. abs(actual - expected_value) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected_value
            error stop 1
        end if
    end subroutine check_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to hold

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_hypothesis_extended
