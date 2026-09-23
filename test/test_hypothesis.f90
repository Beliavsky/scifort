! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_hypothesis
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_stats, only : correlation_test_result, dp, mannwhitneyu, &
        mannwhitneyu_result, pearsonr, spearmanr, ttest_1samp, ttest_ind, &
        ttest_rel, ttest_result
    implicit none

    real(dp), parameter :: x(6) = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]
    real(dp), parameter :: y(6) = [3.0_dp, -1.0_dp, 5.0_dp, 4.0_dp, 9.0_dp, 12.0_dp]
    real(dp), parameter :: constant_one(3) = [1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: constant_two(3) = [2.0_dp, 2.0_dp, 2.0_dp]

    type(correlation_test_result) :: correlation_result
    type(mannwhitneyu_result) :: mann_whitney_result
    type(ttest_result) :: t_result

    t_result = ttest_1samp(x, 3.0_dp)
    call check_close('one-sample statistic', t_result%statistic, &
        1.0708823421952984_dp, 5.0e-15_dp)
    call check_close('one-sample pvalue', t_result%pvalue, &
        0.3331688866712901_dp, 2.0e-14_dp)
    call check_close('one-sample df', t_result%df, 5.0_dp, 0.0_dp)

    t_result = ttest_1samp(x, 3.0_dp, 'less')
    call check_close('one-sample less pvalue', t_result%pvalue, &
        0.833415556664355_dp, 2.0e-14_dp)
    t_result = ttest_1samp(x, 3.0_dp, 'greater')
    call check_close('one-sample greater pvalue', t_result%pvalue, &
        0.16658444333564504_dp, 2.0e-14_dp)

    t_result = ttest_ind(x, y)
    call check_close('pooled statistic', t_result%statistic, &
        0.05567597337415254_dp, 5.0e-15_dp)
    call check_close('pooled pvalue', t_result%pvalue, &
        0.9566966323293252_dp, 3.0e-14_dp)
    call check_close('pooled df', t_result%df, 10.0_dp, 0.0_dp)

    t_result = ttest_ind(x, y, .false.)
    call check_close('welch statistic', t_result%statistic, &
        0.05567597337415254_dp, 5.0e-15_dp)
    call check_close('welch pvalue', t_result%pvalue, &
        0.9567471281556149_dp, 3.0e-14_dp)
    call check_close('welch df', t_result%df, &
        9.552789170096602_dp, 2.0e-14_dp)

    t_result = ttest_rel(x, y)
    call check_close('paired statistic', t_result%statistic, &
        0.14648968382726188_dp, 7.0e-15_dp)
    call check_close('paired pvalue', t_result%pvalue, &
        0.889257952096569_dp, 3.0e-14_dp)
    call check_close('paired df', t_result%df, 5.0_dp, 0.0_dp)

    mann_whitney_result = mannwhitneyu(x, y)
    call check_close('mann-whitney statistic', mann_whitney_result%statistic, &
        15.5_dp, 0.0_dp)
    call check_close('mann-whitney pvalue', mann_whitney_result%pvalue, &
        0.747920927964895_dp, 3.0e-14_dp)
    mann_whitney_result = mannwhitneyu(x, y, 'less')
    call check_close('mann-whitney less pvalue', mann_whitney_result%pvalue, &
        0.3739604639824475_dp, 2.0e-14_dp)
    mann_whitney_result = mannwhitneyu(x, y, 'greater')
    call check_close('mann-whitney greater pvalue', mann_whitney_result%pvalue, &
        0.6851229512270187_dp, 2.0e-14_dp)
    mann_whitney_result = mannwhitneyu(x, y, use_continuity=.false.)
    call check_close('mann-whitney no continuity', mann_whitney_result%pvalue, &
        0.6878845906625766_dp, 3.0e-14_dp)

    correlation_result = pearsonr(x, y)
    call check_close('pearson statistic', correlation_result%statistic, &
        0.8763065299100673_dp, 3.0e-15_dp)
    call check_close('pearson pvalue', correlation_result%pvalue, &
        0.02200385215791128_dp, 2.0e-14_dp)
    correlation_result = pearsonr(x, y, 'less')
    call check_close('pearson less pvalue', correlation_result%pvalue, &
        0.9889980739210443_dp, 2.0e-14_dp)
    correlation_result = pearsonr(x, y, 'greater')
    call check_close('pearson greater pvalue', correlation_result%pvalue, &
        0.01100192607895564_dp, 2.0e-14_dp)

    correlation_result = spearmanr(x, y)
    call check_close('spearman statistic', correlation_result%statistic, &
        0.8116794499134278_dp, 3.0e-15_dp)
    call check_close('spearman pvalue', correlation_result%pvalue, &
        0.049857585101340404_dp, 3.0e-14_dp)
    correlation_result = spearmanr(x, y, 'greater')
    call check_close('spearman greater pvalue', correlation_result%pvalue, &
        0.024928792550670202_dp, 2.0e-14_dp)

    t_result = ttest_1samp(constant_one, 1.0_dp)
    call check_true('constant equal mean statistic nan', ieee_is_nan(t_result%statistic))
    call check_true('constant equal mean pvalue nan', ieee_is_nan(t_result%pvalue))
    t_result = ttest_1samp(constant_one, 0.0_dp)
    call check_true('constant shifted mean statistic inf', t_result%statistic > huge(1.0_dp))
    call check_close('constant shifted mean pvalue', t_result%pvalue, 0.0_dp, 0.0_dp)

    t_result = ttest_ind(constant_one, constant_two, .false.)
    call check_true('constant welch statistic negative inf', t_result%statistic < -huge(1.0_dp))
    call check_close('constant welch pvalue', t_result%pvalue, 0.0_dp, 0.0_dp)
    call check_close('constant welch df', t_result%df, 2.0_dp, 0.0_dp)

    correlation_result = pearsonr([1.0_dp, 2.0_dp], [2.0_dp, 4.0_dp])
    call check_close('pearson n=2 statistic', correlation_result%statistic, 1.0_dp, 0.0_dp)
    call check_close('pearson n=2 pvalue', correlation_result%pvalue, 1.0_dp, 0.0_dp)

    correlation_result = pearsonr(constant_one, constant_two)
    call check_true('pearson constant statistic nan', ieee_is_nan(correlation_result%statistic))
    call check_true('pearson constant pvalue nan', ieee_is_nan(correlation_result%pvalue))

    t_result = ttest_1samp(x, 0.0_dp, 'bad')
    call check_true('invalid t alternative', ieee_is_nan(t_result%pvalue))
    correlation_result = pearsonr(x, y(1:5))
    call check_true('mismatched pearson sizes', ieee_is_nan(correlation_result%statistic))
    correlation_result = spearmanr(x(1:2), y(1:2))
    call check_true('spearman too short', ieee_is_nan(correlation_result%pvalue))

    print *, 'hypothesis tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_hypothesis
