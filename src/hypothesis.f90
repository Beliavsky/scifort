! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_hypothesis
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_beta, only : beta_cdf, beta_sf
    use scifort_descriptive, only : mean, pearson_correlation, rankdata, variance
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_normal, only : normal_sf
    use scifort_student_t, only : t_cdf, t_sf
    implicit none
    private

    type, public :: correlation_test_result
        real(dp) :: statistic !! sample correlation coefficient
        real(dp) :: pvalue !! p-value for the requested alternative
    end type correlation_test_result

    type, public :: mannwhitneyu_result
        real(dp) :: statistic !! U statistic for the first sample
        real(dp) :: pvalue !! asymptotic p-value for the requested alternative
    end type mannwhitneyu_result

    type, public :: ttest_result
        real(dp) :: statistic !! Student or Welch t statistic
        real(dp) :: pvalue !! p-value for the requested alternative
        real(dp) :: df !! degrees of freedom
    end type ttest_result

    public :: mannwhitneyu
    public :: pearsonr
    public :: spearmanr
    public :: ttest_1samp
    public :: ttest_ind
    public :: ttest_rel

contains

    pure function ttest_1samp(x, popmean, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! sample observations; at least two finite values
        real(dp), intent(in) :: popmean !! null-hypothesis population mean, finite
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(ttest_result) :: result_value

        character(len=9) :: selected_alternative
        real(dp) :: difference
        real(dp) :: sample_mean
        real(dp) :: sample_variance
        real(dp) :: standard_error

        selected_alternative = choose_alternative(alternative)
        result_value%df = real(size(x) - 1, dp)
        if (size(x) < 2 .or. .not. ieee_is_finite(popmean) .or. &
                .not. valid_alternative(selected_alternative)) then
            call set_invalid_ttest(result_value)
            return
        end if

        sample_mean = mean(x)
        sample_variance = variance(x)
        if (.not. ieee_is_finite(sample_mean) .or. ieee_is_nan(sample_variance)) then
            call set_invalid_ttest(result_value)
            return
        end if

        difference = sample_mean - popmean
        standard_error = sqrt(sample_variance / real(size(x), dp))
        result_value%statistic = ratio_statistic(difference, standard_error)
        result_value%pvalue = student_t_pvalue(&
            result_value%statistic, result_value%df, selected_alternative)
    end function ttest_1samp

    pure function ttest_ind(x, y, equal_var, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first independent sample, at least two finite values
        real(dp), intent(in) :: y(:) !! second independent sample, at least two finite values
        logical, intent(in), optional :: equal_var !! use pooled variance when true, default true
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(ttest_result) :: result_value

        character(len=9) :: selected_alternative
        logical :: use_equal_variance
        real(dp) :: denominator
        real(dp) :: difference
        real(dp) :: mean_x
        real(dp) :: mean_y
        real(dp) :: pooled_variance
        real(dp) :: standard_error
        real(dp) :: term_x
        real(dp) :: term_y
        real(dp) :: variance_x
        real(dp) :: variance_y

        selected_alternative = choose_alternative(alternative)
        use_equal_variance = .true.
        if (present(equal_var)) use_equal_variance = equal_var

        if (size(x) < 2 .or. size(y) < 2 .or. &
                .not. valid_alternative(selected_alternative)) then
            call set_invalid_ttest(result_value)
            return
        end if

        mean_x = mean(x)
        mean_y = mean(y)
        variance_x = variance(x)
        variance_y = variance(y)
        if (.not. ieee_is_finite(mean_x) .or. .not. ieee_is_finite(mean_y) .or. &
                ieee_is_nan(variance_x) .or. ieee_is_nan(variance_y)) then
            call set_invalid_ttest(result_value)
            return
        end if

        difference = mean_x - mean_y
        if (use_equal_variance) then
            result_value%df = real(size(x) + size(y) - 2, dp)
            pooled_variance = (real(size(x) - 1, dp) * variance_x + &
                real(size(y) - 1, dp) * variance_y) / result_value%df
            standard_error = sqrt(pooled_variance * &
                (1.0_dp / real(size(x), dp) + 1.0_dp / real(size(y), dp)))
        else
            term_x = variance_x / real(size(x), dp)
            term_y = variance_y / real(size(y), dp)
            standard_error = sqrt(term_x + term_y)
            denominator = term_x * term_x / real(size(x) - 1, dp) + &
                term_y * term_y / real(size(y) - 1, dp)
            if (denominator > 0.0_dp) then
                result_value%df = (term_x + term_y) ** 2 / denominator
            else
                result_value%df = real(min(size(x), size(y)) - 1, dp)
            end if
        end if

        result_value%statistic = ratio_statistic(difference, standard_error)
        result_value%pvalue = student_t_pvalue(&
            result_value%statistic, result_value%df, selected_alternative)
    end function ttest_ind

    pure function ttest_rel(x, y, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first paired observations
        real(dp), intent(in) :: y(:) !! second paired observations, same size as x
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(ttest_result) :: result_value

        real(dp), allocatable :: differences(:)

        if (size(x) /= size(y)) then
            call set_invalid_ttest(result_value)
            return
        end if
        allocate(differences(size(x)))
        differences = x - y
        if (present(alternative)) then
            result_value = ttest_1samp(differences, 0.0_dp, alternative)
        else
            result_value = ttest_1samp(differences, 0.0_dp)
        end if
    end function ttest_rel

    pure function mannwhitneyu(x, y, alternative, use_continuity) result(result_value)
        real(dp), intent(in) :: x(:) !! first independent sample, nonempty and without NaN
        real(dp), intent(in) :: y(:) !! second independent sample, nonempty and without NaN
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        logical, intent(in), optional :: use_continuity !! apply 0.5 correction, default true
        type(mannwhitneyu_result) :: result_value

        character(len=9) :: selected_alternative
        integer, allocatable :: tie_counts(:)
        integer :: group
        integer :: i
        integer :: n
        integer :: n_groups
        logical :: continuity
        real(dp), allocatable :: dense_ranks(:)
        real(dp) :: mean_u
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: ranks(:)
        real(dp) :: sigma_u
        real(dp) :: tie_term
        real(dp) :: u1
        real(dp) :: u2
        real(dp) :: u_for_tail
        real(dp) :: variance_u
        real(dp) :: z

        selected_alternative = choose_alternative(alternative)
        continuity = .true.
        if (present(use_continuity)) continuity = use_continuity
        if (size(x) == 0 .or. size(y) == 0 .or. &
                .not. valid_alternative(selected_alternative)) then
            call set_invalid_mannwhitneyu(result_value)
            return
        end if

        n = size(x) + size(y)
        allocate(pooled(n))
        pooled(1:size(x)) = x
        pooled(size(x) + 1:) = y
        ranks = rankdata(pooled)
        if (any(ieee_is_nan(ranks))) then
            call set_invalid_mannwhitneyu(result_value)
            return
        end if

        u1 = sum(ranks(1:size(x))) - &
            0.5_dp * real(size(x), dp) * real(size(x) + 1, dp)
        u2 = real(size(x), dp) * real(size(y), dp) - u1
        result_value%statistic = u1

        dense_ranks = rankdata(pooled, 'dense')
        n_groups = int(maxval(dense_ranks))
        allocate(tie_counts(n_groups))
        tie_counts = 0
        do i = 1, n
            group = int(dense_ranks(i))
            tie_counts(group) = tie_counts(group) + 1
        end do
        tie_term = 0.0_dp
        do group = 1, n_groups
            tie_term = tie_term + real(tie_counts(group), dp) ** 3 - &
                real(tie_counts(group), dp)
        end do

        mean_u = 0.5_dp * real(size(x), dp) * real(size(y), dp)
        if (n <= 1) then
            call set_invalid_mannwhitneyu(result_value)
            return
        end if
        variance_u = real(size(x), dp) * real(size(y), dp) / 12.0_dp * &
            (real(n + 1, dp) - tie_term / (real(n, dp) * real(n - 1, dp)))
        if (variance_u <= 0.0_dp) then
            result_value%pvalue = quiet_nan(0.0_dp)
            return
        end if
        sigma_u = sqrt(variance_u)

        u_for_tail = u1
        select case (selected_alternative)
        case ('greater')
            u_for_tail = u1
        case ('less')
            u_for_tail = u2
        case ('two-sided')
            u_for_tail = max(u1, u2)
        end select
        z = u_for_tail - mean_u
        if (continuity) z = z - 0.5_dp
        z = z / sigma_u
        result_value%pvalue = normal_sf(z)
        if (selected_alternative == 'two-sided') then
            result_value%pvalue = min(1.0_dp, 2.0_dp * result_value%pvalue)
        end if
    end function mannwhitneyu

    pure function pearsonr(x, y, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first observations; at least two finite values
        real(dp), intent(in) :: y(:) !! second observations, same size as x
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(correlation_test_result) :: result_value

        character(len=9) :: selected_alternative
        real(dp) :: beta_shape
        real(dp) :: transformed

        selected_alternative = choose_alternative(alternative)
        if (size(x) /= size(y) .or. size(x) < 2 .or. &
                .not. valid_alternative(selected_alternative)) then
            call set_invalid_correlation(result_value)
            return
        end if

        result_value%statistic = pearson_correlation(x, y)
        if (ieee_is_nan(result_value%statistic)) then
            result_value%pvalue = quiet_nan(0.0_dp)
            return
        end if

        if (size(x) == 2) then
            result_value%pvalue = 1.0_dp
            return
        end if

        beta_shape = 0.5_dp * real(size(x) - 2, dp)
        select case (selected_alternative)
        case ('two-sided')
            transformed = 0.5_dp * (1.0_dp - abs(result_value%statistic))
            result_value%pvalue = min(1.0_dp, 2.0_dp * &
                beta_cdf(transformed, beta_shape, beta_shape))
        case ('less')
            transformed = 0.5_dp * (result_value%statistic + 1.0_dp)
            result_value%pvalue = beta_cdf(transformed, beta_shape, beta_shape)
        case ('greater')
            transformed = 0.5_dp * (result_value%statistic + 1.0_dp)
            result_value%pvalue = beta_sf(transformed, beta_shape, beta_shape)
        end select
    end function pearsonr

    pure function spearmanr(x, y, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first observations; at least three finite values
        real(dp), intent(in) :: y(:) !! second observations, same size as x
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(correlation_test_result) :: result_value

        character(len=9) :: selected_alternative
        real(dp) :: degrees_of_freedom
        real(dp) :: statistic_t
        real(dp), allocatable :: ranks_x(:)
        real(dp), allocatable :: ranks_y(:)

        selected_alternative = choose_alternative(alternative)
        if (size(x) /= size(y) .or. size(x) < 3 .or. &
                .not. valid_alternative(selected_alternative)) then
            call set_invalid_correlation(result_value)
            return
        end if

        ranks_x = rankdata(x)
        ranks_y = rankdata(y)
        result_value%statistic = pearson_correlation(ranks_x, ranks_y)
        if (ieee_is_nan(result_value%statistic)) then
            result_value%pvalue = quiet_nan(0.0_dp)
            return
        end if

        degrees_of_freedom = real(size(x) - 2, dp)
        if (abs(result_value%statistic) >= 1.0_dp) then
            if (result_value%statistic > 0.0_dp) then
                statistic_t = positive_infinity(0.0_dp)
            else
                statistic_t = negative_infinity(0.0_dp)
            end if
        else
            statistic_t = result_value%statistic * sqrt(degrees_of_freedom / &
                ((result_value%statistic + 1.0_dp) * &
                (1.0_dp - result_value%statistic)))
        end if
        result_value%pvalue = student_t_pvalue(&
            statistic_t, degrees_of_freedom, selected_alternative)
    end function spearmanr

    pure function student_t_pvalue(statistic, df, alternative) result(pvalue)
        real(dp), intent(in) :: statistic !! t statistic, possibly infinite
        real(dp), intent(in) :: df !! positive degrees of freedom
        character(len=*), intent(in) :: alternative !! validated alternative hypothesis
        real(dp) :: pvalue

        if (ieee_is_nan(statistic) .or. .not. ieee_is_finite(df) .or. df <= 0.0_dp) then
            pvalue = quiet_nan(0.0_dp)
            return
        end if

        select case (alternative)
        case ('two-sided')
            pvalue = min(1.0_dp, 2.0_dp * t_sf(abs(statistic), df))
        case ('less')
            pvalue = t_cdf(statistic, df)
        case ('greater')
            pvalue = t_sf(statistic, df)
        case default
            pvalue = quiet_nan(0.0_dp)
        end select
    end function student_t_pvalue

    pure function ratio_statistic(numerator, denominator) result(value)
        real(dp), intent(in) :: numerator !! finite numerator of a standardized statistic
        real(dp), intent(in) :: denominator !! nonnegative standard error
        real(dp) :: value

        if (ieee_is_nan(numerator) .or. ieee_is_nan(denominator) .or. denominator < 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else if (denominator > 0.0_dp) then
            value = numerator / denominator
        else if (numerator > 0.0_dp) then
            value = positive_infinity(0.0_dp)
        else if (numerator < 0.0_dp) then
            value = negative_infinity(0.0_dp)
        else
            value = quiet_nan(0.0_dp)
        end if
    end function ratio_statistic

    pure function choose_alternative(alternative) result(selected)
        character(len=*), intent(in), optional :: alternative !! requested alternative or absent
        character(len=9) :: selected

        selected = 'two-sided'
        if (present(alternative)) selected = alternative
    end function choose_alternative

    pure logical function valid_alternative(alternative) result(valid)
        character(len=*), intent(in) :: alternative !! alternative-hypothesis name

        valid = alternative == 'two-sided' .or. alternative == 'less' .or. &
            alternative == 'greater'
    end function valid_alternative

    pure subroutine set_invalid_ttest(result_value)
        type(ttest_result), intent(out) :: result_value !! result populated entirely with NaNs

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%df = quiet_nan(0.0_dp)
    end subroutine set_invalid_ttest

    pure subroutine set_invalid_mannwhitneyu(result_value)
        type(mannwhitneyu_result), intent(out) :: result_value !! result populated with NaNs

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
    end subroutine set_invalid_mannwhitneyu

    pure subroutine set_invalid_correlation(result_value)
        type(correlation_test_result), intent(out) :: result_value !! result populated with NaNs

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
    end subroutine set_invalid_correlation

end module scifort_hypothesis
