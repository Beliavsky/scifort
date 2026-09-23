! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Additional scipy.stats-style hypothesis tests. The scalar formulas and
! default semantics are checked against SciPy 1.17.0. Multi-sample routines
! use sample_group so groups may have different lengths without padding.

module scifort_hypothesis_extended
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_chi2, only : chi2_sf
    use scifort_descriptive, only : mean, median, rankdata, variance
    use scifort_f_distribution, only : f_sf
    use scifort_hypergeom, only : hypergeom_cdf, hypergeom_pmf, hypergeom_sf
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_normal, only : normal_cdf, normal_ppf, normal_sf
    implicit none
    private

    type, public :: sample_group
        real(dp), allocatable :: values(:) !! observations in one independent sample
    end type sample_group

    type, public :: ranksums_result
        real(dp) :: statistic !! standardized Wilcoxon rank-sum statistic
        real(dp) :: pvalue !! p-value for the requested alternative
    end type ranksums_result

    type, public :: kruskal_result
        real(dp) :: statistic !! tie-corrected Kruskal-Wallis H statistic
        real(dp) :: pvalue !! chi-square approximation p-value
    end type kruskal_result

    type, public :: friedmanchisquare_result
        real(dp) :: statistic !! tie-corrected Friedman statistic
        real(dp) :: pvalue !! chi-square approximation p-value
    end type friedmanchisquare_result

    type, public :: wilcoxon_result
        real(dp) :: statistic !! signed-rank statistic reported by SciPy semantics
        real(dp) :: pvalue !! exact, permutation, or asymptotic p-value
        real(dp) :: zstatistic !! asymptotic z statistic; NaN otherwise
    end type wilcoxon_result

    type, public :: power_divergence_result
        real(dp) :: statistic !! Cressie-Read power-divergence statistic
        real(dp) :: pvalue !! chi-square reference p-value
    end type power_divergence_result

    type, public :: fisher_exact_result
        real(dp) :: statistic !! sample odds ratio for a 2 by 2 table
        real(dp) :: pvalue !! exact hypergeometric p-value
    end type fisher_exact_result

    type, public :: chi2_contingency_result
        real(dp) :: statistic !! contingency-table power-divergence statistic
        real(dp) :: pvalue !! chi-square reference p-value
        integer :: dof !! degrees of freedom
        real(dp), allocatable :: expected_freq(:, :) !! expected frequencies under independence
    end type chi2_contingency_result

    type, public :: f_oneway_result
        real(dp) :: statistic !! classical or Welch one-way ANOVA F statistic
        real(dp) :: pvalue !! upper-tail F p-value
    end type f_oneway_result

    type, public :: bartlett_result
        real(dp) :: statistic !! Bartlett chi-square statistic
        real(dp) :: pvalue !! upper-tail chi-square p-value
    end type bartlett_result

    type, public :: levene_result
        real(dp) :: statistic !! Levene/Brown-Forsythe F statistic
        real(dp) :: pvalue !! upper-tail F p-value
    end type levene_result

    type, public :: fligner_result
        real(dp) :: statistic !! Fligner-Killeen statistic
        real(dp) :: pvalue !! upper-tail chi-square p-value
    end type fligner_result

    public :: bartlett
    public :: chi2_contingency
    public :: chisquare
    public :: contingency_expected_freq
    public :: fisher_exact
    public :: fligner
    public :: friedmanchisquare
    public :: f_oneway
    public :: kruskal
    public :: levene
    public :: make_sample_group
    public :: power_divergence
    public :: power_divergence_named
    public :: ranksums
    public :: wilcoxon

contains

    function make_sample_group(values) result(group)
        real(dp), intent(in) :: values(:) !! observations copied into the returned group
        type(sample_group) :: group

        group%values = values
    end function make_sample_group

    function ranksums(x, y, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first independent sample
        real(dp), intent(in) :: y(:) !! second independent sample
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(ranksums_result) :: result_value

        character(len=9) :: selected_alternative
        integer :: n1
        integer :: n2
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: ranks(:)
        real(dp) :: expected
        real(dp) :: rank_sum
        real(dp) :: se

        selected_alternative = choose_alternative(alternative)
        if (size(x) < 1 .or. size(y) < 1 .or. .not. valid_alternative(selected_alternative)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if

        n1 = size(x)
        n2 = size(y)
        allocate(pooled(n1 + n2))
        pooled(1:n1) = x
        pooled(n1 + 1:) = y
        ranks = rankdata(pooled)
        if (any(ieee_is_nan(ranks))) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if

        rank_sum = sum(ranks(1:n1))
        expected = real(n1 * (n1 + n2 + 1), dp) / 2.0_dp
        se = sqrt(real(n1 * n2 * (n1 + n2 + 1), dp) / 12.0_dp)
        result_value%statistic = (rank_sum - expected) / se
        result_value%pvalue = normal_pvalue(result_value%statistic, selected_alternative)
    end function ranksums

    function kruskal(groups) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples
        type(kruskal_result) :: result_value

        integer :: i
        integer :: first
        integer :: last
        integer :: n_total
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: ranks(:)
        real(dp) :: correction
        real(dp) :: h
        real(dp) :: ssbn
        real(dp) :: tie_sum

        if (.not. valid_groups(groups, 2, .false.)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if

        n_total = total_group_size(groups)
        allocate(pooled(n_total))
        first = 1
        do i = 1, size(groups)
            last = first + size(groups(i)%values) - 1
            pooled(first:last) = groups(i)%values
            first = last + 1
        end do
        ranks = rankdata(pooled)
        if (any(ieee_is_nan(ranks))) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if

        tie_sum = tie_cubic_sum(pooled)
        if (n_total <= 1) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        correction = 1.0_dp - tie_sum / &
            (real(n_total, dp) ** 3 - real(n_total, dp))
        if (correction <= 0.0_dp) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if

        ssbn = 0.0_dp
        first = 1
        do i = 1, size(groups)
            last = first + size(groups(i)%values) - 1
            ssbn = ssbn + sum(ranks(first:last)) ** 2 / real(size(groups(i)%values), dp)
            first = last + 1
        end do
        h = 12.0_dp * ssbn / (real(n_total, dp) * real(n_total + 1, dp)) - &
            3.0_dp * real(n_total + 1, dp)
        h = h / correction
        result_value%statistic = h
        result_value%pvalue = chi2_sf(h, real(size(groups) - 1, dp))
    end function kruskal

    function friedmanchisquare(samples) result(result_value)
        real(dp), intent(in) :: samples(:, :) !! rows are repeated samples; columns are subjects/blocks
        type(friedmanchisquare_result) :: result_value

        integer :: i
        integer :: j
        integer :: k
        integer :: n
        real(dp) :: correction
        real(dp), allocatable :: ranks(:, :)
        real(dp) :: ssbn
        real(dp) :: tie_sum

        k = size(samples, 1)
        n = size(samples, 2)
        if (k < 3 .or. n < 1 .or. any(ieee_is_nan(samples))) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if

        allocate(ranks(k, n))
        tie_sum = 0.0_dp
        do j = 1, n
            ranks(:, j) = rankdata(samples(:, j))
            tie_sum = tie_sum + tie_cubic_sum(samples(:, j))
        end do
        correction = 1.0_dp - tie_sum / &
            (real(k * (k * k - 1) * n, dp))
        if (correction <= 0.0_dp) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if

        ssbn = 0.0_dp
        do i = 1, k
            ssbn = ssbn + sum(ranks(i, :)) ** 2
        end do
        result_value%statistic = &
            (12.0_dp * ssbn / real(k * n * (k + 1), dp) - &
            3.0_dp * real(n * (k + 1), dp)) / correction
        result_value%pvalue = chi2_sf(result_value%statistic, real(k - 1, dp))
    end function friedmanchisquare

    function wilcoxon(x, y, zero_method, correction, alternative, method) result(result_value)
        real(dp), intent(in) :: x(:) !! differences, or first paired sample when y is present
        real(dp), intent(in), optional :: y(:) !! optional second paired sample
        character(len=*), intent(in), optional :: zero_method !! wilcox, pratt, or zsplit
        logical, intent(in), optional :: correction !! apply asymptotic continuity correction
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        character(len=*), intent(in), optional :: method !! auto, exact, asymptotic, or permutation
        type(wilcoxon_result) :: result_value

        character(len=11) :: selected_method
        character(len=9) :: selected_alternative
        character(len=6) :: selected_zero
        logical :: apply_correction
        logical :: has_ties
        integer :: rank_count
        integer :: n_zero
        real(dp), allocatable :: d(:)
        real(dp) :: mn
        real(dp) :: r_minus
        real(dp) :: r_plus
        real(dp) :: se
        real(dp) :: z

        result_value%zstatistic = quiet_nan(0.0_dp)
        selected_zero = 'wilcox'
        if (present(zero_method)) selected_zero = zero_method
        selected_alternative = choose_alternative(alternative)
        selected_method = 'auto'
        if (present(method)) selected_method = method
        apply_correction = .false.
        if (present(correction)) apply_correction = correction

        if (.not. valid_zero_method(selected_zero) .or. &
                .not. valid_alternative(selected_alternative) .or. &
                .not. valid_wilcoxon_method(selected_method)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        if (present(y)) then
            if (size(y) /= size(x)) then
                call invalidate_pair(result_value%statistic, result_value%pvalue)
                return
            end if
            allocate(d(size(x)))
            d = x - y
        else
            allocate(d(size(x)))
            d = x
        end if
        if (size(d) == 0 .or. any(ieee_is_nan(d))) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if

        call wilcoxon_components(d, selected_zero, r_plus, r_minus, mn, se, &
            rank_count, n_zero, has_ties)
        if (selected_alternative == 'two-sided') then
            result_value%statistic = min(r_plus, r_minus)
        else
            result_value%statistic = r_plus
        end if

        if (selected_method == 'auto') then
            if (size(d) > 50) then
                selected_method = 'asymptotic'
            else if (.not. has_ties .and. n_zero == 0) then
                selected_method = 'exact'
            else if (size(d) <= 13) then
                selected_method = 'permutation'
            else
                selected_method = 'asymptotic'
            end if
        end if

        select case (trim(selected_method))
        case ('asymptotic')
            if (se <= 0.0_dp .or. ieee_is_nan(se)) then
                result_value%pvalue = quiet_nan(0.0_dp)
                result_value%zstatistic = quiet_nan(0.0_dp)
                return
            end if
            z = (r_plus - mn) / se
            if (apply_correction) z = z - correction_sign(z, selected_alternative) * 0.5_dp / se
            result_value%pvalue = normal_pvalue(z, selected_alternative)
            if (selected_alternative == 'two-sided') then
                result_value%zstatistic = -abs(z)
            else
                result_value%zstatistic = z
            end if
        case ('exact')
            result_value%pvalue = wilcoxon_exact_pvalue(r_plus, rank_count, selected_alternative)
        case ('permutation')
            if (size(d) > 22) then
                result_value%pvalue = quiet_nan(0.0_dp)
            else
                result_value%pvalue = wilcoxon_permutation_pvalue(&
                    d, selected_zero, r_plus, selected_alternative)
            end if
        end select
    end function wilcoxon

    function power_divergence(f_obs, f_exp, ddof, lambda_) result(result_value)
        real(dp), intent(in) :: f_obs(:) !! observed category frequencies
        real(dp), intent(in), optional :: f_exp(:) !! expected frequencies; uniform when absent
        integer, intent(in), optional :: ddof !! delta degrees of freedom, default zero
        real(dp), intent(in), optional :: lambda_ !! Cressie-Read lambda, default one
        type(power_divergence_result) :: result_value

        integer :: delta_df
        integer :: df
        integer :: i
        real(dp), allocatable :: expected(:)
        real(dp) :: lambda_value
        real(dp) :: mean_expected
        real(dp) :: obs_sum
        real(dp) :: exp_sum
        real(dp) :: ratio
        real(dp) :: rtol
        real(dp) :: term

        delta_df = 0
        if (present(ddof)) delta_df = ddof
        lambda_value = 1.0_dp
        if (present(lambda_)) lambda_value = lambda_

        if (size(f_obs) < 1 .or. any(ieee_is_nan(f_obs)) .or. any(f_obs < 0.0_dp) .or. &
                .not. ieee_is_finite(lambda_value)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        allocate(expected(size(f_obs)))
        if (present(f_exp)) then
            if (size(f_exp) /= size(f_obs) .or. any(ieee_is_nan(f_exp)) .or. any(f_exp < 0.0_dp)) then
                call invalidate_pair(result_value%statistic, result_value%pvalue)
                return
            end if
            expected = f_exp
            obs_sum = sum(f_obs)
            exp_sum = sum(expected)
            rtol = sqrt(epsilon(1.0_dp))
            if (min(obs_sum, exp_sum) <= 0.0_dp) then
                if (obs_sum /= exp_sum) then
                    call invalidate_pair(result_value%statistic, result_value%pvalue)
                    return
                end if
            else if (abs(obs_sum - exp_sum) / min(obs_sum, exp_sum) > rtol) then
                call invalidate_pair(result_value%statistic, result_value%pvalue)
                return
            end if
        else
            mean_expected = sum(f_obs) / real(size(f_obs), dp)
            expected = mean_expected
        end if

        result_value%statistic = 0.0_dp
        do i = 1, size(f_obs)
            if (lambda_value == 1.0_dp) then
                if (expected(i) <= 0.0_dp) then
                    if (f_obs(i) == 0.0_dp) then
                        term = quiet_nan(0.0_dp)
                    else
                        term = positive_infinity(0.0_dp)
                    end if
                else
                    term = (f_obs(i) - expected(i)) ** 2 / expected(i)
                end if
            else if (lambda_value == 0.0_dp) then
                if (f_obs(i) == 0.0_dp) then
                    term = 0.0_dp
                else if (expected(i) <= 0.0_dp) then
                    term = positive_infinity(0.0_dp)
                else
                    term = 2.0_dp * f_obs(i) * log(f_obs(i) / expected(i))
                end if
            else if (lambda_value == -1.0_dp) then
                if (expected(i) == 0.0_dp) then
                    term = 0.0_dp
                else if (f_obs(i) <= 0.0_dp) then
                    term = positive_infinity(0.0_dp)
                else
                    term = 2.0_dp * expected(i) * log(expected(i) / f_obs(i))
                end if
            else
                if (f_obs(i) == 0.0_dp .or. expected(i) == 0.0_dp) then
                    term = power_divergence_zero_term(f_obs(i), expected(i), lambda_value)
                else
                    ratio = f_obs(i) / expected(i)
                    term = f_obs(i) * (ratio ** lambda_value - 1.0_dp) / &
                        (0.5_dp * lambda_value * (lambda_value + 1.0_dp))
                end if
            end if
            result_value%statistic = result_value%statistic + term
        end do

        df = size(f_obs) - 1 - delta_df
        if (df <= 0 .or. ieee_is_nan(result_value%statistic)) then
            result_value%pvalue = quiet_nan(0.0_dp)
        else
            result_value%pvalue = chi2_sf(result_value%statistic, real(df, dp))
        end if
    end function power_divergence

    function power_divergence_named(f_obs, f_exp, ddof, lambda_name) result(result_value)
        real(dp), intent(in) :: f_obs(:) !! observed category frequencies
        real(dp), intent(in), optional :: f_exp(:) !! expected frequencies; uniform when absent
        integer, intent(in), optional :: ddof !! delta degrees of freedom, default zero
        character(len=*), intent(in) :: lambda_name !! SciPy Cressie-Read alias
        type(power_divergence_result) :: result_value

        real(dp) :: lambda_value

        select case (trim(lambda_name))
        case ('pearson')
            lambda_value = 1.0_dp
        case ('log-likelihood')
            lambda_value = 0.0_dp
        case ('freeman-tukey')
            lambda_value = -0.5_dp
        case ('mod-log-likelihood')
            lambda_value = -1.0_dp
        case ('neyman')
            lambda_value = -2.0_dp
        case ('cressie-read')
            lambda_value = 2.0_dp / 3.0_dp
        case default
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end select

        if (present(f_exp)) then
            if (present(ddof)) then
                result_value = power_divergence(f_obs, f_exp, ddof, lambda_value)
            else
                result_value = power_divergence(f_obs, f_exp, lambda_=lambda_value)
            end if
        else
            if (present(ddof)) then
                result_value = power_divergence(f_obs, ddof=ddof, lambda_=lambda_value)
            else
                result_value = power_divergence(f_obs, lambda_=lambda_value)
            end if
        end if
    end function power_divergence_named

    function chisquare(f_obs, f_exp, ddof) result(result_value)
        real(dp), intent(in) :: f_obs(:) !! observed category frequencies
        real(dp), intent(in), optional :: f_exp(:) !! expected frequencies; uniform when absent
        integer, intent(in), optional :: ddof !! delta degrees of freedom, default zero
        type(power_divergence_result) :: result_value

        if (present(f_exp)) then
            if (present(ddof)) then
                result_value = power_divergence(f_obs, f_exp, ddof, 1.0_dp)
            else
                result_value = power_divergence(f_obs, f_exp, lambda_=1.0_dp)
            end if
        else
            if (present(ddof)) then
                result_value = power_divergence(f_obs, ddof=ddof, lambda_=1.0_dp)
            else
                result_value = power_divergence(f_obs, lambda_=1.0_dp)
            end if
        end if
    end function chisquare

    function fisher_exact(table, alternative) result(result_value)
        integer, intent(in) :: table(2, 2) !! nonnegative two-by-two contingency table
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(fisher_exact_result) :: result_value

        integer :: a
        integer :: b
        integer :: c
        integer :: d
        integer :: lower
        integer :: m
        integer :: n_col
        integer :: n_row
        integer :: upper
        integer :: x
        character(len=9) :: selected_alternative
        real(dp) :: p_exact
        real(dp) :: p_x

        selected_alternative = choose_alternative(alternative)
        if (any(table < 0) .or. .not. valid_alternative(selected_alternative)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        a = table(1, 1)
        b = table(1, 2)
        c = table(2, 1)
        d = table(2, 2)
        if (a + b == 0 .or. c + d == 0 .or. a + c == 0 .or. b + d == 0) then
            result_value%statistic = quiet_nan(0.0_dp)
            result_value%pvalue = 1.0_dp
            return
        end if
        if (c > 0 .and. b > 0) then
            result_value%statistic = real(a, dp) * real(d, dp) / (real(c, dp) * real(b, dp))
        else
            result_value%statistic = positive_infinity(0.0_dp)
        end if

        m = a + b + c + d
        n_row = a + b
        n_col = a + c
        select case (selected_alternative)
        case ('less')
            result_value%pvalue = hypergeom_cdf(a, m, n_row, n_col)
        case ('greater')
            result_value%pvalue = hypergeom_sf(a - 1, m, n_row, n_col)
        case ('two-sided')
            p_exact = hypergeom_pmf(a, m, n_row, n_col)
            lower = max(0, n_row + n_col - m)
            upper = min(n_row, n_col)
            result_value%pvalue = 0.0_dp
            do x = lower, upper
                p_x = hypergeom_pmf(x, m, n_row, n_col)
                if (p_x <= p_exact * (1.0_dp + 1.0e-14_dp)) then
                    result_value%pvalue = result_value%pvalue + p_x
                end if
            end do
            result_value%pvalue = min(1.0_dp, result_value%pvalue)
        end select
    end function fisher_exact

    function contingency_expected_freq(observed) result(expected)
        real(dp), intent(in) :: observed(:, :) !! nonnegative contingency-table frequencies
        real(dp), allocatable :: expected(:, :)

        integer :: i
        integer :: j
        real(dp), allocatable :: col_sum(:)
        real(dp), allocatable :: row_sum(:)
        real(dp) :: total

        allocate(expected(size(observed, 1), size(observed, 2)))
        if (size(observed, 1) < 1 .or. size(observed, 2) < 1 .or. &
                any(ieee_is_nan(observed)) .or. any(observed < 0.0_dp)) then
            expected = quiet_nan(0.0_dp)
            return
        end if
        total = sum(observed)
        if (total <= 0.0_dp) then
            expected = quiet_nan(0.0_dp)
            return
        end if
        allocate(row_sum(size(observed, 1)), col_sum(size(observed, 2)))
        do i = 1, size(observed, 1)
            row_sum(i) = sum(observed(i, :))
        end do
        do j = 1, size(observed, 2)
            col_sum(j) = sum(observed(:, j))
        end do
        do j = 1, size(observed, 2)
            do i = 1, size(observed, 1)
                expected(i, j) = row_sum(i) * col_sum(j) / total
            end do
        end do
    end function contingency_expected_freq

    function chi2_contingency(observed, correction, lambda_) result(result_value)
        real(dp), intent(in) :: observed(:, :) !! nonnegative two-dimensional table
        logical, intent(in), optional :: correction !! use Yates correction when dof is one
        real(dp), intent(in), optional :: lambda_ !! Cressie-Read lambda, default Pearson one
        type(chi2_contingency_result) :: result_value

        logical :: use_correction
        integer :: ddof
        integer :: i
        integer :: j
        real(dp) :: difference
        real(dp) :: lambda_value
        real(dp) :: magnitude
        real(dp), allocatable :: adjusted(:, :)
        type(power_divergence_result) :: pd

        result_value%expected_freq = contingency_expected_freq(observed)
        result_value%dof = (size(observed, 1) - 1) * (size(observed, 2) - 1)
        if (any(ieee_is_nan(result_value%expected_freq)) .or. &
                any(result_value%expected_freq <= 0.0_dp)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        if (result_value%dof == 0) then
            result_value%statistic = 0.0_dp
            result_value%pvalue = 1.0_dp
            return
        end if

        use_correction = .true.
        if (present(correction)) use_correction = correction
        lambda_value = 1.0_dp
        if (present(lambda_)) lambda_value = lambda_
        allocate(adjusted(size(observed, 1), size(observed, 2)))
        adjusted = observed
        if (result_value%dof == 1 .and. use_correction) then
            do j = 1, size(observed, 2)
                do i = 1, size(observed, 1)
                    difference = result_value%expected_freq(i, j) - adjusted(i, j)
                    magnitude = min(0.5_dp, abs(difference))
                    adjusted(i, j) = adjusted(i, j) + sign(magnitude, difference)
                end do
            end do
        end if
        ddof = size(observed) - 1 - result_value%dof
        pd = power_divergence(reshape(adjusted, [size(adjusted)]), &
            reshape(result_value%expected_freq, [size(result_value%expected_freq)]), &
            ddof, lambda_value)
        result_value%statistic = pd%statistic
        result_value%pvalue = pd%pvalue
    end function chi2_contingency

    function f_oneway(groups, equal_var) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples
        logical, intent(in), optional :: equal_var !! pooled ANOVA when true, Welch when false
        type(f_oneway_result) :: result_value

        logical :: pooled
        integer :: i
        integer :: k
        integer :: n_total
        real(dp) :: denominator
        real(dp) :: grand_mean
        real(dp), allocatable :: means(:)
        real(dp) :: ms_between
        real(dp) :: ms_within
        real(dp) :: numerator
        real(dp) :: ss_between
        real(dp) :: ss_within
        real(dp) :: sum_weights
        real(dp) :: term
        real(dp), allocatable :: variances(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: weighted_mean
        real(dp) :: welch_df2

        pooled = .true.
        if (present(equal_var)) pooled = equal_var
        if (.not. valid_groups(groups, 2, .true.)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        k = size(groups)
        n_total = total_group_size(groups)
        if (all([(size(groups(i)%values) == 1, i=1,k)])) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        allocate(means(k), variances(k))
        do i = 1, k
            means(i) = mean(groups(i)%values)
            variances(i) = variance(groups(i)%values)
        end do
        if (any(ieee_is_nan(means))) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        if (all_groups_constant(groups)) then
            if (all_group_constants_equal(groups)) then
                call invalidate_pair(result_value%statistic, result_value%pvalue)
            else
                result_value%statistic = positive_infinity(0.0_dp)
                if (pooled) then
                    result_value%pvalue = 0.0_dp
                else
                    result_value%pvalue = quiet_nan(0.0_dp)
                end if
            end if
            return
        end if

        if (pooled) then
            grand_mean = 0.0_dp
            do i = 1, k
                grand_mean = grand_mean + real(size(groups(i)%values), dp) * means(i)
            end do
            grand_mean = grand_mean / real(n_total, dp)
            ss_between = 0.0_dp
            ss_within = 0.0_dp
            do i = 1, k
                ss_between = ss_between + real(size(groups(i)%values), dp) * &
                    (means(i) - grand_mean) ** 2
                if (size(groups(i)%values) > 1) then
                    ss_within = ss_within + real(size(groups(i)%values) - 1, dp) * variances(i)
                end if
            end do
            ms_between = ss_between / real(k - 1, dp)
            ms_within = ss_within / real(n_total - k, dp)
            if (ms_within == 0.0_dp) then
                result_value%statistic = positive_infinity(0.0_dp)
                result_value%pvalue = 0.0_dp
            else
                result_value%statistic = ms_between / ms_within
                result_value%pvalue = f_sf(result_value%statistic, &
                    real(k - 1, dp), real(n_total - k, dp))
            end if
        else
            if (any([(size(groups(i)%values) < 2, i=1,k)]) .or. &
                    any(.not. ieee_is_finite(variances)) .or. any(variances <= 0.0_dp)) then
                call invalidate_pair(result_value%statistic, result_value%pvalue)
                return
            end if
            allocate(weights(k))
            do i = 1, k
                weights(i) = real(size(groups(i)%values), dp) / variances(i)
            end do
            sum_weights = sum(weights)
            weighted_mean = sum(weights * means) / sum_weights
            numerator = sum(weights * (means - weighted_mean) ** 2) / real(k - 1, dp)
            term = 0.0_dp
            do i = 1, k
                term = term + (1.0_dp - weights(i) / sum_weights) ** 2 / &
                    real(size(groups(i)%values) - 1, dp)
            end do
            denominator = 1.0_dp + 2.0_dp * real(k - 2, dp) * term / real(k * k - 1, dp)
            result_value%statistic = numerator / denominator
            welch_df2 = real(k * k - 1, dp) / (3.0_dp * term)
            result_value%pvalue = f_sf(result_value%statistic, real(k - 1, dp), welch_df2)
        end if
    end function f_oneway

    function bartlett(groups) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples
        type(bartlett_result) :: result_value

        integer :: i
        integer :: k
        integer :: n_total
        real(dp) :: denominator
        real(dp) :: numerator
        real(dp) :: pooled_variance
        real(dp) :: sum_inverse_df
        real(dp), allocatable :: variances(:)

        if (.not. valid_groups(groups, 2, .true.)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        k = size(groups)
        if (any([(size(groups(i)%values) < 2, i=1,k)])) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        n_total = total_group_size(groups)
        allocate(variances(k))
        do i = 1, k
            variances(i) = variance(groups(i)%values)
        end do
        if (any(.not. ieee_is_finite(variances)) .or. any(variances <= 0.0_dp)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        pooled_variance = 0.0_dp
        sum_inverse_df = 0.0_dp
        numerator = 0.0_dp
        do i = 1, k
            pooled_variance = pooled_variance + real(size(groups(i)%values) - 1, dp) * variances(i)
            sum_inverse_df = sum_inverse_df + 1.0_dp / real(size(groups(i)%values) - 1, dp)
            numerator = numerator - real(size(groups(i)%values) - 1, dp) * log(variances(i))
        end do
        pooled_variance = pooled_variance / real(n_total - k, dp)
        numerator = numerator + real(n_total - k, dp) * log(pooled_variance)
        denominator = 1.0_dp + (sum_inverse_df - 1.0_dp / real(n_total - k, dp)) / &
            (3.0_dp * real(k - 1, dp))
        result_value%statistic = max(0.0_dp, numerator / denominator)
        result_value%pvalue = chi2_sf(result_value%statistic, real(k - 1, dp))
    end function bartlett

    function levene(groups, center, proportiontocut) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples
        character(len=*), intent(in), optional :: center !! median, mean, or trimmed
        real(dp), intent(in), optional :: proportiontocut !! trim fraction per tail, default 0.05
        type(levene_result) :: result_value

        character(len=7) :: selected_center
        integer :: i
        integer :: k
        integer :: n_total
        real(dp), allocatable :: centers(:)
        real(dp) :: denominator
        real(dp) :: grand_z
        real(dp) :: numerator
        real(dp) :: trim_fraction
        real(dp), allocatable :: zbar(:)

        selected_center = 'median'
        if (present(center)) selected_center = center
        trim_fraction = 0.05_dp
        if (present(proportiontocut)) trim_fraction = proportiontocut
        if (.not. valid_center(selected_center, trim_fraction) .or. &
                .not. valid_groups(groups, 2, .false.)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        k = size(groups)
        n_total = total_group_size(groups)
        if (n_total <= k) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        allocate(centers(k), zbar(k))
        do i = 1, k
            centers(i) = group_center(groups(i)%values, selected_center, trim_fraction)
            zbar(i) = mean(abs(groups(i)%values - centers(i)))
        end do
        grand_z = 0.0_dp
        do i = 1, k
            grand_z = grand_z + real(size(groups(i)%values), dp) * zbar(i)
        end do
        grand_z = grand_z / real(n_total, dp)
        numerator = 0.0_dp
        denominator = 0.0_dp
        do i = 1, k
            numerator = numerator + real(size(groups(i)%values), dp) * (zbar(i) - grand_z) ** 2
            denominator = denominator + sum((abs(groups(i)%values - centers(i)) - zbar(i)) ** 2)
        end do
        numerator = real(n_total - k, dp) * numerator
        denominator = real(k - 1, dp) * denominator
        if (denominator <= 0.0_dp) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        result_value%statistic = numerator / denominator
        result_value%pvalue = f_sf(result_value%statistic, real(k - 1, dp), real(n_total - k, dp))
    end function levene

    function fligner(groups, center, proportiontocut) result(result_value)
        type(sample_group), intent(in) :: groups(:) !! two or more independent samples
        character(len=*), intent(in), optional :: center !! median, mean, or trimmed
        real(dp), intent(in), optional :: proportiontocut !! trim fraction per tail, default 0.05
        type(fligner_result) :: result_value

        character(len=7) :: selected_center
        integer :: first
        integer :: i
        integer :: k
        integer :: last
        integer :: n_total
        real(dp), allocatable :: centers(:)
        real(dp), allocatable :: deviations(:)
        real(dp) :: grand_score
        real(dp), allocatable :: ranks(:)
        real(dp), allocatable :: scores(:)
        real(dp), allocatable :: score_means(:)
        real(dp) :: score_variance
        real(dp) :: trim_fraction

        selected_center = 'median'
        if (present(center)) selected_center = center
        trim_fraction = 0.05_dp
        if (present(proportiontocut)) trim_fraction = proportiontocut
        if (.not. valid_center(selected_center, trim_fraction) .or. &
                .not. valid_groups(groups, 2, .false.)) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        k = size(groups)
        n_total = total_group_size(groups)
        if (n_total < 2) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        allocate(centers(k), deviations(n_total), score_means(k))
        first = 1
        do i = 1, k
            centers(i) = group_center(groups(i)%values, selected_center, trim_fraction)
            last = first + size(groups(i)%values) - 1
            deviations(first:last) = abs(groups(i)%values - centers(i))
            first = last + 1
        end do
        ranks = rankdata(deviations)
        allocate(scores(n_total))
        do i = 1, n_total
            scores(i) = normal_ppf(0.5_dp + ranks(i) / (2.0_dp * real(n_total + 1, dp)))
        end do
        grand_score = mean(scores)
        score_variance = variance(scores)
        if (.not. ieee_is_finite(score_variance) .or. score_variance <= 0.0_dp) then
            call invalidate_pair(result_value%statistic, result_value%pvalue)
            return
        end if
        first = 1
        do i = 1, k
            last = first + size(groups(i)%values) - 1
            score_means(i) = mean(scores(first:last))
            first = last + 1
        end do
        result_value%statistic = 0.0_dp
        do i = 1, k
            result_value%statistic = result_value%statistic + &
                real(size(groups(i)%values), dp) * (score_means(i) - grand_score) ** 2
        end do
        result_value%statistic = result_value%statistic / score_variance
        result_value%pvalue = chi2_sf(result_value%statistic, real(k - 1, dp))
    end function fligner

    subroutine wilcoxon_components(d, zero_method, r_plus, r_minus, mn, se, rank_count, n_zero, has_ties)
        real(dp), intent(in) :: d(:) !! paired differences
        character(len=*), intent(in) :: zero_method !! validated zero method
        real(dp), intent(out) :: r_plus !! sum of positive ranks
        real(dp), intent(out) :: r_minus !! sum of negative ranks
        real(dp), intent(out) :: mn !! null mean of positive rank sum
        real(dp), intent(out) :: se !! tie/zero-corrected standard error
        integer, intent(out) :: rank_count !! number used in null mean/variance
        integer, intent(out) :: n_zero !! number of exact zero differences
        logical, intent(out) :: has_ties !! whether ranked absolute differences contain ties

        integer :: i
        integer :: j
        integer :: m
        real(dp), allocatable :: abs_values(:)
        real(dp), allocatable :: ranks(:)
        real(dp) :: tie_sum
        real(dp) :: zero_rank_sum

        n_zero = count(d == 0.0_dp)
        if (zero_method == 'wilcox') then
            m = size(d) - n_zero
            allocate(abs_values(m))
            j = 0
            do i = 1, size(d)
                if (d(i) /= 0.0_dp) then
                    j = j + 1
                    abs_values(j) = abs(d(i))
                end if
            end do
            rank_count = m
            if (m == 0) then
                r_plus = 0.0_dp
                r_minus = 0.0_dp
                mn = 0.0_dp
                se = 0.0_dp
                has_ties = .false.
                return
            end if
            ranks = rankdata(abs_values)
            r_plus = 0.0_dp
            r_minus = 0.0_dp
            j = 0
            do i = 1, size(d)
                if (d(i) /= 0.0_dp) then
                    j = j + 1
                    if (d(i) > 0.0_dp) r_plus = r_plus + ranks(j)
                    if (d(i) < 0.0_dp) r_minus = r_minus + ranks(j)
                end if
            end do
            tie_sum = tie_cubic_sum(abs_values)
            has_ties = has_duplicate(abs_values)
        else
            rank_count = size(d)
            allocate(abs_values(rank_count))
            abs_values = abs(d)
            ranks = rankdata(abs_values)
            r_plus = 0.0_dp
            r_minus = 0.0_dp
            zero_rank_sum = 0.0_dp
            do i = 1, size(d)
                if (d(i) > 0.0_dp) r_plus = r_plus + ranks(i)
                if (d(i) < 0.0_dp) r_minus = r_minus + ranks(i)
                if (d(i) == 0.0_dp) zero_rank_sum = zero_rank_sum + ranks(i)
            end do
            if (zero_method == 'zsplit') then
                r_plus = r_plus + 0.5_dp * zero_rank_sum
                r_minus = r_minus + 0.5_dp * zero_rank_sum
            end if
            tie_sum = tie_cubic_sum(abs_values)
            if (zero_method == 'pratt' .and. n_zero > 1) then
                tie_sum = tie_sum - real(n_zero ** 3 - n_zero, dp)
            end if
            has_ties = has_duplicate_nonzero(abs_values)
        end if

        mn = real(rank_count * (rank_count + 1), dp) / 4.0_dp
        se = real(rank_count * (rank_count + 1) * (2 * rank_count + 1), dp)
        if (zero_method == 'pratt') then
            mn = mn - real(n_zero * (n_zero + 1), dp) / 4.0_dp
            se = se - real(n_zero * (n_zero + 1) * (2 * n_zero + 1), dp)
        end if
        se = sqrt(max(0.0_dp, (se - 0.5_dp * tie_sum) / 24.0_dp))
    end subroutine wilcoxon_components

    function wilcoxon_exact_pvalue(r_plus, count, alternative) result(pvalue)
        real(dp), intent(in) :: r_plus !! observed positive-rank sum
        integer, intent(in) :: count !! number of ranks in exact null distribution
        character(len=*), intent(in) :: alternative !! validated alternative
        real(dp) :: pvalue

        integer :: i
        integer :: max_sum
        integer :: upper_index
        integer :: lower_index
        real(dp), allocatable :: next_prob(:)
        real(dp), allocatable :: prob(:)

        if (count < 0 .or. count > 2000) then
            pvalue = quiet_nan(0.0_dp)
            return
        end if
        max_sum = count * (count + 1) / 2
        allocate(prob(0:max_sum), next_prob(0:max_sum))
        prob = 0.0_dp
        prob(0) = 1.0_dp
        do i = 1, count
            next_prob = 0.5_dp * prob
            next_prob(i:max_sum) = next_prob(i:max_sum) + 0.5_dp * prob(0:max_sum-i)
            prob = next_prob
        end do
        lower_index = max(0, min(max_sum, int(floor(r_plus))))
        upper_index = max(0, min(max_sum, int(ceiling(r_plus))))
        select case (alternative)
        case ('less')
            pvalue = sum(prob(0:upper_index))
        case ('greater')
            pvalue = sum(prob(lower_index:max_sum))
        case ('two-sided')
            pvalue = 2.0_dp * min(sum(prob(0:upper_index)), sum(prob(lower_index:max_sum)))
            pvalue = min(1.0_dp, pvalue)
        end select
    end function wilcoxon_exact_pvalue

    function wilcoxon_permutation_pvalue(d, zero_method, observed_r_plus, alternative) result(pvalue)
        real(dp), intent(in) :: d(:) !! paired differences, at most 22 entries
        character(len=*), intent(in) :: zero_method !! validated zero treatment
        real(dp), intent(in) :: observed_r_plus !! observed positive rank sum
        character(len=*), intent(in) :: alternative !! validated alternative
        real(dp) :: pvalue

        integer :: bits
        integer :: i
        integer :: n_outcomes
        integer :: n_less
        integer :: n_greater
        real(dp), allocatable :: d_perm(:)
        real(dp) :: dummy_mn
        real(dp) :: dummy_minus
        real(dp) :: dummy_se
        real(dp) :: r_plus
        integer :: dummy_count
        integer :: dummy_zero
        logical :: dummy_ties
        real(dp), parameter :: tol = 64.0_dp * epsilon(1.0_dp)

        n_outcomes = 2 ** size(d)
        allocate(d_perm(size(d)))
        n_less = 0
        n_greater = 0
        do bits = 0, n_outcomes - 1
            do i = 1, size(d)
                if (btest(bits, i - 1)) then
                    d_perm(i) = abs(d(i))
                else
                    d_perm(i) = -abs(d(i))
                end if
                if (d(i) == 0.0_dp) d_perm(i) = 0.0_dp
            end do
            call wilcoxon_components(d_perm, zero_method, r_plus, dummy_minus, dummy_mn, &
                dummy_se, dummy_count, dummy_zero, dummy_ties)
            if (r_plus <= observed_r_plus + tol) n_less = n_less + 1
            if (r_plus >= observed_r_plus - tol) n_greater = n_greater + 1
        end do
        select case (alternative)
        case ('less')
            pvalue = real(n_less, dp) / real(n_outcomes, dp)
        case ('greater')
            pvalue = real(n_greater, dp) / real(n_outcomes, dp)
        case ('two-sided')
            pvalue = min(1.0_dp, 2.0_dp * real(min(n_less, n_greater), dp) / real(n_outcomes, dp))
        end select
    end function wilcoxon_permutation_pvalue

    pure function normal_pvalue(statistic, alternative) result(pvalue)
        real(dp), intent(in) :: statistic !! normal-reference statistic
        character(len=*), intent(in) :: alternative !! validated alternative
        real(dp) :: pvalue

        select case (alternative)
        case ('two-sided')
            pvalue = min(1.0_dp, 2.0_dp * normal_sf(abs(statistic)))
        case ('less')
            pvalue = normal_cdf(statistic)
        case ('greater')
            pvalue = normal_sf(statistic)
        end select
    end function normal_pvalue

    pure function correction_sign(z, alternative) result(value)
        real(dp), intent(in) :: z !! uncorrected asymptotic z statistic
        character(len=*), intent(in) :: alternative !! validated alternative
        real(dp) :: value

        select case (alternative)
        case ('greater')
            value = 1.0_dp
        case ('less')
            value = -1.0_dp
        case default
            value = sign(1.0_dp, z)
            if (z == 0.0_dp) value = 0.0_dp
        end select
    end function correction_sign

    function tie_cubic_sum(values) result(value)
        real(dp), intent(in) :: values(:) !! values whose exact ties are counted
        real(dp) :: value

        integer :: group
        integer :: i
        integer :: n_groups
        integer, allocatable :: counts(:)
        real(dp), allocatable :: dense(:)

        value = 0.0_dp
        if (size(values) < 2) return
        dense = rankdata(values, 'dense')
        if (any(ieee_is_nan(dense))) then
            value = quiet_nan(0.0_dp)
            return
        end if
        n_groups = int(maxval(dense))
        allocate(counts(n_groups))
        counts = 0
        do i = 1, size(values)
            counts(int(dense(i))) = counts(int(dense(i))) + 1
        end do
        do group = 1, n_groups
            value = value + real(counts(group) ** 3 - counts(group), dp)
        end do
    end function tie_cubic_sum

    function has_duplicate(values) result(found)
        real(dp), intent(in) :: values(:) !! values checked for exact ties
        logical :: found

        integer :: i
        integer :: j

        found = .false.
        do i = 1, size(values) - 1
            do j = i + 1, size(values)
                if (values(i) == values(j)) then
                    found = .true.
                    return
                end if
            end do
        end do
    end function has_duplicate

    function has_duplicate_nonzero(values) result(found)
        real(dp), intent(in) :: values(:) !! absolute values checked for nonzero exact ties
        logical :: found

        integer :: i
        integer :: j

        found = .false.
        do i = 1, size(values) - 1
            if (values(i) == 0.0_dp) cycle
            do j = i + 1, size(values)
                if (values(i) == values(j)) then
                    found = .true.
                    return
                end if
            end do
        end do
    end function has_duplicate_nonzero

    function group_center(values, center, proportiontocut) result(value)
        real(dp), intent(in) :: values(:) !! sample observations
        character(len=*), intent(in) :: center !! mean, median, or trimmed
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail for trimmed center
        real(dp) :: value

        select case (center)
        case ('mean')
            value = mean(values)
        case ('median')
            value = median(values)
        case ('trimmed')
            value = trimmed_mean(values, proportiontocut)
        end select
    end function group_center

    function trimmed_mean(values, proportiontocut) result(value)
        real(dp), intent(in) :: values(:) !! sample observations
        real(dp), intent(in) :: proportiontocut !! fraction discarded from each tail
        real(dp) :: value

        integer :: cut
        real(dp), allocatable :: sorted(:)

        cut = int(proportiontocut * real(size(values), dp))
        if (2 * cut >= size(values)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        sorted = values
        call sort_real(sorted)
        value = mean(sorted(cut + 1:size(sorted) - cut))
    end function trimmed_mean

    subroutine sort_real(values)
        real(dp), intent(inout) :: values(:) !! sorted ascending in place

        integer :: i
        integer :: j
        real(dp) :: key

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort_real

    function power_divergence_zero_term(obs, expected, lambda_value) result(term)
        real(dp), intent(in) :: obs !! observed frequency, zero or positive
        real(dp), intent(in) :: expected !! expected frequency, zero or positive
        real(dp), intent(in) :: lambda_value !! finite Cressie-Read exponent excluding 0, 1, and -1
        real(dp) :: term

        if (obs == 0.0_dp .and. expected == 0.0_dp) then
            term = quiet_nan(0.0_dp)
        else if (obs == 0.0_dp) then
            if (lambda_value > -1.0_dp) then
                term = 0.0_dp
            else
                term = positive_infinity(0.0_dp)
            end if
        else
            if (lambda_value > 0.0_dp) then
                term = positive_infinity(0.0_dp)
            else
                term = -obs / (0.5_dp * lambda_value * (lambda_value + 1.0_dp))
            end if
        end if
    end function power_divergence_zero_term

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

    pure logical function valid_zero_method(zero_method) result(valid)
        character(len=*), intent(in) :: zero_method !! signed-rank zero treatment

        valid = zero_method == 'wilcox' .or. zero_method == 'pratt' .or. zero_method == 'zsplit'
    end function valid_zero_method

    pure logical function valid_wilcoxon_method(method) result(valid)
        character(len=*), intent(in) :: method !! signed-rank p-value method

        valid = trim(method) == 'auto' .or. trim(method) == 'exact' .or. &
            trim(method) == 'asymptotic' .or. trim(method) == 'permutation'
    end function valid_wilcoxon_method

    pure logical function valid_center(center, proportiontocut) result(valid)
        character(len=*), intent(in) :: center !! center definition
        real(dp), intent(in) :: proportiontocut !! trim fraction

        valid = (center == 'mean' .or. center == 'median' .or. center == 'trimmed') .and. &
            proportiontocut >= 0.0_dp .and. proportiontocut < 0.5_dp
    end function valid_center

    function valid_groups(groups, min_groups, allow_singletons) result(valid)
        type(sample_group), intent(in) :: groups(:) !! group collection to validate
        integer, intent(in) :: min_groups !! minimum number of groups
        logical, intent(in) :: allow_singletons !! whether one-observation groups are valid
        logical :: valid

        integer :: i

        valid = size(groups) >= min_groups
        if (.not. valid) return
        do i = 1, size(groups)
            if (.not. allocated(groups(i)%values) .or. size(groups(i)%values) == 0) then
                valid = .false.
                return
            end if
            if (.not. allow_singletons .and. size(groups(i)%values) < 1) then
                valid = .false.
                return
            end if
            if (any(ieee_is_nan(groups(i)%values))) then
                valid = .false.
                return
            end if
        end do
    end function valid_groups

    pure function total_group_size(groups) result(total)
        type(sample_group), intent(in) :: groups(:) !! allocated sample groups
        integer :: total

        integer :: i

        total = 0
        do i = 1, size(groups)
            total = total + size(groups(i)%values)
        end do
    end function total_group_size

    function all_groups_constant(groups) result(all_constant)
        type(sample_group), intent(in) :: groups(:) !! allocated sample groups
        logical :: all_constant

        integer :: i

        all_constant = .true.
        do i = 1, size(groups)
            if (any(groups(i)%values /= groups(i)%values(1))) then
                all_constant = .false.
                return
            end if
        end do
    end function all_groups_constant

    function all_group_constants_equal(groups) result(all_equal)
        type(sample_group), intent(in) :: groups(:) !! constant sample groups
        logical :: all_equal

        integer :: i

        all_equal = .true.
        do i = 2, size(groups)
            if (groups(i)%values(1) /= groups(1)%values(1)) then
                all_equal = .false.
                return
            end if
        end do
    end function all_group_constants_equal

    pure subroutine invalidate_pair(statistic, pvalue)
        real(dp), intent(out) :: statistic !! statistic set to NaN
        real(dp), intent(out) :: pvalue !! p-value set to NaN

        statistic = quiet_nan(0.0_dp)
        pvalue = quiet_nan(0.0_dp)
    end subroutine invalidate_pair

end module scifort_hypothesis_extended
