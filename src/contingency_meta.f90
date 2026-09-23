! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Binomial exact tests, unconditional 2x2 exact tests, epidemiological
! effect measures, contingency-table association measures, p-value
! combination, and false-discovery-rate adjustment.  The public semantics
! and formulas are checked against SciPy 1.17.0.

module scifort_contingency_meta
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_beta, only : beta_ppf
    use scifort_binomial, only : binomial_cdf, binomial_pmf, binomial_sf
    use scifort_chi2, only : chi2_cdf, chi2_sf
    use scifort_hypothesis_extended, only : chi2_contingency, chi2_contingency_result
    use scifort_hypergeom, only : hypergeom_cdf
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_nchypergeom_fisher, only : nchypergeom_fisher_cdf, &
        nchypergeom_fisher_pmf, nchypergeom_fisher_sf
    use scifort_normal, only : normal_isf, normal_ppf, normal_sf
    use scifort_student_t, only : t_sf
    implicit none
    private

    real(dp), parameter :: pi_dp = acos(-1.0_dp)

    type, public :: confidence_interval
        real(dp) :: low !! lower confidence limit
        real(dp) :: high !! upper confidence limit
    end type confidence_interval

    type, public :: binomtest_result
        integer :: k !! observed number of successes
        integer :: n !! number of Bernoulli trials
        character(len=9) :: alternative !! two-sided, less, or greater
        real(dp) :: statistic !! observed proportion k/n
        real(dp) :: pvalue !! exact binomial-test p-value
    contains
        procedure :: proportion_ci => binomtest_proportion_ci
    end type binomtest_result

    type, public :: exact_2x2_result
        real(dp) :: statistic !! Wald statistic or Fisher p-value statistic
        real(dp) :: pvalue !! unconditional exact-test p-value
    end type exact_2x2_result

    type, public :: odds_ratio_result
        real(dp) :: statistic !! sample or conditional odds-ratio estimate
        integer :: table(2, 2) !! original 2x2 table
        character(len=11) :: kind !! conditional or sample
    contains
        procedure :: confidence_interval => odds_ratio_confidence_interval
    end type odds_ratio_result

    type, public :: relative_risk_result
        real(dp) :: relative_risk !! exposed risk divided by control risk
        integer :: exposed_cases !! exposed cases
        integer :: exposed_total !! exposed sample size
        integer :: control_cases !! control cases
        integer :: control_total !! control sample size
    contains
        procedure :: confidence_interval => relative_risk_confidence_interval
    end type relative_risk_result

    type, public :: combined_pvalue_result
        real(dp) :: statistic !! combined-test statistic
        real(dp) :: pvalue !! combined p-value
    end type combined_pvalue_result

    public :: association
    public :: barnard_exact
    public :: binomtest
    public :: boschloo_exact
    public :: combine_pvalues
    public :: false_discovery_control
    public :: margins
    public :: odds_ratio
    public :: relative_risk

    interface margins
        module procedure contingency_margins_2d
        module procedure contingency_margins_2d_int
    end interface margins

contains

    function binomtest(k, n, p, alternative) result(result_value)
        integer, intent(in) :: k !! observed number of successes, 0 <= k <= n
        integer, intent(in) :: n !! number of trials, n >= 1
        real(dp), intent(in), optional :: p !! null success probability, default 0.5
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(binomtest_result) :: result_value

        character(len=9) :: alt
        integer :: j
        real(dp) :: d
        real(dp) :: p0
        real(dp) :: threshold

        p0 = 0.5_dp
        if (present(p)) p0 = p
        alt = choose_alternative(alternative)
        result_value%k = k
        result_value%n = n
        result_value%alternative = alt
        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        if (n < 1 .or. k < 0 .or. k > n .or. p0 < 0.0_dp .or. p0 > 1.0_dp .or. &
                .not. valid_alternative(alt)) return

        result_value%statistic = real(k, dp) / real(n, dp)
        select case (trim(alt))
        case ('less')
            result_value%pvalue = binomial_cdf(k, n, p0)
        case ('greater')
            result_value%pvalue = binomial_sf(k - 1, n, p0)
        case default
            d = binomial_pmf(k, n, p0)
            if (real(k, dp) == p0 * real(n, dp)) then
                result_value%pvalue = 1.0_dp
            else
                threshold = d * (1.0_dp + 1.0e-7_dp)
                result_value%pvalue = 0.0_dp
                do j = 0, n
                    if (binomial_pmf(j, n, p0) <= threshold) then
                        result_value%pvalue = result_value%pvalue + binomial_pmf(j, n, p0)
                    end if
                end do
                result_value%pvalue = min(1.0_dp, result_value%pvalue)
            end if
        end select
    end function binomtest

    function binomtest_proportion_ci(self, confidence_level, method) result(ci)
        class(binomtest_result), intent(in) :: self !! binomial-test result
        real(dp), intent(in), optional :: confidence_level !! confidence level in [0,1], default 0.95
        character(len=*), intent(in), optional :: method !! exact, wilson, or wilsoncc
        type(confidence_interval) :: ci

        character(len=8) :: selected_method
        real(dp) :: level

        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        selected_method = 'exact'
        if (present(method)) selected_method = adjustl(method)
        if (level < 0.0_dp .or. level > 1.0_dp) then
            ci%low = quiet_nan(level)
            ci%high = quiet_nan(level)
            return
        end if
        select case (trim(selected_method))
        case ('exact')
            ci = binom_exact_confidence_interval(self%k, self%n, level, self%alternative)
        case ('wilson')
            ci = binom_wilson_confidence_interval(self%k, self%n, level, self%alternative, .false.)
        case ('wilsoncc')
            ci = binom_wilson_confidence_interval(self%k, self%n, level, self%alternative, .true.)
        case default
            ci%low = quiet_nan(level)
            ci%high = quiet_nan(level)
        end select
    end function binomtest_proportion_ci

    function binom_exact_confidence_interval(k, n, level, alternative) result(ci)
        integer, intent(in) :: k !! observed successes
        integer, intent(in) :: n !! number of trials
        real(dp), intent(in) :: level !! confidence level
        character(len=*), intent(in) :: alternative !! validated alternative
        type(confidence_interval) :: ci

        real(dp) :: alpha

        if (trim(alternative) == 'two-sided') then
            alpha = 0.5_dp * (1.0_dp - level)
            if (k == 0) then
                ci%low = 0.0_dp
            else
                ci%low = beta_ppf(alpha, real(k, dp), real(n - k + 1, dp))
            end if
            if (k == n) then
                ci%high = 1.0_dp
            else
                ci%high = beta_ppf(1.0_dp - alpha, real(k + 1, dp), real(n - k, dp))
            end if
        else if (trim(alternative) == 'less') then
            alpha = 1.0_dp - level
            ci%low = 0.0_dp
            if (k == n) then
                ci%high = 1.0_dp
            else
                ci%high = beta_ppf(1.0_dp - alpha, real(k + 1, dp), real(n - k, dp))
            end if
        else
            alpha = 1.0_dp - level
            if (k == 0) then
                ci%low = 0.0_dp
            else
                ci%low = beta_ppf(alpha, real(k, dp), real(n - k + 1, dp))
            end if
            ci%high = 1.0_dp
        end if
    end function binom_exact_confidence_interval

    function binom_wilson_confidence_interval(k, n, level, alternative, correction) result(ci)
        integer, intent(in) :: k !! observed successes
        integer, intent(in) :: n !! number of trials
        real(dp), intent(in) :: level !! confidence level
        character(len=*), intent(in) :: alternative !! validated alternative
        logical, intent(in) :: correction !! include continuity correction
        type(confidence_interval) :: ci

        real(dp) :: center
        real(dp) :: delta
        real(dp) :: denom
        real(dp) :: dhi
        real(dp) :: dlo
        real(dp) :: p
        real(dp) :: q
        real(dp) :: z

        p = real(k, dp) / real(n, dp)
        q = 1.0_dp - p
        if (trim(alternative) == 'two-sided') then
            z = normal_ppf(0.5_dp + 0.5_dp * level)
        else
            z = normal_ppf(level)
        end if
        denom = 2.0_dp * (real(n, dp) + z * z)
        center = (2.0_dp * real(n, dp) * p + z * z) / denom
        if (correction) then
            if (trim(alternative) == 'less' .or. k == 0) then
                ci%low = 0.0_dp
            else
                dlo = (1.0_dp + z * sqrt(max(0.0_dp, z * z - 2.0_dp - 1.0_dp / real(n, dp) + &
                    4.0_dp * p * (real(n, dp) * q + 1.0_dp)))) / denom
                ci%low = center - dlo
            end if
            if (trim(alternative) == 'greater' .or. k == n) then
                ci%high = 1.0_dp
            else
                dhi = (1.0_dp + z * sqrt(max(0.0_dp, z * z + 2.0_dp - 1.0_dp / real(n, dp) + &
                    4.0_dp * p * (real(n, dp) * q - 1.0_dp)))) / denom
                ci%high = center + dhi
            end if
        else
            delta = z / denom * sqrt(max(0.0_dp, 4.0_dp * real(n, dp) * p * q + z * z))
            if (trim(alternative) == 'less' .or. k == 0) then
                ci%low = 0.0_dp
            else
                ci%low = center - delta
            end if
            if (trim(alternative) == 'greater' .or. k == n) then
                ci%high = 1.0_dp
            else
                ci%high = center + delta
            end if
        end if
        ci%low = max(0.0_dp, ci%low)
        ci%high = min(1.0_dp, ci%high)
    end function binom_wilson_confidence_interval

    function barnard_exact(table, alternative, pooled, n) result(result_value)
        integer, intent(in) :: table(2, 2) !! nonnegative 2x2 table, columns are binomial experiments
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        logical, intent(in), optional :: pooled !! use pooled Wald variance, default true
        integer, intent(in), optional :: n !! optimization sampling control, default 32
        type(exact_2x2_result) :: result_value

        character(len=9) :: alt
        logical :: use_pooled
        integer :: n_points
        integer :: c1
        integer :: c2
        integer :: i
        integer :: j
        real(dp), allocatable :: log_weight(:, :)
        logical, allocatable :: selected(:, :)
        real(dp) :: obs
        real(dp) :: stat

        alt = choose_alternative(alternative)
        use_pooled = .true.
        if (present(pooled)) use_pooled = pooled
        n_points = 32
        if (present(n)) n_points = n
        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        if (n_points <= 0 .or. any(table < 0) .or. .not. valid_alternative(alt)) return
        c1 = sum(table(:, 1))
        c2 = sum(table(:, 2))
        if (c1 == 0 .or. c2 == 0) then
            result_value%pvalue = 1.0_dp
            return
        end if
        obs = wald_statistic(table(1, 1), table(1, 2), c1, c2, use_pooled)
        result_value%statistic = obs
        allocate(selected(0:c1, 0:c2), log_weight(0:c1, 0:c2))
        do i = 0, c1
            do j = 0, c2
                stat = wald_statistic(i, j, c1, c2, use_pooled)
                select case (trim(alt))
                case ('less')
                    selected(i, j) = stat <= obs
                case ('greater')
                    selected(i, j) = stat >= obs
                case default
                    selected(i, j) = abs(stat) >= abs(obs)
                end select
                log_weight(i, j) = log_combination(c1, i) + log_combination(c2, j)
            end do
        end do
        result_value%pvalue = maximize_nuisance_probability(selected, log_weight, c1, c2, n_points)
    end function barnard_exact

    recursive function boschloo_exact(table, alternative, n) result(result_value)
        integer, intent(in) :: table(2, 2) !! nonnegative 2x2 table, columns are binomial experiments
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        integer, intent(in), optional :: n !! optimization sampling control, default 32
        type(exact_2x2_result) :: result_value

        character(len=9) :: alt
        integer :: n_points
        integer :: c1
        integer :: c2
        integer :: i
        integer :: j
        integer :: total
        real(dp), allocatable :: fisher(:, :)
        real(dp), allocatable :: log_weight(:, :)
        logical, allocatable :: selected(:, :)
        type(exact_2x2_result) :: greater_result
        type(exact_2x2_result) :: less_result
        real(dp) :: fisher_stat

        alt = choose_alternative(alternative)
        n_points = 32
        if (present(n)) n_points = n
        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        if (n_points <= 0 .or. any(table < 0) .or. .not. valid_alternative(alt)) return
        c1 = sum(table(:, 1))
        c2 = sum(table(:, 2))
        if (c1 == 0 .or. c2 == 0) return
        if (trim(alt) == 'two-sided') then
            less_result = boschloo_exact(table, 'less', n_points)
            greater_result = boschloo_exact(table, 'greater', n_points)
            if (less_result%pvalue < greater_result%pvalue) then
                result_value%statistic = less_result%statistic
                result_value%pvalue = min(1.0_dp, 2.0_dp * less_result%pvalue)
            else
                result_value%statistic = greater_result%statistic
                result_value%pvalue = min(1.0_dp, 2.0_dp * greater_result%pvalue)
            end if
            return
        end if
        total = c1 + c2
        allocate(fisher(0:c1, 0:c2), selected(0:c1, 0:c2), log_weight(0:c1, 0:c2))
        do i = 0, c1
            do j = 0, c2
                if (trim(alt) == 'less') then
                    fisher(i, j) = hypergeom_cdf(i, total, i + j, c1)
                else
                    fisher(i, j) = hypergeom_cdf(j, total, i + j, c2)
                end if
                log_weight(i, j) = log_combination(c1, i) + log_combination(c2, j)
            end do
        end do
        fisher_stat = fisher(table(1, 1), table(1, 2))
        selected = fisher <= fisher_stat * (1.0_dp + 1.0e-13_dp)
        result_value%statistic = fisher_stat
        result_value%pvalue = maximize_nuisance_probability(selected, log_weight, c1, c2, n_points)
    end function boschloo_exact

    pure function wald_statistic(x1, x2, c1, c2, pooled) result(value)
        integer, intent(in) :: x1 !! successes in first binomial experiment
        integer, intent(in) :: x2 !! successes in second binomial experiment
        integer, intent(in) :: c1 !! first binomial sample size
        integer, intent(in) :: c2 !! second binomial sample size
        logical, intent(in) :: pooled !! pooled-variance flag
        real(dp) :: value

        real(dp) :: p
        real(dp) :: p1
        real(dp) :: p2
        real(dp) :: variance_value

        p1 = real(x1, dp) / real(c1, dp)
        p2 = real(x2, dp) / real(c2, dp)
        if (p1 == p2) then
            value = 0.0_dp
            return
        end if
        if (pooled) then
            p = real(x1 + x2, dp) / real(c1 + c2, dp)
            variance_value = p * (1.0_dp - p) * (1.0_dp / real(c1, dp) + 1.0_dp / real(c2, dp))
        else
            variance_value = p1 * (1.0_dp - p1) / real(c1, dp) + &
                p2 * (1.0_dp - p2) / real(c2, dp)
        end if
        if (variance_value <= 0.0_dp) then
            value = sign(positive_infinity(0.0_dp), p1 - p2)
        else
            value = (p1 - p2) / sqrt(variance_value)
        end if
    end function wald_statistic

    function maximize_nuisance_probability(selected, log_weight, c1, c2, n) result(maximum)
        logical, intent(in) :: selected(0:, 0:) !! selected extreme tables
        real(dp), intent(in) :: log_weight(0:, 0:) !! log binomial coefficients for each table
        integer, intent(in) :: c1 !! first column total
        integer, intent(in) :: c2 !! second column total
        integer, intent(in) :: n !! requested sampling control
        real(dp) :: maximum

        integer :: best_index
        integer :: i
        integer :: n_grid
        real(dp) :: a
        real(dp) :: b
        real(dp) :: p
        real(dp) :: value

        n_grid = max(2048, 64 * next_power_of_two(max(1, n)))
        n_grid = min(n_grid, 65536)
        maximum = -1.0_dp
        best_index = 0
        do i = 0, n_grid
            p = real(i, dp) / real(n_grid, dp)
            value = nuisance_probability(p, selected, log_weight, c1, c2)
            if (value > maximum) then
                maximum = value
                best_index = i
            end if
        end do
        if (best_index > 0 .and. best_index < n_grid) then
            a = real(best_index - 1, dp) / real(n_grid, dp)
            b = real(best_index + 1, dp) / real(n_grid, dp)
            value = golden_maximum(a, b, selected, log_weight, c1, c2)
            maximum = max(maximum, value)
        end if
        maximum = min(1.0_dp, max(0.0_dp, maximum))
    end function maximize_nuisance_probability

    function golden_maximum(a, b, selected, log_weight, c1, c2) result(maximum)
        real(dp), intent(in) :: a !! lower refinement bound
        real(dp), intent(in) :: b !! upper refinement bound
        logical, intent(in) :: selected(0:, 0:) !! selected extreme tables
        real(dp), intent(in) :: log_weight(0:, 0:) !! log binomial coefficients
        integer, intent(in) :: c1 !! first column total
        integer, intent(in) :: c2 !! second column total
        real(dp) :: maximum

        integer :: iter
        real(dp), parameter :: ratio = 0.6180339887498948482_dp
        real(dp) :: left
        real(dp) :: right
        real(dp) :: x1
        real(dp) :: x2
        real(dp) :: f1
        real(dp) :: f2

        left = a
        right = b
        x1 = right - ratio * (right - left)
        x2 = left + ratio * (right - left)
        f1 = nuisance_probability(x1, selected, log_weight, c1, c2)
        f2 = nuisance_probability(x2, selected, log_weight, c1, c2)
        do iter = 1, 100
            if (right - left <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(left) + abs(right))) exit
            if (f1 < f2) then
                left = x1
                x1 = x2
                f1 = f2
                x2 = left + ratio * (right - left)
                f2 = nuisance_probability(x2, selected, log_weight, c1, c2)
            else
                right = x2
                x2 = x1
                f2 = f1
                x1 = right - ratio * (right - left)
                f1 = nuisance_probability(x1, selected, log_weight, c1, c2)
            end if
        end do
        maximum = max(f1, f2, nuisance_probability(0.5_dp * (left + right), selected, log_weight, c1, c2))
    end function golden_maximum

    function nuisance_probability(p, selected, log_weight, c1, c2) result(probability)
        real(dp), intent(in) :: p !! nuisance Bernoulli probability in [0,1]
        logical, intent(in) :: selected(0:, 0:) !! selected extreme tables
        real(dp), intent(in) :: log_weight(0:, 0:) !! log binomial coefficients
        integer, intent(in) :: c1 !! first column total
        integer, intent(in) :: c2 !! second column total
        real(dp) :: probability

        integer :: i
        integer :: j
        integer :: k
        integer :: total
        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: term

        total = c1 + c2
        probability = 0.0_dp
        if (p <= 0.0_dp) then
            if (selected(0, 0)) probability = 1.0_dp
            return
        else if (p >= 1.0_dp) then
            if (selected(c1, c2)) probability = 1.0_dp
            return
        end if
        logp = log(p)
        logq = log(1.0_dp - p)
        do i = 0, c1
            do j = 0, c2
                if (.not. selected(i, j)) cycle
                k = i + j
                term = exp(log_weight(i, j) + real(k, dp) * logp + real(total - k, dp) * logq)
                probability = probability + term
            end do
        end do
    end function nuisance_probability

    function odds_ratio(table, kind) result(result_value)
        integer, intent(in) :: table(2, 2) !! nonnegative 2x2 contingency table
        character(len=*), intent(in), optional :: kind !! conditional or sample
        type(odds_ratio_result) :: result_value

        character(len=11) :: selected_kind

        selected_kind = 'conditional'
        if (present(kind)) selected_kind = adjustl(kind)
        result_value%table = table
        result_value%kind = selected_kind
        result_value%statistic = quiet_nan(0.0_dp)
        if (any(table < 0)) return
        if (trim(selected_kind) /= 'conditional' .and. trim(selected_kind) /= 'sample') return
        if (any(sum(table, dim=1) == 0) .or. any(sum(table, dim=2) == 0)) return
        if (trim(selected_kind) == 'sample') then
            result_value%statistic = sample_odds_ratio(table)
        else
            result_value%statistic = conditional_odds_ratio(table)
        end if
    end function odds_ratio

    function odds_ratio_confidence_interval(self, confidence_level, alternative) result(ci)
        class(odds_ratio_result), intent(in) :: self !! odds-ratio result
        real(dp), intent(in), optional :: confidence_level !! confidence level in [0,1], default 0.95
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(confidence_interval) :: ci

        character(len=9) :: alt
        real(dp) :: level

        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        alt = choose_alternative(alternative)
        if (level < 0.0_dp .or. level > 1.0_dp .or. .not. valid_alternative(alt)) then
            ci%low = quiet_nan(level)
            ci%high = quiet_nan(level)
            return
        end if
        if (any(sum(self%table, dim=1) == 0) .or. any(sum(self%table, dim=2) == 0)) then
            ci%low = 0.0_dp
            ci%high = positive_infinity(0.0_dp)
            return
        end if
        if (trim(self%kind) == 'sample') then
            ci = sample_odds_ratio_ci(self%table, level, alt)
        else
            ci = conditional_odds_ratio_ci(self%table, level, alt)
        end if
    end function odds_ratio_confidence_interval

    pure function sample_odds_ratio(table) result(value)
        integer, intent(in) :: table(2, 2) !! validated 2x2 table
        real(dp) :: value

        if (table(2, 1) > 0 .and. table(1, 2) > 0) then
            value = real(table(1, 1), dp) * real(table(2, 2), dp) / &
                (real(table(2, 1), dp) * real(table(1, 2), dp))
        else if (table(1, 1) == 0 .or. table(2, 2) == 0) then
            value = quiet_nan(0.0_dp)
        else
            value = positive_infinity(0.0_dp)
        end if
    end function sample_odds_ratio

    function conditional_odds_ratio(table) result(value)
        integer, intent(in) :: table(2, 2) !! validated nondegenerate 2x2 table
        real(dp) :: value

        integer :: x
        integer :: m
        integer :: nrow
        integer :: ncol
        integer :: lo
        integer :: hi

        call table_hypergeom_parameters(table, x, m, nrow, ncol, lo, hi)
        if (x == lo) then
            value = 0.0_dp
        else if (x == hi) then
            value = positive_infinity(0.0_dp)
        else
            value = solve_odds_mean(x, m, nrow, ncol)
        end if
    end function conditional_odds_ratio

    function conditional_odds_ratio_ci(table, level, alternative) result(ci)
        integer, intent(in) :: table(2, 2) !! validated nondegenerate 2x2 table
        real(dp), intent(in) :: level !! confidence level
        character(len=*), intent(in) :: alternative !! validated alternative
        type(confidence_interval) :: ci

        real(dp) :: alpha

        if (trim(alternative) == 'two-sided') then
            alpha = 0.5_dp * (1.0_dp - level)
            ci%low = conditional_ci_lower(table, alpha)
            ci%high = conditional_ci_upper(table, alpha)
        else if (trim(alternative) == 'less') then
            ci%low = 0.0_dp
            ci%high = conditional_ci_upper(table, 1.0_dp - level)
        else
            ci%low = conditional_ci_lower(table, 1.0_dp - level)
            ci%high = positive_infinity(0.0_dp)
        end if
    end function conditional_odds_ratio_ci

    function sample_odds_ratio_ci(table, level, alternative) result(ci)
        integer, intent(in) :: table(2, 2) !! validated nondegenerate 2x2 table
        real(dp), intent(in) :: level !! confidence level
        character(len=*), intent(in) :: alternative !! validated alternative
        type(confidence_interval) :: ci

        integer :: i
        integer :: j
        real(dp) :: log_or
        real(dp) :: odds
        real(dp) :: se
        real(dp) :: z

        odds = sample_odds_ratio(table)
        if (.not. ieee_is_finite(odds) .or. odds <= 0.0_dp) then
            ! Reproduce limiting log-Wald behavior without trapping on log(0).
            if (ieee_is_nan(odds)) then
                ci%low = quiet_nan(0.0_dp)
                ci%high = quiet_nan(0.0_dp)
            else if (odds == 0.0_dp) then
                ci%low = 0.0_dp
                ci%high = quiet_nan(0.0_dp)
            else
                ci%low = quiet_nan(0.0_dp)
                ci%high = positive_infinity(0.0_dp)
            end if
            return
        end if
        se = 0.0_dp
        do j = 1, 2
            do i = 1, 2
                if (table(i, j) == 0) then
                    ci%low = 0.0_dp
                    ci%high = positive_infinity(0.0_dp)
                    return
                end if
                se = se + 1.0_dp / real(table(i, j), dp)
            end do
        end do
        se = sqrt(se)
        log_or = log(odds)
        if (trim(alternative) == 'less') then
            z = normal_ppf(level)
            ci%low = 0.0_dp
            ci%high = exp(log_or + z * se)
        else if (trim(alternative) == 'greater') then
            z = normal_ppf(level)
            ci%low = exp(log_or - z * se)
            ci%high = positive_infinity(0.0_dp)
        else
            z = normal_ppf(0.5_dp + 0.5_dp * level)
            ci%low = exp(log_or - z * se)
            ci%high = exp(log_or + z * se)
        end if
    end function sample_odds_ratio_ci

    function conditional_ci_lower(table, alpha) result(value)
        integer, intent(in) :: table(2, 2) !! validated table
        real(dp), intent(in) :: alpha !! one-tail error probability
        real(dp) :: value

        integer :: x
        integer :: m
        integer :: nrow
        integer :: ncol
        integer :: lo
        integer :: hi

        if (sample_odds_ratio(table) == 0.0_dp) then
            value = 0.0_dp
            return
        end if
        call table_hypergeom_parameters(table, x, m, nrow, ncol, lo, hi)
        value = solve_odds_tail(x, m, nrow, ncol, alpha, .true.)
    end function conditional_ci_lower

    function conditional_ci_upper(table, alpha) result(value)
        integer, intent(in) :: table(2, 2) !! validated table
        real(dp), intent(in) :: alpha !! one-tail error probability
        real(dp) :: value

        integer :: x
        integer :: m
        integer :: nrow
        integer :: ncol
        integer :: lo
        integer :: hi
        real(dp) :: sample_or

        sample_or = sample_odds_ratio(table)
        if (.not. ieee_is_finite(sample_or) .and. .not. ieee_is_nan(sample_or)) then
            value = positive_infinity(0.0_dp)
            return
        end if
        call table_hypergeom_parameters(table, x, m, nrow, ncol, lo, hi)
        value = solve_odds_tail(x, m, nrow, ncol, alpha, .false.)
    end function conditional_ci_upper

    subroutine table_hypergeom_parameters(table, x, m, nrow, ncol, lo, hi)
        integer, intent(in) :: table(2, 2) !! 2x2 table
        integer, intent(out) :: x !! upper-left cell
        integer, intent(out) :: m !! grand total
        integer, intent(out) :: nrow !! first-row total
        integer, intent(out) :: ncol !! first-column total
        integer, intent(out) :: lo !! support lower endpoint
        integer, intent(out) :: hi !! support upper endpoint

        x = table(1, 1)
        m = sum(table)
        nrow = sum(table(1, :))
        ncol = sum(table(:, 1))
        lo = max(0, nrow + ncol - m)
        hi = min(nrow, ncol)
    end subroutine table_hypergeom_parameters

    function solve_odds_mean(x, m, nrow, ncol) result(odds)
        integer, intent(in) :: x !! target conditional mean
        integer, intent(in) :: m !! population size
        integer, intent(in) :: nrow !! Type-I population count
        integer, intent(in) :: ncol !! draws
        real(dp) :: odds

        integer :: iter
        real(dp) :: eta_hi
        real(dp) :: eta_lo
        real(dp) :: eta_mid
        real(dp) :: f_hi
        real(dp) :: f_lo
        real(dp) :: f_mid

        eta_lo = 0.0_dp
        f_lo = nchg_mean(m, nrow, ncol, exp(eta_lo)) - real(x, dp)
        eta_hi = eta_lo
        f_hi = f_lo
        if (f_lo > 0.0_dp) then
            do while (f_lo > 0.0_dp .and. eta_lo > -700.0_dp)
                eta_hi = eta_lo
                f_hi = f_lo
                eta_lo = eta_lo - log(2.0_dp)
                f_lo = nchg_mean(m, nrow, ncol, exp(eta_lo)) - real(x, dp)
            end do
        else
            do while (f_hi < 0.0_dp .and. eta_hi < 700.0_dp)
                eta_lo = eta_hi
                f_lo = f_hi
                eta_hi = eta_hi + log(2.0_dp)
                f_hi = nchg_mean(m, nrow, ncol, exp(eta_hi)) - real(x, dp)
            end do
        end if
        do iter = 1, 120
            eta_mid = 0.5_dp * (eta_lo + eta_hi)
            f_mid = nchg_mean(m, nrow, ncol, exp(eta_mid)) - real(x, dp)
            if (f_mid < 0.0_dp) then
                eta_lo = eta_mid
            else
                eta_hi = eta_mid
            end if
            if (eta_hi - eta_lo < 2.0e-13_dp) exit
        end do
        odds = exp(0.5_dp * (eta_lo + eta_hi))
    end function solve_odds_mean

    function solve_odds_tail(x, m, nrow, ncol, alpha, lower_limit) result(odds)
        integer, intent(in) :: x !! observed upper-left count
        integer, intent(in) :: m !! population size
        integer, intent(in) :: nrow !! Type-I population count
        integer, intent(in) :: ncol !! draws
        real(dp), intent(in) :: alpha !! target one-tail probability
        logical, intent(in) :: lower_limit !! true solves SF, false solves CDF
        real(dp) :: odds

        integer :: iter
        real(dp) :: eta_hi
        real(dp) :: eta_lo
        real(dp) :: eta_mid
        real(dp) :: f_hi
        real(dp) :: f_lo
        real(dp) :: f_mid

        eta_lo = 0.0_dp
        f_lo = odds_tail_equation(eta_lo, x, m, nrow, ncol, alpha, lower_limit)
        eta_hi = eta_lo
        f_hi = f_lo
        if (f_lo > 0.0_dp) then
            do while (f_lo > 0.0_dp .and. eta_lo > -700.0_dp)
                eta_hi = eta_lo
                f_hi = f_lo
                eta_lo = eta_lo - log(2.0_dp)
                f_lo = odds_tail_equation(eta_lo, x, m, nrow, ncol, alpha, lower_limit)
            end do
        else
            do while (f_hi < 0.0_dp .and. eta_hi < 700.0_dp)
                eta_lo = eta_hi
                f_lo = f_hi
                eta_hi = eta_hi + log(2.0_dp)
                f_hi = odds_tail_equation(eta_hi, x, m, nrow, ncol, alpha, lower_limit)
            end do
        end if
        do iter = 1, 140
            eta_mid = 0.5_dp * (eta_lo + eta_hi)
            f_mid = odds_tail_equation(eta_mid, x, m, nrow, ncol, alpha, lower_limit)
            if (f_mid < 0.0_dp) then
                eta_lo = eta_mid
            else
                eta_hi = eta_mid
            end if
            if (eta_hi - eta_lo < 1.0e-13_dp) exit
        end do
        odds = exp(0.5_dp * (eta_lo + eta_hi))
    end function solve_odds_tail

    function odds_tail_equation(eta, x, m, nrow, ncol, alpha, lower_limit) result(value)
        real(dp), intent(in) :: eta !! logarithm of odds ratio
        integer, intent(in) :: x !! observed upper-left count
        integer, intent(in) :: m !! population size
        integer, intent(in) :: nrow !! Type-I population count
        integer, intent(in) :: ncol !! draws
        real(dp), intent(in) :: alpha !! target one-tail probability
        logical, intent(in) :: lower_limit !! true uses SF; false uses negative CDF
        real(dp) :: value

        real(dp) :: odds

        odds = exp(eta)
        if (lower_limit) then
            value = nchypergeom_fisher_sf(x - 1, m, nrow, ncol, odds) - alpha
        else
            value = -nchypergeom_fisher_cdf(x, m, nrow, ncol, odds) + alpha
        end if
    end function odds_tail_equation

    function nchg_mean(m, nrow, ncol, odds) result(value)
        integer, intent(in) :: m !! population size
        integer, intent(in) :: nrow !! Type-I population count
        integer, intent(in) :: ncol !! draws
        real(dp), intent(in) :: odds !! noncentral odds ratio
        real(dp) :: value

        integer :: k
        integer :: lo
        integer :: hi

        lo = max(0, nrow + ncol - m)
        hi = min(nrow, ncol)
        value = 0.0_dp
        do k = lo, hi
            value = value + real(k, dp) * nchypergeom_fisher_pmf(k, m, nrow, ncol, odds)
        end do
    end function nchg_mean

    function relative_risk(exposed_cases, exposed_total, control_cases, control_total) result(result_value)
        integer, intent(in) :: exposed_cases !! cases in exposed group
        integer, intent(in) :: exposed_total !! total exposed sample size
        integer, intent(in) :: control_cases !! cases in control group
        integer, intent(in) :: control_total !! total control sample size
        type(relative_risk_result) :: result_value

        result_value%exposed_cases = exposed_cases
        result_value%exposed_total = exposed_total
        result_value%control_cases = control_cases
        result_value%control_total = control_total
        result_value%relative_risk = quiet_nan(0.0_dp)
        if (exposed_cases < 0 .or. exposed_total < 1 .or. control_cases < 0 .or. control_total < 1) return
        if (exposed_cases > exposed_total .or. control_cases > control_total) return
        if (exposed_cases == 0 .and. control_cases == 0) then
            result_value%relative_risk = quiet_nan(0.0_dp)
        else if (exposed_cases == 0) then
            result_value%relative_risk = 0.0_dp
        else if (control_cases == 0) then
            result_value%relative_risk = positive_infinity(0.0_dp)
        else
            result_value%relative_risk = &
                (real(exposed_cases, dp) / real(exposed_total, dp)) / &
                (real(control_cases, dp) / real(control_total, dp))
        end if
    end function relative_risk

    function relative_risk_confidence_interval(self, confidence_level) result(ci)
        class(relative_risk_result), intent(in) :: self !! relative-risk result
        real(dp), intent(in), optional :: confidence_level !! confidence level in [0,1], default 0.95
        type(confidence_interval) :: ci

        real(dp) :: alpha
        real(dp) :: delta
        real(dp) :: level
        real(dp) :: se
        real(dp) :: z

        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        if (level < 0.0_dp .or. level > 1.0_dp) then
            ci%low = quiet_nan(level)
            ci%high = quiet_nan(level)
            return
        end if
        if (self%exposed_cases == 0 .and. self%control_cases == 0) then
            ci%low = quiet_nan(0.0_dp)
            ci%high = quiet_nan(0.0_dp)
        else if (self%exposed_cases == 0) then
            ci%low = 0.0_dp
            ci%high = quiet_nan(0.0_dp)
        else if (self%control_cases == 0) then
            ci%low = quiet_nan(0.0_dp)
            ci%high = positive_infinity(0.0_dp)
        else
            alpha = 1.0_dp - level
            z = normal_ppf(1.0_dp - 0.5_dp * alpha)
            se = sqrt(1.0_dp / real(self%exposed_cases, dp) - 1.0_dp / real(self%exposed_total, dp) + &
                1.0_dp / real(self%control_cases, dp) - 1.0_dp / real(self%control_total, dp))
            delta = z * se
            ci%low = self%relative_risk * exp(-delta)
            ci%high = self%relative_risk * exp(delta)
        end if
    end function relative_risk_confidence_interval

    function association(observed, method, correction, lambda_) result(value)
        integer, intent(in) :: observed(:, :) !! nonnegative integer contingency table
        character(len=*), intent(in), optional :: method !! cramer, tschuprow, or pearson
        logical, intent(in), optional :: correction !! Yates correction passed to chi2_contingency
        real(dp), intent(in), optional :: lambda_ !! Cressie-Read lambda passed to chi2_contingency
        real(dp) :: value

        character(len=9) :: selected_method
        logical :: use_correction
        real(dp) :: lambda_value
        real(dp) :: phi2
        real(dp), allocatable :: table(:, :)
        type(chi2_contingency_result) :: chi2_result

        selected_method = 'cramer'
        if (present(method)) selected_method = adjustl(method)
        value = quiet_nan(0.0_dp)
        if (any(observed < 0) .or. size(observed, 1) < 2 .or. size(observed, 2) < 2) return
        use_correction = .false.
        if (present(correction)) use_correction = correction
        lambda_value = 1.0_dp
        if (present(lambda_)) lambda_value = lambda_
        table = real(observed, dp)
        chi2_result = chi2_contingency(table, use_correction, lambda_value)
        if (ieee_is_nan(chi2_result%statistic) .or. sum(observed) <= 0) return
        phi2 = chi2_result%statistic / real(sum(observed), dp)
        select case (trim(selected_method))
        case ('cramer')
            value = sqrt(phi2 / real(min(size(observed, 1) - 1, size(observed, 2) - 1), dp))
        case ('tschuprow')
            value = sqrt(phi2 / sqrt(real((size(observed, 1) - 1) * (size(observed, 2) - 1), dp)))
        case ('pearson')
            value = sqrt(phi2 / (1.0_dp + phi2))
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function association

    subroutine contingency_margins_2d(a, margin0, margin1)
        real(dp), intent(in) :: a(:, :) !! two-dimensional array
        real(dp), allocatable, intent(out) :: margin0(:, :) !! row sums as shape (rows,1)
        real(dp), allocatable, intent(out) :: margin1(:, :) !! column sums as shape (1,columns)

        integer :: i
        integer :: j

        allocate(margin0(size(a, 1), 1), margin1(1, size(a, 2)))
        do i = 1, size(a, 1)
            margin0(i, 1) = sum(a(i, :))
        end do
        do j = 1, size(a, 2)
            margin1(1, j) = sum(a(:, j))
        end do
    end subroutine contingency_margins_2d


    subroutine contingency_margins_2d_int(a, margin0, margin1)
        integer, intent(in) :: a(:, :) !! two-dimensional integer array
        integer, allocatable, intent(out) :: margin0(:, :) !! row sums as shape (rows,1)
        integer, allocatable, intent(out) :: margin1(:, :) !! column sums as shape (1,columns)

        integer :: i
        integer :: j

        allocate(margin0(size(a, 1), 1), margin1(1, size(a, 2)))
        do i = 1, size(a, 1)
            margin0(i, 1) = sum(a(i, :))
        end do
        do j = 1, size(a, 2)
            margin1(1, j) = sum(a(:, j))
        end do
    end subroutine contingency_margins_2d_int

    function combine_pvalues(pvalues, method, weights) result(result_value)
        real(dp), intent(in) :: pvalues(:) !! independent p-values in [0,1]
        character(len=*), intent(in), optional :: method !! fisher, pearson, mudholkar_george, tippett, or stouffer
        real(dp), intent(in), optional :: weights(:) !! Stouffer weights, same length as pvalues
        type(combined_pvalue_result) :: result_value

        character(len=18) :: selected_method
        integer :: i
        integer :: n
        real(dp) :: approx_factor
        real(dp) :: normalizing_factor
        real(dp) :: nu
        real(dp) :: sum_log_p
        real(dp) :: sum_log_q
        real(dp) :: sum_w2
        real(dp) :: weighted_z

        selected_method = 'fisher'
        if (present(method)) selected_method = adjustl(method)
        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        n = size(pvalues)
        if (n < 1 .or. any(pvalues < 0.0_dp) .or. any(pvalues > 1.0_dp) .or. any(ieee_is_nan(pvalues))) return
        select case (trim(selected_method))
        case ('fisher')
            sum_log_p = 0.0_dp
            do i = 1, n
                if (pvalues(i) == 0.0_dp) then
                    result_value%statistic = positive_infinity(0.0_dp)
                    result_value%pvalue = 0.0_dp
                    return
                end if
                sum_log_p = sum_log_p + log(pvalues(i))
            end do
            result_value%statistic = -2.0_dp * sum_log_p
            result_value%pvalue = chi2_sf(result_value%statistic, 2.0_dp * real(n, dp))
        case ('pearson')
            if (any(pvalues == 1.0_dp)) then
                result_value%statistic = negative_infinity(0.0_dp)
                result_value%pvalue = 1.0_dp
                return
            end if
            sum_log_q = 0.0_dp
            do i = 1, n
                sum_log_q = sum_log_q + log(1.0_dp - pvalues(i))
            end do
            result_value%statistic = 2.0_dp * sum_log_q
            result_value%pvalue = chi2_cdf(-result_value%statistic, 2.0_dp * real(n, dp))
        case ('mudholkar_george')
            if (any(pvalues == 0.0_dp) .and. any(pvalues == 1.0_dp)) then
                result_value%statistic = quiet_nan(0.0_dp)
                result_value%pvalue = quiet_nan(0.0_dp)
                return
            else if (any(pvalues == 0.0_dp)) then
                result_value%statistic = positive_infinity(0.0_dp)
                result_value%pvalue = 0.0_dp
                return
            else if (any(pvalues == 1.0_dp)) then
                result_value%statistic = negative_infinity(0.0_dp)
                result_value%pvalue = 1.0_dp
                return
            end if
            sum_log_p = 0.0_dp
            sum_log_q = 0.0_dp
            do i = 1, n
                sum_log_p = sum_log_p + log(pvalues(i))
                sum_log_q = sum_log_q + log(1.0_dp - pvalues(i))
            end do
            result_value%statistic = -sum_log_p + sum_log_q
            normalizing_factor = sqrt(3.0_dp / real(n, dp)) / pi_dp
            nu = 5.0_dp * real(n, dp) + 4.0_dp
            approx_factor = sqrt(nu / (nu - 2.0_dp))
            result_value%pvalue = t_sf(result_value%statistic * normalizing_factor * approx_factor, nu)
        case ('tippett')
            result_value%statistic = minval(pvalues)
            result_value%pvalue = 1.0_dp - (1.0_dp - result_value%statistic) ** n
        case ('stouffer')
            if (present(weights)) then
                if (size(weights) /= n .or. any(ieee_is_nan(weights))) return
                sum_w2 = sum(weights * weights)
                if (sum_w2 <= 0.0_dp) return
                weighted_z = 0.0_dp
                do i = 1, n
                    weighted_z = weighted_z + weights(i) * normal_isf(pvalues(i))
                end do
                result_value%statistic = weighted_z / sqrt(sum_w2)
            else
                weighted_z = 0.0_dp
                do i = 1, n
                    weighted_z = weighted_z + normal_isf(pvalues(i))
                end do
                result_value%statistic = weighted_z / sqrt(real(n, dp))
            end if
            result_value%pvalue = normal_sf(result_value%statistic)
        case default
            return
        end select
        result_value%pvalue = min(1.0_dp, max(0.0_dp, result_value%pvalue))
    end function combine_pvalues

    function false_discovery_control(ps, method) result(adjusted)
        real(dp), intent(in) :: ps(:) !! p-values in [0,1]
        character(len=*), intent(in), optional :: method !! bh or by
        real(dp), allocatable :: adjusted(:)

        character(len=2) :: selected_method
        integer :: i
        integer :: j
        integer :: m
        integer :: temp_index
        integer, allocatable :: order(:)
        real(dp), allocatable :: sorted(:)
        real(dp) :: factor
        real(dp) :: harmonic

        selected_method = 'bh'
        if (present(method)) selected_method = adjustl(method)
        m = size(ps)
        allocate(adjusted(m))
        if (m == 0) return
        if (any(ps < 0.0_dp) .or. any(ps > 1.0_dp) .or. any(ieee_is_nan(ps)) .or. &
                (trim(selected_method) /= 'bh' .and. trim(selected_method) /= 'by')) then
            adjusted = quiet_nan(0.0_dp)
            return
        end if
        if (m == 1) then
            adjusted = ps
            return
        end if
        allocate(order(m), sorted(m))
        do i = 1, m
            order(i) = i
        end do
        do i = 2, m
            temp_index = order(i)
            j = i - 1
            do while (j >= 1)
                if (ps(order(j)) <= ps(temp_index)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = temp_index
        end do
        harmonic = 1.0_dp
        if (trim(selected_method) == 'by') then
            harmonic = 0.0_dp
            do i = 1, m
                harmonic = harmonic + 1.0_dp / real(i, dp)
            end do
        end if
        do i = 1, m
            factor = real(m, dp) / real(i, dp) * harmonic
            sorted(i) = ps(order(i)) * factor
        end do
        do i = m - 1, 1, -1
            sorted(i) = min(sorted(i), sorted(i + 1))
        end do
        do i = 1, m
            adjusted(order(i)) = min(1.0_dp, max(0.0_dp, sorted(i)))
        end do
    end function false_discovery_control

    pure function log_combination(n, k) result(value)
        integer, intent(in) :: n !! nonnegative total
        integer, intent(in) :: k !! selected count, 0 <= k <= n
        real(dp) :: value

        value = log_gamma(real(n + 1, dp)) - log_gamma(real(k + 1, dp)) - &
            log_gamma(real(n - k + 1, dp))
    end function log_combination

    pure integer function next_power_of_two(n) result(value)
        integer, intent(in) :: n !! positive integer

        value = 1
        do while (value < n)
            if (value > huge(value) - value) exit
            value = 2 * value
        end do
    end function next_power_of_two

    pure function choose_alternative(alternative) result(selected)
        character(len=*), intent(in), optional :: alternative !! optional alternative name
        character(len=9) :: selected

        selected = 'two-sided'
        if (present(alternative)) selected = adjustl(alternative)
    end function choose_alternative

    pure logical function valid_alternative(alternative) result(valid)
        character(len=*), intent(in) :: alternative !! candidate alternative name

        valid = trim(alternative) == 'two-sided' .or. trim(alternative) == 'less' .or. &
            trim(alternative) == 'greater'
    end function valid_alternative

end module scifort_contingency_meta
