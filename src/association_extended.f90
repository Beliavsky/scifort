! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Adapted portions: SciPy 1.17.0 scipy/stats/_mstats_basic.py (Kendall exact recurrence),
! scipy/stats/_stats.pyx (weighted Kendall tau), and _stats_pythran.py / _hypotests.py
! (Somers D contingency formulas). SciPy 1.17.0 spearmanrho public semantics were used
! as an API and numerical-parity reference; its implementation was not copied.
! Retained BSD terms: THIRD_PARTY_LICENSES.md; details: CODE_PROVENANCE.md.
!
! Association, robust-regression, and ordered-trend procedures following
! scipy.stats 1.17.0 scalar semantics where practical.

module scifort_association_extended
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_descriptive, only : covariance, mean, median, pearson_correlation, &
        rankdata, variance
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_cdf, normal_ppf, normal_sf
    use scifort_random, only : rng_state, rng_uniform
    use scifort_student_t, only : t_cdf, t_sf
    implicit none
    private

    integer(int64), parameter :: integer_grid_size = 4503599627370496_int64

    abstract interface
        function weightedtau_weigher(rank) result(weight)
            import :: dp
            integer, intent(in) :: rank !! zero-based nonnegative importance rank
            real(dp) :: weight !! nonnegative exchange weight contribution
        end function weightedtau_weigher
    end interface

    type, public :: association_result
        real(dp) :: statistic !! association or test statistic
        real(dp) :: pvalue !! p-value for the requested alternative
    end type association_result

    type, public :: somersd_result
        real(dp) :: statistic !! Somers D(Y|X) statistic
        real(dp) :: pvalue !! asymptotic p-value for the requested alternative
        real(dp) :: correlation !! alias of statistic, matching SciPy result semantics
        real(dp), allocatable :: table(:, :) !! contingency table used for the calculation
    end type somersd_result

    interface somersd
        module procedure somersd_samples
        module procedure somersd_table_real
        module procedure somersd_table_integer
    end interface somersd

    interface weightedtau
        module procedure weightedtau_default
        module procedure weightedtau_rank_integer
        module procedure weightedtau_rank_real
    end interface weightedtau

    type, public :: linregress_result
        real(dp) :: slope !! least-squares slope
        real(dp) :: intercept !! least-squares intercept
        real(dp) :: rvalue !! Pearson correlation coefficient
        real(dp) :: pvalue !! p-value for a zero-slope test
        real(dp) :: stderr !! standard error of the slope
        real(dp) :: intercept_stderr !! standard error of the intercept
    end type linregress_result

    type, public :: theilslopes_result
        real(dp) :: slope !! Theil-Sen median slope
        real(dp) :: intercept !! selected robust intercept estimate
        real(dp) :: low_slope !! lower confidence bound for the slope
        real(dp) :: high_slope !! upper confidence bound for the slope
    end type theilslopes_result

    type, public :: siegelslopes_result
        real(dp) :: slope !! Siegel repeated-median slope
        real(dp) :: intercept !! selected robust intercept estimate
    end type siegelslopes_result

    type, public :: page_trend_result
        real(dp) :: statistic !! Page L statistic
        real(dp) :: pvalue !! exact or asymptotic one-sided p-value
        character(len=10) :: method !! exact or asymptotic
    end type page_trend_result

    public :: brunnermunzel
    public :: chatterjeexi
    public :: kendalltau
    public :: linregress
    public :: page_trend_test
    public :: pointbiserialr
    public :: siegelslopes
    public :: somersd
    public :: spearmanrho
    public :: theilslopes
    public :: weightedtau
    public :: weightedtau_weigher

contains

    function somersd_samples(x, y, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! independent row-variable rankings
        real(dp), intent(in) :: y(:) !! dependent column-variable rankings
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(somersd_result) :: result_value

        character(len=9) :: selected_alternative

        call invalidate_somersd(result_value)
        if (size(x) /= size(y)) return
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return

        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        if (.not. valid_alternative(selected_alternative)) return

        call build_somersd_table(x, y, result_value%table)
        call somersd_from_table(result_value%table, selected_alternative, &
            result_value%statistic, result_value%pvalue)
        result_value%correlation = result_value%statistic
    end function somersd_samples

    function somersd_table_real(table, alternative) result(result_value)
        real(dp), intent(in) :: table(:, :) !! nonnegative integer-valued contingency table
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(somersd_result) :: result_value

        character(len=9) :: selected_alternative

        call invalidate_somersd(result_value)
        if (any(.not. ieee_is_finite(table))) return
        if (any(table < 0.0_dp)) return
        if (any(table /= floor(table))) return
        if (count(table /= 0.0_dp) < 2) return

        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        if (.not. valid_alternative(selected_alternative)) return

        allocate(result_value%table(size(table, 1), size(table, 2)))
        result_value%table = table
        call somersd_from_table(result_value%table, selected_alternative, &
            result_value%statistic, result_value%pvalue)
        result_value%correlation = result_value%statistic
    end function somersd_table_real

    function somersd_table_integer(table, alternative) result(result_value)
        integer, intent(in) :: table(:, :) !! nonnegative integer contingency table
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(somersd_result) :: result_value

        character(len=9) :: selected_alternative

        call invalidate_somersd(result_value)
        if (any(table < 0)) return
        if (count(table /= 0) < 2) return

        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        if (.not. valid_alternative(selected_alternative)) return

        allocate(result_value%table(size(table, 1), size(table, 2)))
        result_value%table = real(table, dp)
        call somersd_from_table(result_value%table, selected_alternative, &
            result_value%statistic, result_value%pvalue)
        result_value%correlation = result_value%statistic
    end function somersd_table_integer

    subroutine build_somersd_table(x, y, table)
        real(dp), intent(in) :: x(:) !! row-variable rankings
        real(dp), intent(in) :: y(:) !! column-variable rankings
        real(dp), allocatable, intent(out) :: table(:, :) !! sorted-category contingency table

        integer :: i
        integer :: ix
        integer :: iy
        real(dp), allocatable :: unique_x(:)
        real(dp), allocatable :: unique_y(:)

        call unique_sorted_values(x, unique_x)
        call unique_sorted_values(y, unique_y)
        allocate(table(size(unique_x), size(unique_y)))
        table = 0.0_dp
        do i = 1, size(x)
            ix = sorted_value_index(unique_x, x(i))
            iy = sorted_value_index(unique_y, y(i))
            if (ix > 0 .and. iy > 0) table(ix, iy) = table(ix, iy) + 1.0_dp
        end do
    end subroutine build_somersd_table

    subroutine unique_sorted_values(values, unique_values)
        real(dp), intent(in) :: values(:) !! finite values to sort and deduplicate
        real(dp), allocatable, intent(out) :: unique_values(:) !! sorted distinct values

        integer :: count_unique
        integer :: i
        integer :: j
        real(dp), allocatable :: work(:)

        if (size(values) == 0) then
            allocate(unique_values(0))
            return
        end if

        work = values
        call sort_real(work)
        count_unique = 1
        do i = 2, size(work)
            if (.not. same_value(work(i), work(i - 1))) count_unique = count_unique + 1
        end do

        allocate(unique_values(count_unique))
        unique_values(1) = work(1)
        j = 1
        do i = 2, size(work)
            if (.not. same_value(work(i), work(i - 1))) then
                j = j + 1
                unique_values(j) = work(i)
            end if
        end do
    end subroutine unique_sorted_values

    integer function sorted_value_index(values, target) result(index_value)
        real(dp), intent(in) :: values(:) !! sorted distinct finite values
        real(dp), intent(in) :: target !! value known to occur in values

        integer :: high
        integer :: low
        integer :: middle

        index_value = 0
        low = 1
        high = size(values)
        do while (low <= high)
            middle = low + (high - low) / 2
            if (same_value(values(middle), target)) then
                index_value = middle
                return
            else if (values(middle) < target) then
                low = middle + 1
            else
                high = middle - 1
            end if
        end do
    end function sorted_value_index

    subroutine somersd_from_table(table, alternative, statistic, pvalue)
        real(dp), intent(in) :: table(:, :) !! nonnegative contingency table
        character(len=*), intent(in) :: alternative !! valid alternative
        real(dp), intent(out) :: statistic !! Somers D(Y|X)
        real(dp), intent(out) :: pvalue !! asymptotic p-value

        real(dp) :: aij
        real(dp) :: dij
        real(dp) :: difference
        real(dp) :: denominator
        integer :: i
        integer :: j
        integer :: ncol
        integer :: nrow
        real(dp) :: n_total
        real(dp) :: p_concordant
        real(dp), allocatable :: prefix(:, :)
        real(dp) :: q_discordant
        real(dp) :: row_square_sum
        real(dp) :: s_term
        real(dp) :: z

        statistic = quiet_nan(0.0_dp)
        pvalue = quiet_nan(0.0_dp)
        nrow = size(table, 1)
        ncol = size(table, 2)
        if (nrow <= 1 .or. ncol <= 1) return

        allocate(prefix(0:nrow, 0:ncol))
        prefix = 0.0_dp
        do i = 1, nrow
            do j = 1, ncol
                prefix(i, j) = table(i, j) + prefix(i - 1, j) + &
                    prefix(i, j - 1) - prefix(i - 1, j - 1)
            end do
        end do

        n_total = prefix(nrow, ncol)
        if (n_total <= 0.0_dp) return
        p_concordant = 0.0_dp
        q_discordant = 0.0_dp
        s_term = 0.0_dp
        do i = 1, nrow
            do j = 1, ncol
                aij = prefix(i - 1, j - 1) + prefix(nrow, ncol) - &
                    prefix(i, ncol) - prefix(nrow, j) + prefix(i, j)
                dij = prefix(nrow, j - 1) - prefix(i, j - 1) + &
                    prefix(i - 1, ncol) - prefix(i - 1, j)
                p_concordant = p_concordant + table(i, j) * aij
                q_discordant = q_discordant + table(i, j) * dij
                s_term = s_term + table(i, j) * (aij - dij) ** 2
            end do
        end do

        row_square_sum = 0.0_dp
        do i = 1, nrow
            row_square_sum = row_square_sum + sum(table(i, :)) ** 2
        end do
        denominator = n_total ** 2 - row_square_sum
        if (denominator <= 0.0_dp) return

        difference = p_concordant - q_discordant
        statistic = difference / denominator
        s_term = s_term - difference ** 2 / n_total
        if (s_term < 0.0_dp) then
            pvalue = quiet_nan(0.0_dp)
        else if (s_term == 0.0_dp) then
            if (difference > 0.0_dp) then
                select case (alternative)
                case ('less')
                    pvalue = 1.0_dp
                case default
                    pvalue = 0.0_dp
                end select
            else if (difference < 0.0_dp) then
                select case (alternative)
                case ('greater')
                    pvalue = 1.0_dp
                case default
                    pvalue = 0.0_dp
                end select
            else
                pvalue = quiet_nan(0.0_dp)
            end if
        else
            z = difference / sqrt(4.0_dp * s_term)
            pvalue = normal_pvalue(z, alternative)
        end if
    end subroutine somersd_from_table

    function spearmanrho(x, y, alternative, method, state, n_resamples) result(result_value)
        real(dp), intent(in) :: x(:) !! first observations; at least two non-NaN values
        real(dp), intent(in) :: y(:) !! second observations, same size as x
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        character(len=*), intent(in), optional :: method !! asymptotic, permutation, or monte_carlo
        type(rng_state), intent(inout), optional :: state !! RNG state for randomized inference
        integer, intent(in), optional :: n_resamples !! maximum random draws, default 9999
        type(association_result) :: result_value

        integer :: draws
        real(dp) :: degrees_of_freedom
        real(dp), allocatable :: ranks_x(:)
        real(dp), allocatable :: ranks_y(:)
        character(len=11) :: selected_method
        character(len=9) :: selected_alternative
        real(dp) :: tstat

        call invalidate_association(result_value)
        if (size(x) /= size(y) .or. size(x) < 2) return
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return

        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        selected_method = 'asymptotic'
        if (present(method)) selected_method = adjustl(method)
        if (.not. valid_alternative(selected_alternative)) return
        if (selected_method /= 'asymptotic' .and. selected_method /= 'permutation' .and. &
                selected_method /= 'monte_carlo') return

        ranks_x = rankdata(x)
        ranks_y = rankdata(y)
        if (any(ieee_is_nan(ranks_x)) .or. any(ieee_is_nan(ranks_y))) return
        result_value%statistic = pearson_correlation(ranks_x, ranks_y)
        if (ieee_is_nan(result_value%statistic)) return

        if (selected_method == 'asymptotic') then
            if (size(x) == 2) then
                result_value%pvalue = 1.0_dp
                return
            end if
            if (result_value%statistic >= 1.0_dp) then
                select case (selected_alternative)
                case ('less')
                    result_value%pvalue = 1.0_dp
                case default
                    result_value%pvalue = 0.0_dp
                end select
                return
            else if (result_value%statistic <= -1.0_dp) then
                select case (selected_alternative)
                case ('greater')
                    result_value%pvalue = 1.0_dp
                case default
                    result_value%pvalue = 0.0_dp
                end select
                return
            end if
            degrees_of_freedom = real(size(x) - 2, dp)
            tstat = result_value%statistic * sqrt(degrees_of_freedom / &
                ((1.0_dp + result_value%statistic) * (1.0_dp - result_value%statistic)))
            result_value%pvalue = student_t_pvalue(tstat, degrees_of_freedom, selected_alternative)
            return
        end if

        draws = 9999
        if (present(n_resamples)) draws = n_resamples
        if (draws < 1) then
            result_value%pvalue = quiet_nan(0.0_dp)
            return
        end if

        if (selected_method == 'permutation') then
            if (present(state)) then
                result_value%pvalue = spearmanrho_permutation_pvalue(&
                    ranks_x, ranks_y, result_value%statistic, selected_alternative, draws, state)
            else
                result_value%pvalue = spearmanrho_permutation_pvalue(&
                    ranks_x, ranks_y, result_value%statistic, selected_alternative, draws)
            end if
        else
            if (.not. present(state)) then
                result_value%pvalue = quiet_nan(0.0_dp)
                return
            end if
            result_value%pvalue = spearmanrho_monte_carlo_pvalue(&
                size(x), result_value%statistic, selected_alternative, draws, state)
        end if
    end function spearmanrho

    function spearmanrho_permutation_pvalue(ranks_x, ranks_y, observed, alternative, draws, state) &
            result(pvalue)
        real(dp), intent(in) :: ranks_x(:) !! fixed ranks of the first sample
        real(dp), intent(in) :: ranks_y(:) !! ranks of the second sample to re-pair
        real(dp), intent(in) :: observed !! observed Spearman correlation
        character(len=*), intent(in) :: alternative !! valid alternative hypothesis
        integer, intent(in) :: draws !! requested maximum number of permutations
        type(rng_state), intent(inout), optional :: state !! RNG state for randomized pairings
        real(dp) :: pvalue

        integer :: b
        integer :: count_greater
        integer :: count_less
        real(dp) :: gamma
        integer :: i
        integer(int64) :: n_exact
        integer, allocatable :: permutation(:)
        real(dp) :: p_greater
        real(dp) :: p_less
        integer :: selected
        real(dp) :: statistic
        real(dp) :: temporary
        real(dp), allocatable :: work_y(:)

        pvalue = quiet_nan(0.0_dp)
        n_exact = factorial_count_capped(size(ranks_y), int(draws, int64) + 1_int64)
        gamma = abs(100.0_dp * epsilon(1.0_dp) * observed)
        count_less = 0
        count_greater = 0
        allocate(work_y(size(ranks_y)))

        if (n_exact <= int(draws, int64)) then
            allocate(permutation(size(ranks_y)))
            do i = 1, size(ranks_y)
                permutation(i) = i
            end do
            do
                work_y = ranks_y(permutation)
                statistic = pearson_correlation(ranks_x, work_y)
                if (ieee_is_nan(statistic)) return
                if (statistic <= observed + gamma) count_less = count_less + 1
                if (statistic >= observed - gamma) count_greater = count_greater + 1
                if (.not. next_integer_permutation(permutation)) exit
            end do
            p_less = real(count_less, dp) / real(n_exact, dp)
            p_greater = real(count_greater, dp) / real(n_exact, dp)
            pvalue = association_tail_pvalue(p_less, p_greater, alternative)
            return
        end if

        if (.not. present(state)) return
        do b = 1, draws
            work_y = ranks_y
            do i = size(work_y), 2, -1
                selected = unbiased_index(state, i)
                temporary = work_y(i)
                work_y(i) = work_y(selected)
                work_y(selected) = temporary
            end do
            statistic = pearson_correlation(ranks_x, work_y)
            if (ieee_is_nan(statistic)) return
            if (statistic <= observed + gamma) count_less = count_less + 1
            if (statistic >= observed - gamma) count_greater = count_greater + 1
        end do
        p_less = real(count_less + 1, dp) / real(draws + 1, dp)
        p_greater = real(count_greater + 1, dp) / real(draws + 1, dp)
        pvalue = association_tail_pvalue(p_less, p_greater, alternative)
    end function spearmanrho_permutation_pvalue

    function spearmanrho_monte_carlo_pvalue(n, observed, alternative, draws, state) &
            result(pvalue)
        integer, intent(in) :: n !! sample size, at least two
        real(dp), intent(in) :: observed !! observed Spearman correlation
        character(len=*), intent(in) :: alternative !! valid alternative hypothesis
        integer, intent(in) :: draws !! number of null simulations
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by simulation
        real(dp) :: pvalue

        integer :: b
        integer :: count_greater
        integer :: count_less
        real(dp), allocatable :: fixed_ranks(:)
        real(dp) :: gamma
        integer :: i
        real(dp) :: p_greater
        real(dp) :: p_less
        integer :: selected
        real(dp) :: statistic
        real(dp) :: temporary
        real(dp), allocatable :: work_ranks(:)

        pvalue = quiet_nan(0.0_dp)
        if (n < 2 .or. draws < 1) return
        allocate(fixed_ranks(n), work_ranks(n))
        do i = 1, n
            fixed_ranks(i) = real(i, dp)
        end do
        gamma = abs(100.0_dp * epsilon(1.0_dp) * observed)
        count_less = 0
        count_greater = 0

        do b = 1, draws
            work_ranks = fixed_ranks
            do i = n, 2, -1
                selected = unbiased_index(state, i)
                temporary = work_ranks(i)
                work_ranks(i) = work_ranks(selected)
                work_ranks(selected) = temporary
            end do
            statistic = pearson_correlation(fixed_ranks, work_ranks)
            if (ieee_is_nan(statistic)) return
            if (statistic <= observed + gamma) count_less = count_less + 1
            if (statistic >= observed - gamma) count_greater = count_greater + 1
        end do
        p_less = real(count_less + 1, dp) / real(draws + 1, dp)
        p_greater = real(count_greater + 1, dp) / real(draws + 1, dp)
        pvalue = association_tail_pvalue(p_less, p_greater, alternative)
    end function spearmanrho_monte_carlo_pvalue

    pure function association_tail_pvalue(p_less, p_greater, alternative) result(pvalue)
        real(dp), intent(in) :: p_less !! lower-tail probability
        real(dp), intent(in) :: p_greater !! upper-tail probability
        character(len=*), intent(in) :: alternative !! valid alternative hypothesis
        real(dp) :: pvalue

        select case (alternative)
        case ('less')
            pvalue = p_less
        case ('greater')
            pvalue = p_greater
        case default
            pvalue = min(1.0_dp, 2.0_dp * min(p_less, p_greater))
        end select
    end function association_tail_pvalue

    function chatterjeexi(x, y, y_continuous, method, state, n_resamples) result(result_value)
        real(dp), intent(in) :: x(:) !! independent-variable observations
        real(dp), intent(in) :: y(:) !! dependent-variable observations, same size as x
        logical, intent(in), optional :: y_continuous !! assume the distribution of y is continuous
        character(len=*), intent(in), optional :: method !! asymptotic or permutation
        type(rng_state), intent(inout), optional :: state !! RNG state required for permutation inference
        integer, intent(in), optional :: n_resamples !! maximum permutation draws, default 9999
        type(association_result) :: result_value

        integer :: draws
        character(len=11) :: selected_method
        logical :: continuous
        real(dp), allocatable :: l(:)
        real(dp), allocatable :: r(:)
        real(dp) :: standard_error

        call invalidate_association(result_value)
        if (size(x) /= size(y) .or. size(x) < 2) return
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return

        continuous = .false.
        if (present(y_continuous)) continuous = y_continuous
        selected_method = 'asymptotic'
        if (present(method)) selected_method = adjustl(method)
        if (selected_method /= 'asymptotic' .and. selected_method /= 'permutation') return

        call chatterjee_statistic(x, y, continuous, result_value%statistic, r, l)
        if (ieee_is_nan(result_value%statistic)) return

        if (selected_method == 'asymptotic') then
            standard_error = chatterjee_standard_error(r, l, continuous)
            if (.not. ieee_is_finite(standard_error) .or. standard_error <= 0.0_dp) return
            result_value%pvalue = normal_sf(result_value%statistic / standard_error)
            return
        end if

        draws = 9999
        if (present(n_resamples)) draws = n_resamples
        if (draws < 1 .or. .not. present(state)) then
            result_value%pvalue = quiet_nan(0.0_dp)
            return
        end if
        result_value%pvalue = chatterjee_permutation_pvalue(&
            x, y, continuous, result_value%statistic, state, draws)
    end function chatterjeexi

    subroutine chatterjee_statistic(x, y, y_continuous, statistic, r, l)
        real(dp), intent(in) :: x(:) !! independent observations, same length as y
        real(dp), intent(in) :: y(:) !! dependent observations
        logical, intent(in) :: y_continuous !! use the continuous-y normalization
        real(dp), intent(out) :: statistic !! Chatterjee xi statistic
        real(dp), allocatable, intent(out) :: r(:) !! max ranks of y after sorting by x
        real(dp), allocatable, intent(out) :: l(:) !! max ranks of -y after sorting by x

        real(dp) :: denominator
        integer :: i
        integer :: n
        real(dp) :: numerator
        integer, allocatable :: order(:)
        real(dp), allocatable :: sorted_y(:)

        n = size(x)
        statistic = quiet_nan(0.0_dp)
        allocate(r(n), l(n))
        if (n < 2 .or. size(y) /= n) then
            r = quiet_nan(0.0_dp)
            l = quiet_nan(0.0_dp)
            return
        end if

        allocate(order(n), sorted_y(n))
        do i = 1, n
            order(i) = i
        end do
        call stable_sort_indices(x, order)
        sorted_y = y(order)
        r = rankdata(sorted_y, 'max')
        l = rankdata(-sorted_y, 'max')
        if (any(ieee_is_nan(r)) .or. any(ieee_is_nan(l))) return

        numerator = 0.0_dp
        do i = 1, n - 1
            numerator = numerator + abs(r(i + 1) - r(i))
        end do

        if (y_continuous) then
            denominator = real(n, dp) ** 2 - 1.0_dp
            statistic = 1.0_dp - 3.0_dp * numerator / denominator
        else
            denominator = 2.0_dp * sum((real(n, dp) - l) * l)
            if (denominator <= 0.0_dp) return
            statistic = 1.0_dp - real(n, dp) * numerator / denominator
        end if
    end subroutine chatterjee_statistic

    function chatterjee_standard_error(r, l, y_continuous) result(standard_error)
        real(dp), intent(in) :: r(:) !! max ranks of y after sorting by x
        real(dp), intent(in) :: l(:) !! max ranks of -y after sorting by x
        logical, intent(in) :: y_continuous !! use the continuous-y null variance
        real(dp) :: standard_error

        real(dp) :: an
        real(dp) :: bn
        real(dp) :: cn
        real(dp) :: cumulative
        real(dp) :: dn
        integer :: i
        integer :: n
        real(dp) :: nreal
        real(dp) :: tau2
        real(dp), allocatable :: u(:)

        n = size(r)
        standard_error = quiet_nan(0.0_dp)
        if (n < 2 .or. size(l) /= n) return
        nreal = real(n, dp)
        if (y_continuous) then
            standard_error = sqrt(2.0_dp / 5.0_dp) / sqrt(nreal)
            return
        end if

        u = r
        call sort_real(u)
        an = 0.0_dp
        bn = 0.0_dp
        cn = 0.0_dp
        cumulative = 0.0_dp
        do i = 1, n
            cumulative = cumulative + u(i)
            an = an + real(2 * n - 2 * i + 1, dp) * u(i) ** 2
            bn = bn + (cumulative + real(n - i, dp) * u(i)) ** 2
            cn = cn + real(2 * n - 2 * i + 1, dp) * u(i)
        end do
        an = an / nreal ** 4
        bn = bn / nreal ** 5
        cn = cn / nreal ** 3
        dn = sum(l * (nreal - l)) / nreal ** 3
        if (dn <= 0.0_dp) return
        tau2 = (an - 2.0_dp * bn + cn ** 2) / dn ** 2
        if (.not. ieee_is_finite(tau2) .or. tau2 <= 0.0_dp) return
        standard_error = sqrt(tau2) / sqrt(nreal)
    end function chatterjee_standard_error

    function chatterjee_permutation_pvalue(x, y, y_continuous, observed, state, draws) &
            result(pvalue)
        real(dp), intent(in) :: x(:) !! fixed independent observations
        real(dp), intent(in) :: y(:) !! dependent observations to re-pair
        logical, intent(in) :: y_continuous !! statistic normalization
        real(dp), intent(in) :: observed !! observed xi statistic
        type(rng_state), intent(inout) :: state !! explicit RNG state for randomized pairings
        integer, intent(in) :: draws !! requested maximum number of permutations
        real(dp) :: pvalue

        integer :: b
        integer :: count_greater
        real(dp) :: gamma
        integer :: i
        integer(int64) :: n_exact
        integer, allocatable :: permutation(:)
        integer :: selected
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: l(:)
        real(dp) :: statistic
        real(dp) :: temporary
        real(dp), allocatable :: work_y(:)

        pvalue = quiet_nan(0.0_dp)
        n_exact = factorial_count_capped(size(y), int(draws, int64) + 1_int64)
        gamma = abs(100.0_dp * epsilon(1.0_dp) * observed)
        count_greater = 0

        allocate(work_y(size(y)))
        if (n_exact <= int(draws, int64)) then
            allocate(permutation(size(y)))
            do i = 1, size(y)
                permutation(i) = i
            end do
            do
                work_y = y(permutation)
                call chatterjee_statistic(x, work_y, y_continuous, statistic, r, l)
                if (ieee_is_nan(statistic)) return
                if (statistic >= observed - gamma) count_greater = count_greater + 1
                if (.not. next_integer_permutation(permutation)) exit
            end do
            pvalue = real(count_greater, dp) / real(n_exact, dp)
            return
        end if

        do b = 1, draws
            work_y = y
            do i = size(work_y), 2, -1
                selected = unbiased_index(state, i)
                temporary = work_y(i)
                work_y(i) = work_y(selected)
                work_y(selected) = temporary
            end do
            call chatterjee_statistic(x, work_y, y_continuous, statistic, r, l)
            if (ieee_is_nan(statistic)) return
            if (statistic >= observed - gamma) count_greater = count_greater + 1
        end do
        pvalue = real(count_greater + 1, dp) / real(draws + 1, dp)
    end function chatterjee_permutation_pvalue

    function weightedtau_default(x, y, weigher, additive, rank_mode) result(result_value)
        real(dp), intent(in) :: x(:) !! first score array
        real(dp), intent(in) :: y(:) !! second score array, same size as x
        procedure(weightedtau_weigher), optional :: weigher !! rank-to-weight callback
        logical, intent(in), optional :: additive !! add rather than multiply pair weights
        character(len=*), intent(in), optional :: rank_mode !! average, lexicographic, or index
        type(association_result) :: result_value

        logical :: use_additive
        character(len=13) :: selected_rank_mode
        integer :: i
        integer, allocatable :: rank(:)
        real(dp) :: tau_first
        real(dp) :: tau_second

        call invalidate_association(result_value)
        if (size(x) /= size(y) .or. size(x) == 0) return

        use_additive = .true.
        if (present(additive)) use_additive = additive
        selected_rank_mode = 'average'
        if (present(rank_mode)) selected_rank_mode = adjustl(rank_mode)

        select case (selected_rank_mode)
        case ('average')
            call decreasing_lexicographic_rank(x, y, rank)
            tau_first = weighted_tau_core(x, y, rank, use_additive, weigher)
            call decreasing_lexicographic_rank(y, x, rank)
            tau_second = weighted_tau_core(y, x, rank, use_additive, weigher)
            if (ieee_is_nan(tau_first) .or. ieee_is_nan(tau_second)) return
            result_value%statistic = 0.5_dp * (tau_first + tau_second)
        case ('lexicographic')
            call decreasing_lexicographic_rank(x, y, rank)
            result_value%statistic = weighted_tau_core(x, y, rank, use_additive, weigher)
        case ('index')
            allocate(rank(size(x)))
            rank = [(i - 1, i = 1, size(x))]
            result_value%statistic = weighted_tau_core(x, y, rank, use_additive, weigher)
        case default
            return
        end select
        if (.not. ieee_is_nan(result_value%statistic)) then
            result_value%statistic = max(-1.0_dp, min(1.0_dp, result_value%statistic))
        end if
    end function weightedtau_default

    function weightedtau_rank_integer(x, y, rank, weigher, additive) result(result_value)
        real(dp), intent(in) :: x(:) !! first score array
        real(dp), intent(in) :: y(:) !! second score array, same size as x
        integer, intent(in) :: rank(:) !! external importance scores, smaller values more important
        procedure(weightedtau_weigher), optional :: weigher !! rank-to-weight callback
        logical, intent(in), optional :: additive !! add rather than multiply pair weights
        type(association_result) :: result_value

        logical :: use_additive
        integer, allocatable :: dense_rank(:)

        call invalidate_association(result_value)
        if (size(x) /= size(y) .or. size(x) /= size(rank) .or. size(x) == 0) return
        use_additive = .true.
        if (present(additive)) use_additive = additive
        call dense_rank_integer(rank, dense_rank)
        result_value%statistic = weighted_tau_core(x, y, dense_rank, use_additive, weigher)
        if (.not. ieee_is_nan(result_value%statistic)) then
            result_value%statistic = max(-1.0_dp, min(1.0_dp, result_value%statistic))
        end if
    end function weightedtau_rank_integer

    function weightedtau_rank_real(x, y, rank, weigher, additive) result(result_value)
        real(dp), intent(in) :: x(:) !! first score array
        real(dp), intent(in) :: y(:) !! second score array, same size as x
        real(dp), intent(in) :: rank(:) !! external importance scores, smaller values more important
        procedure(weightedtau_weigher), optional :: weigher !! rank-to-weight callback
        logical, intent(in), optional :: additive !! add rather than multiply pair weights
        type(association_result) :: result_value

        logical :: use_additive
        integer, allocatable :: dense_rank(:)

        call invalidate_association(result_value)
        if (size(x) /= size(y) .or. size(x) /= size(rank) .or. size(x) == 0) return
        use_additive = .true.
        if (present(additive)) use_additive = additive
        call dense_rank_real(rank, dense_rank)
        result_value%statistic = weighted_tau_core(x, y, dense_rank, use_additive, weigher)
        if (.not. ieee_is_nan(result_value%statistic)) then
            result_value%statistic = max(-1.0_dp, min(1.0_dp, result_value%statistic))
        end if
    end function weightedtau_rank_real

    function weighted_tau_core(x, y, rank, additive, weigher) result(tau)
        real(dp), intent(in) :: x(:) !! first score array
        real(dp), intent(in) :: y(:) !! second score array
        integer, intent(in) :: rank(:) !! zero-based nonnegative importance ranks
        logical, intent(in) :: additive !! additive or multiplicative exchange weights
        procedure(weightedtau_weigher), optional :: weigher !! rank-to-weight callback
        real(dp) :: tau

        real(dp) :: discordant_weight
        real(dp) :: joint_tie_weight
        integer :: i
        integer :: j
        real(dp) :: non_tied_weight
        real(dp) :: pair_weight
        real(dp) :: total_weight
        real(dp) :: tie_x_weight
        real(dp) :: tie_y_weight
        real(dp), allocatable :: weights(:)

        tau = quiet_nan(0.0_dp)
        if (size(x) /= size(y) .or. size(x) /= size(rank) .or. size(x) == 0) return
        if (any(rank < 0)) return

        allocate(weights(size(rank)))
        do i = 1, size(rank)
            if (present(weigher)) then
                weights(i) = weigher(rank(i))
            else
                weights(i) = 1.0_dp / real(rank(i) + 1, dp)
            end if
            if (.not. ieee_is_finite(weights(i)) .or. weights(i) < 0.0_dp) return
        end do

        total_weight = 0.0_dp
        tie_x_weight = 0.0_dp
        tie_y_weight = 0.0_dp
        joint_tie_weight = 0.0_dp
        discordant_weight = 0.0_dp
        do i = 1, size(x) - 1
            do j = i + 1, size(x)
                if (additive) then
                    if (weights(i) > huge(1.0_dp) - weights(j)) return
                    pair_weight = weights(i) + weights(j)
                else
                    if (weights(i) > 1.0_dp) then
                        if (weights(j) > huge(1.0_dp) / weights(i)) return
                    end if
                    pair_weight = weights(i) * weights(j)
                end if
                if (pair_weight > huge(1.0_dp) - total_weight) return
                total_weight = total_weight + pair_weight

                if (score_equal(x(i), x(j))) then
                    if (pair_weight > huge(1.0_dp) - tie_x_weight) return
                    tie_x_weight = tie_x_weight + pair_weight
                end if
                if (score_equal(y(i), y(j))) then
                    if (pair_weight > huge(1.0_dp) - tie_y_weight) return
                    tie_y_weight = tie_y_weight + pair_weight
                end if
                if (score_equal(x(i), x(j)) .and. score_equal(y(i), y(j))) then
                    if (pair_weight > huge(1.0_dp) - joint_tie_weight) return
                    joint_tie_weight = joint_tie_weight + pair_weight
                else if (.not. score_equal(x(i), x(j)) .and. .not. score_equal(y(i), y(j))) then
                    if (score_less(x(i), x(j)) .neqv. score_less(y(i), y(j))) then
                        if (pair_weight > huge(1.0_dp) - discordant_weight) return
                        discordant_weight = discordant_weight + pair_weight
                    end if
                end if
            end do
        end do

        if (total_weight <= tie_x_weight .or. total_weight <= tie_y_weight) return
        non_tied_weight = total_weight - tie_x_weight - tie_y_weight + joint_tie_weight
        tau = ((non_tied_weight - discordant_weight) - discordant_weight) / &
            sqrt(total_weight - tie_x_weight) / sqrt(total_weight - tie_y_weight)
    end function weighted_tau_core

    subroutine decreasing_lexicographic_rank(x, y, rank)
        real(dp), intent(in) :: x(:) !! primary score, sorted ascending internally
        real(dp), intent(in) :: y(:) !! tie-breaking score, sorted ascending internally
        integer, allocatable, intent(out) :: rank(:) !! zero-based decreasing lexicographic ranks

        integer :: i
        integer, allocatable :: order(:)

        allocate(order(size(x)), rank(size(x)))
        order = [(i, i = 1, size(x))]
        call stable_lexicographic_index_sort(x, y, order)
        do i = 1, size(x)
            rank(order(size(x) - i + 1)) = i - 1
        end do
    end subroutine decreasing_lexicographic_rank

    subroutine dense_rank_integer(values, rank)
        integer, intent(in) :: values(:) !! arbitrary integer importance scores
        integer, allocatable, intent(out) :: rank(:) !! zero-based dense ranks

        integer :: i
        integer :: j
        integer :: key
        integer, allocatable :: order(:)

        allocate(order(size(values)), rank(size(values)))
        order = [(i, i = 1, size(values))]
        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (values(order(j)) <= values(key)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
        if (size(values) == 0) return
        rank(order(1)) = 0
        do i = 2, size(values)
            rank(order(i)) = rank(order(i - 1))
            if (values(order(i)) /= values(order(i - 1))) rank(order(i)) = rank(order(i)) + 1
        end do
    end subroutine dense_rank_integer

    subroutine dense_rank_real(values, rank)
        real(dp), intent(in) :: values(:) !! arbitrary real importance scores; NaN is smallest
        integer, allocatable, intent(out) :: rank(:) !! zero-based dense ranks

        integer :: i
        integer :: j
        integer :: key
        integer, allocatable :: order(:)

        allocate(order(size(values)), rank(size(values)))
        order = [(i, i = 1, size(values))]
        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (.not. score_less(values(key), values(order(j)))) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
        if (size(values) == 0) return
        rank(order(1)) = 0
        do i = 2, size(values)
            rank(order(i)) = rank(order(i - 1))
            if (.not. score_equal(values(order(i)), values(order(i - 1)))) then
                rank(order(i)) = rank(order(i)) + 1
            end if
        end do
    end subroutine dense_rank_real

    subroutine stable_lexicographic_index_sort(x, y, order)
        real(dp), intent(in) :: x(:) !! primary sort keys
        real(dp), intent(in) :: y(:) !! secondary sort keys
        integer, intent(inout) :: order(:) !! indices sorted stably in ascending lexicographic order

        integer :: i
        integer :: j
        integer :: key

        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (.not. lexicographic_less(key, order(j), x, y)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
    end subroutine stable_lexicographic_index_sort

    pure logical function lexicographic_less(i, j, x, y) result(less)
        integer, intent(in) :: i !! first element index
        integer, intent(in) :: j !! second element index
        real(dp), intent(in) :: x(:) !! primary keys
        real(dp), intent(in) :: y(:) !! secondary keys

        if (score_less(x(i), x(j))) then
            less = .true.
        else if (score_less(x(j), x(i))) then
            less = .false.
        else
            less = score_less(y(i), y(j))
        end if
    end function lexicographic_less

    pure logical function score_equal(x, y) result(equal)
        real(dp), intent(in) :: x !! first score
        real(dp), intent(in) :: y !! second score

        if (ieee_is_nan(x)) then
            equal = ieee_is_nan(y)
        else if (ieee_is_nan(y)) then
            equal = .false.
        else
            equal = x <= y .and. y <= x
        end if
    end function score_equal

    pure logical function score_less(x, y) result(less)
        real(dp), intent(in) :: x !! first score
        real(dp), intent(in) :: y !! second score

        if (ieee_is_nan(x)) then
            less = .not. ieee_is_nan(y)
        else if (ieee_is_nan(y)) then
            less = .false.
        else
            less = x < y
        end if
    end function score_less

    function pointbiserialr(x, y) result(result_value)
        real(dp), intent(in) :: x(:) !! binary-coded observations; SciPy treats these as numeric
        real(dp), intent(in) :: y(:) !! continuous observations, same size as x
        type(association_result) :: result_value

        real(dp) :: r
        real(dp) :: tstat
        real(dp) :: df

        call invalidate_association(result_value)
        if (size(x) /= size(y) .or. size(x) < 2) return
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return

        r = pearson_correlation(x, y)
        result_value%statistic = r
        if (ieee_is_nan(r)) return
        if (size(x) == 2) then
            if (abs(r) >= 1.0_dp) then
                result_value%pvalue = 1.0_dp
            else
                result_value%pvalue = quiet_nan(0.0_dp)
            end if
            return
        end if

        df = real(size(x) - 2, dp)
        if (abs(r) >= 1.0_dp) then
            result_value%pvalue = 0.0_dp
        else
            tstat = r * sqrt(df / max(tiny(1.0_dp), (1.0_dp - r) * (1.0_dp + r)))
            result_value%pvalue = min(1.0_dp, 2.0_dp * t_sf(abs(tstat), df))
        end if
    end function pointbiserialr

    function kendalltau(x, y, method, variant, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first ordinal observations
        real(dp), intent(in) :: y(:) !! second ordinal observations, same size as x
        character(len=*), intent(in), optional :: method !! auto, exact, or asymptotic
        character(len=*), intent(in), optional :: variant !! b or c
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(association_result) :: result_value

        character(len=10) :: selected_method
        character(len=9) :: selected_alternative
        character(len=1) :: selected_variant
        integer(int64) :: concordant
        integer(int64) :: discordant
        integer :: i
        integer :: j
        integer :: minclasses
        integer :: n
        integer(int64) :: ntie
        integer(int64) :: tot
        integer(int64) :: x0
        integer(int64) :: x1
        integer(int64) :: xtie
        integer(int64) :: y0
        integer(int64) :: y1
        integer(int64) :: ytie
        real(dp) :: con_minus_dis
        real(dp) :: m
        real(dp) :: var_s
        real(dp) :: z

        call invalidate_association(result_value)
        if (size(x) /= size(y) .or. size(x) < 2) return
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return

        selected_method = 'auto'
        if (present(method)) selected_method = adjustl(method)
        selected_variant = 'b'
        if (present(variant)) selected_variant = variant(1:1)
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        if (.not. valid_alternative(selected_alternative)) return
        if (selected_method /= 'auto' .and. selected_method /= 'exact' .and. &
                selected_method /= 'asymptotic') return
        if (selected_variant /= 'b' .and. selected_variant /= 'c') return

        n = size(x)
        concordant = 0
        discordant = 0
        ntie = 0
        do i = 1, n - 1
            do j = i + 1, n
                if (same_value(x(i), x(j)) .and. same_value(y(i), y(j))) then
                    ntie = ntie + 1
                else if (.not. same_value(x(i), x(j)) .and. &
                        .not. same_value(y(i), y(j))) then
                    if ((x(i) - x(j)) * (y(i) - y(j)) > 0.0_dp) then
                        concordant = concordant + 1
                    else
                        discordant = discordant + 1
                    end if
                end if
            end do
        end do

        call tie_statistics(x, xtie, x0, x1, minclasses)
        call tie_statistics(y, ytie, y0, y1, i)
        minclasses = min(minclasses, i)
        tot = int(n, int64) * int(n - 1, int64) / 2_int64
        if (xtie == tot .or. ytie == tot) return

        con_minus_dis = real(concordant - discordant, dp)
        if (selected_variant == 'b') then
            result_value%statistic = con_minus_dis / &
                sqrt(real(tot - xtie, dp) * real(tot - ytie, dp))
        else
            if (minclasses <= 1) return
            result_value%statistic = 2.0_dp * con_minus_dis / &
                (real(n, dp) ** 2 * real(minclasses - 1, dp) / real(minclasses, dp))
        end if
        result_value%statistic = max(-1.0_dp, min(1.0_dp, result_value%statistic))

        if (selected_method == 'exact' .and. (xtie /= 0 .or. ytie /= 0)) return
        if (selected_method == 'auto') then
            if (xtie == 0 .and. ytie == 0 .and. &
                    (n <= 33 .or. min(discordant, tot - discordant) <= 1_int64)) then
                selected_method = 'exact'
            else
                selected_method = 'asymptotic'
            end if
        end if

        if (selected_method == 'exact') then
            if (concordant > int(huge(0), int64)) return
            result_value%pvalue = kendall_p_exact(n, int(concordant), selected_alternative)
        else
            m = real(n, dp) * real(n - 1, dp)
            var_s = (m * real(2 * n + 5, dp) - real(x1 + y1, dp)) / 18.0_dp
            var_s = var_s + 2.0_dp * real(xtie, dp) * real(ytie, dp) / m
            if (n > 2) then
                var_s = var_s + real(x0, dp) * real(y0, dp) / &
                    (9.0_dp * m * real(n - 2, dp))
            end if
            if (var_s <= 0.0_dp) then
                result_value%pvalue = quiet_nan(0.0_dp)
            else
                z = con_minus_dis / sqrt(var_s)
                result_value%pvalue = normal_pvalue(z, selected_alternative)
            end if
        end if
    end function kendalltau

    function linregress(x, y, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! independent variable
        real(dp), intent(in) :: y(:) !! dependent variable, same size as x
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(linregress_result) :: result_value

        character(len=9) :: selected_alternative
        integer :: n
        real(dp), parameter :: tiny_r = 1.0e-20_dp
        real(dp) :: df
        real(dp) :: r
        real(dp) :: ssxm
        real(dp) :: ssxym
        real(dp) :: ssym
        real(dp) :: tstat
        real(dp) :: xmean
        real(dp) :: ymean

        call invalidate_linregress(result_value)
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        if (.not. valid_alternative(selected_alternative)) return
        if (size(x) /= size(y) .or. size(x) < 2) return
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return
        if (maxval(x) == minval(x)) return

        n = size(x)
        xmean = mean(x)
        ymean = mean(y)
        ssxm = variance(x, 0)
        ssym = variance(y, 0)
        ssxym = covariance(x, y, 0)
        if (ssxm <= 0.0_dp) return

        if (ssxm == 0.0_dp .or. ssym == 0.0_dp) then
            if (ssxym == 0.0_dp) then
                r = quiet_nan(0.0_dp)
            else
                r = 0.0_dp
            end if
        else
            r = ssxym / sqrt(ssxm * ssym)
            r = max(-1.0_dp, min(1.0_dp, r))
        end if

        result_value%slope = ssxym / ssxm
        result_value%intercept = ymean - result_value%slope * xmean
        result_value%rvalue = r

        if (n == 2) then
            if (y(1) == y(2)) then
                result_value%pvalue = 1.0_dp
            else
                result_value%pvalue = 0.0_dp
            end if
            result_value%stderr = 0.0_dp
            result_value%intercept_stderr = 0.0_dp
            return
        end if

        df = real(n - 2, dp)
        if (ieee_is_nan(r)) then
            result_value%pvalue = quiet_nan(0.0_dp)
            result_value%stderr = quiet_nan(0.0_dp)
            result_value%intercept_stderr = quiet_nan(0.0_dp)
            return
        end if
        tstat = r * sqrt(df / ((1.0_dp - r + tiny_r) * (1.0_dp + r + tiny_r)))
        result_value%pvalue = student_t_pvalue(tstat, df, selected_alternative)
        result_value%stderr = sqrt(max(0.0_dp, (1.0_dp - r ** 2) * ssym / ssxm / df))
        result_value%intercept_stderr = result_value%stderr * sqrt(ssxm + xmean ** 2)
    end function linregress

    function theilslopes(y, x, alpha, method) result(result_value)
        real(dp), intent(in) :: y(:) !! dependent variable
        real(dp), intent(in), optional :: x(:) !! independent variable; default 0,1,...
        real(dp), intent(in), optional :: alpha !! confidence degree, default 0.95
        character(len=*), intent(in), optional :: method !! separate or joint
        type(theilslopes_result) :: result_value

        character(len=8) :: selected_method
        integer :: i
        integer :: j
        integer :: low_index
        integer :: n
        integer :: nslopes
        integer :: position
        integer :: upper_index
        real(dp) :: a
        real(dp), allocatable :: intercept_data(:)
        real(dp), allocatable :: slopes(:)
        real(dp), allocatable :: xv(:)
        real(dp) :: sigsq
        real(dp) :: sigma
        real(dp) :: z

        call invalidate_theil(result_value)
        n = size(y)
        if (n < 2) return
        if (any(ieee_is_nan(y))) return
        selected_method = 'separate'
        if (present(method)) selected_method = adjustl(method)
        if (selected_method /= 'separate' .and. selected_method /= 'joint') return

        allocate(xv(n))
        if (present(x)) then
            if (size(x) /= n .or. any(ieee_is_nan(x))) return
            xv = x
        else
            do i = 1, n
                xv(i) = real(i - 1, dp)
            end do
        end if

        nslopes = 0
        do i = 1, n
            do j = 1, n
                if (xv(i) > xv(j)) nslopes = nslopes + 1
            end do
        end do
        if (nslopes <= 0) return
        allocate(slopes(nslopes))
        position = 0
        do i = 1, n
            do j = 1, n
                if (xv(i) > xv(j)) then
                    position = position + 1
                    slopes(position) = (y(i) - y(j)) / (xv(i) - xv(j))
                end if
            end do
        end do
        call sort_real(slopes)
        result_value%slope = median_sorted(slopes)

        if (selected_method == 'joint') then
            allocate(intercept_data(n))
            intercept_data = y - result_value%slope * xv
            result_value%intercept = median(intercept_data)
        else
            result_value%intercept = median(y) - result_value%slope * median(xv)
        end if

        a = 0.95_dp
        if (present(alpha)) a = alpha
        if (.not. ieee_is_finite(a) .or. a <= 0.0_dp .or. a >= 1.0_dp) return
        if (a > 0.5_dp) a = 1.0_dp - a
        z = normal_ppf(a / 2.0_dp)
        sigsq = real(n * (n - 1) * (2 * n + 5), dp)
        sigsq = sigsq - tie_cubic5_sum(xv) - tie_cubic5_sum(y)
        sigsq = sigsq / 18.0_dp
        if (sigsq < 0.0_dp) return
        sigma = sqrt(sigsq)
        upper_index = min(round_nearest_even((real(nslopes, dp) - z * sigma) / 2.0_dp), &
            nslopes - 1) + 1
        low_index = max(round_nearest_even((real(nslopes, dp) + z * sigma) / 2.0_dp) - 1, 0) + 1
        result_value%low_slope = slopes(low_index)
        result_value%high_slope = slopes(upper_index)
    end function theilslopes

    function siegelslopes(y, x, method) result(result_value)
        real(dp), intent(in) :: y(:) !! dependent variable
        real(dp), intent(in), optional :: x(:) !! independent variable; default 0,1,...
        character(len=*), intent(in), optional :: method !! hierarchical or separate
        type(siegelslopes_result) :: result_value

        character(len=12) :: selected_method
        integer :: count_nonzero
        integer :: i
        integer :: j
        integer :: n
        integer :: position
        real(dp), allocatable :: intercepts(:)
        real(dp), allocatable :: local_values(:)
        real(dp), allocatable :: slopes(:)
        real(dp), allocatable :: xv(:)

        call invalidate_siegel(result_value)
        n = size(y)
        if (n < 2 .or. any(ieee_is_nan(y))) return
        selected_method = 'hierarchical'
        if (present(method)) selected_method = adjustl(method)
        if (selected_method /= 'hierarchical' .and. selected_method /= 'separate') return

        allocate(xv(n))
        if (present(x)) then
            if (size(x) /= n .or. any(ieee_is_nan(x))) return
            xv = x
        else
            do i = 1, n
                xv(i) = real(i - 1, dp)
            end do
        end if

        allocate(slopes(n))
        if (selected_method == 'separate') allocate(intercepts(n))
        do i = 1, n
            count_nonzero = count(xv(i) /= xv)
            if (count_nonzero <= 0) return
            allocate(local_values(count_nonzero))
            position = 0
            do j = 1, n
                if (xv(i) /= xv(j)) then
                    position = position + 1
                    local_values(position) = (y(i) - y(j)) / (xv(i) - xv(j))
                end if
            end do
            slopes(i) = median(local_values)
            if (selected_method == 'separate') then
                position = 0
                do j = 1, n
                    if (xv(i) /= xv(j)) then
                        position = position + 1
                        local_values(position) = &
                            (y(j) * xv(i) - y(i) * xv(j)) / (xv(i) - xv(j))
                    end if
                end do
                intercepts(i) = median(local_values)
            end if
            deallocate(local_values)
        end do

        result_value%slope = median(slopes)
        if (selected_method == 'separate') then
            result_value%intercept = median(intercepts)
        else
            allocate(local_values(n))
            local_values = y - result_value%slope * xv
            result_value%intercept = median(local_values)
        end if
    end function siegelslopes

    function brunnermunzel(x, y, alternative, distribution) result(result_value)
        real(dp), intent(in) :: x(:) !! first independent sample
        real(dp), intent(in) :: y(:) !! second independent sample
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        character(len=*), intent(in), optional :: distribution !! t or normal
        type(association_result) :: result_value

        character(len=9) :: selected_alternative
        character(len=6) :: selected_distribution
        integer :: nx
        integer :: ny
        real(dp) :: df
        real(dp) :: df_denom
        real(dp) :: df_numer
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: rankc(:)
        real(dp), allocatable :: rankcx(:)
        real(dp), allocatable :: rankcy(:)
        real(dp), allocatable :: rankx(:)
        real(dp), allocatable :: ranky(:)
        real(dp) :: rankcx_mean
        real(dp) :: rankcy_mean
        real(dp) :: rankx_mean
        real(dp) :: ranky_mean
        real(dp) :: sx
        real(dp) :: sy
        real(dp) :: w

        call invalidate_association(result_value)
        nx = size(x)
        ny = size(y)
        if (nx < 2 .or. ny < 2 .or. any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = adjustl(alternative)
        selected_distribution = 't'
        if (present(distribution)) selected_distribution = adjustl(distribution)
        if (.not. valid_alternative(selected_alternative)) return
        if (selected_distribution /= 't' .and. selected_distribution /= 'normal') return

        allocate(pooled(nx + ny))
        pooled(1:nx) = x
        pooled(nx + 1:) = y
        rankc = rankdata(pooled)
        rankcx = rankc(1:nx)
        rankcy = rankc(nx + 1:)
        rankx = rankdata(x)
        ranky = rankdata(y)
        rankcx_mean = mean(rankcx)
        rankcy_mean = mean(rankcy)
        rankx_mean = mean(rankx)
        ranky_mean = mean(ranky)
        sx = sum((rankcx - rankx - rankcx_mean + rankx_mean) ** 2) / real(nx - 1, dp)
        sy = sum((rankcy - ranky - rankcy_mean + ranky_mean) ** 2) / real(ny - 1, dp)
        if (real(nx, dp) * sx + real(ny, dp) * sy <= 0.0_dp) return

        w = real(nx * ny, dp) * (rankcy_mean - rankcx_mean) / &
            (real(nx + ny, dp) * sqrt(real(nx, dp) * sx + real(ny, dp) * sy))
        result_value%statistic = w
        if (selected_distribution == 't') then
            df_numer = (real(nx, dp) * sx + real(ny, dp) * sy) ** 2
            df_denom = (real(nx, dp) * sx) ** 2 / real(nx - 1, dp) + &
                (real(ny, dp) * sy) ** 2 / real(ny - 1, dp)
            if (df_denom <= 0.0_dp) return
            df = df_numer / df_denom
            result_value%pvalue = student_t_pvalue(-w, df, selected_alternative)
        else
            result_value%pvalue = normal_pvalue(-w, selected_alternative)
        end if
    end function brunnermunzel

    function page_trend_test(data, ranked, predicted_ranks, method) result(result_value)
        real(dp), intent(in) :: data(:, :) !! subjects by treatments
        logical, intent(in), optional :: ranked !! data already contain row ranks
        integer, intent(in), optional :: predicted_ranks(:) !! permutation of 1..n treatments
        character(len=*), intent(in), optional :: method !! auto, exact, or asymptotic
        type(page_trend_result) :: result_value

        logical :: already_ranked
        character(len=10) :: selected_method
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        integer, allocatable :: expected_order(:)
        real(dp) :: e0
        real(dp) :: lambda
        real(dp), allocatable :: ranks(:, :)
        real(dp) :: v0

        call invalidate_page(result_value)
        m = size(data, 1)
        n = size(data, 2)
        if (m < 2 .or. n < 3 .or. any(ieee_is_nan(data))) return
        already_ranked = .false.
        if (present(ranked)) already_ranked = ranked
        selected_method = 'auto'
        if (present(method)) selected_method = adjustl(method)
        if (selected_method /= 'auto' .and. selected_method /= 'exact' .and. &
                selected_method /= 'asymptotic') return

        allocate(expected_order(n))
        if (present(predicted_ranks)) then
            if (size(predicted_ranks) /= n) return
            expected_order = predicted_ranks
            if (.not. valid_permutation(expected_order)) return
        else
            do j = 1, n
                expected_order(j) = j
            end do
        end if

        allocate(ranks(m, n))
        if (already_ranked) then
            if (minval(data) < 1.0_dp .or. maxval(data) > real(n, dp)) return
            ranks = data
        else
            do i = 1, m
                ranks(i, :) = rankdata(data(i, :))
            end do
        end if

        result_value%statistic = 0.0_dp
        do j = 1, n
            result_value%statistic = result_value%statistic + &
                real(expected_order(j), dp) * sum(ranks(:, j))
        end do

        if (selected_method == 'auto') then
            if (n > 8 .or. (m > 12 .and. n > 3) .or. m > 20) then
                selected_method = 'asymptotic'
            else
                selected_method = 'exact'
            end if
        end if
        result_value%method = selected_method

        if (selected_method == 'asymptotic') then
            e0 = real(m * n * (n + 1) ** 2, dp) / 4.0_dp
            v0 = real(m * n ** 2 * (n + 1) * (n ** 2 - 1), dp) / 144.0_dp
            if (v0 <= 0.0_dp) return
            lambda = (result_value%statistic - e0) / sqrt(v0)
            result_value%pvalue = normal_sf(lambda)
        else
            result_value%pvalue = page_exact_sf(int(result_value%statistic), m, n)
        end if
    end function page_trend_test

    function kendall_p_exact(n, c_input, alternative) result(probability)
        integer, intent(in) :: n !! sample size
        integer, intent(in) :: c_input !! number of concordant pairs
        character(len=*), intent(in) :: alternative !! requested alternative
        real(dp) :: probability

        integer :: c
        integer :: i
        integer :: j
        integer :: total
        logical :: in_right_tail
        logical :: alternative_greater
        real(dp) :: fact
        real(dp) :: p_mass
        real(dp), allocatable :: new(:)
        real(dp), allocatable :: temp(:)
        real(dp), allocatable :: cumulative(:)

        probability = quiet_nan(0.0_dp)
        if (n <= 0) return
        total = n * (n - 1) / 2
        c = min(c_input, total - c_input)
        if (c < 0 .or. 2 * c > total) return
        in_right_tail = c_input >= total - c_input
        alternative_greater = alternative == 'greater'

        if (n == 1) then
            probability = 1.0_dp
            p_mass = 1.0_dp
        else if (n == 2) then
            probability = 1.0_dp
            p_mass = 0.5_dp
        else if (c == 0) then
            if (n < 171) then
                fact = factorial_real(n)
                probability = 2.0_dp / fact
                p_mass = probability / 2.0_dp
            else
                probability = 0.0_dp
                p_mass = 0.0_dp
            end if
        else if (c == 1) then
            if (n < 172) then
                fact = factorial_real(n)
                probability = 2.0_dp * real(n, dp) / fact
                p_mass = real(n - 1, dp) / fact
            else
                probability = 0.0_dp
                p_mass = 0.0_dp
            end if
        else if (2 * c == total .and. alternative == 'two-sided') then
            probability = 1.0_dp
            return
        else
            allocate(new(0:c), temp(0:c), cumulative(0:c))
            new = 0.0_dp
            new(0:min(1, c)) = 1.0_dp
            if (n < 171) then
                do j = 3, n
                    temp = 0.0_dp
                    temp(0) = new(0)
                    do i = 1, c
                        temp(i) = temp(i - 1) + new(i)
                    end do
                    if (j <= c) then
                        cumulative = temp
                        do i = j, c
                            temp(i) = temp(i) - cumulative(i - j)
                        end do
                    end if
                    new = temp
                end do
                fact = factorial_real(n)
                probability = 2.0_dp * sum(new) / fact
                p_mass = new(c) / fact
            else
                do j = 3, n
                    temp = 0.0_dp
                    temp(0) = new(0) / real(j, dp)
                    do i = 1, c
                        temp(i) = (temp(i - 1) * real(j, dp) + new(i)) / real(j, dp)
                    end do
                    if (j <= c) then
                        cumulative = temp
                        do i = j, c
                            temp(i) = temp(i) - cumulative(i - j)
                        end do
                    end if
                    new = temp
                end do
                probability = sum(new)
                p_mass = new(c) / 2.0_dp
            end if
        end if

        if (alternative /= 'two-sided') then
            if (in_right_tail .eqv. alternative_greater) then
                probability = probability / 2.0_dp
            else
                probability = 1.0_dp - probability / 2.0_dp + p_mass
            end if
        end if
        probability = max(0.0_dp, min(1.0_dp, probability))
    end function kendall_p_exact

    subroutine tie_statistics(values, pair_ties, cubic0, cubic1, nclasses)
        real(dp), intent(in) :: values(:) !! values without NaNs
        integer(int64), intent(out) :: pair_ties !! sum choose(group_size, 2)
        integer(int64), intent(out) :: cubic0 !! sum t(t-1)(t-2)
        integer(int64), intent(out) :: cubic1 !! sum t(t-1)(2t+5)
        integer, intent(out) :: nclasses !! number of distinct values

        integer :: group_end
        integer :: group_start
        integer(int64) :: t
        real(dp), allocatable :: sorted(:)

        pair_ties = 0
        cubic0 = 0
        cubic1 = 0
        nclasses = 0
        if (size(values) == 0) return
        sorted = values
        call sort_real(sorted)
        group_start = 1
        do while (group_start <= size(sorted))
            group_end = group_start
            do while (group_end < size(sorted))
                if (.not. same_value(sorted(group_end), sorted(group_end + 1))) exit
                group_end = group_end + 1
            end do
            t = int(group_end - group_start + 1, int64)
            nclasses = nclasses + 1
            if (t > 1) then
                pair_ties = pair_ties + t * (t - 1_int64) / 2_int64
                cubic0 = cubic0 + t * (t - 1_int64) * (t - 2_int64)
                cubic1 = cubic1 + t * (t - 1_int64) * (2_int64 * t + 5_int64)
            end if
            group_start = group_end + 1
        end do
    end subroutine tie_statistics

    function tie_cubic5_sum(values) result(value)
        real(dp), intent(in) :: values(:) !! values without NaNs
        real(dp) :: value

        integer(int64) :: dummy0
        integer(int64) :: dummy_pairs
        integer :: dummy_classes
        integer(int64) :: cubic1

        call tie_statistics(values, dummy_pairs, dummy0, cubic1, dummy_classes)
        value = real(cubic1, dp)
    end function tie_cubic5_sum

    function page_exact_sf(l_observed, m, k) result(probability)
        integer, intent(in) :: l_observed !! observed L truncated toward zero as in SciPy
        integer, intent(in) :: m !! number of rows/subjects
        integer, intent(in) :: k !! number of treatments
        real(dp) :: probability

        integer :: a
        integer :: b
        integer :: i
        integer :: l
        integer :: max_sum
        integer :: min_sum
        integer :: row
        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next(:)
        real(dp), allocatable :: single(:)

        probability = quiet_nan(0.0_dp)
        if (m < 1 .or. k < 1 .or. k > 10) return
        a = k * (k + 1) * (k + 2) / 6
        b = k * (k + 1) * (2 * k + 1) / 6
        allocate(single(a:b))
        single = 0.0_dp
        call page_single_row_pmf(k, a, b, single)

        min_sum = a
        max_sum = b
        allocate(current(min_sum:max_sum))
        current = single
        do row = 2, m
            allocate(next(row * a:row * b))
            next = 0.0_dp
            do i = min_sum, max_sum
                if (current(i) == 0.0_dp) cycle
                do l = a, b
                    if (single(l) == 0.0_dp) cycle
                    next(i + l) = next(i + l) + current(i) * single(l)
                end do
            end do
            call move_alloc(next, current)
            min_sum = row * a
            max_sum = row * b
        end do

        if (l_observed <= min_sum) then
            probability = 1.0_dp
        else if (l_observed > max_sum) then
            probability = 0.0_dp
        else
            probability = sum(current(l_observed:max_sum))
        end if
        probability = max(0.0_dp, min(1.0_dp, probability))
    end function page_exact_sf

    subroutine page_single_row_pmf(k, a, b, pmf)
        integer, intent(in) :: k !! number of treatments
        integer, intent(in) :: a !! minimum one-row L
        integer, intent(in) :: b !! maximum one-row L
        real(dp), intent(out) :: pmf(a:b) !! one-row null PMF

        integer :: perm(k)
        logical :: used(k)
        real(dp) :: fact

        pmf = 0.0_dp
        used = .false.
        call enumerate_page_permutations(1, k, a, perm, used, pmf)
        fact = factorial_real(k)
        pmf = pmf / fact
    end subroutine page_single_row_pmf

    recursive subroutine enumerate_page_permutations(position, k, a, perm, used, pmf)
        integer, intent(in) :: position !! current permutation position
        integer, intent(in) :: k !! permutation length
        integer, intent(in) :: a !! minimum one-row L used as PMF offset
        integer, intent(inout) :: perm(k) !! current permutation
        logical, intent(inout) :: used(k) !! value-used flags
        real(dp), intent(inout) :: pmf(:) !! offset-indexed by declared lower bound

        integer :: j
        integer :: lvalue

        if (position > k) then
            lvalue = 0
            do j = 1, k
                lvalue = lvalue + j * perm(j)
            end do
            pmf(lvalue - a + 1) = pmf(lvalue - a + 1) + 1.0_dp
            return
        end if
        do j = 1, k
            if (.not. used(j)) then
                used(j) = .true.
                perm(position) = j
                call enumerate_page_permutations(position + 1, k, a, perm, used, pmf)
                used(j) = .false.
            end if
        end do
    end subroutine enumerate_page_permutations

    function factorial_real(n) result(value)
        integer, intent(in) :: n !! nonnegative integer
        real(dp) :: value
        integer :: i

        value = 1.0_dp
        do i = 2, n
            value = value * real(i, dp)
        end do
    end function factorial_real

    subroutine stable_sort_indices(values, order)
        real(dp), intent(in) :: values(:) !! source values used as sort keys
        integer, intent(inout) :: order(:) !! indices sorted by value, preserving tie order

        integer, allocatable :: workspace(:)

        if (size(order) <= 1) return
        allocate(workspace(size(order)))
        call merge_sort_indices(values, order, workspace, 1, size(order))
    end subroutine stable_sort_indices

    recursive subroutine merge_sort_indices(values, order, workspace, left, right)
        real(dp), intent(in) :: values(:) !! source values used as sort keys
        integer, intent(inout) :: order(:) !! indices being sorted
        integer, intent(inout) :: workspace(:) !! merge workspace
        integer, intent(in) :: left !! first index of inclusive range
        integer, intent(in) :: right !! last index of inclusive range

        integer :: i
        integer :: j
        integer :: k
        integer :: middle

        if (left >= right) return
        middle = left + (right - left) / 2
        call merge_sort_indices(values, order, workspace, left, middle)
        call merge_sort_indices(values, order, workspace, middle + 1, right)
        i = left
        j = middle + 1
        do k = left, right
            if (i > middle) then
                workspace(k) = order(j)
                j = j + 1
            else if (j > right) then
                workspace(k) = order(i)
                i = i + 1
            else if (values(order(i)) <= values(order(j))) then
                workspace(k) = order(i)
                i = i + 1
            else
                workspace(k) = order(j)
                j = j + 1
            end if
        end do
        order(left:right) = workspace(left:right)
    end subroutine merge_sort_indices

    function factorial_count_capped(n, cap) result(value)
        integer, intent(in) :: n !! nonnegative factorial argument
        integer(int64), intent(in) :: cap !! positive value returned when factorial reaches cap
        integer(int64) :: value

        integer :: i

        value = 1_int64
        do i = 2, n
            if (value > (cap - 1_int64) / int(i, int64)) then
                value = cap
                return
            end if
            value = value * int(i, int64)
            if (value >= cap) then
                value = cap
                return
            end if
        end do
    end function factorial_count_capped

    logical function next_integer_permutation(values) result(has_next)
        integer, intent(inout) :: values(:) !! current permutation, replaced by lexicographic successor

        integer :: i
        integer :: j
        integer :: left
        integer :: right
        integer :: temporary

        i = size(values) - 1
        do while (i >= 1)
            if (values(i) < values(i + 1)) exit
            i = i - 1
        end do
        if (i < 1) then
            has_next = .false.
            return
        end if
        j = size(values)
        do while (values(j) <= values(i))
            j = j - 1
        end do
        temporary = values(i)
        values(i) = values(j)
        values(j) = temporary
        left = i + 1
        right = size(values)
        do while (left < right)
            temporary = values(left)
            values(left) = values(right)
            values(right) = temporary
            left = left + 1
            right = right - 1
        end do
        has_next = .true.
    end function next_integer_permutation

    function unbiased_index(state, upper) result(index_value)
        type(rng_state), intent(inout) :: state !! explicit generator state advanced until accepted
        integer, intent(in) :: upper !! inclusive upper bound, at least one
        integer :: index_value

        integer(int64) :: bucket
        integer(int64) :: limit
        integer(int64) :: upper64
        real(dp) :: u

        upper64 = int(upper, int64)
        limit = (integer_grid_size / upper64) * upper64
        do
            u = rng_uniform(state)
            bucket = int(u * real(integer_grid_size, dp), int64)
            if (bucket < limit) exit
        end do
        index_value = int(modulo(bucket, upper64)) + 1
    end function unbiased_index

    function median_sorted(values) result(value)
        real(dp), intent(in) :: values(:) !! nonempty sorted values
        real(dp) :: value
        integer :: n

        n = size(values)
        if (n <= 0) then
            value = quiet_nan(0.0_dp)
        else if (mod(n, 2) == 1) then
            value = values((n + 1) / 2)
        else
            value = 0.5_dp * (values(n / 2) + values(n / 2 + 1))
        end if
    end function median_sorted

    subroutine sort_real(values)
        real(dp), intent(inout) :: values(:) !! sorted ascending in place
        real(dp), allocatable :: workspace(:)

        if (size(values) <= 1) return
        allocate(workspace(size(values)))
        call merge_sort_real(values, workspace, 1, size(values))
    end subroutine sort_real

    recursive subroutine merge_sort_real(values, workspace, left, right)
        real(dp), intent(inout) :: values(:) !! values being sorted
        real(dp), intent(inout) :: workspace(:) !! merge workspace
        integer, intent(in) :: left !! first index
        integer, intent(in) :: right !! last index

        integer :: i
        integer :: j
        integer :: k
        integer :: mid

        if (left >= right) return
        mid = left + (right - left) / 2
        call merge_sort_real(values, workspace, left, mid)
        call merge_sort_real(values, workspace, mid + 1, right)
        i = left
        j = mid + 1
        do k = left, right
            if (i > mid) then
                workspace(k) = values(j)
                j = j + 1
            else if (j > right) then
                workspace(k) = values(i)
                i = i + 1
            else if (values(i) <= values(j)) then
                workspace(k) = values(i)
                i = i + 1
            else
                workspace(k) = values(j)
                j = j + 1
            end if
        end do
        values(left:right) = workspace(left:right)
    end subroutine merge_sort_real

    integer function round_nearest_even(value) result(rounded)
        real(dp), intent(in) :: value !! finite nonnegative value
        integer :: lower
        real(dp) :: fraction

        lower = int(floor(value))
        fraction = value - real(lower, dp)
        if (fraction < 0.5_dp) then
            rounded = lower
        else if (fraction > 0.5_dp) then
            rounded = lower + 1
        else if (mod(lower, 2) == 0) then
            rounded = lower
        else
            rounded = lower + 1
        end if
    end function round_nearest_even

    pure logical function valid_permutation(values) result(valid)
        integer, intent(in) :: values(:) !! candidate permutation
        logical, allocatable :: seen(:)
        integer :: i

        valid = .false.
        if (size(values) == 0) return
        allocate(seen(size(values)))
        seen = .false.
        do i = 1, size(values)
            if (values(i) < 1 .or. values(i) > size(values)) return
            if (seen(values(i))) return
            seen(values(i)) = .true.
        end do
        valid = .true.
    end function valid_permutation

    pure logical function same_value(x, y) result(same)
        real(dp), intent(in) :: x !! first non-NaN value
        real(dp), intent(in) :: y !! second non-NaN value

        same = x <= y .and. y <= x
    end function same_value

    pure function normal_pvalue(statistic, alternative) result(pvalue)
        real(dp), intent(in) :: statistic !! standard-normal test statistic
        character(len=*), intent(in) :: alternative !! valid alternative
        real(dp) :: pvalue

        select case (alternative)
        case ('less')
            pvalue = normal_cdf(statistic)
        case ('greater')
            pvalue = normal_sf(statistic)
        case default
            pvalue = min(1.0_dp, 2.0_dp * normal_sf(abs(statistic)))
        end select
    end function normal_pvalue

    pure function student_t_pvalue(statistic, df, alternative) result(pvalue)
        real(dp), intent(in) :: statistic !! Student t statistic
        real(dp), intent(in) :: df !! positive degrees of freedom
        character(len=*), intent(in) :: alternative !! valid alternative
        real(dp) :: pvalue

        select case (alternative)
        case ('less')
            pvalue = t_cdf(statistic, df)
        case ('greater')
            pvalue = t_sf(statistic, df)
        case default
            pvalue = min(1.0_dp, 2.0_dp * t_sf(abs(statistic), df))
        end select
    end function student_t_pvalue

    pure logical function valid_alternative(alternative) result(valid)
        character(len=*), intent(in) :: alternative !! alternative name

        valid = alternative == 'two-sided' .or. alternative == 'less' .or. &
            alternative == 'greater'
    end function valid_alternative

    subroutine invalidate_association(result_value)
        type(association_result), intent(out) :: result_value !! result set to NaN

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
    end subroutine invalidate_association

    subroutine invalidate_somersd(result_value)
        type(somersd_result), intent(out) :: result_value !! result set to NaN and no table

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%correlation = quiet_nan(0.0_dp)
        if (allocated(result_value%table)) deallocate(result_value%table)
    end subroutine invalidate_somersd

    subroutine invalidate_linregress(result_value)
        type(linregress_result), intent(out) :: result_value !! result set to NaN

        result_value%slope = quiet_nan(0.0_dp)
        result_value%intercept = quiet_nan(0.0_dp)
        result_value%rvalue = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%stderr = quiet_nan(0.0_dp)
        result_value%intercept_stderr = quiet_nan(0.0_dp)
    end subroutine invalidate_linregress

    subroutine invalidate_theil(result_value)
        type(theilslopes_result), intent(out) :: result_value !! result set to NaN

        result_value%slope = quiet_nan(0.0_dp)
        result_value%intercept = quiet_nan(0.0_dp)
        result_value%low_slope = quiet_nan(0.0_dp)
        result_value%high_slope = quiet_nan(0.0_dp)
    end subroutine invalidate_theil

    subroutine invalidate_siegel(result_value)
        type(siegelslopes_result), intent(out) :: result_value !! result set to NaN

        result_value%slope = quiet_nan(0.0_dp)
        result_value%intercept = quiet_nan(0.0_dp)
    end subroutine invalidate_siegel

    subroutine invalidate_page(result_value)
        type(page_trend_result), intent(out) :: result_value !! result set to NaN

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%method = ''
    end subroutine invalidate_page

end module scifort_association_extended
