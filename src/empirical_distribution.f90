! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Empirical-distribution and exact quantile-inference helpers matching the
! public numerical semantics of scipy.stats.ecdf, scipy.stats.quantile_test,
! and scipy.stats.order_statistic in SciPy 1.17.0.  Right-censored ECDFs use
! the Kaplan-Meier product-limit estimator; confidence intervals use the
! Greenwood and exponential-Greenwood formulas.

module scifort_empirical_distribution
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_beta, only : beta_cdf, beta_isf, beta_pdf, beta_ppf, beta_sf
    use scifort_binomial, only : binomial_cdf, binomial_isf, binomial_ppf, binomial_sf
    use scifort_contingency_meta, only : confidence_interval
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_normal, only : normal_cdf, normal_ppf, normal_sf
    implicit none
    private

    abstract interface
        pure function order_statistic_callback(x, parameters) result(y)
            import dp
            real(dp), intent(in) :: x !! scalar argument or probability
            real(dp), intent(in) :: parameters(:) !! user-defined base-distribution parameters
            real(dp) :: y
        end function order_statistic_callback
    end interface

    type, public :: censored_data
        real(dp), allocatable :: uncensored(:) !! observed event/failure times
        real(dp), allocatable :: right(:) !! right-censoring/loss times
        real(dp), allocatable :: left(:) !! finite upper limits; used by parametric fitting only
        real(dp), allocatable :: interval(:,:) !! (n,2) lower/upper limits; parametric fitting only
    contains
        procedure :: n_observations => censored_data_n_observations
    end type censored_data

    type, public :: logrank_result
        real(dp) :: statistic = 0.0_dp !! signed standard-normal logrank statistic
        real(dp) :: pvalue = 1.0_dp !! asymptotic normal p-value
    end type logrank_result

    type, public :: empirical_distribution_function
        real(dp), allocatable :: quantiles(:) !! unique observed times/values
        real(dp), allocatable :: probabilities(:) !! CDF or SF values at quantiles
        integer, allocatable :: at_risk(:) !! risk-set sizes used by Greenwood variance
        integer, allocatable :: events(:) !! event/death counts at each quantile
        logical :: is_survival = .false. !! true for an empirical survival function
        logical :: has_variance_data = .false. !! false for confidence-interval bounds
    contains
        procedure :: evaluate => edf_evaluate_scalar
        procedure :: confidence_interval => edf_confidence_interval_method
    end type empirical_distribution_function

    type, public :: ecdf_confidence_interval
        type(empirical_distribution_function) :: low !! lower confidence-bound function
        type(empirical_distribution_function) :: high !! upper confidence-bound function
    end type ecdf_confidence_interval

    type, public :: ecdf_result
        type(empirical_distribution_function) :: cdf !! empirical cumulative distribution function
        type(empirical_distribution_function) :: sf !! empirical survival function
    end type ecdf_result

    type, public :: quantile_test_result
        integer :: statistic = 0 !! T1 or T2 count used for the p-value
        integer :: statistic_type = 1 !! 1 for T1 (<= q), 2 for T2 (< q)
        real(dp) :: pvalue = 1.0_dp !! exact binomial p-value
        character(len=9) :: alternative = 'two-sided' !! two-sided, less, or greater
        real(dp) :: probability = 0.5_dp !! target population quantile probability
        real(dp), allocatable :: sample(:) !! original one-dimensional sample
    contains
        procedure :: confidence_interval => quantile_test_confidence_interval
    end type quantile_test_result

    type, public :: order_statistic_distribution
        integer :: rank = 1 !! order-statistic rank r
        integer :: sample_size = 1 !! sample size n
        real(dp), allocatable :: parameters(:) !! base-distribution callback parameters
        procedure(order_statistic_callback), pointer, nopass :: base_cdf => null()
        procedure(order_statistic_callback), pointer, nopass :: base_pdf => null()
        procedure(order_statistic_callback), pointer, nopass :: base_ppf => null()
    contains
        procedure :: cdf => order_statistic_cdf
        procedure :: sf => order_statistic_sf
        procedure :: pdf => order_statistic_pdf
        procedure :: ppf => order_statistic_ppf
        procedure :: isf => order_statistic_isf
    end type order_statistic_distribution

    interface ecdf
        module procedure ecdf_uncensored
        module procedure ecdf_right_censored
        module procedure ecdf_censored
    end interface ecdf

    interface logrank
        module procedure logrank_arrays
        module procedure logrank_censored
        module procedure logrank_censored_array
        module procedure logrank_array_censored
    end interface logrank

    public :: ecdf
    public :: logrank
    public :: order_statistic
    public :: order_statistic_callback
    public :: quantile_test

contains

    pure function censored_data_n_observations(self) result(n)
        implicit none
        class(censored_data), intent(in) :: self !! exact and censored observations
        integer :: n

        n = 0
        if (allocated(self%uncensored)) n = n + size(self%uncensored)
        if (allocated(self%right)) n = n + size(self%right)
        if (allocated(self%left)) n = n + size(self%left)
        if (allocated(self%interval)) n = n + size(self%interval, 1)
    end function censored_data_n_observations

    function ecdf_censored(sample) result(res)
        type(censored_data), intent(in) :: sample
        type(ecdf_result) :: res
        real(dp), allocatable :: right(:)
        real(dp), allocatable :: uncensored(:)

        if (.not. right_censoring_only(sample)) then
            call set_edf(res%cdf, [quiet_nan(0.0_dp)], [quiet_nan(0.0_dp)], &
                [0], [0], .false., .false.)
            call set_edf(res%sf, [quiet_nan(0.0_dp)], [quiet_nan(0.0_dp)], &
                [0], [0], .true., .false.)
            return
        end if
        call unpack_censored_data(sample, uncensored, right)
        res = ecdf_right_censored(uncensored, right)
    end function ecdf_censored

    function logrank_arrays(x, y, alternative) result(res)
        real(dp), intent(in) :: x(:) !! uncensored observations in group X
        real(dp), intent(in) :: y(:) !! uncensored observations in group Y
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(logrank_result) :: res
        real(dp), allocatable :: empty(:)

        allocate(empty(0))
        res = logrank_core(x, empty, y, empty, alternative)
    end function logrank_arrays

    function logrank_censored(x, y, alternative) result(res)
        type(censored_data), intent(in) :: x !! group X right-censored sample
        type(censored_data), intent(in) :: y !! group Y right-censored sample
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(logrank_result) :: res
        real(dp), allocatable :: xr(:)
        real(dp), allocatable :: xu(:)
        real(dp), allocatable :: yr(:)
        real(dp), allocatable :: yu(:)

        if (.not. right_censoring_only(x) .or. .not. right_censoring_only(y)) then
            res%statistic = quiet_nan(0.0_dp)
            res%pvalue = quiet_nan(0.0_dp)
            return
        end if
        call unpack_censored_data(x, xu, xr)
        call unpack_censored_data(y, yu, yr)
        res = logrank_core(xu, xr, yu, yr, alternative)
    end function logrank_censored

    function logrank_censored_array(x, y, alternative) result(res)
        type(censored_data), intent(in) :: x !! group X right-censored sample
        real(dp), intent(in) :: y(:) !! uncensored observations in group Y
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(logrank_result) :: res
        real(dp), allocatable :: empty(:)
        real(dp), allocatable :: xr(:)
        real(dp), allocatable :: xu(:)

        if (.not. right_censoring_only(x)) then
            res%statistic = quiet_nan(0.0_dp)
            res%pvalue = quiet_nan(0.0_dp)
            return
        end if
        allocate(empty(0))
        call unpack_censored_data(x, xu, xr)
        res = logrank_core(xu, xr, y, empty, alternative)
    end function logrank_censored_array

    function logrank_array_censored(x, y, alternative) result(res)
        real(dp), intent(in) :: x(:) !! uncensored observations in group X
        type(censored_data), intent(in) :: y !! group Y right-censored sample
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(logrank_result) :: res
        real(dp), allocatable :: empty(:)
        real(dp), allocatable :: yr(:)
        real(dp), allocatable :: yu(:)

        if (.not. right_censoring_only(y)) then
            res%statistic = quiet_nan(0.0_dp)
            res%pvalue = quiet_nan(0.0_dp)
            return
        end if
        allocate(empty(0))
        call unpack_censored_data(y, yu, yr)
        res = logrank_core(x, empty, yu, yr, alternative)
    end function logrank_array_censored

    function logrank_core(x_uncensored, x_right, y_uncensored, y_right, alternative) result(res)
        real(dp), intent(in) :: x_uncensored(:)
        real(dp), intent(in) :: x_right(:)
        real(dp), intent(in) :: y_uncensored(:)
        real(dp), intent(in) :: y_right(:)
        character(len=*), intent(in), optional :: alternative
        type(logrank_result) :: res

        character(len=9) :: alt
        integer :: deaths
        integer :: deaths_x
        integer :: i
        integer :: nx
        integer :: ny
        integer :: nxy
        integer :: nevents
        real(dp) :: expected_x
        real(dp) :: numerator
        real(dp) :: denominator
        real(dp) :: sum_var
        real(dp) :: time
        real(dp), allocatable :: event_times(:)

        alt = 'two-sided'
        if (present(alternative)) alt = normalize_alternative(alternative)
        if (.not. valid_alternative(alt)) then
            error stop 'logrank: alternative must be two-sided, less, or greater'
        end if
        if (.not. all_finite(x_uncensored) .or. .not. all_finite(x_right) .or. &
            .not. all_finite(y_uncensored) .or. .not. all_finite(y_right)) then
            error stop 'logrank: samples must contain only finite values'
        end if

        nevents = size(x_uncensored) + size(y_uncensored)
        allocate(event_times(nevents))
        if (size(x_uncensored) > 0) event_times(1:size(x_uncensored)) = x_uncensored
        if (size(y_uncensored) > 0) then
            event_times(size(x_uncensored) + 1:nevents) = y_uncensored
        end if
        if (nevents > 0) call sort_real(event_times)

        expected_x = 0.0_dp
        sum_var = 0.0_dp
        i = 1
        do while (i <= nevents)
            time = event_times(i)
            deaths_x = count_equal_value(x_uncensored, time)
            deaths = deaths_x + count_equal_value(y_uncensored, time)
            nx = count_at_risk(x_uncensored, x_right, time)
            ny = count_at_risk(y_uncensored, y_right, time)
            nxy = nx + ny
            if (nxy > 0) then
                expected_x = expected_x + real(nx, dp) * real(deaths, dp) / real(nxy, dp)
                if (nxy > 1) then
                    numerator = real(nx, dp) * real(ny, dp) * real(deaths, dp) * &
                        real(nxy - deaths, dp)
                    denominator = real(nxy, dp) * real(nxy, dp) * real(nxy - 1, dp)
                    sum_var = sum_var + numerator / denominator
                end if
            end if
            do while (i <= nevents)
                if (.not. real_equal(event_times(i), time)) exit
                i = i + 1
            end do
        end do

        if (sum_var > 0.0_dp .and. ieee_is_finite(sum_var)) then
            res%statistic = (real(size(x_uncensored), dp) - expected_x) / sqrt(sum_var)
        else
            res%statistic = quiet_nan(0.0_dp)
        end if

        if (ieee_is_nan(res%statistic)) then
            res%pvalue = quiet_nan(0.0_dp)
        else
            select case (trim(alt))
            case ('less')
                res%pvalue = normal_cdf(res%statistic)
            case ('greater')
                res%pvalue = normal_sf(res%statistic)
            case default
                res%pvalue = min(1.0_dp, 2.0_dp * normal_sf(abs(res%statistic)))
            end select
        end if
    end function logrank_core

    pure function right_censoring_only(sample) result(valid)
        implicit none
        type(censored_data), intent(in) :: sample !! sample for ECDF or logrank
        logical :: valid
        valid = .true.
        if (allocated(sample%left)) then
            if (size(sample%left) > 0) valid = .false.
        end if
        if (allocated(sample%interval)) then
            if (size(sample%interval, 1) > 0 .or. size(sample%interval, 2) /= 2) valid = .false.
        end if
    end function right_censoring_only

    subroutine unpack_censored_data(sample, uncensored, right)
        type(censored_data), intent(in) :: sample
        real(dp), allocatable, intent(out) :: uncensored(:)
        real(dp), allocatable, intent(out) :: right(:)

        if (allocated(sample%uncensored)) then
            uncensored = sample%uncensored
        else
            allocate(uncensored(0))
        end if
        if (allocated(sample%right)) then
            right = sample%right
        else
            allocate(right(0))
        end if
    end subroutine unpack_censored_data

    pure integer function count_equal_value(values, target) result(n)
        real(dp), intent(in) :: values(:)
        real(dp), intent(in) :: target
        integer :: i

        n = 0
        do i = 1, size(values)
            if (real_equal(values(i), target)) n = n + 1
        end do
    end function count_equal_value

    pure integer function count_at_risk(uncensored, right, time) result(n)
        real(dp), intent(in) :: uncensored(:)
        real(dp), intent(in) :: right(:)
        real(dp), intent(in) :: time
        integer :: i

        n = 0
        do i = 1, size(uncensored)
            if (uncensored(i) >= time) n = n + 1
        end do
        do i = 1, size(right)
            if (right(i) >= time) n = n + 1
        end do
    end function count_at_risk

    function ecdf_uncensored(sample) result(res)
        real(dp), intent(in) :: sample(:) !! one-dimensional uncensored sample
        type(ecdf_result) :: res

        integer :: cumulative
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        real(dp), allocatable :: cdf(:)
        integer, allocatable :: counts(:)
        integer, allocatable :: events(:)
        integer, allocatable :: n_risk(:)
        real(dp), allocatable :: sorted(:)
        real(dp), allocatable :: sf(:)
        real(dp), allocatable :: unique_values(:)

        if (.not. all_finite(sample)) then
            error stop 'ecdf: sample must contain only finite values'
        end if

        n = size(sample)
        if (n == 0) then
            call set_edf(res%cdf, [real(dp) ::], [real(dp) ::], [integer ::], &
                [integer ::], .false., .true.)
            call set_edf(res%sf, [real(dp) ::], [real(dp) ::], [integer ::], &
                [integer ::], .true., .true.)
            return
        end if

        allocate(sorted(n))
        sorted = sample
        call sort_real(sorted)

        m = 1
        do i = 2, n
            if (.not. real_equal(sorted(i), sorted(i - 1))) m = m + 1
        end do
        allocate(unique_values(m), counts(m), cdf(m), sf(m), n_risk(m), events(m))

        j = 1
        unique_values(1) = sorted(1)
        counts(1) = 1
        do i = 2, n
            if (real_equal(sorted(i), sorted(i - 1))) then
                counts(j) = counts(j) + 1
            else
                j = j + 1
                unique_values(j) = sorted(i)
                counts(j) = 1
            end if
        end do

        cumulative = 0
        do i = 1, m
            n_risk(i) = n - cumulative
            cumulative = cumulative + counts(i)
            events(i) = counts(i)
            cdf(i) = real(cumulative, dp) / real(n, dp)
            sf(i) = 1.0_dp - cdf(i)
        end do

        call set_edf(res%cdf, unique_values, cdf, n_risk, events, .false., .true.)
        call set_edf(res%sf, unique_values, sf, n_risk, events, .true., .true.)
    end function ecdf_uncensored

    function ecdf_right_censored(uncensored, right) result(res)
        real(dp), intent(in) :: uncensored(:) !! observed event/death times
        real(dp), intent(in) :: right(:) !! right-censoring/loss times
        type(ecdf_result) :: res

        real(dp), allocatable :: cdf(:)
        integer :: cumulative_events
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        integer :: nu
        integer, allocatable :: died(:)
        integer, allocatable :: died_sorted(:)
        integer, allocatable :: events(:)
        integer, allocatable :: n_risk(:)
        integer, allocatable :: order(:)
        real(dp), allocatable :: sf(:)
        real(dp) :: survival
        real(dp), allocatable :: times(:)
        real(dp), allocatable :: times_sorted(:)
        real(dp), allocatable :: unique_times(:)

        if (.not. all_finite(uncensored) .or. .not. all_finite(right)) then
            error stop 'ecdf: sample must contain only finite values'
        end if

        nu = size(uncensored)
        n = nu + size(right)
        if (n == 0) then
            call set_edf(res%cdf, [real(dp) ::], [real(dp) ::], [integer ::], &
                [integer ::], .false., .true.)
            call set_edf(res%sf, [real(dp) ::], [real(dp) ::], [integer ::], &
                [integer ::], .true., .true.)
            return
        end if

        allocate(times(n), died(n), order(n), times_sorted(n), died_sorted(n))
        if (nu > 0) then
            times(1:nu) = uncensored
            died(1:nu) = 1
        end if
        if (size(right) > 0) then
            times(nu + 1:n) = right
            died(nu + 1:n) = 0
        end if
        call stable_sort_indices(times, order)
        do i = 1, n
            times_sorted(i) = times(order(i))
            died_sorted(i) = died(order(i))
        end do

        m = 1
        do i = 2, n
            if (.not. real_equal(times_sorted(i), times_sorted(i - 1))) m = m + 1
        end do
        allocate(unique_times(m), n_risk(m), events(m), cdf(m), sf(m))

        j = 1
        unique_times(1) = times_sorted(1)
        n_risk(1) = n
        events(1) = died_sorted(1)
        do i = 2, n
            if (real_equal(times_sorted(i), times_sorted(i - 1))) then
                events(j) = events(j) + died_sorted(i)
            else
                j = j + 1
                unique_times(j) = times_sorted(i)
                n_risk(j) = n - i + 1
                events(j) = died_sorted(i)
            end if
        end do

        survival = 1.0_dp
        cumulative_events = 0
        do i = 1, m
            if (n_risk(i) > 0) then
                survival = survival * real(n_risk(i) - events(i), dp) / real(n_risk(i), dp)
            end if
            sf(i) = survival
            cdf(i) = 1.0_dp - survival
            cumulative_events = cumulative_events + events(i)
        end do

        call set_edf(res%cdf, unique_times, cdf, n_risk, events, .false., .true.)
        call set_edf(res%sf, unique_times, sf, n_risk, events, .true., .true.)
    end function ecdf_right_censored

    function quantile_test(x, q, p, alternative) result(res)
        real(dp), intent(in) :: x(:) !! one-dimensional sample
        real(dp), intent(in), optional :: q !! hypothesized quantile value, default zero
        real(dp), intent(in), optional :: p !! target probability, default 0.5
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(quantile_test_result) :: res

        character(len=9) :: alt
        integer :: i
        integer :: n
        real(dp) :: p0
        real(dp) :: p_greater
        real(dp) :: p_less
        real(dp) :: q0
        integer :: t1
        integer :: t2

        q0 = 0.0_dp
        if (present(q)) q0 = q
        p0 = 0.5_dp
        if (present(p)) p0 = p
        alt = 'two-sided'
        if (present(alternative)) alt = normalize_alternative(alternative)

        if (ieee_is_nan(p0)) then
            ! Match SciPy: NaN p propagates to the p-value rather than raising.
        else if (.not. (p0 > 0.0_dp .and. p0 < 1.0_dp)) then
            error stop 'quantile_test: p must be strictly between zero and one'
        end if
        if (.not. valid_alternative(alt)) then
            error stop 'quantile_test: alternative must be two-sided, less, or greater'
        end if

        n = size(x)
        t1 = 0
        t2 = 0
        do i = 1, n
            if (x(i) <= q0) t1 = t1 + 1
            if (x(i) < q0) t2 = t2 + 1
        end do

        res%alternative = alt
        res%probability = p0
        res%sample = x

        if (ieee_is_nan(p0)) then
            if (trim(alt) == 'less') then
                res%statistic = t2
                res%statistic_type = 2
            else
                res%statistic = t1
                res%statistic_type = 1
            end if
            res%pvalue = quiet_nan(p0)
            return
        end if

        select case (trim(alt))
        case ('less')
            res%statistic = t2
            res%statistic_type = 2
            res%pvalue = binomial_sf(t2 - 1, n, p0)
        case ('greater')
            res%statistic = t1
            res%statistic_type = 1
            res%pvalue = binomial_cdf(t1, n, p0)
        case default
            p_greater = binomial_cdf(t1, n, p0)
            p_less = binomial_sf(t2 - 1, n, p0)
            if (p_less < p_greater) then
                res%statistic = t2
                res%statistic_type = 2
                res%pvalue = min(1.0_dp, max(0.0_dp, 2.0_dp * p_less))
            else
                res%statistic = t1
                res%statistic_type = 1
                res%pvalue = min(1.0_dp, max(0.0_dp, 2.0_dp * p_greater))
            end if
        end select
    end function quantile_test

    function quantile_test_confidence_interval(self, confidence_level) result(ci)
        class(quantile_test_result), intent(in) :: self !! quantile-test result
        real(dp), intent(in), optional :: confidence_level !! confidence level in (0,1), default 0.95
        type(confidence_interval) :: ci

        real(dp) :: alpha
        real(dp) :: level
        integer :: high_index
        integer :: low_index
        integer :: n
        real(dp), allocatable :: x(:)

        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        if (.not. (level > 0.0_dp .and. level < 1.0_dp)) then
            error stop 'quantile_test confidence_interval: confidence level must be in (0,1)'
        end if

        if (ieee_is_nan(self%probability)) then
            error stop 'quantile_test confidence_interval: quantile probability is NaN'
        end if

        n = size(self%sample)
        allocate(x(n))
        x = self%sample
        call sort_nan_last(x)

        select case (trim(self%alternative))
        case ('less')
            alpha = 1.0_dp - level
            ci%low = negative_infinity(level)
            high_index = int(binomial_isf(alpha, real(n, dp), self%probability)) + 1
            if (high_index >= 1 .and. high_index <= n) then
                ci%high = x(high_index)
            else
                ci%high = quiet_nan(level)
            end if
        case ('greater')
            alpha = 1.0_dp - level
            low_index = int(binomial_ppf(alpha, real(n, dp), self%probability))
            if (low_index >= 1 .and. low_index <= n) then
                ci%low = x(low_index)
            else
                ci%low = quiet_nan(level)
            end if
            ci%high = positive_infinity(level)
        case default
            alpha = 0.5_dp * (1.0_dp - level)
            low_index = int(binomial_ppf(alpha, real(n, dp), self%probability))
            high_index = int(binomial_isf(alpha, real(n, dp), self%probability)) + 1
            if (low_index >= 1 .and. low_index <= n) then
                ci%low = x(low_index)
            else
                ci%low = quiet_nan(level)
            end if
            if (high_index >= 1 .and. high_index <= n) then
                ci%high = x(high_index)
            else
                ci%high = quiet_nan(level)
            end if
        end select
    end function quantile_test_confidence_interval

    function order_statistic(base_cdf, base_pdf, base_ppf, r, n, parameters) result(dist)
        procedure(order_statistic_callback) :: base_cdf !! base-distribution CDF callback
        procedure(order_statistic_callback) :: base_pdf !! base-distribution PDF callback
        procedure(order_statistic_callback) :: base_ppf !! base-distribution PPF callback
        integer, intent(in) :: r !! order-statistic rank
        integer, intent(in) :: n !! sample size
        real(dp), intent(in), optional :: parameters(:) !! callback parameter vector
        type(order_statistic_distribution) :: dist

        if (r < 0 .or. n < 0) then
            error stop 'order_statistic: r and n must be nonnegative integers'
        end if
        dist%rank = r
        dist%sample_size = n
        dist%base_cdf => base_cdf
        dist%base_pdf => base_pdf
        dist%base_ppf => base_ppf
        if (present(parameters)) then
            allocate(dist%parameters(size(parameters)))
            dist%parameters = parameters
        else
            allocate(dist%parameters(0))
        end if
    end function order_statistic

    pure elemental function order_statistic_cdf(self, x) result(y)
        class(order_statistic_distribution), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: u

        if (.not. valid_order_statistic(self)) then
            y = quiet_nan(x)
            return
        end if
        u = self%base_cdf(x, self%parameters)
        y = beta_cdf(u, real(self%rank, dp), real(self%sample_size - self%rank + 1, dp))
    end function order_statistic_cdf

    pure elemental function order_statistic_sf(self, x) result(y)
        class(order_statistic_distribution), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: u

        if (.not. valid_order_statistic(self)) then
            y = quiet_nan(x)
            return
        end if
        u = self%base_cdf(x, self%parameters)
        y = beta_sf(u, real(self%rank, dp), real(self%sample_size - self%rank + 1, dp))
    end function order_statistic_sf

    pure elemental function order_statistic_pdf(self, x) result(y)
        class(order_statistic_distribution), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: u

        if (.not. valid_order_statistic(self)) then
            y = quiet_nan(x)
            return
        end if
        u = self%base_cdf(x, self%parameters)
        y = beta_pdf(u, real(self%rank, dp), real(self%sample_size - self%rank + 1, dp)) * &
            self%base_pdf(x, self%parameters)
    end function order_statistic_pdf

    pure elemental function order_statistic_ppf(self, probability) result(x)
        class(order_statistic_distribution), intent(in) :: self
        real(dp), intent(in) :: probability
        real(dp) :: x
        real(dp) :: u

        if (.not. valid_order_statistic(self)) then
            x = quiet_nan(probability)
            return
        end if
        u = beta_ppf(probability, real(self%rank, dp), real(self%sample_size - self%rank + 1, dp))
        x = self%base_ppf(u, self%parameters)
    end function order_statistic_ppf

    pure elemental function order_statistic_isf(self, probability) result(x)
        class(order_statistic_distribution), intent(in) :: self
        real(dp), intent(in) :: probability
        real(dp) :: x
        real(dp) :: u

        if (.not. valid_order_statistic(self)) then
            x = quiet_nan(probability)
            return
        end if
        u = beta_isf(probability, real(self%rank, dp), real(self%sample_size - self%rank + 1, dp))
        x = self%base_ppf(u, self%parameters)
    end function order_statistic_isf

    pure logical function valid_order_statistic(self) result(valid)
        class(order_statistic_distribution), intent(in) :: self

        valid = self%rank >= 1 .and. self%sample_size >= 1 .and. self%rank <= self%sample_size .and. &
            associated(self%base_cdf) .and. associated(self%base_pdf) .and. associated(self%base_ppf)
    end function valid_order_statistic

    pure elemental function edf_evaluate_scalar(self, x) result(y)
        class(empirical_distribution_function), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp) :: y

        integer :: i
        integer :: n

        n = size(self%quantiles)
        if (ieee_is_nan(x)) then
            if (self%is_survival) then
                y = 0.0_dp
            else
                y = 1.0_dp
            end if
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                if (self%is_survival) then
                    y = 0.0_dp
                else
                    y = 1.0_dp
                end if
            else
                if (self%is_survival) then
                    y = 1.0_dp
                else
                    y = 0.0_dp
                end if
            end if
            return
        end if
        if (n == 0) then
            if (self%is_survival) then
                y = 1.0_dp
            else
                y = 0.0_dp
            end if
            return
        end if
        if (x < self%quantiles(1)) then
            if (self%is_survival) then
                y = 1.0_dp
            else
                y = 0.0_dp
            end if
            return
        end if

        y = self%probabilities(n)
        do i = 1, n - 1
            if (x < self%quantiles(i + 1)) then
                y = self%probabilities(i)
                return
            end if
        end do
    end function edf_evaluate_scalar

    function edf_confidence_interval_method(self, confidence_level, method) result(ci)
        class(empirical_distribution_function), intent(in) :: self
        real(dp), intent(in), optional :: confidence_level !! confidence level in [0,1], default 0.95
        character(len=*), intent(in), optional :: method !! linear or log-log
        type(ecdf_confidence_interval) :: ci

        real(dp), allocatable :: high(:)
        real(dp) :: level
        real(dp), allocatable :: low(:)
        character(len=7) :: selected

        if (.not. self%has_variance_data) then
            error stop 'ecdf confidence_interval: bounds do not support nested confidence intervals'
        end if

        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        if (.not. (level >= 0.0_dp .and. level <= 1.0_dp)) then
            error stop 'ecdf confidence_interval: confidence level must be in [0,1]'
        end if
        selected = 'linear'
        if (present(method)) selected = lower_ascii(adjustl(method))
        if (trim(selected) /= 'linear' .and. trim(selected) /= 'log-log') then
            error stop 'ecdf confidence_interval: method must be linear or log-log'
        end if

        allocate(low(size(self%probabilities)), high(size(self%probabilities)))
        if (trim(selected) == 'linear') then
            call linear_greenwood(self, level, low, high)
        else
            call loglog_greenwood(self, level, low, high)
        end if
        call clip_unit_interval(low)
        call clip_unit_interval(high)
        call set_edf(ci%low, self%quantiles, low, [integer ::], [integer ::], &
            self%is_survival, .false.)
        call set_edf(ci%high, self%quantiles, high, [integer ::], [integer ::], &
            self%is_survival, .false.)
    end function edf_confidence_interval_method

    subroutine linear_greenwood(self, level, low, high)
        class(empirical_distribution_function), intent(in) :: self
        real(dp), intent(in) :: level
        real(dp), intent(out) :: low(:)
        real(dp), intent(out) :: high(:)

        real(dp) :: cumulative
        integer :: d
        integer :: i
        integer :: n
        real(dp) :: se
        real(dp) :: sf
        real(dp) :: var
        real(dp) :: z
        real(dp) :: zse

        z = normal_ppf(0.5_dp + 0.5_dp * level)
        cumulative = 0.0_dp
        do i = 1, size(self%probabilities)
            d = self%events(i)
            n = self%at_risk(i)
            call accumulate_greenwood(cumulative, n, d)
            if (self%is_survival) then
                sf = self%probabilities(i)
            else
                sf = 1.0_dp - self%probabilities(i)
            end if
            var = product_with_infinity(sf * sf, cumulative)
            if (ieee_is_nan(var)) then
                se = quiet_nan(level)
            else if (var < 0.0_dp) then
                se = quiet_nan(level)
            else
                se = sqrt(var)
            end if
            zse = product_with_infinity(z, se)
            low(i) = safe_subtract(self%probabilities(i), zse)
            high(i) = safe_add(self%probabilities(i), zse)
        end do
    end subroutine linear_greenwood

    subroutine loglog_greenwood(self, level, low, high)
        class(empirical_distribution_function), intent(in) :: self
        real(dp), intent(in) :: level
        real(dp), intent(out) :: low(:)
        real(dp), intent(out) :: high(:)

        real(dp) :: cumulative
        integer :: d
        real(dp) :: denom
        real(dp) :: high_sf
        integer :: i
        real(dp) :: lnl
        real(dp) :: low_sf
        integer :: n
        real(dp) :: se
        real(dp) :: sf
        real(dp) :: var
        real(dp) :: z
        real(dp) :: zse

        z = normal_ppf(0.5_dp + 0.5_dp * level)
        cumulative = 0.0_dp
        do i = 1, size(self%probabilities)
            d = self%events(i)
            n = self%at_risk(i)
            call accumulate_greenwood(cumulative, n, d)
            if (self%is_survival) then
                sf = self%probabilities(i)
            else
                sf = 1.0_dp - self%probabilities(i)
            end if

            if (sf > 0.0_dp .and. sf < 1.0_dp) then
                denom = log(sf) * log(sf)
                var = safe_ratio(cumulative, denom)
                se = sqrt(var)
                lnl = log(-log(sf))
                zse = product_with_infinity(z, se)
                low_sf = exp(-exp(safe_add(lnl, zse)))
                high_sf = exp(-exp(safe_subtract(lnl, zse)))
            else
                call loglog_boundary(sf, cumulative, z, low_sf, high_sf)
            end if

            if (self%is_survival) then
                low(i) = low_sf
                high(i) = high_sf
            else
                low(i) = 1.0_dp - high_sf
                high(i) = 1.0_dp - low_sf
            end if
        end do
    end subroutine loglog_greenwood

    subroutine loglog_boundary(sf, cumulative, z, low_sf, high_sf)
        real(dp), intent(in) :: sf
        real(dp), intent(in) :: cumulative
        real(dp), intent(in) :: z
        real(dp), intent(out) :: low_sf
        real(dp), intent(out) :: high_sf

        if (real_equal(sf, 1.0_dp)) then
            if (real_equal(cumulative, 0.0_dp)) then
                low_sf = quiet_nan(sf)
                high_sf = quiet_nan(sf)
            else
                low_sf = 0.0_dp
                high_sf = 1.0_dp
            end if
        else if (real_equal(sf, 0.0_dp)) then
            if (.not. ieee_is_finite(cumulative)) then
                low_sf = quiet_nan(sf)
                high_sf = quiet_nan(sf)
            else
                low_sf = 0.0_dp
                high_sf = 0.0_dp
            end if
        else
            low_sf = quiet_nan(z)
            high_sf = quiet_nan(z)
        end if
    end subroutine loglog_boundary

    subroutine accumulate_greenwood(cumulative, n, d)
        real(dp), intent(inout) :: cumulative
        integer, intent(in) :: n
        integer, intent(in) :: d

        real(dp) :: term

        if (.not. ieee_is_finite(cumulative)) return
        if (d == 0) return
        if (n == d) then
            cumulative = positive_infinity(0.0_dp)
        else if (n > 0 .and. n > d) then
            term = real(d, dp) / (real(n, dp) * real(n - d, dp))
            cumulative = cumulative + term
        else
            cumulative = quiet_nan(0.0_dp)
        end if
    end subroutine accumulate_greenwood

    pure function product_with_infinity(a, b) result(y)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp) :: y

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            y = quiet_nan(a)
            return
        end if
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            if (ieee_is_finite(a)) then
                if (real_equal(a, 0.0_dp)) then
                    y = quiet_nan(a)
                else if ((a < 0.0_dp) .neqv. (b < 0.0_dp)) then
                    y = negative_infinity(0.0_dp)
                else
                    y = positive_infinity(0.0_dp)
                end if
            else if (ieee_is_finite(b)) then
                if (real_equal(b, 0.0_dp)) then
                    y = quiet_nan(b)
                else if ((a < 0.0_dp) .neqv. (b < 0.0_dp)) then
                    y = negative_infinity(0.0_dp)
                else
                    y = positive_infinity(0.0_dp)
                end if
            else if ((a < 0.0_dp) .neqv. (b < 0.0_dp)) then
                y = negative_infinity(0.0_dp)
            else
                y = positive_infinity(0.0_dp)
            end if
            return
        end if
        y = a * b
    end function product_with_infinity

    pure function safe_ratio(a, b) result(y)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp) :: y

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            y = quiet_nan(a)
            return
        end if
        if (real_equal(b, 0.0_dp)) then
            if (real_equal(a, 0.0_dp)) then
                y = quiet_nan(a)
            else if (a < 0.0_dp) then
                y = negative_infinity(a)
            else
                y = positive_infinity(a)
            end if
            return
        end if
        if (.not. ieee_is_finite(a) .and. .not. ieee_is_finite(b)) then
            y = quiet_nan(a)
            return
        end if
        y = a / b
    end function safe_ratio

    pure function safe_add(a, b) result(y)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp) :: y

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            y = quiet_nan(a)
            return
        end if
        if (.not. ieee_is_finite(a) .and. .not. ieee_is_finite(b)) then
            if ((a < 0.0_dp) .neqv. (b < 0.0_dp)) then
                y = quiet_nan(a)
                return
            end if
        end if
        y = a + b
    end function safe_add

    pure function safe_subtract(a, b) result(y)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp) :: y

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            y = quiet_nan(a)
            return
        end if
        if (.not. ieee_is_finite(a) .and. .not. ieee_is_finite(b)) then
            if ((a < 0.0_dp) .eqv. (b < 0.0_dp)) then
                y = quiet_nan(a)
                return
            end if
        end if
        y = a - b
    end function safe_subtract

    subroutine clip_unit_interval(x)
        real(dp), intent(inout) :: x(:)
        integer :: i

        do i = 1, size(x)
            if (.not. ieee_is_nan(x(i))) x(i) = min(1.0_dp, max(0.0_dp, x(i)))
        end do
    end subroutine clip_unit_interval

    subroutine set_edf(edf, quantiles, probabilities, at_risk, events, is_survival, has_variance_data)
        type(empirical_distribution_function), intent(out) :: edf
        real(dp), intent(in) :: quantiles(:)
        real(dp), intent(in) :: probabilities(:)
        integer, intent(in) :: at_risk(:)
        integer, intent(in) :: events(:)
        logical, intent(in) :: is_survival
        logical, intent(in) :: has_variance_data

        edf%quantiles = quantiles
        edf%probabilities = probabilities
        edf%at_risk = at_risk
        edf%events = events
        edf%is_survival = is_survival
        edf%has_variance_data = has_variance_data
    end subroutine set_edf

    pure logical function all_finite(x) result(ok)
        real(dp), intent(in) :: x(:)
        integer :: i

        ok = .true.
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
                ok = .false.
                return
            end if
        end do
    end function all_finite

    pure logical function real_equal(a, b) result(equal)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b

        equal = a <= b .and. b <= a
    end function real_equal

    subroutine sort_real(values)
        real(dp), intent(inout) :: values(:)
        real(dp), allocatable :: workspace(:)

        if (size(values) <= 1) return
        allocate(workspace(size(values)))
        call merge_sort_real(values, workspace, 1, size(values))
    end subroutine sort_real

    recursive subroutine merge_sort_real(values, workspace, left, right)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(inout) :: workspace(:)
        integer, intent(in) :: left
        integer, intent(in) :: right

        integer :: i
        integer :: j
        integer :: k
        integer :: middle

        if (left >= right) return
        middle = (left + right) / 2
        call merge_sort_real(values, workspace, left, middle)
        call merge_sort_real(values, workspace, middle + 1, right)
        i = left
        j = middle + 1
        do k = left, right
            if (i > middle) then
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

    subroutine sort_nan_last(values)
        real(dp), intent(inout) :: values(:)
        real(dp), allocatable :: finite_values(:)
        integer :: finite_count
        integer :: i
        integer :: j

        finite_count = count(.not. ieee_is_nan(values))
        allocate(finite_values(finite_count))
        j = 0
        do i = 1, size(values)
            if (.not. ieee_is_nan(values(i))) then
                j = j + 1
                finite_values(j) = values(i)
            end if
        end do
        call sort_real(finite_values)
        if (finite_count > 0) values(1:finite_count) = finite_values
        do i = finite_count + 1, size(values)
            values(i) = quiet_nan(0.0_dp)
        end do
    end subroutine sort_nan_last

    subroutine stable_sort_indices(values, order)
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: order(size(values))
        integer, allocatable :: workspace(:)
        integer :: i

        do i = 1, size(values)
            order(i) = i
        end do
        if (size(values) <= 1) return
        allocate(workspace(size(values)))
        call merge_sort_indices(values, order, workspace, 1, size(values))
    end subroutine stable_sort_indices

    recursive subroutine merge_sort_indices(values, order, workspace, left, right)
        real(dp), intent(in) :: values(:)
        integer, intent(inout) :: order(:)
        integer, intent(inout) :: workspace(:)
        integer, intent(in) :: left
        integer, intent(in) :: right

        integer :: i
        integer :: j
        integer :: k
        integer :: middle

        if (left >= right) return
        middle = (left + right) / 2
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

    pure function normalize_alternative(alternative) result(value)
        character(len=*), intent(in) :: alternative
        character(len=9) :: value
        character(len=len(alternative)) :: lowered

        lowered = lower_ascii(alternative)
        value = adjustl(lowered)
    end function normalize_alternative

    pure logical function valid_alternative(alternative) result(valid)
        character(len=*), intent(in) :: alternative

        valid = trim(alternative) == 'two-sided' .or. trim(alternative) == 'less' .or. &
            trim(alternative) == 'greater'
    end function valid_alternative

    pure function lower_ascii(text) result(lowered)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: lowered
        integer :: code
        integer :: i

        lowered = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lowered(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lower_ascii

end module scifort_empirical_distribution
