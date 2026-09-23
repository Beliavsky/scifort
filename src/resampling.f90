! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_resampling
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_descriptive, only : mean, quantile, standard_deviation
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_normal, only : normal_cdf, normal_ppf
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer(int64), parameter :: integer_grid_size = 4503599627370496_int64

    abstract interface
        pure function one_sample_statistic(sample) result(value)
            import dp
            real(dp), intent(in) :: sample(:) !! resampled observations
            real(dp) :: value
        end function one_sample_statistic

        pure function two_sample_statistic(sample_x, sample_y) result(value)
            import dp
            real(dp), intent(in) :: sample_x(:) !! first resampled group
            real(dp), intent(in) :: sample_y(:) !! second resampled group
            real(dp) :: value
        end function two_sample_statistic

        subroutine sample_generator(state, sample)
            import dp, rng_state
            type(rng_state), intent(inout) :: state !! explicit generator state advanced by sampling
            real(dp), intent(out) :: sample(:) !! generated observations
        end subroutine sample_generator

        pure function one_sample_pvalue(sample) result(value)
            import dp
            real(dp), intent(in) :: sample(:) !! generated sample passed to a hypothesis test
            real(dp) :: value
        end function one_sample_pvalue

        pure function two_sample_pvalue(sample_x, sample_y) result(value)
            import dp
            real(dp), intent(in) :: sample_x(:) !! first generated sample
            real(dp), intent(in) :: sample_y(:) !! second generated sample
            real(dp) :: value
        end function two_sample_pvalue
    end interface

    type, public :: bootstrap_result
        real(dp) :: statistic !! statistic evaluated on the original sample(s)
        real(dp) :: standard_error !! sample standard deviation of bootstrap statistics
        real(dp) :: confidence_low !! lower confidence endpoint
        real(dp) :: confidence_high !! upper confidence endpoint
        integer :: n_resamples !! number of bootstrap resamples generated
        real(dp), allocatable :: bootstrap_distribution(:) !! bootstrap statistics in draw order
    end type bootstrap_result

    type, public :: permutation_test_result
        real(dp) :: statistic !! statistic evaluated on the original samples
        real(dp) :: pvalue !! exact or randomized permutation p-value
        integer :: n_resamples !! number of permutation statistics in the null distribution
        logical :: exact !! true when all distinct permutations were enumerated
        real(dp), allocatable :: null_distribution(:) !! permutation null distribution
    end type permutation_test_result

    type, public :: monte_carlo_test_result
        real(dp) :: statistic !! statistic evaluated on the observed sample(s)
        real(dp) :: pvalue !! Monte Carlo p-value with plus-one adjustment
        integer :: n_resamples !! number of simulated null samples
        real(dp), allocatable :: null_distribution(:) !! simulated null distribution
    end type monte_carlo_test_result

    type, public :: power_result
        real(dp) :: power !! fraction of simulated p-values below the significance level
        real(dp) :: significance !! significance threshold used for the simulation
        integer :: n_resamples !! number of simulated experiments
        real(dp), allocatable :: pvalues(:) !! p-values generated under the alternative
    end type power_result

    interface bootstrap
        module procedure bootstrap_one_sample
        module procedure bootstrap_two_sample
    end interface bootstrap

    interface monte_carlo_test
        module procedure monte_carlo_test_one_sample
        module procedure monte_carlo_test_two_sample
    end interface monte_carlo_test

    interface power
        module procedure power_one_sample
        module procedure power_two_sample
    end interface power

    public :: bootstrap
    public :: bootstrap_percentile
    public :: monte_carlo_test
    public :: permutation_test
    public :: power

contains

    function bootstrap_percentile(x, statistic, state, n_resamples, confidence_level) &
            result(result_value)
        real(dp), intent(in) :: x(:) !! original observations, nonempty
        procedure(one_sample_statistic) :: statistic !! statistic applied to each sample
        type(rng_state), intent(inout) :: state !! explicit generator state advanced on valid calls
        integer, intent(in), optional :: n_resamples !! number of bootstrap draws, default 9999
        real(dp), intent(in), optional :: confidence_level !! level in (0, 1), default 0.95
        type(bootstrap_result) :: result_value

        integer :: draws
        real(dp) :: level

        draws = 9999
        if (present(n_resamples)) draws = n_resamples
        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        result_value = bootstrap_one_sample(x, statistic, state, draws, level, &
            'two-sided', 'percentile', .true.)
    end function bootstrap_percentile

    function bootstrap_one_sample(x, statistic, state, n_resamples, confidence_level, &
            alternative, method, allow_singleton) result(result_value)
        real(dp), intent(in) :: x(:) !! original observations
        procedure(one_sample_statistic) :: statistic !! statistic applied to each resample
        type(rng_state), intent(inout) :: state !! explicit generator state advanced by resampling
        integer, intent(in), optional :: n_resamples !! number of bootstrap resamples, default 9999
        real(dp), intent(in), optional :: confidence_level !! confidence level in (0,1), default 0.95
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        character(len=*), intent(in), optional :: method !! percentile, basic, or bca
        logical, intent(in), optional :: allow_singleton !! compatibility flag used by bootstrap_percentile
        type(bootstrap_result) :: result_value

        character(len=10) :: selected_method
        character(len=9) :: selected_alternative
        integer :: b
        integer :: draws
        integer :: i
        integer :: selected
        integer :: minimum_size
        real(dp) :: level
        real(dp), allocatable :: sample(:)

        draws = 9999
        if (present(n_resamples)) draws = n_resamples
        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = alternative
        selected_method = 'bca'
        if (present(method)) selected_method = lower_method(method)
        minimum_size = 2
        if (present(allow_singleton)) then
            if (allow_singleton) minimum_size = 1
        end if

        if (size(x) < minimum_size .or. draws < 2 .or. .not. valid_level(level) .or. &
                .not. valid_alternative(selected_alternative) .or. &
                .not. valid_bootstrap_method(selected_method)) then
            call set_invalid_bootstrap(result_value)
            return
        end if

        result_value%statistic = statistic(x)
        if (ieee_is_nan(result_value%statistic)) then
            call set_invalid_bootstrap(result_value)
            return
        end if

        allocate(sample(size(x)))
        allocate(result_value%bootstrap_distribution(draws))
        do b = 1, draws
            do i = 1, size(x)
                selected = unbiased_index(state, size(x))
                sample(i) = x(selected)
            end do
            result_value%bootstrap_distribution(b) = statistic(sample)
            if (ieee_is_nan(result_value%bootstrap_distribution(b))) then
                call set_invalid_bootstrap(result_value)
                return
            end if
        end do

        result_value%standard_error = standard_deviation(result_value%bootstrap_distribution)
        call bootstrap_interval_one(x, statistic, result_value%bootstrap_distribution, &
            result_value%statistic, level, selected_alternative, selected_method, &
            result_value%confidence_low, result_value%confidence_high)
        result_value%n_resamples = draws
    end function bootstrap_one_sample

    function bootstrap_two_sample(x, y, statistic, state, n_resamples, confidence_level, &
            alternative, method, paired) result(result_value)
        real(dp), intent(in) :: x(:) !! first original sample
        real(dp), intent(in) :: y(:) !! second original sample
        procedure(two_sample_statistic) :: statistic !! statistic applied to each resample pair
        type(rng_state), intent(inout) :: state !! explicit generator state advanced by resampling
        integer, intent(in), optional :: n_resamples !! number of bootstrap resamples, default 9999
        real(dp), intent(in), optional :: confidence_level !! confidence level in (0,1), default 0.95
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        character(len=*), intent(in), optional :: method !! percentile, basic, or bca
        logical, intent(in), optional :: paired !! resample common indices when true, default false
        type(bootstrap_result) :: result_value

        character(len=10) :: selected_method
        character(len=9) :: selected_alternative
        integer :: b
        integer :: draws
        integer :: i
        integer :: selected
        logical :: use_paired
        real(dp) :: level
        real(dp), allocatable :: sample_x(:)
        real(dp), allocatable :: sample_y(:)

        draws = 9999
        if (present(n_resamples)) draws = n_resamples
        level = 0.95_dp
        if (present(confidence_level)) level = confidence_level
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = alternative
        selected_method = 'bca'
        if (present(method)) selected_method = lower_method(method)
        use_paired = .false.
        if (present(paired)) use_paired = paired

        if (size(x) < 2 .or. size(y) < 2 .or. draws < 2 .or. .not. valid_level(level) .or. &
                .not. valid_alternative(selected_alternative) .or. &
                .not. valid_bootstrap_method(selected_method) .or. &
                (use_paired .and. size(x) /= size(y))) then
            call set_invalid_bootstrap(result_value)
            return
        end if

        result_value%statistic = statistic(x, y)
        if (ieee_is_nan(result_value%statistic)) then
            call set_invalid_bootstrap(result_value)
            return
        end if

        allocate(sample_x(size(x)))
        allocate(sample_y(size(y)))
        allocate(result_value%bootstrap_distribution(draws))
        do b = 1, draws
            if (use_paired) then
                do i = 1, size(x)
                    selected = unbiased_index(state, size(x))
                    sample_x(i) = x(selected)
                    sample_y(i) = y(selected)
                end do
            else
                do i = 1, size(x)
                    selected = unbiased_index(state, size(x))
                    sample_x(i) = x(selected)
                end do
                do i = 1, size(y)
                    selected = unbiased_index(state, size(y))
                    sample_y(i) = y(selected)
                end do
            end if
            result_value%bootstrap_distribution(b) = statistic(sample_x, sample_y)
            if (ieee_is_nan(result_value%bootstrap_distribution(b))) then
                call set_invalid_bootstrap(result_value)
                return
            end if
        end do

        result_value%standard_error = standard_deviation(result_value%bootstrap_distribution)
        call bootstrap_interval_two(x, y, statistic, result_value%bootstrap_distribution, &
            result_value%statistic, level, selected_alternative, selected_method, use_paired, &
            result_value%confidence_low, result_value%confidence_high)
        result_value%n_resamples = draws
    end function bootstrap_two_sample

    subroutine bootstrap_interval_one(x, statistic, distribution, observed, level, &
            alternative, method, ci_low, ci_high)
        real(dp), intent(in) :: x(:) !! original observations
        procedure(one_sample_statistic) :: statistic !! statistic used for the bootstrap
        real(dp), intent(in) :: distribution(:) !! bootstrap statistic distribution
        real(dp), intent(in) :: observed !! statistic of the original sample
        real(dp), intent(in) :: level !! requested confidence level
        character(len=*), intent(in) :: alternative !! requested confidence-interval alternative
        character(len=*), intent(in) :: method !! normalized bootstrap method name
        real(dp), intent(out) :: ci_low !! lower confidence endpoint
        real(dp), intent(out) :: ci_high !! upper confidence endpoint

        real(dp) :: alpha
        real(dp) :: alpha_low
        real(dp) :: alpha_high
        real(dp) :: acceleration

        alpha = bootstrap_alpha(level, alternative)
        if (method == 'bca') then
            acceleration = bca_acceleration_one(x, statistic)
            call bca_probabilities(distribution, observed, alpha, acceleration, &
                alpha_low, alpha_high)
        else
            alpha_low = alpha
            alpha_high = 1.0_dp - alpha
        end if
        call interval_from_probabilities(distribution, observed, alpha_low, alpha_high, &
            alternative, method, ci_low, ci_high)
    end subroutine bootstrap_interval_one

    subroutine bootstrap_interval_two(x, y, statistic, distribution, observed, level, &
            alternative, method, paired, ci_low, ci_high)
        real(dp), intent(in) :: x(:) !! first original sample
        real(dp), intent(in) :: y(:) !! second original sample
        procedure(two_sample_statistic) :: statistic !! statistic used for the bootstrap
        real(dp), intent(in) :: distribution(:) !! bootstrap statistic distribution
        real(dp), intent(in) :: observed !! statistic of the original samples
        real(dp), intent(in) :: level !! requested confidence level
        character(len=*), intent(in) :: alternative !! requested confidence-interval alternative
        character(len=*), intent(in) :: method !! normalized bootstrap method name
        logical, intent(in) :: paired !! whether observations are paired
        real(dp), intent(out) :: ci_low !! lower confidence endpoint
        real(dp), intent(out) :: ci_high !! upper confidence endpoint

        real(dp) :: alpha
        real(dp) :: alpha_low
        real(dp) :: alpha_high
        real(dp) :: acceleration

        alpha = bootstrap_alpha(level, alternative)
        if (method == 'bca') then
            acceleration = bca_acceleration_two(x, y, statistic, paired)
            call bca_probabilities(distribution, observed, alpha, acceleration, &
                alpha_low, alpha_high)
        else
            alpha_low = alpha
            alpha_high = 1.0_dp - alpha
        end if
        call interval_from_probabilities(distribution, observed, alpha_low, alpha_high, &
            alternative, method, ci_low, ci_high)
    end subroutine bootstrap_interval_two

    subroutine interval_from_probabilities(distribution, observed, alpha_low, alpha_high, &
            alternative, method, ci_low, ci_high)
        real(dp), intent(in) :: distribution(:) !! bootstrap statistic distribution
        real(dp), intent(in) :: observed !! statistic of the original data
        real(dp), intent(in) :: alpha_low !! lower quantile probability
        real(dp), intent(in) :: alpha_high !! upper quantile probability
        character(len=*), intent(in) :: alternative !! requested confidence-interval alternative
        character(len=*), intent(in) :: method !! normalized bootstrap method name
        real(dp), intent(out) :: ci_low !! lower confidence endpoint
        real(dp), intent(out) :: ci_high !! upper confidence endpoint

        real(dp) :: percentile_low
        real(dp) :: percentile_high

        if (.not. valid_probability(alpha_low) .or. .not. valid_probability(alpha_high)) then
            ci_low = quiet_nan(0.0_dp)
            ci_high = quiet_nan(0.0_dp)
            return
        end if
        percentile_low = quantile(distribution, alpha_low)
        percentile_high = quantile(distribution, alpha_high)
        if (method == 'basic') then
            ci_low = 2.0_dp * observed - percentile_high
            ci_high = 2.0_dp * observed - percentile_low
        else
            ci_low = percentile_low
            ci_high = percentile_high
        end if
        if (alternative == 'less') ci_low = negative_infinity(0.0_dp)
        if (alternative == 'greater') ci_high = positive_infinity(0.0_dp)
    end subroutine interval_from_probabilities

    subroutine bca_probabilities(distribution, observed, alpha, acceleration, &
            alpha_low, alpha_high)
        real(dp), intent(in) :: distribution(:) !! bootstrap statistic distribution
        real(dp), intent(in) :: observed !! statistic of the original data
        real(dp), intent(in) :: alpha !! lower base tail probability
        real(dp), intent(in) :: acceleration !! BCa acceleration constant
        real(dp), intent(out) :: alpha_low !! BCa lower percentile probability
        real(dp), intent(out) :: alpha_high !! BCa upper percentile probability

        integer :: count_le
        integer :: count_lt
        real(dp) :: denominator
        real(dp) :: percentile
        real(dp) :: term
        real(dp) :: z0
        real(dp) :: z_alpha

        if (.not. ieee_is_finite(acceleration)) then
            alpha_low = quiet_nan(0.0_dp)
            alpha_high = quiet_nan(0.0_dp)
            return
        end if
        count_lt = count(distribution < observed)
        count_le = count(distribution <= observed)
        percentile = real(count_lt + count_le, dp) / real(2 * size(distribution), dp)
        if (percentile <= 0.0_dp .or. percentile >= 1.0_dp) then
            alpha_low = quiet_nan(0.0_dp)
            alpha_high = quiet_nan(0.0_dp)
            return
        end if
        z0 = normal_ppf(percentile)
        z_alpha = normal_ppf(alpha)

        term = z0 + z_alpha
        denominator = 1.0_dp - acceleration * term
        if (denominator == 0.0_dp) then
            alpha_low = quiet_nan(0.0_dp)
        else
            alpha_low = normal_cdf(z0 + term / denominator)
        end if

        term = z0 - z_alpha
        denominator = 1.0_dp - acceleration * term
        if (denominator == 0.0_dp) then
            alpha_high = quiet_nan(0.0_dp)
        else
            alpha_high = normal_cdf(z0 + term / denominator)
        end if
    end subroutine bca_probabilities

    function bca_acceleration_one(x, statistic) result(acceleration)
        real(dp), intent(in) :: x(:) !! original observations
        procedure(one_sample_statistic) :: statistic !! statistic used for the bootstrap
        real(dp) :: acceleration

        integer :: i
        integer :: j
        integer :: k
        real(dp) :: denominator
        real(dp) :: jackknife_mean
        real(dp) :: numerator
        real(dp), allocatable :: jackknife(:)
        real(dp), allocatable :: leave_one_out(:)
        real(dp), allocatable :: influence(:)

        allocate(jackknife(size(x)))
        allocate(leave_one_out(size(x) - 1))
        do i = 1, size(x)
            k = 0
            do j = 1, size(x)
                if (j == i) cycle
                k = k + 1
                leave_one_out(k) = x(j)
            end do
            jackknife(i) = statistic(leave_one_out)
        end do
        if (any(ieee_is_nan(jackknife))) then
            acceleration = quiet_nan(0.0_dp)
            return
        end if
        jackknife_mean = mean(jackknife)
        allocate(influence(size(x)))
        influence = real(size(x) - 1, dp) * (jackknife_mean - jackknife)
        numerator = sum(influence**3)
        denominator = sum(influence**2)
        if (denominator <= 0.0_dp) then
            acceleration = quiet_nan(0.0_dp)
        else
            acceleration = numerator / (6.0_dp * denominator**1.5_dp)
        end if
    end function bca_acceleration_one

    function bca_acceleration_two(x, y, statistic, paired) result(acceleration)
        real(dp), intent(in) :: x(:) !! first original sample
        real(dp), intent(in) :: y(:) !! second original sample
        procedure(two_sample_statistic) :: statistic !! statistic used for the bootstrap
        logical, intent(in) :: paired !! whether observations are paired
        real(dp) :: acceleration

        real(dp) :: denominator
        real(dp) :: numerator

        if (paired) then
            call bca_terms_paired(x, y, statistic, numerator, denominator)
        else
            call bca_terms_unpaired(x, y, statistic, numerator, denominator)
        end if
        if (.not. ieee_is_finite(numerator) .or. denominator <= 0.0_dp .or. &
                .not. ieee_is_finite(denominator)) then
            acceleration = quiet_nan(0.0_dp)
        else
            acceleration = numerator / (6.0_dp * denominator**1.5_dp)
        end if
    end function bca_acceleration_two

    subroutine bca_terms_paired(x, y, statistic, numerator, denominator)
        real(dp), intent(in) :: x(:) !! first paired sample
        real(dp), intent(in) :: y(:) !! second paired sample
        procedure(two_sample_statistic) :: statistic !! statistic used for the bootstrap
        real(dp), intent(out) :: numerator !! scaled third-moment numerator
        real(dp), intent(out) :: denominator !! scaled second-moment denominator

        integer :: i
        integer :: j
        integer :: k
        real(dp) :: jackknife_mean
        real(dp), allocatable :: jackknife(:)
        real(dp), allocatable :: leave_x(:)
        real(dp), allocatable :: leave_y(:)
        real(dp), allocatable :: influence(:)

        allocate(jackknife(size(x)))
        allocate(leave_x(size(x) - 1))
        allocate(leave_y(size(y) - 1))
        do i = 1, size(x)
            k = 0
            do j = 1, size(x)
                if (j == i) cycle
                k = k + 1
                leave_x(k) = x(j)
                leave_y(k) = y(j)
            end do
            jackknife(i) = statistic(leave_x, leave_y)
        end do
        if (any(ieee_is_nan(jackknife))) then
            numerator = quiet_nan(0.0_dp)
            denominator = quiet_nan(0.0_dp)
            return
        end if
        jackknife_mean = mean(jackknife)
        allocate(influence(size(x)))
        influence = real(size(x) - 1, dp) * (jackknife_mean - jackknife)
        numerator = sum(influence**3) / real(size(x), dp)**3
        denominator = sum(influence**2) / real(size(x), dp)**2
    end subroutine bca_terms_paired

    subroutine bca_terms_unpaired(x, y, statistic, numerator, denominator)
        real(dp), intent(in) :: x(:) !! first independent sample
        real(dp), intent(in) :: y(:) !! second independent sample
        procedure(two_sample_statistic) :: statistic !! statistic used for the bootstrap
        real(dp), intent(out) :: numerator !! sum of sample-wise scaled third moments
        real(dp), intent(out) :: denominator !! sum of sample-wise scaled second moments

        real(dp) :: denominator_part
        real(dp) :: numerator_part

        call jackknife_terms_first(x, y, statistic, numerator_part, denominator_part)
        numerator = numerator_part
        denominator = denominator_part
        call jackknife_terms_second(x, y, statistic, numerator_part, denominator_part)
        numerator = numerator + numerator_part
        denominator = denominator + denominator_part
    end subroutine bca_terms_unpaired

    subroutine jackknife_terms_first(x, y, statistic, numerator, denominator)
        real(dp), intent(in) :: x(:) !! sample from which observations are deleted
        real(dp), intent(in) :: y(:) !! sample held fixed
        procedure(two_sample_statistic) :: statistic !! statistic used for the bootstrap
        real(dp), intent(out) :: numerator !! scaled third moment of jackknife influence values
        real(dp), intent(out) :: denominator !! scaled second moment of jackknife influence values

        integer :: i
        integer :: j
        integer :: k
        real(dp) :: jackknife_mean
        real(dp), allocatable :: jackknife(:)
        real(dp), allocatable :: leave_x(:)
        real(dp), allocatable :: influence(:)

        allocate(jackknife(size(x)))
        allocate(leave_x(size(x) - 1))
        do i = 1, size(x)
            k = 0
            do j = 1, size(x)
                if (j == i) cycle
                k = k + 1
                leave_x(k) = x(j)
            end do
            jackknife(i) = statistic(leave_x, y)
        end do
        if (any(ieee_is_nan(jackknife))) then
            numerator = quiet_nan(0.0_dp)
            denominator = quiet_nan(0.0_dp)
            return
        end if
        jackknife_mean = mean(jackknife)
        allocate(influence(size(x)))
        influence = real(size(x) - 1, dp) * (jackknife_mean - jackknife)
        numerator = sum(influence**3) / real(size(x), dp)**3
        denominator = sum(influence**2) / real(size(x), dp)**2
    end subroutine jackknife_terms_first

    subroutine jackknife_terms_second(x, y, statistic, numerator, denominator)
        real(dp), intent(in) :: x(:) !! sample held fixed
        real(dp), intent(in) :: y(:) !! sample from which observations are deleted
        procedure(two_sample_statistic) :: statistic !! statistic used for the bootstrap
        real(dp), intent(out) :: numerator !! scaled third moment of jackknife influence values
        real(dp), intent(out) :: denominator !! scaled second moment of jackknife influence values

        integer :: i
        integer :: j
        integer :: k
        real(dp) :: jackknife_mean
        real(dp), allocatable :: jackknife(:)
        real(dp), allocatable :: leave_y(:)
        real(dp), allocatable :: influence(:)

        allocate(jackknife(size(y)))
        allocate(leave_y(size(y) - 1))
        do i = 1, size(y)
            k = 0
            do j = 1, size(y)
                if (j == i) cycle
                k = k + 1
                leave_y(k) = y(j)
            end do
            jackknife(i) = statistic(x, leave_y)
        end do
        if (any(ieee_is_nan(jackknife))) then
            numerator = quiet_nan(0.0_dp)
            denominator = quiet_nan(0.0_dp)
            return
        end if
        jackknife_mean = mean(jackknife)
        allocate(influence(size(y)))
        influence = real(size(y) - 1, dp) * (jackknife_mean - jackknife)
        numerator = sum(influence**3) / real(size(y), dp)**3
        denominator = sum(influence**2) / real(size(y), dp)**2
    end subroutine jackknife_terms_second

    function permutation_test(x, y, statistic, state, n_resamples, alternative, &
            permutation_type) result(result_value)
        real(dp), intent(in) :: x(:) !! first original sample, nonempty
        real(dp), intent(in) :: y(:) !! second original sample, nonempty
        procedure(two_sample_statistic) :: statistic !! statistic applied to each group pair
        type(rng_state), intent(inout) :: state !! explicit generator state advanced for random tests
        integer, intent(in), optional :: n_resamples !! requested resamples, default 9999
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        character(len=*), intent(in), optional :: permutation_type !! independent, samples, or pairings
        type(permutation_test_result) :: result_value

        character(len=11) :: selected_type
        character(len=9) :: selected_alternative
        integer :: draws

        draws = 9999
        if (present(n_resamples)) draws = n_resamples
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = alternative
        selected_type = 'independent'
        if (present(permutation_type)) selected_type = permutation_type

        if (size(x) == 0 .or. size(y) == 0 .or. draws < 1 .or. &
                .not. valid_alternative(selected_alternative) .or. &
                .not. valid_permutation_type(selected_type) .or. &
                (selected_type /= 'independent' .and. size(x) /= size(y))) then
            call set_invalid_permutation(result_value)
            return
        end if

        result_value%statistic = statistic(x, y)
        if (ieee_is_nan(result_value%statistic)) then
            call set_invalid_permutation(result_value)
            return
        end if

        select case (selected_type)
        case ('independent')
            call permutation_independent(x, y, statistic, state, draws, result_value)
        case ('samples')
            call permutation_samples(x, y, statistic, state, draws, result_value)
        case ('pairings')
            call permutation_pairings(x, y, statistic, state, draws, result_value)
        end select
        call permutation_pvalue(result_value, selected_alternative)
    end function permutation_test

    subroutine permutation_independent(x, y, statistic, state, draws, result_value)
        real(dp), intent(in) :: x(:) !! first sample
        real(dp), intent(in) :: y(:) !! second sample
        procedure(two_sample_statistic) :: statistic !! statistic applied to each partition
        type(rng_state), intent(inout) :: state !! explicit RNG state for randomized partitions
        integer, intent(in) :: draws !! requested number of permutation statistics
        type(permutation_test_result), intent(inout) :: result_value !! result under construction

        integer :: b
        integer :: j
        integer :: selected
        integer(int64) :: n_exact
        real(dp) :: swap_value
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: work(:)

        n_exact = combination_count_capped(size(x) + size(y), size(x), int(draws, int64) + 1_int64)
        if (n_exact <= int(draws, int64)) then
            result_value%exact = .true.
            result_value%n_resamples = int(n_exact)
            allocate(result_value%null_distribution(result_value%n_resamples))
            call enumerate_independent_exact(x, y, statistic, result_value%null_distribution)
            return
        end if

        result_value%exact = .false.
        result_value%n_resamples = draws
        allocate(result_value%null_distribution(draws))
        allocate(pooled(size(x) + size(y)))
        allocate(work(size(pooled)))
        pooled(1:size(x)) = x
        pooled(size(x) + 1:) = y
        do b = 1, draws
            work = pooled
            do j = size(work), 2, -1
                selected = unbiased_index(state, j)
                swap_value = work(j)
                work(j) = work(selected)
                work(selected) = swap_value
            end do
            result_value%null_distribution(b) = statistic(&
                work(1:size(x)), work(size(x) + 1:))
        end do
    end subroutine permutation_independent

    subroutine enumerate_independent_exact(x, y, statistic, distribution)
        real(dp), intent(in) :: x(:) !! first original sample
        real(dp), intent(in) :: y(:) !! second original sample
        procedure(two_sample_statistic) :: statistic !! statistic applied to each partition
        real(dp), intent(out) :: distribution(:) !! complete exact null distribution

        integer :: b
        integer :: i
        integer :: j
        integer :: next_index
        integer, allocatable :: combination(:)
        logical, allocatable :: selected(:)
        real(dp), allocatable :: pooled(:)
        real(dp), allocatable :: work_x(:)
        real(dp), allocatable :: work_y(:)

        allocate(pooled(size(x) + size(y)))
        pooled(1:size(x)) = x
        pooled(size(x) + 1:) = y
        allocate(combination(size(x)))
        allocate(selected(size(pooled)))
        allocate(work_x(size(x)))
        allocate(work_y(size(y)))
        do i = 1, size(x)
            combination(i) = i
        end do

        b = 0
        do
            b = b + 1
            selected = .false.
            selected(combination) = .true.
            i = 0
            j = 0
            do next_index = 1, size(pooled)
                if (selected(next_index)) then
                    i = i + 1
                    work_x(i) = pooled(next_index)
                else
                    j = j + 1
                    work_y(j) = pooled(next_index)
                end if
            end do
            distribution(b) = statistic(work_x, work_y)
            if (.not. next_combination(combination, size(pooled))) exit
        end do
    end subroutine enumerate_independent_exact

    subroutine permutation_samples(x, y, statistic, state, draws, result_value)
        real(dp), intent(in) :: x(:) !! first paired sample
        real(dp), intent(in) :: y(:) !! second paired sample
        procedure(two_sample_statistic) :: statistic !! statistic applied after within-pair swaps
        type(rng_state), intent(inout) :: state !! explicit RNG state for randomized swaps
        integer, intent(in) :: draws !! requested number of permutation statistics
        type(permutation_test_result), intent(inout) :: result_value !! result under construction

        integer :: b
        integer :: i
        integer(int64) :: mask
        integer(int64) :: n_exact
        real(dp) :: temporary
        real(dp), allocatable :: work_x(:)
        real(dp), allocatable :: work_y(:)

        n_exact = power_two_count_capped(size(x), int(draws, int64) + 1_int64)
        allocate(work_x(size(x)))
        allocate(work_y(size(y)))
        if (n_exact <= int(draws, int64)) then
            result_value%exact = .true.
            result_value%n_resamples = int(n_exact)
            allocate(result_value%null_distribution(result_value%n_resamples))
            do mask = 0_int64, n_exact - 1_int64
                work_x = x
                work_y = y
                do i = 1, size(x)
                    if (btest(mask, i - 1)) then
                        temporary = work_x(i)
                        work_x(i) = work_y(i)
                        work_y(i) = temporary
                    end if
                end do
                result_value%null_distribution(int(mask) + 1) = statistic(work_x, work_y)
            end do
            return
        end if

        result_value%exact = .false.
        result_value%n_resamples = draws
        allocate(result_value%null_distribution(draws))
        do b = 1, draws
            work_x = x
            work_y = y
            do i = 1, size(x)
                if (rng_uniform(state) < 0.5_dp) then
                    temporary = work_x(i)
                    work_x(i) = work_y(i)
                    work_y(i) = temporary
                end if
            end do
            result_value%null_distribution(b) = statistic(work_x, work_y)
        end do
    end subroutine permutation_samples

    subroutine permutation_pairings(x, y, statistic, state, draws, result_value)
        real(dp), intent(in) :: x(:) !! first paired sample
        real(dp), intent(in) :: y(:) !! second paired sample
        procedure(two_sample_statistic) :: statistic !! association statistic after re-pairing
        type(rng_state), intent(inout) :: state !! explicit RNG state for randomized pairings
        integer, intent(in) :: draws !! requested number of permutation statistics
        type(permutation_test_result), intent(inout) :: result_value !! result under construction

        integer :: b
        integer :: j
        integer :: selected
        integer(int64) :: factorial_n
        integer(int64) :: n_exact
        real(dp) :: temporary
        real(dp), allocatable :: work_x(:)
        real(dp), allocatable :: work_y(:)

        factorial_n = factorial_count_capped(size(x), int(draws, int64) + 1_int64)
        if (factorial_n > int(draws, int64)) then
            n_exact = int(draws, int64) + 1_int64
        else if (factorial_n > 0_int64 .and. factorial_n > int(draws, int64) / factorial_n) then
            n_exact = int(draws, int64) + 1_int64
        else
            n_exact = factorial_n * factorial_n
        end if
        allocate(work_x(size(x)))
        allocate(work_y(size(y)))
        if (n_exact <= int(draws, int64)) then
            result_value%exact = .true.
            result_value%n_resamples = int(n_exact)
            allocate(result_value%null_distribution(result_value%n_resamples))
            call enumerate_pairings_exact(x, y, statistic, result_value%null_distribution)
            return
        end if

        result_value%exact = .false.
        result_value%n_resamples = draws
        allocate(result_value%null_distribution(draws))
        do b = 1, draws
            work_x = x
            work_y = y
            do j = size(work_x), 2, -1
                selected = unbiased_index(state, j)
                temporary = work_x(j)
                work_x(j) = work_x(selected)
                work_x(selected) = temporary
            end do
            do j = size(work_y), 2, -1
                selected = unbiased_index(state, j)
                temporary = work_y(j)
                work_y(j) = work_y(selected)
                work_y(selected) = temporary
            end do
            result_value%null_distribution(b) = statistic(work_x, work_y)
        end do
    end subroutine permutation_pairings

    subroutine enumerate_pairings_exact(x, y, statistic, distribution)
        real(dp), intent(in) :: x(:) !! first paired sample
        real(dp), intent(in) :: y(:) !! second paired sample
        procedure(two_sample_statistic) :: statistic !! association statistic after re-pairing
        real(dp), intent(out) :: distribution(:) !! complete exact null distribution

        integer :: b
        integer :: i
        integer, allocatable :: permutation_x(:)
        integer, allocatable :: permutation_y(:)
        real(dp), allocatable :: work_x(:)
        real(dp), allocatable :: work_y(:)

        allocate(permutation_x(size(x)))
        allocate(permutation_y(size(y)))
        allocate(work_x(size(x)))
        allocate(work_y(size(y)))
        do i = 1, size(x)
            permutation_x(i) = i
            permutation_y(i) = i
        end do
        b = 0
        do
            permutation_y = [(i, i = 1, size(y))]
            do
                b = b + 1
                work_x = x(permutation_x)
                work_y = y(permutation_y)
                distribution(b) = statistic(work_x, work_y)
                if (.not. next_permutation(permutation_y)) exit
            end do
            if (.not. next_permutation(permutation_x)) exit
        end do
    end subroutine enumerate_pairings_exact

    subroutine permutation_pvalue(result_value, alternative)
        type(permutation_test_result), intent(inout) :: result_value !! result containing null distribution
        character(len=*), intent(in) :: alternative !! two-sided, less, or greater

        integer :: adjustment
        integer :: count_greater
        integer :: count_less
        real(dp) :: gamma
        real(dp) :: p_greater
        real(dp) :: p_less

        if (any(ieee_is_nan(result_value%null_distribution))) then
            call set_invalid_permutation(result_value)
            return
        end if
        adjustment = 1
        if (result_value%exact) adjustment = 0
        gamma = abs(100.0_dp * epsilon(1.0_dp) * result_value%statistic)
        count_less = count(result_value%null_distribution <= result_value%statistic + gamma)
        count_greater = count(result_value%null_distribution >= result_value%statistic - gamma)
        p_less = real(count_less + adjustment, dp) / &
            real(result_value%n_resamples + adjustment, dp)
        p_greater = real(count_greater + adjustment, dp) / &
            real(result_value%n_resamples + adjustment, dp)
        select case (alternative)
        case ('less')
            result_value%pvalue = p_less
        case ('greater')
            result_value%pvalue = p_greater
        case ('two-sided')
            result_value%pvalue = min(1.0_dp, 2.0_dp * min(p_less, p_greater))
        end select
    end subroutine permutation_pvalue

    function monte_carlo_test_one_sample(x, rvs, statistic, state, n_resamples, &
            alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! observed sample
        procedure(sample_generator) :: rvs !! null generator producing samples of size size(x)
        procedure(one_sample_statistic) :: statistic !! statistic applied to observed/generated samples
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by null simulation
        integer, intent(in), optional :: n_resamples !! number of null simulations, default 9999
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(monte_carlo_test_result) :: result_value

        character(len=9) :: selected_alternative
        integer :: b
        integer :: draws
        real(dp), allocatable :: sample(:)

        draws = 9999
        if (present(n_resamples)) draws = n_resamples
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = alternative
        if (size(x) == 0 .or. draws < 1 .or. .not. valid_alternative(selected_alternative)) then
            call set_invalid_monte_carlo(result_value)
            return
        end if
        result_value%statistic = statistic(x)
        if (ieee_is_nan(result_value%statistic)) then
            call set_invalid_monte_carlo(result_value)
            return
        end if
        result_value%n_resamples = draws
        allocate(result_value%null_distribution(draws))
        allocate(sample(size(x)))
        do b = 1, draws
            call rvs(state, sample)
            result_value%null_distribution(b) = statistic(sample)
        end do
        call monte_carlo_pvalue(result_value, selected_alternative)
    end function monte_carlo_test_one_sample

    function monte_carlo_test_two_sample(x, y, rvs_x, rvs_y, statistic, state, &
            n_resamples, alternative) result(result_value)
        real(dp), intent(in) :: x(:) !! first observed sample
        real(dp), intent(in) :: y(:) !! second observed sample
        procedure(sample_generator) :: rvs_x !! null generator for the first sample
        procedure(sample_generator) :: rvs_y !! null generator for the second sample
        procedure(two_sample_statistic) :: statistic !! statistic applied to observed/generated samples
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by null simulation
        integer, intent(in), optional :: n_resamples !! number of null simulations, default 9999
        character(len=*), intent(in), optional :: alternative !! two-sided, less, or greater
        type(monte_carlo_test_result) :: result_value

        character(len=9) :: selected_alternative
        integer :: b
        integer :: draws
        real(dp), allocatable :: sample_x(:)
        real(dp), allocatable :: sample_y(:)

        draws = 9999
        if (present(n_resamples)) draws = n_resamples
        selected_alternative = 'two-sided'
        if (present(alternative)) selected_alternative = alternative
        if (size(x) == 0 .or. size(y) == 0 .or. draws < 1 .or. &
                .not. valid_alternative(selected_alternative)) then
            call set_invalid_monte_carlo(result_value)
            return
        end if
        result_value%statistic = statistic(x, y)
        if (ieee_is_nan(result_value%statistic)) then
            call set_invalid_monte_carlo(result_value)
            return
        end if
        result_value%n_resamples = draws
        allocate(result_value%null_distribution(draws))
        allocate(sample_x(size(x)))
        allocate(sample_y(size(y)))
        do b = 1, draws
            call rvs_x(state, sample_x)
            call rvs_y(state, sample_y)
            result_value%null_distribution(b) = statistic(sample_x, sample_y)
        end do
        call monte_carlo_pvalue(result_value, selected_alternative)
    end function monte_carlo_test_two_sample

    subroutine monte_carlo_pvalue(result_value, alternative)
        type(monte_carlo_test_result), intent(inout) :: result_value !! Monte Carlo result under construction
        character(len=*), intent(in) :: alternative !! two-sided, less, or greater

        integer :: count_greater
        integer :: count_less
        real(dp) :: gamma
        real(dp) :: p_greater
        real(dp) :: p_less

        if (any(ieee_is_nan(result_value%null_distribution))) then
            call set_invalid_monte_carlo(result_value)
            return
        end if
        gamma = abs(100.0_dp * epsilon(1.0_dp) * result_value%statistic)
        count_less = count(result_value%null_distribution <= result_value%statistic + gamma)
        count_greater = count(result_value%null_distribution >= result_value%statistic - gamma)
        p_less = real(count_less + 1, dp) / real(result_value%n_resamples + 1, dp)
        p_greater = real(count_greater + 1, dp) / real(result_value%n_resamples + 1, dp)
        select case (alternative)
        case ('less')
            result_value%pvalue = p_less
        case ('greater')
            result_value%pvalue = p_greater
        case ('two-sided')
            result_value%pvalue = min(1.0_dp, 2.0_dp * min(p_less, p_greater))
        end select
    end subroutine monte_carlo_pvalue

    function power_one_sample(rvs, test, n_observations, state, significance, &
            n_resamples) result(result_value)
        procedure(sample_generator) :: rvs !! alternative generator for one sample
        procedure(one_sample_pvalue) :: test !! hypothesis test returning a p-value
        integer, intent(in) :: n_observations !! generated sample size, > 0
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by simulation
        real(dp), intent(in), optional :: significance !! significance level in [0,1], default 0.01
        integer, intent(in), optional :: n_resamples !! number of simulated experiments, default 10000
        type(power_result) :: result_value

        integer :: b
        integer :: draws
        integer :: significant
        real(dp) :: alpha
        real(dp), allocatable :: sample(:)

        draws = 10000
        if (present(n_resamples)) draws = n_resamples
        alpha = 0.01_dp
        if (present(significance)) alpha = significance
        if (n_observations <= 0 .or. draws <= 0 .or. ieee_is_nan(alpha) .or. &
                alpha < 0.0_dp .or. alpha > 1.0_dp) then
            call set_invalid_power(result_value)
            return
        end if
        result_value%significance = alpha
        result_value%n_resamples = draws
        allocate(result_value%pvalues(draws))
        allocate(sample(n_observations))
        significant = 0
        do b = 1, draws
            call rvs(state, sample)
            result_value%pvalues(b) = test(sample)
            if (ieee_is_nan(result_value%pvalues(b))) then
                call set_invalid_power(result_value)
                return
            end if
            if (result_value%pvalues(b) < alpha) significant = significant + 1
        end do
        result_value%power = real(significant, dp) / real(draws, dp)
    end function power_one_sample

    function power_two_sample(rvs_x, rvs_y, test, n_observations_x, n_observations_y, &
            state, significance, n_resamples) result(result_value)
        procedure(sample_generator) :: rvs_x !! alternative generator for the first sample
        procedure(sample_generator) :: rvs_y !! alternative generator for the second sample
        procedure(two_sample_pvalue) :: test !! hypothesis test returning a p-value
        integer, intent(in) :: n_observations_x !! first generated sample size, > 0
        integer, intent(in) :: n_observations_y !! second generated sample size, > 0
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by simulation
        real(dp), intent(in), optional :: significance !! significance level in [0,1], default 0.01
        integer, intent(in), optional :: n_resamples !! number of simulated experiments, default 10000
        type(power_result) :: result_value

        integer :: b
        integer :: draws
        integer :: significant
        real(dp) :: alpha
        real(dp), allocatable :: sample_x(:)
        real(dp), allocatable :: sample_y(:)

        draws = 10000
        if (present(n_resamples)) draws = n_resamples
        alpha = 0.01_dp
        if (present(significance)) alpha = significance
        if (n_observations_x <= 0 .or. n_observations_y <= 0 .or. draws <= 0 .or. &
                ieee_is_nan(alpha) .or. alpha < 0.0_dp .or. alpha > 1.0_dp) then
            call set_invalid_power(result_value)
            return
        end if
        result_value%significance = alpha
        result_value%n_resamples = draws
        allocate(result_value%pvalues(draws))
        allocate(sample_x(n_observations_x))
        allocate(sample_y(n_observations_y))
        significant = 0
        do b = 1, draws
            call rvs_x(state, sample_x)
            call rvs_y(state, sample_y)
            result_value%pvalues(b) = test(sample_x, sample_y)
            if (ieee_is_nan(result_value%pvalues(b))) then
                call set_invalid_power(result_value)
                return
            end if
            if (result_value%pvalues(b) < alpha) significant = significant + 1
        end do
        result_value%power = real(significant, dp) / real(draws, dp)
    end function power_two_sample

    function unbiased_index(state, upper) result(index_value)
        type(rng_state), intent(inout) :: state !! explicit generator state advanced until accepted
        integer, intent(in) :: upper !! inclusive upper bound, >= 1
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

    logical function next_combination(combination, n) result(has_next)
        integer, intent(inout) :: combination(:) !! increasing selected indices, modified to next combination
        integer, intent(in) :: n !! size of the full index set

        integer :: i
        integer :: j

        has_next = .false.
        do i = size(combination), 1, -1
            if (combination(i) < n - size(combination) + i) then
                combination(i) = combination(i) + 1
                do j = i + 1, size(combination)
                    combination(j) = combination(j - 1) + 1
                end do
                has_next = .true.
                return
            end if
        end do
    end function next_combination

    logical function next_permutation(values) result(has_next)
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
    end function next_permutation

    pure function combination_count_capped(n, k, cap) result(value)
        integer, intent(in) :: n !! size of the full set
        integer, intent(in) :: k !! number selected
        integer(int64), intent(in) :: cap !! first value at which exact count is no longer required
        integer(int64) :: value

        integer :: i
        integer :: kk
        integer(int64) :: numerator

        if (k < 0 .or. k > n) then
            value = 0_int64
            return
        end if
        kk = min(k, n - k)
        value = 1_int64
        do i = 1, kk
            numerator = int(n - kk + i, int64)
            if (value > (cap * int(i, int64)) / max(1_int64, numerator)) then
                value = cap
                return
            end if
            value = (value * numerator) / int(i, int64)
            if (value >= cap) then
                value = cap
                return
            end if
        end do
    end function combination_count_capped

    pure function factorial_count_capped(n, cap) result(value)
        integer, intent(in) :: n !! factorial argument, >= 0
        integer(int64), intent(in) :: cap !! first value at which exact count is no longer required
        integer(int64) :: value

        integer :: i

        value = 1_int64
        do i = 2, n
            if (value > cap / int(i, int64)) then
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

    pure function power_two_count_capped(n, cap) result(value)
        integer, intent(in) :: n !! nonnegative exponent
        integer(int64), intent(in) :: cap !! first value at which exact count is no longer required
        integer(int64) :: value

        integer :: i

        value = 1_int64
        do i = 1, n
            if (value > cap / 2_int64) then
                value = cap
                return
            end if
            value = 2_int64 * value
            if (value >= cap) then
                value = cap
                return
            end if
        end do
    end function power_two_count_capped

    pure real(dp) function bootstrap_alpha(level, alternative) result(alpha)
        real(dp), intent(in) :: level !! confidence level in (0,1)
        character(len=*), intent(in) :: alternative !! two-sided, less, or greater

        if (alternative == 'two-sided') then
            alpha = 0.5_dp * (1.0_dp - level)
        else
            alpha = 1.0_dp - level
        end if
    end function bootstrap_alpha

    pure character(len=10) function lower_method(method) result(normalized)
        character(len=*), intent(in) :: method !! bootstrap method spelling

        integer :: code
        integer :: i
        integer :: n

        normalized = ' '
        n = min(len_trim(method), len(normalized))
        do i = 1, n
            code = iachar(method(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                normalized(i:i) = achar(code + iachar('a') - iachar('A'))
            else
                normalized(i:i) = method(i:i)
            end if
        end do
    end function lower_method

    pure logical function valid_level(level) result(valid)
        real(dp), intent(in) :: level !! candidate confidence level

        valid = .not. ieee_is_nan(level) .and. level > 0.0_dp .and. level < 1.0_dp
    end function valid_level

    pure logical function valid_probability(probability) result(valid)
        real(dp), intent(in) :: probability !! candidate probability

        valid = .not. ieee_is_nan(probability) .and. probability >= 0.0_dp .and. &
            probability <= 1.0_dp
    end function valid_probability

    pure logical function valid_alternative(alternative) result(valid)
        character(len=*), intent(in) :: alternative !! alternative-hypothesis name

        valid = alternative == 'two-sided' .or. alternative == 'less' .or. &
            alternative == 'greater'
    end function valid_alternative

    pure logical function valid_bootstrap_method(method) result(valid)
        character(len=*), intent(in) :: method !! normalized bootstrap method name

        valid = method == 'percentile' .or. method == 'basic' .or. method == 'bca'
    end function valid_bootstrap_method

    pure logical function valid_permutation_type(permutation_type) result(valid)
        character(len=*), intent(in) :: permutation_type !! permutation strategy name

        valid = permutation_type == 'independent' .or. permutation_type == 'samples' .or. &
            permutation_type == 'pairings'
    end function valid_permutation_type

    pure subroutine set_invalid_bootstrap(result_value)
        type(bootstrap_result), intent(out) :: result_value !! NaNs and zero draws on failure

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%standard_error = quiet_nan(0.0_dp)
        result_value%confidence_low = quiet_nan(0.0_dp)
        result_value%confidence_high = quiet_nan(0.0_dp)
        result_value%n_resamples = 0
    end subroutine set_invalid_bootstrap

    pure subroutine set_invalid_permutation(result_value)
        type(permutation_test_result), intent(out) :: result_value !! NaNs and zero draws on failure

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%n_resamples = 0
        result_value%exact = .false.
    end subroutine set_invalid_permutation

    pure subroutine set_invalid_monte_carlo(result_value)
        type(monte_carlo_test_result), intent(out) :: result_value !! NaNs and zero draws on failure

        result_value%statistic = quiet_nan(0.0_dp)
        result_value%pvalue = quiet_nan(0.0_dp)
        result_value%n_resamples = 0
    end subroutine set_invalid_monte_carlo

    pure subroutine set_invalid_power(result_value)
        type(power_result), intent(out) :: result_value !! NaNs and zero draws on failure

        result_value%power = quiet_nan(0.0_dp)
        result_value%significance = quiet_nan(0.0_dp)
        result_value%n_resamples = 0
    end subroutine set_invalid_power

end module scifort_resampling
