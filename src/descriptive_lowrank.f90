! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_lowrank
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, negative_infinity, quiet_nan
    implicit none
    private

    interface quantile
        module procedure quantile_scalar
        module procedure quantile_array
    end interface quantile

    interface gmean
        module procedure gmean_vector
        module procedure gmean_matrix
        module procedure gmean_matrix_weights2
        module procedure gmean_vector_nan_policy
        module procedure gmean_matrix_nan_policy
        module procedure gmean_matrix_weights2_nan_policy
        module procedure gmean_tensor
        module procedure gmean_tensor_weights1
        module procedure gmean_tensor_weights2
        module procedure gmean_tensor_weights3
        module procedure gmean_tensor_nan_policy
        module procedure gmean_tensor_weights1_nan_policy
        module procedure gmean_tensor_weights2_nan_policy
        module procedure gmean_tensor_weights3_nan_policy
        module procedure gmean_rank4
        module procedure gmean_rank4_nan_policy
    end interface gmean

    interface hmean
        module procedure hmean_vector
        module procedure hmean_matrix
        module procedure hmean_matrix_weights2
        module procedure hmean_vector_nan_policy
        module procedure hmean_matrix_nan_policy
        module procedure hmean_matrix_weights2_nan_policy
        module procedure hmean_tensor
        module procedure hmean_tensor_weights1
        module procedure hmean_tensor_weights2
        module procedure hmean_tensor_weights3
        module procedure hmean_tensor_nan_policy
        module procedure hmean_tensor_weights1_nan_policy
        module procedure hmean_tensor_weights2_nan_policy
        module procedure hmean_tensor_weights3_nan_policy
        module procedure hmean_rank4
        module procedure hmean_rank4_nan_policy
    end interface hmean

    interface pmean
        module procedure pmean_vector
        module procedure pmean_matrix
        module procedure pmean_matrix_weights2
        module procedure pmean_vector_nan_policy
        module procedure pmean_matrix_nan_policy
        module procedure pmean_matrix_weights2_nan_policy
        module procedure pmean_tensor
        module procedure pmean_tensor_weights1
        module procedure pmean_tensor_weights2
        module procedure pmean_tensor_weights3
        module procedure pmean_tensor_nan_policy
        module procedure pmean_tensor_weights1_nan_policy
        module procedure pmean_tensor_weights2_nan_policy
        module procedure pmean_tensor_weights3_nan_policy
        module procedure pmean_rank4
        module procedure pmean_rank4_nan_policy
    end interface pmean

    interface sem
        module procedure sem_vector
        module procedure sem_matrix
        module procedure sem_vector_nan_policy
        module procedure sem_matrix_nan_policy
        module procedure sem_tensor
        module procedure sem_tensor_nan_policy
        module procedure sem_rank4
        module procedure sem_rank4_nan_policy
    end interface sem

    interface variation
        module procedure variation_vector
        module procedure variation_matrix
        module procedure variation_vector_nan_policy
        module procedure variation_matrix_nan_policy
        module procedure variation_tensor
        module procedure variation_tensor_nan_policy
        module procedure variation_rank4
        module procedure variation_rank4_nan_policy
    end interface variation

    interface iqr
        module procedure iqr_vector
        module procedure iqr_matrix_all
        module procedure iqr_matrix_axis
        module procedure iqr_vector_nan_policy
        module procedure iqr_matrix_all_nan_policy
        module procedure iqr_matrix_axis_nan_policy
        module procedure iqr_tensor_all
        module procedure iqr_tensor_axis
        module procedure iqr_tensor_all_nan_policy
        module procedure iqr_tensor_axis_nan_policy
        module procedure iqr_rank4_all
        module procedure iqr_rank4_axis
        module procedure iqr_rank4_all_nan_policy
        module procedure iqr_rank4_axis_nan_policy
    end interface iqr

    interface median_abs_deviation
        module procedure median_abs_deviation_vector
        module procedure median_abs_deviation_matrix
        module procedure median_abs_deviation_vector_nan_policy
        module procedure median_abs_deviation_matrix_nan_policy
        module procedure median_abs_deviation_tensor
        module procedure median_abs_deviation_tensor_nan_policy
        module procedure median_abs_deviation_rank4
        module procedure median_abs_deviation_rank4_nan_policy
    end interface median_abs_deviation

    interface moment
        module procedure moment_vector
        module procedure moment_matrix
        module procedure moment_vector_nan_policy
        module procedure moment_matrix_nan_policy
        module procedure moment_tensor
        module procedure moment_tensor_nan_policy
        module procedure moment_rank4
        module procedure moment_rank4_nan_policy
    end interface moment

    interface skew
        module procedure skew_vector
        module procedure skew_matrix
        module procedure skew_vector_nan_policy
        module procedure skew_matrix_nan_policy
        module procedure skew_tensor
        module procedure skew_tensor_nan_policy
        module procedure skew_rank4
        module procedure skew_rank4_nan_policy
    end interface skew

    interface kurtosis
        module procedure kurtosis_vector
        module procedure kurtosis_matrix
        module procedure kurtosis_vector_nan_policy
        module procedure kurtosis_matrix_nan_policy
        module procedure kurtosis_tensor
        module procedure kurtosis_tensor_nan_policy
        module procedure kurtosis_rank4
        module procedure kurtosis_rank4_nan_policy
    end interface kurtosis

    interface gmean_keepdims
        module procedure gmean_vector_keepdims
        module procedure gmean_matrix_keepdims
        module procedure gmean_matrix_keepdims_weights2
        module procedure gmean_vector_keepdims_nan_policy
        module procedure gmean_matrix_keepdims_nan_policy
        module procedure gmean_matrix_keepdims_weights2_nan_policy
        module procedure gmean_tensor_keepdims
        module procedure gmean_tensor_keepdims_weights1
        module procedure gmean_tensor_keepdims_weights2
        module procedure gmean_tensor_keepdims_weights3
        module procedure gmean_tensor_keepdims_nan_policy
        module procedure gmean_tensor_keepdims_weights1_nan_policy
        module procedure gmean_tensor_keepdims_weights2_nan_policy
        module procedure gmean_tensor_keepdims_weights3_nan_policy
        module procedure gmean_rank4_keepdims
        module procedure gmean_rank4_keepdims_nan_policy
    end interface gmean_keepdims

    interface hmean_keepdims
        module procedure hmean_vector_keepdims
        module procedure hmean_matrix_keepdims
        module procedure hmean_matrix_keepdims_weights2
        module procedure hmean_vector_keepdims_nan_policy
        module procedure hmean_matrix_keepdims_nan_policy
        module procedure hmean_matrix_keepdims_weights2_nan_policy
        module procedure hmean_tensor_keepdims
        module procedure hmean_tensor_keepdims_weights1
        module procedure hmean_tensor_keepdims_weights2
        module procedure hmean_tensor_keepdims_weights3
        module procedure hmean_tensor_keepdims_nan_policy
        module procedure hmean_tensor_keepdims_weights1_nan_policy
        module procedure hmean_tensor_keepdims_weights2_nan_policy
        module procedure hmean_tensor_keepdims_weights3_nan_policy
        module procedure hmean_rank4_keepdims
        module procedure hmean_rank4_keepdims_nan_policy
    end interface hmean_keepdims

    interface pmean_keepdims
        module procedure pmean_vector_keepdims
        module procedure pmean_matrix_keepdims
        module procedure pmean_matrix_keepdims_weights2
        module procedure pmean_vector_keepdims_nan_policy
        module procedure pmean_matrix_keepdims_nan_policy
        module procedure pmean_matrix_keepdims_weights2_nan_policy
        module procedure pmean_tensor_keepdims
        module procedure pmean_tensor_keepdims_weights1
        module procedure pmean_tensor_keepdims_weights2
        module procedure pmean_tensor_keepdims_weights3
        module procedure pmean_tensor_keepdims_nan_policy
        module procedure pmean_tensor_keepdims_weights1_nan_policy
        module procedure pmean_tensor_keepdims_weights2_nan_policy
        module procedure pmean_tensor_keepdims_weights3_nan_policy
        module procedure pmean_rank4_keepdims
        module procedure pmean_rank4_keepdims_nan_policy
    end interface pmean_keepdims

    interface sem_keepdims
        module procedure sem_vector_keepdims
        module procedure sem_matrix_keepdims
        module procedure sem_vector_keepdims_nan_policy
        module procedure sem_matrix_keepdims_nan_policy
        module procedure sem_tensor_keepdims
        module procedure sem_tensor_keepdims_nan_policy
        module procedure sem_rank4_keepdims
        module procedure sem_rank4_keepdims_nan_policy
    end interface sem_keepdims

    interface variation_keepdims
        module procedure variation_vector_keepdims
        module procedure variation_matrix_keepdims
        module procedure variation_vector_keepdims_nan_policy
        module procedure variation_matrix_keepdims_nan_policy
        module procedure variation_tensor_keepdims
        module procedure variation_tensor_keepdims_nan_policy
        module procedure variation_rank4_keepdims
        module procedure variation_rank4_keepdims_nan_policy
    end interface variation_keepdims

    interface iqr_keepdims
        module procedure iqr_vector_keepdims
        module procedure iqr_matrix_keepdims
        module procedure iqr_vector_keepdims_nan_policy
        module procedure iqr_matrix_keepdims_nan_policy
        module procedure iqr_tensor_keepdims
        module procedure iqr_tensor_keepdims_nan_policy
        module procedure iqr_rank4_keepdims
        module procedure iqr_rank4_keepdims_nan_policy
    end interface iqr_keepdims

    interface median_abs_deviation_keepdims
        module procedure median_abs_deviation_vector_keepdims
        module procedure median_abs_deviation_matrix_keepdims
        module procedure median_abs_deviation_vector_keepdims_nan_policy
        module procedure median_abs_deviation_matrix_keepdims_nan_policy
        module procedure median_abs_deviation_tensor_keepdims
        module procedure median_abs_deviation_tensor_keepdims_nan_policy
        module procedure median_abs_deviation_rank4_keepdims
        module procedure median_abs_deviation_rank4_keepdims_nan_policy
    end interface median_abs_deviation_keepdims

    interface moment_keepdims
        module procedure moment_vector_keepdims
        module procedure moment_matrix_keepdims
        module procedure moment_vector_keepdims_nan_policy
        module procedure moment_matrix_keepdims_nan_policy
        module procedure moment_tensor_keepdims
        module procedure moment_tensor_keepdims_nan_policy
        module procedure moment_rank4_keepdims
        module procedure moment_rank4_keepdims_nan_policy
    end interface moment_keepdims

    interface skew_keepdims
        module procedure skew_vector_keepdims
        module procedure skew_matrix_keepdims
        module procedure skew_vector_keepdims_nan_policy
        module procedure skew_matrix_keepdims_nan_policy
        module procedure skew_tensor_keepdims
        module procedure skew_tensor_keepdims_nan_policy
        module procedure skew_rank4_keepdims
        module procedure skew_rank4_keepdims_nan_policy
    end interface skew_keepdims

    interface kurtosis_keepdims
        module procedure kurtosis_vector_keepdims
        module procedure kurtosis_matrix_keepdims
        module procedure kurtosis_vector_keepdims_nan_policy
        module procedure kurtosis_matrix_keepdims_nan_policy
        module procedure kurtosis_tensor_keepdims
        module procedure kurtosis_tensor_keepdims_nan_policy
        module procedure kurtosis_rank4_keepdims
        module procedure kurtosis_rank4_keepdims_nan_policy
    end interface kurtosis_keepdims

    public :: central_moment
    public :: covariance
    public :: gmean
    public :: hmean
    public :: iqr
    public :: gmean_keepdims
    public :: hmean_keepdims
    public :: iqr_keepdims
    public :: kurtosis
    public :: kurtosis_keepdims
    public :: mean
    public :: median
    public :: median_abs_deviation
    public :: median_abs_deviation_keepdims
    public :: moment
    public :: moment_keepdims
    public :: pearson_correlation
    public :: pmean
    public :: pmean_keepdims
    public :: quantile
    public :: sem
    public :: sem_keepdims
    public :: skew
    public :: skew_keepdims
    public :: standard_deviation
    public :: variation
    public :: variation_keepdims
    public :: variance

contains

    pure function mean(x) result(value)
        real(dp), intent(in) :: x(:) !! observations; NaN propagates
        real(dp) :: value

        integer :: i
        integer :: n_negative_inf
        integer :: n_positive_inf
        real(dp) :: correction
        real(dp) :: max_abs
        real(dp) :: scaled
        real(dp) :: sum_scaled
        real(dp) :: term

        if (size(x) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if

        n_negative_inf = 0
        n_positive_inf = 0
        max_abs = 0.0_dp
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            else if (.not. ieee_is_finite(x(i))) then
                if (x(i) > 0.0_dp) then
                    n_positive_inf = n_positive_inf + 1
                else
                    n_negative_inf = n_negative_inf + 1
                end if
            else
                max_abs = max(max_abs, abs(x(i)))
            end if
        end do

        if (n_positive_inf > 0 .and. n_negative_inf > 0) then
            value = quiet_nan(0.0_dp)
            return
        else if (n_positive_inf > 0) then
            value = positive_infinity(0.0_dp)
            return
        else if (n_negative_inf > 0) then
            value = negative_infinity(0.0_dp)
            return
        end if

        if (max_abs <= 0.0_dp) then
            value = 0.0_dp
            return
        end if

        sum_scaled = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(x)
            scaled = x(i) / max_abs
            term = scaled - correction
            value = sum_scaled + term
            correction = (value - sum_scaled) - term
            sum_scaled = value
        end do
        value = max_abs * (sum_scaled / real(size(x), dp))
    end function mean

    pure function standard_deviation(x, ddof) result(value)
        real(dp), intent(in) :: x(:) !! observations; all values must be finite
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp) :: value

        integer :: divisor_adjustment
        integer :: i
        real(dp) :: correction
        real(dp) :: data_scale
        real(dp) :: factor
        real(dp) :: mu
        real(dp) :: mu_scaled
        real(dp) :: scaled
        real(dp) :: sum_squares
        real(dp) :: term
        real(dp) :: updated

        divisor_adjustment = 1
        if (present(ddof)) divisor_adjustment = ddof

        if (size(x) == 0 .or. divisor_adjustment < 0 .or. &
                size(x) <= divisor_adjustment) then
            value = quiet_nan(0.0_dp)
            return
        end if

        data_scale = 0.0_dp
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            end if
            data_scale = max(data_scale, abs(x(i)))
        end do
        if (data_scale <= 0.0_dp) then
            value = 0.0_dp
            return
        end if

        mu = mean(x)
        mu_scaled = mu / data_scale
        sum_squares = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(x)
            scaled = x(i) / data_scale - mu_scaled
            term = scaled * scaled - correction
            updated = sum_squares + term
            correction = (updated - sum_squares) - term
            sum_squares = updated
        end do
        factor = sqrt(sum_squares / real(size(x) - divisor_adjustment, dp))
        value = safe_product(data_scale, factor)
    end function standard_deviation

    pure function variance(x, ddof) result(value)
        real(dp), intent(in) :: x(:) !! observations; all values must be finite
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp) :: value

        real(dp) :: sd

        if (present(ddof)) then
            sd = standard_deviation(x, ddof)
        else
            sd = standard_deviation(x)
        end if
        value = safe_product(sd, sd)
    end function variance

    pure function central_moment(x, order) result(value)
        real(dp), intent(in) :: x(:) !! observations; all values must be finite
        integer, intent(in) :: order !! nonnegative central-moment order
        real(dp) :: value

        integer :: i
        integer :: k
        real(dp) :: correction
        real(dp) :: data_scale
        real(dp) :: mu
        real(dp) :: mu_scaled
        real(dp) :: scaled
        real(dp) :: sum_power
        real(dp) :: term
        real(dp) :: updated

        if (size(x) == 0 .or. order < 0) then
            value = quiet_nan(0.0_dp)
            return
        else if (order == 0) then
            value = 1.0_dp
            return
        end if

        data_scale = 0.0_dp
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            end if
            data_scale = max(data_scale, abs(x(i)))
        end do
        if (order == 1 .or. data_scale <= 0.0_dp) then
            value = 0.0_dp
            return
        end if

        mu = mean(x)
        mu_scaled = mu / data_scale
        sum_power = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(x)
            scaled = x(i) / data_scale - mu_scaled
            term = scaled ** order - correction
            updated = sum_power + term
            correction = (updated - sum_power) - term
            sum_power = updated
        end do
        value = sum_power / real(size(x), dp)
        do k = 1, order
            value = safe_product(value, data_scale)
            if (.not. ieee_is_finite(value)) exit
        end do
    end function central_moment

    pure function gmean_vector(x, weights) result(value)
        real(dp), intent(in) :: x(:) !! nonnegative observations
        real(dp), intent(in), optional :: weights(:) !! finite observation weights
        real(dp) :: value

        integer :: i
        logical :: negative_infinite_sum
        logical :: positive_infinite_sum
        real(dp) :: correction
        real(dp) :: log_sum
        real(dp) :: term
        real(dp) :: updated
        real(dp) :: weight
        real(dp) :: weight_sum

        if (.not. valid_weight_shape(x, weights)) then
            value = quiet_nan(0.0_dp)
            return
        end if

        call weight_total(x, weights, weight_sum, value)
        if (ieee_is_nan(value)) return

        log_sum = 0.0_dp
        correction = 0.0_dp
        positive_infinite_sum = .false.
        negative_infinite_sum = .false.
        do i = 1, size(x)
            weight = observation_weight(weights, i)
            if (ieee_is_nan(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            else if (x(i) < 0.0_dp) then
                value = quiet_nan(0.0_dp)
                return
            else if (x(i) == 0.0_dp) then
                if (weight == 0.0_dp) then
                    value = quiet_nan(0.0_dp)
                    return
                else if (weight > 0.0_dp) then
                    negative_infinite_sum = .true.
                else
                    positive_infinite_sum = .true.
                end if
            else if (.not. ieee_is_finite(x(i))) then
                if (x(i) < 0.0_dp .or. weight == 0.0_dp) then
                    value = quiet_nan(0.0_dp)
                    return
                else if (weight > 0.0_dp) then
                    positive_infinite_sum = .true.
                else
                    negative_infinite_sum = .true.
                end if
            else
                term = weight * log(x(i))
                if (.not. ieee_is_finite(term)) then
                    if (term > 0.0_dp) then
                        positive_infinite_sum = .true.
                    else
                        negative_infinite_sum = .true.
                    end if
                else
                    term = term - correction
                    updated = log_sum + term
                    correction = (updated - log_sum) - term
                    log_sum = updated
                end if
            end if
        end do

        if (positive_infinite_sum .and. negative_infinite_sum) then
            value = quiet_nan(0.0_dp)
        else if (positive_infinite_sum) then
            if (weight_sum > 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else
                value = 0.0_dp
            end if
        else if (negative_infinite_sum) then
            if (weight_sum > 0.0_dp) then
                value = 0.0_dp
            else
                value = positive_infinity(0.0_dp)
            end if
        else
            value = exp(log_sum / weight_sum)
        end if
    end function gmean_vector

    pure function hmean_vector(x, weights) result(value)
        real(dp), intent(in) :: x(:) !! nonnegative observations
        real(dp), intent(in), optional :: weights(:) !! finite observation weights
        real(dp) :: value

        integer :: i
        logical :: negative_infinite_sum
        logical :: positive_infinite_sum
        real(dp) :: correction
        real(dp) :: denominator
        real(dp) :: term
        real(dp) :: updated
        real(dp) :: weight
        real(dp) :: weight_sum

        if (.not. valid_weight_shape(x, weights)) then
            value = quiet_nan(0.0_dp)
            return
        end if

        call weight_total(x, weights, weight_sum, value)
        if (ieee_is_nan(value)) return

        denominator = 0.0_dp
        correction = 0.0_dp
        positive_infinite_sum = .false.
        negative_infinite_sum = .false.
        do i = 1, size(x)
            weight = observation_weight(weights, i)
            if (ieee_is_nan(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            else if (x(i) < 0.0_dp) then
                value = quiet_nan(0.0_dp)
                return
            else if (x(i) == 0.0_dp) then
                if (weight == 0.0_dp) then
                    value = quiet_nan(0.0_dp)
                    return
                else if (weight > 0.0_dp) then
                    positive_infinite_sum = .true.
                else
                    negative_infinite_sum = .true.
                end if
            else if (.not. ieee_is_finite(x(i))) then
                cycle
            else
                term = weight / x(i)
                if (.not. ieee_is_finite(term)) then
                    if (term > 0.0_dp) then
                        positive_infinite_sum = .true.
                    else
                        negative_infinite_sum = .true.
                    end if
                else
                    term = term - correction
                    updated = denominator + term
                    correction = (updated - denominator) - term
                    denominator = updated
                end if
            end if
        end do

        if (positive_infinite_sum .and. negative_infinite_sum) then
            value = quiet_nan(0.0_dp)
        else if (positive_infinite_sum .or. negative_infinite_sum) then
            value = 0.0_dp
        else if (denominator == 0.0_dp) then
            if (weight_sum > 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else
                value = negative_infinity(0.0_dp)
            end if
        else
            value = weight_sum / denominator
        end if
    end function hmean_vector

    pure function pmean_vector(x, p, weights) result(value)
        real(dp), intent(in) :: x(:) !! nonnegative observations
        real(dp), intent(in) :: p !! finite power; zero selects the geometric mean
        real(dp), intent(in), optional :: weights(:) !! finite observation weights
        real(dp) :: value

        integer :: i
        logical :: negative_infinite_sum
        logical :: positive_infinite_sum
        real(dp) :: correction
        real(dp) :: mean_power
        real(dp) :: power_sum
        real(dp) :: powered
        real(dp) :: term
        real(dp) :: updated
        real(dp) :: weight
        real(dp) :: weight_sum

        if (.not. ieee_is_finite(p) .or. .not. valid_weight_shape(x, weights)) then
            value = quiet_nan(0.0_dp)
            return
        else if (p == 0.0_dp) then
            if (present(weights)) then
                value = gmean_vector(x, weights)
            else
                value = gmean_vector(x)
            end if
            return
        end if

        call weight_total(x, weights, weight_sum, value)
        if (ieee_is_nan(value)) return

        power_sum = 0.0_dp
        correction = 0.0_dp
        positive_infinite_sum = .false.
        negative_infinite_sum = .false.
        do i = 1, size(x)
            weight = observation_weight(weights, i)
            if (ieee_is_nan(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            else if (x(i) < 0.0_dp) then
                value = quiet_nan(0.0_dp)
                return
            else if (x(i) == 0.0_dp .and. p < 0.0_dp) then
                if (weight == 0.0_dp) then
                    value = quiet_nan(0.0_dp)
                    return
                else if (weight > 0.0_dp) then
                    positive_infinite_sum = .true.
                else
                    negative_infinite_sum = .true.
                end if
                cycle
            else if (.not. ieee_is_finite(x(i))) then
                if (p < 0.0_dp) cycle
                if (weight == 0.0_dp) then
                    value = quiet_nan(0.0_dp)
                    return
                else if (weight > 0.0_dp) then
                    positive_infinite_sum = .true.
                else
                    negative_infinite_sum = .true.
                end if
                cycle
            end if

            powered = nonnegative_power(x(i), p)
            if (.not. ieee_is_finite(powered)) then
                if (weight == 0.0_dp) then
                    value = quiet_nan(0.0_dp)
                    return
                else if (weight > 0.0_dp) then
                    positive_infinite_sum = .true.
                else
                    negative_infinite_sum = .true.
                end if
                cycle
            end if
            term = safe_product(weight, powered)
            if (.not. ieee_is_finite(term)) then
                if (term > 0.0_dp) then
                    positive_infinite_sum = .true.
                else
                    negative_infinite_sum = .true.
                end if
            else
                term = term - correction
                updated = power_sum + term
                correction = (updated - power_sum) - term
                power_sum = updated
            end if
        end do

        if (positive_infinite_sum .and. negative_infinite_sum) then
            value = quiet_nan(0.0_dp)
            return
        else if (positive_infinite_sum) then
            if (weight_sum > 0.0_dp) then
                mean_power = positive_infinity(0.0_dp)
            else
                mean_power = negative_infinity(0.0_dp)
            end if
        else if (negative_infinite_sum) then
            if (weight_sum > 0.0_dp) then
                mean_power = negative_infinity(0.0_dp)
            else
                mean_power = positive_infinity(0.0_dp)
            end if
        else
            mean_power = power_sum / weight_sum
        end if

        if (p == 1.0_dp) then
            value = mean_power
        else if (p == -1.0_dp) then
            if (mean_power == 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else
                value = 1.0_dp / mean_power
            end if
        else if (mean_power < 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else if (mean_power == 0.0_dp .and. p < 0.0_dp) then
            value = positive_infinity(0.0_dp)
        else
            value = mean_power ** (1.0_dp / p)
        end if
    end function pmean_vector

    pure function sem_vector(x, ddof) result(value)
        real(dp), intent(in) :: x(:) !! observations; all values must be finite
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp) :: value

        integer :: divisor_adjustment
        real(dp) :: sd

        divisor_adjustment = 1
        if (present(ddof)) divisor_adjustment = ddof
        sd = standard_deviation(x, divisor_adjustment)
        if (ieee_is_nan(sd)) then
            value = sd
        else
            value = sd / sqrt(real(size(x), dp))
        end if
    end function sem_vector

    pure function variation_vector(x, ddof) result(value)
        real(dp), intent(in) :: x(:) !! observations; all values must be finite
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        real(dp) :: value

        integer :: divisor_adjustment
        real(dp) :: mu
        real(dp) :: sd

        divisor_adjustment = 0
        if (present(ddof)) divisor_adjustment = ddof
        sd = standard_deviation(x, divisor_adjustment)
        mu = mean(x)
        if (ieee_is_nan(sd) .or. ieee_is_nan(mu)) then
            value = quiet_nan(0.0_dp)
        else if (mu == 0.0_dp) then
            if (sd == 0.0_dp) then
                value = quiet_nan(0.0_dp)
            else
                value = positive_infinity(0.0_dp)
            end if
        else
            value = sd / mu
        end if
    end function variation_vector

    pure function iqr_vector(x, rng, scale) result(value)
        real(dp), intent(in) :: x(:) !! observations; NaN propagates
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to the percentile range, default 1
        real(dp) :: value

        real(dp) :: high
        real(dp) :: low
        real(dp) :: scale_value

        low = 25.0_dp
        high = 75.0_dp
        if (present(rng)) then
            low = min(rng(1), rng(2))
            high = max(rng(1), rng(2))
        end if
        scale_value = 1.0_dp
        if (present(scale)) scale_value = scale

        if (low < 0.0_dp .or. high > 100.0_dp .or. ieee_is_nan(scale_value)) then
            value = quiet_nan(0.0_dp)
            return
        end if

        value = quantile(x, high / 100.0_dp) - quantile(x, low / 100.0_dp)
        value = divide_by_scale(value, scale_value)
    end function iqr_vector

    pure function median_abs_deviation_vector(x, scale) result(value)
        real(dp), intent(in) :: x(:) !! observations; NaN propagates
        real(dp), intent(in), optional :: scale !! divisor applied to the MAD, default 1
        real(dp) :: value

        real(dp), allocatable :: deviations(:)
        real(dp) :: center
        real(dp) :: scale_value

        scale_value = 1.0_dp
        if (present(scale)) scale_value = scale
        if (size(x) == 0 .or. ieee_is_nan(scale_value)) then
            value = quiet_nan(0.0_dp)
            return
        end if

        center = median(x)
        if (ieee_is_nan(center)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        deviations = abs(x - center)
        value = divide_by_scale(median(deviations), scale_value)
    end function median_abs_deviation_vector

    pure function moment_vector(x, order, center) result(value)
        real(dp), intent(in) :: x(:) !! observations; all values must be finite
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        real(dp), intent(in), optional :: center !! finite center; default sample mean
        real(dp) :: value

        integer :: moment_order
        real(dp) :: center_value

        moment_order = 1
        if (present(order)) moment_order = order
        if (size(x) == 0 .or. moment_order < 0) then
            value = quiet_nan(0.0_dp)
            return
        end if

        if (present(center)) then
            center_value = center
            value = centered_moment(x, moment_order, center_value)
        else
            value = central_moment(x, moment_order)
        end if
    end function moment_vector

    pure function skew_vector(x, bias) result(value)
        real(dp), intent(in) :: x(:) !! observations; all values must be finite
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp) :: value

        logical :: use_bias
        real(dp) :: m2
        real(dp) :: m3
        real(dp) :: n

        use_bias = .true.
        if (present(bias)) use_bias = bias
        if (size(x) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        m2 = central_moment(x, 2)
        m3 = central_moment(x, 3)
        if (.not. ieee_is_finite(m2)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (.not. ieee_is_finite(m3)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (m2 <= 0.0_dp) then
            value = quiet_nan(0.0_dp)
            return
        end if

        value = m3 / (m2 * sqrt(m2))
        if (.not. use_bias .and. size(x) > 2) then
            n = real(size(x), dp)
            value = sqrt(n * (n - 1.0_dp)) / (n - 2.0_dp) * value
        end if
    end function skew_vector

    pure function kurtosis_vector(x, fisher, bias) result(value)
        real(dp), intent(in) :: x(:) !! observations; all values must be finite
        logical, intent(in), optional :: fisher !! subtract 3 for Fisher kurtosis, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp) :: value

        logical :: fisher_value
        logical :: use_bias
        real(dp) :: m2
        real(dp) :: m4
        real(dp) :: n
        real(dp) :: pearson

        fisher_value = .true.
        if (present(fisher)) fisher_value = fisher
        use_bias = .true.
        if (present(bias)) use_bias = bias
        if (size(x) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        m2 = central_moment(x, 2)
        m4 = central_moment(x, 4)
        if (.not. ieee_is_finite(m2)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (.not. ieee_is_finite(m4)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (m2 <= 0.0_dp) then
            value = quiet_nan(0.0_dp)
            return
        end if

        pearson = m4 / (m2 * m2)
        if (.not. use_bias .and. size(x) > 3) then
            n = real(size(x), dp)
            pearson = ((n * n - 1.0_dp) * pearson - 3.0_dp * (n - 1.0_dp) ** 2) / &
                ((n - 2.0_dp) * (n - 3.0_dp)) + 3.0_dp
        end if
        if (fisher_value) then
            value = pearson - 3.0_dp
        else
            value = pearson
        end if
    end function kurtosis_vector


    pure function gmean_matrix(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        logical :: valid
        real(dp), allocatable :: work_weights(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (.not. present(weights)) then
            if (axis_value == 1) then
                do j = 1, size(x, 2)
                    values(j) = gmean_vector(x(:, j))
                end do
            else
                do j = 1, size(x, 1)
                    values(j) = gmean_vector(x(j, :))
                end do
            end if
            return
        end if

        valid = matrix_weight_vector_valid(x, weights)
        if (.not. valid) then
            values = quiet_nan(0.0_dp)
            return
        end if
        do j = 1, size(values)
            call matrix_weight_vector_slice(x, axis_value, j, weights, work_weights)
            if (axis_value == 1) then
                values(j) = gmean_vector(x(:, j), work_weights)
            else
                values(j) = gmean_vector(x(j, :), work_weights)
            end if
        end do
    end function gmean_matrix

    pure function gmean_matrix_weights2(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (.not. matrix_weight_matrix_valid(x, weights)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        do j = 1, size(values)
            call matrix_weight_matrix_slice(x, axis_value, j, weights, work_weights)
            if (axis_value == 1) then
                values(j) = gmean_vector(x(:, j), work_weights)
            else
                values(j) = gmean_vector(x(j, :), work_weights)
            end if
        end do
    end function gmean_matrix_weights2

    pure function hmean_matrix(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (present(weights)) then
            if (.not. matrix_weight_vector_valid(x, weights)) then
                values = quiet_nan(0.0_dp)
                return
            end if
        end if
        do j = 1, size(values)
            if (present(weights)) then
                call matrix_weight_vector_slice(x, axis_value, j, weights, work_weights)
                if (axis_value == 1) then
                    values(j) = hmean_vector(x(:, j), work_weights)
                else
                    values(j) = hmean_vector(x(j, :), work_weights)
                end if
            else if (axis_value == 1) then
                values(j) = hmean_vector(x(:, j))
            else
                values(j) = hmean_vector(x(j, :))
            end if
        end do
    end function hmean_matrix

    pure function hmean_matrix_weights2(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (.not. matrix_weight_matrix_valid(x, weights)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        do j = 1, size(values)
            call matrix_weight_matrix_slice(x, axis_value, j, weights, work_weights)
            if (axis_value == 1) then
                values(j) = hmean_vector(x(:, j), work_weights)
            else
                values(j) = hmean_vector(x(j, :), work_weights)
            end if
        end do
    end function hmean_matrix_weights2

    pure function pmean_matrix(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power; zero selects the geometric mean
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (present(weights)) then
            if (.not. matrix_weight_vector_valid(x, weights)) then
                values = quiet_nan(0.0_dp)
                return
            end if
        end if
        do j = 1, size(values)
            if (present(weights)) then
                call matrix_weight_vector_slice(x, axis_value, j, weights, work_weights)
                if (axis_value == 1) then
                    values(j) = pmean_vector(x(:, j), p, work_weights)
                else
                    values(j) = pmean_vector(x(j, :), p, work_weights)
                end if
            else if (axis_value == 1) then
                values(j) = pmean_vector(x(:, j), p)
            else
                values(j) = pmean_vector(x(j, :), p)
            end if
        end do
    end function pmean_matrix

    pure function pmean_matrix_weights2(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power; zero selects the geometric mean
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (.not. matrix_weight_matrix_valid(x, weights)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        do j = 1, size(values)
            call matrix_weight_matrix_slice(x, axis_value, j, weights, work_weights)
            if (axis_value == 1) then
                values(j) = pmean_vector(x(:, j), p, work_weights)
            else
                values(j) = pmean_vector(x(j, :), p, work_weights)
            end if
        end do
    end function pmean_matrix_weights2

    pure function sem_matrix(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = sem_vector(x(:, j), ddof)
            else
                values(j) = sem_vector(x(j, :), ddof)
            end if
        end do
    end function sem_matrix

    pure function variation_matrix(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = variation_vector(x(:, j), ddof)
            else
                values(j) = variation_vector(x(j, :), ddof)
            end if
        end do
    end function variation_matrix

    pure function iqr_matrix_all(x, rng, scale) result(value)
        real(dp), intent(in) :: x(:, :) !! observations reduced over both axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to the range, default 1
        real(dp) :: value
        real(dp), allocatable :: flattened(:)

        flattened = reshape(x, [size(x)])
        value = iqr_vector(flattened, rng, scale)
    end function iqr_matrix_all

    pure function iqr_matrix_axis(x, axis, rng, scale) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 or 2
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to the range, default 1
        real(dp), allocatable :: values(:)

        integer :: j

        if (axis /= 1 .and. axis /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis == 1)))
        do j = 1, size(values)
            if (axis == 1) then
                values(j) = iqr_vector(x(:, j), rng, scale)
            else
                values(j) = iqr_vector(x(j, :), rng, scale)
            end if
        end do
    end function iqr_matrix_axis

    pure function median_abs_deviation_matrix(x, axis, scale) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = median_abs_deviation_vector(x(:, j), scale)
            else
                values(j) = median_abs_deviation_vector(x(j, :), scale)
            end if
        end do
    end function median_abs_deviation_matrix

    pure function moment_matrix(x, order, axis, center) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = moment_vector(x(:, j), order, center)
            else
                values(j) = moment_vector(x(j, :), order, center)
            end if
        end do
    end function moment_matrix

    pure function skew_matrix(x, axis, bias) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = skew_vector(x(:, j), bias)
            else
                values(j) = skew_vector(x(j, :), bias)
            end if
        end do
    end function skew_matrix

    pure function kurtosis_matrix(x, axis, fisher, bias) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = kurtosis_vector(x(:, j), fisher, bias)
            else
                values(j) = kurtosis_vector(x(j, :), fisher, bias)
            end if
        end do
    end function kurtosis_matrix

    pure function gmean_vector_keepdims(x, weights) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in), optional :: weights(:) !! observation weights
        real(dp) :: values(1)

        values(1) = gmean_vector(x, weights)
    end function gmean_vector_keepdims

    pure function hmean_vector_keepdims(x, weights) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in), optional :: weights(:) !! observation weights
        real(dp) :: values(1)

        values(1) = hmean_vector(x, weights)
    end function hmean_vector_keepdims

    pure function pmean_vector_keepdims(x, p, weights) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in) :: p !! finite power
        real(dp), intent(in), optional :: weights(:) !! observation weights
        real(dp) :: values(1)

        values(1) = pmean_vector(x, p, weights)
    end function pmean_vector_keepdims

    pure function sem_vector_keepdims(x, ddof) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp) :: values(1)

        values(1) = sem_vector(x, ddof)
    end function sem_vector_keepdims

    pure function variation_vector_keepdims(x, ddof) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        real(dp) :: values(1)

        values(1) = variation_vector(x, ddof)
    end function variation_vector_keepdims

    pure function iqr_vector_keepdims(x, rng, scale) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        real(dp) :: values(1)

        values(1) = iqr_vector(x, rng, scale)
    end function iqr_vector_keepdims

    pure function median_abs_deviation_vector_keepdims(x, scale) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        real(dp) :: values(1)

        values(1) = median_abs_deviation_vector(x, scale)
    end function median_abs_deviation_vector_keepdims

    pure function moment_vector_keepdims(x, order, center) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        real(dp), intent(in), optional :: center !! finite center; default sample mean
        real(dp) :: values(1)

        values(1) = moment_vector(x, order, center)
    end function moment_vector_keepdims

    pure function skew_vector_keepdims(x, bias) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp) :: values(1)

        values(1) = skew_vector(x, bias)
    end function skew_vector_keepdims

    pure function kurtosis_vector_keepdims(x, fisher, bias) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp) :: values(1)

        values(1) = kurtosis_vector(x, fisher, bias)
    end function kurtosis_vector_keepdims

    pure function gmean_matrix_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = gmean_matrix(x, axis, weights)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function gmean_matrix_keepdims

    pure function gmean_matrix_keepdims_weights2(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = gmean_matrix_weights2(x, axis, weights)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function gmean_matrix_keepdims_weights2

    pure function hmean_matrix_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = hmean_matrix(x, axis, weights)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function hmean_matrix_keepdims

    pure function hmean_matrix_keepdims_weights2(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = hmean_matrix_weights2(x, axis, weights)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function hmean_matrix_keepdims_weights2

    pure function pmean_matrix_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = pmean_matrix(x, p, axis, weights)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function pmean_matrix_keepdims

    pure function pmean_matrix_keepdims_weights2(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = pmean_matrix_weights2(x, p, axis, weights)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function pmean_matrix_keepdims_weights2

    pure function sem_matrix_keepdims(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = sem_matrix(x, axis, ddof)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function sem_matrix_keepdims

    pure function variation_matrix_keepdims(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = variation_matrix(x, axis, ddof)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function variation_matrix_keepdims

    pure function iqr_matrix_keepdims(x, axis, rng, scale) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimensions retained
        integer, intent(in), optional :: axis !! Fortran axis; absent reduces all dimensions
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: flattened(:)
        real(dp), allocatable :: reduced(:)

        if (.not. present(axis)) then
            allocate(values(1, 1))
            flattened = reshape(x, [size(x)])
            values(1, 1) = iqr_vector(flattened, rng, scale)
            return
        end if
        reduced = iqr_matrix_axis(x, axis, rng, scale)
        values = retain_matrix_axis(x, axis, reduced)
    end function iqr_matrix_keepdims

    pure function median_abs_deviation_matrix_keepdims(x, axis, scale) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = median_abs_deviation_matrix(x, axis, scale)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function median_abs_deviation_matrix_keepdims

    pure function moment_matrix_keepdims(x, order, axis, center) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = moment_matrix(x, order, axis, center)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function moment_matrix_keepdims

    pure function skew_matrix_keepdims(x, axis, bias) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = skew_matrix(x, axis, bias)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function skew_matrix_keepdims

    pure function kurtosis_matrix_keepdims(x, axis, fisher, bias) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = kurtosis_matrix(x, axis, fisher, bias)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function kurtosis_matrix_keepdims

    pure function gmean_vector_nan_policy(x, weights, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! nonnegative observations
        real(dp), intent(in), optional :: weights(:) !! observation weights; NaN pairs may be omitted
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: mask(size(x))
        real(dp), allocatable :: clean_weights(:)
        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (present(weights)) then
            if (size(weights) /= size(x)) then
                value = quiet_nan(0.0_dp)
                return
            end if
        end if
        if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x, weights)
            value = gmean_vector(x, weights)
            return
        end if
        if (.not. nan_policy_omits(nan_policy)) then
            value = gmean_vector(x, weights)
            return
        end if
        if (present(weights)) then
            mask = .not. ieee_is_nan(x) .and. .not. ieee_is_nan(weights)
            clean_x = pack(x, mask)
            clean_weights = pack(weights, mask)
            value = gmean_vector(clean_x, clean_weights)
        else
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = gmean_vector(clean_x)
        end if
    end function gmean_vector_nan_policy

    pure function hmean_vector_nan_policy(x, weights, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! nonnegative observations
        real(dp), intent(in), optional :: weights(:) !! observation weights; NaN pairs may be omitted
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: mask(size(x))
        real(dp), allocatable :: clean_weights(:)
        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (present(weights)) then
            if (size(weights) /= size(x)) then
                value = quiet_nan(0.0_dp)
                return
            end if
        end if
        if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x, weights)
            value = hmean_vector(x, weights)
            return
        end if
        if (.not. nan_policy_omits(nan_policy)) then
            value = hmean_vector(x, weights)
            return
        end if
        if (present(weights)) then
            mask = .not. ieee_is_nan(x) .and. .not. ieee_is_nan(weights)
            clean_x = pack(x, mask)
            clean_weights = pack(weights, mask)
            value = hmean_vector(clean_x, clean_weights)
        else
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = hmean_vector(clean_x)
        end if
    end function hmean_vector_nan_policy

    pure function pmean_vector_nan_policy(x, p, weights, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! nonnegative observations
        real(dp), intent(in) :: p !! finite power; zero selects the geometric mean
        real(dp), intent(in), optional :: weights(:) !! observation weights; NaN pairs may be omitted
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: mask(size(x))
        real(dp), allocatable :: clean_weights(:)
        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (present(weights)) then
            if (size(weights) /= size(x)) then
                value = quiet_nan(0.0_dp)
                return
            end if
        end if
        if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x, weights)
            value = pmean_vector(x, p, weights)
            return
        end if
        if (.not. nan_policy_omits(nan_policy)) then
            value = pmean_vector(x, p, weights)
            return
        end if
        if (present(weights)) then
            mask = .not. ieee_is_nan(x) .and. .not. ieee_is_nan(weights)
            clean_x = pack(x, mask)
            clean_weights = pack(weights, mask)
            value = pmean_vector(clean_x, p, clean_weights)
        else
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = pmean_vector(clean_x, p)
        end if
    end function pmean_vector_nan_policy

    pure function sem_vector_nan_policy(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
        else if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x)
            value = sem_vector(x, ddof)
        else if (nan_policy_omits(nan_policy)) then
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = sem_vector(clean_x, ddof)
        else
            value = sem_vector(x, ddof)
        end if
    end function sem_vector_nan_policy

    pure function variation_vector_nan_policy(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
        else if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x)
            value = variation_vector(x, ddof)
        else if (nan_policy_omits(nan_policy)) then
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = variation_vector(clean_x, ddof)
        else
            value = variation_vector(x, ddof)
        end if
    end function variation_vector_nan_policy

    pure function iqr_vector_nan_policy(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to the percentile range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
        else if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x)
            value = iqr_vector(x, rng, scale)
        else if (nan_policy_omits(nan_policy)) then
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = iqr_vector(clean_x, rng, scale)
        else
            value = iqr_vector(x, rng, scale)
        end if
    end function iqr_vector_nan_policy

    pure function median_abs_deviation_vector_nan_policy(x, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
        else if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x)
            value = median_abs_deviation_vector(x, scale)
        else if (nan_policy_omits(nan_policy)) then
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = median_abs_deviation_vector(clean_x, scale)
        else
            value = median_abs_deviation_vector(x, scale)
        end if
    end function median_abs_deviation_vector_nan_policy

    pure function moment_vector_nan_policy(x, order, center, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        real(dp), intent(in), optional :: center !! finite center; default sample mean
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
        else if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x)
            value = moment_vector(x, order, center)
        else if (nan_policy_omits(nan_policy)) then
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = moment_vector(clean_x, order, center)
        else
            value = moment_vector(x, order, center)
        end if
    end function moment_vector_nan_policy

    pure function skew_vector_nan_policy(x, bias, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
        else if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x)
            value = skew_vector(x, bias)
        else if (nan_policy_omits(nan_policy)) then
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = skew_vector(clean_x, bias)
        else
            value = skew_vector(x, bias)
        end if
    end function skew_vector_nan_policy

    pure function kurtosis_vector_nan_policy(x, fisher, bias, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        logical, intent(in), optional :: fisher !! subtract 3 for Fisher kurtosis, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        real(dp), allocatable :: clean_x(:)

        if (.not. valid_nan_policy(nan_policy)) then
            value = quiet_nan(0.0_dp)
        else if (nan_policy_raises(nan_policy)) then
            call raise_if_nan_vector(x)
            value = kurtosis_vector(x, fisher, bias)
        else if (nan_policy_omits(nan_policy)) then
            clean_x = pack(x, .not. ieee_is_nan(x))
            value = kurtosis_vector(clean_x, fisher, bias)
        else
            value = kurtosis_vector(x, fisher, bias)
        end if
    end function kurtosis_vector_nan_policy

    pure function gmean_matrix_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        if (.not. valid_nan_policy(nan_policy)) then
            allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
            values = quiet_nan(0.0_dp)
            return
        end if
        if (.not. nan_policy_omits(nan_policy) .and. .not. nan_policy_raises(nan_policy)) then
            values = gmean_matrix(x, axis, weights)
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (present(weights)) then
            if (.not. matrix_weight_vector_valid(x, weights)) then
                values = quiet_nan(0.0_dp)
                return
            end if
        end if
        do j = 1, size(values)
            if (present(weights)) then
                call matrix_weight_vector_slice(x, axis_value, j, weights, work_weights)
                if (axis_value == 1) then
                    values(j) = gmean_vector_nan_policy(x(:, j), work_weights, nan_policy)
                else
                    values(j) = gmean_vector_nan_policy(x(j, :), work_weights, nan_policy)
                end if
            else if (axis_value == 1) then
                values(j) = gmean_vector_nan_policy(x(:, j), nan_policy=nan_policy)
            else
                values(j) = gmean_vector_nan_policy(x(j, :), nan_policy=nan_policy)
            end if
        end do
    end function gmean_matrix_nan_policy

    pure function gmean_matrix_weights2_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (.not. valid_nan_policy(nan_policy)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        if (.not. nan_policy_omits(nan_policy) .and. .not. nan_policy_raises(nan_policy)) then
            values = gmean_matrix_weights2(x, axis, weights)
            return
        end if
        if (.not. matrix_weight_matrix_valid(x, weights)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        do j = 1, size(values)
            call matrix_weight_matrix_slice(x, axis_value, j, weights, work_weights)
            if (axis_value == 1) then
                values(j) = gmean_vector_nan_policy(x(:, j), work_weights, nan_policy)
            else
                values(j) = gmean_vector_nan_policy(x(j, :), work_weights, nan_policy)
            end if
        end do
    end function gmean_matrix_weights2_nan_policy

    pure function hmean_matrix_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (.not. valid_nan_policy(nan_policy)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        if (.not. nan_policy_omits(nan_policy) .and. .not. nan_policy_raises(nan_policy)) then
            values = hmean_matrix(x, axis, weights)
            return
        end if
        if (present(weights)) then
            if (.not. matrix_weight_vector_valid(x, weights)) then
                values = quiet_nan(0.0_dp)
                return
            end if
        end if
        do j = 1, size(values)
            if (present(weights)) then
                call matrix_weight_vector_slice(x, axis_value, j, weights, work_weights)
                if (axis_value == 1) then
                    values(j) = hmean_vector_nan_policy(x(:, j), work_weights, nan_policy)
                else
                    values(j) = hmean_vector_nan_policy(x(j, :), work_weights, nan_policy)
                end if
            else if (axis_value == 1) then
                values(j) = hmean_vector_nan_policy(x(:, j), nan_policy=nan_policy)
            else
                values(j) = hmean_vector_nan_policy(x(j, :), nan_policy=nan_policy)
            end if
        end do
    end function hmean_matrix_nan_policy

    pure function hmean_matrix_weights2_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (.not. valid_nan_policy(nan_policy)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        if (.not. nan_policy_omits(nan_policy) .and. .not. nan_policy_raises(nan_policy)) then
            values = hmean_matrix_weights2(x, axis, weights)
            return
        end if
        if (.not. matrix_weight_matrix_valid(x, weights)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        do j = 1, size(values)
            call matrix_weight_matrix_slice(x, axis_value, j, weights, work_weights)
            if (axis_value == 1) then
                values(j) = hmean_vector_nan_policy(x(:, j), work_weights, nan_policy)
            else
                values(j) = hmean_vector_nan_policy(x(j, :), work_weights, nan_policy)
            end if
        end do
    end function hmean_matrix_weights2_nan_policy

    pure function pmean_matrix_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power; zero selects the geometric mean
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (.not. valid_nan_policy(nan_policy)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        if (.not. nan_policy_omits(nan_policy) .and. .not. nan_policy_raises(nan_policy)) then
            values = pmean_matrix(x, p, axis, weights)
            return
        end if
        if (present(weights)) then
            if (.not. matrix_weight_vector_valid(x, weights)) then
                values = quiet_nan(0.0_dp)
                return
            end if
        end if
        do j = 1, size(values)
            if (present(weights)) then
                call matrix_weight_vector_slice(x, axis_value, j, weights, work_weights)
                if (axis_value == 1) then
                    values(j) = pmean_vector_nan_policy(x(:, j), p, work_weights, nan_policy)
                else
                    values(j) = pmean_vector_nan_policy(x(j, :), p, work_weights, nan_policy)
                end if
            else if (axis_value == 1) then
                values(j) = pmean_vector_nan_policy(x(:, j), p, nan_policy=nan_policy)
            else
                values(j) = pmean_vector_nan_policy(x(j, :), p, nan_policy=nan_policy)
            end if
        end do
    end function pmean_matrix_nan_policy

    pure function pmean_matrix_weights2_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power; zero selects the geometric mean
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j
        real(dp), allocatable :: work_weights(:)

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        if (.not. valid_nan_policy(nan_policy)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        if (.not. nan_policy_omits(nan_policy) .and. .not. nan_policy_raises(nan_policy)) then
            values = pmean_matrix_weights2(x, p, axis, weights)
            return
        end if
        if (.not. matrix_weight_matrix_valid(x, weights)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        do j = 1, size(values)
            call matrix_weight_matrix_slice(x, axis_value, j, weights, work_weights)
            if (axis_value == 1) then
                values(j) = pmean_vector_nan_policy(x(:, j), p, work_weights, nan_policy)
            else
                values(j) = pmean_vector_nan_policy(x(j, :), p, work_weights, nan_policy)
            end if
        end do
    end function pmean_matrix_weights2_nan_policy

    pure function sem_matrix_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = sem_vector_nan_policy(x(:, j), ddof, nan_policy)
            else
                values(j) = sem_vector_nan_policy(x(j, :), ddof, nan_policy)
            end if
        end do
    end function sem_matrix_nan_policy

    pure function variation_matrix_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = variation_vector_nan_policy(x(:, j), ddof, nan_policy)
            else
                values(j) = variation_vector_nan_policy(x(j, :), ddof, nan_policy)
            end if
        end do
    end function variation_matrix_nan_policy

    pure function iqr_matrix_all_nan_policy(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations reduced over both axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to the range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp), allocatable :: flattened(:)

        flattened = reshape(x, [size(x)])
        value = iqr_vector_nan_policy(flattened, rng, scale, nan_policy)
    end function iqr_matrix_all_nan_policy

    pure function iqr_matrix_axis_nan_policy(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 or 2
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to the range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: j

        if (axis /= 1 .and. axis /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis == 1)))
        do j = 1, size(values)
            if (axis == 1) then
                values(j) = iqr_vector_nan_policy(x(:, j), rng, scale, nan_policy)
            else
                values(j) = iqr_vector_nan_policy(x(j, :), rng, scale, nan_policy)
            end if
        end do
    end function iqr_matrix_axis_nan_policy

    pure function median_abs_deviation_matrix_nan_policy(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = median_abs_deviation_vector_nan_policy(x(:, j), scale, nan_policy)
            else
                values(j) = median_abs_deviation_vector_nan_policy(x(j, :), scale, nan_policy)
            end if
        end do
    end function median_abs_deviation_matrix_nan_policy

    pure function moment_matrix_nan_policy(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = moment_vector_nan_policy(x(:, j), order, center, nan_policy)
            else
                values(j) = moment_vector_nan_policy(x(j, :), order, center, nan_policy)
            end if
        end do
    end function moment_matrix_nan_policy

    pure function skew_matrix_nan_policy(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = skew_vector_nan_policy(x(:, j), bias, nan_policy)
            else
                values(j) = skew_vector_nan_policy(x(j, :), bias, nan_policy)
            end if
        end do
    end function skew_matrix_nan_policy

    pure function kurtosis_matrix_nan_policy(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: j

        axis_value = default_matrix_axis(axis)
        if (axis_value /= 1 .and. axis_value /= 2) then
            allocate(values(0))
            return
        end if
        allocate(values(merge(size(x, 2), size(x, 1), axis_value == 1)))
        do j = 1, size(values)
            if (axis_value == 1) then
                values(j) = kurtosis_vector_nan_policy(x(:, j), fisher, bias, nan_policy)
            else
                values(j) = kurtosis_vector_nan_policy(x(j, :), fisher, bias, nan_policy)
            end if
        end do
    end function kurtosis_matrix_nan_policy

    pure function gmean_vector_keepdims_nan_policy(x, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in), optional :: weights(:) !! observation weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = gmean_vector_nan_policy(x, weights, nan_policy)
    end function gmean_vector_keepdims_nan_policy

    pure function hmean_vector_keepdims_nan_policy(x, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in), optional :: weights(:) !! observation weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = hmean_vector_nan_policy(x, weights, nan_policy)
    end function hmean_vector_keepdims_nan_policy

    pure function pmean_vector_keepdims_nan_policy(x, p, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in) :: p !! finite power
        real(dp), intent(in), optional :: weights(:) !! observation weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = pmean_vector_nan_policy(x, p, weights, nan_policy)
    end function pmean_vector_keepdims_nan_policy

    pure function sem_vector_keepdims_nan_policy(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = sem_vector_nan_policy(x, ddof, nan_policy)
    end function sem_vector_keepdims_nan_policy

    pure function variation_vector_keepdims_nan_policy(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = variation_vector_nan_policy(x, ddof, nan_policy)
    end function variation_vector_keepdims_nan_policy

    pure function iqr_vector_keepdims_nan_policy(x, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = iqr_vector_nan_policy(x, rng, scale, nan_policy)
    end function iqr_vector_keepdims_nan_policy

    pure function median_abs_deviation_vector_keepdims_nan_policy(x, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = median_abs_deviation_vector_nan_policy(x, scale, nan_policy)
    end function median_abs_deviation_vector_keepdims_nan_policy

    pure function moment_vector_keepdims_nan_policy(x, order, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        real(dp), intent(in), optional :: center !! finite center; default sample mean
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = moment_vector_nan_policy(x, order, center, nan_policy)
    end function moment_vector_keepdims_nan_policy

    pure function skew_vector_keepdims_nan_policy(x, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = skew_vector_nan_policy(x, bias, nan_policy)
    end function skew_vector_keepdims_nan_policy

    pure function kurtosis_vector_keepdims_nan_policy(x, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced to one retained dimension
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: values(1)

        values(1) = kurtosis_vector_nan_policy(x, fisher, bias, nan_policy)
    end function kurtosis_vector_keepdims_nan_policy

    pure function gmean_matrix_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = gmean_matrix_nan_policy(x, axis, weights, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function gmean_matrix_keepdims_nan_policy

    pure function gmean_matrix_keepdims_weights2_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = gmean_matrix_weights2_nan_policy(x, axis, weights, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function gmean_matrix_keepdims_weights2_nan_policy

    pure function hmean_matrix_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = hmean_matrix_nan_policy(x, axis, weights, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function hmean_matrix_keepdims_nan_policy

    pure function hmean_matrix_keepdims_weights2_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = hmean_matrix_weights2_nan_policy(x, axis, weights, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function hmean_matrix_keepdims_weights2_nan_policy

    pure function pmean_matrix_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = pmean_matrix_nan_policy(x, p, axis, weights, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function pmean_matrix_keepdims_nan_policy

    pure function pmean_matrix_keepdims_weights2_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = pmean_matrix_weights2_nan_policy(x, p, axis, weights, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function pmean_matrix_keepdims_weights2_nan_policy

    pure function sem_matrix_keepdims_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = sem_matrix_nan_policy(x, axis, ddof, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function sem_matrix_keepdims_nan_policy

    pure function variation_matrix_keepdims_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = variation_matrix_nan_policy(x, axis, ddof, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function variation_matrix_keepdims_nan_policy

    pure function iqr_matrix_keepdims_nan_policy(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimensions retained
        integer, intent(in), optional :: axis !! Fortran axis; absent reduces all dimensions
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: flattened(:)
        real(dp), allocatable :: reduced(:)

        if (.not. present(axis)) then
            allocate(values(1, 1))
            flattened = reshape(x, [size(x)])
            values(1, 1) = iqr_vector_nan_policy(flattened, rng, scale, nan_policy)
            return
        end if
        reduced = iqr_matrix_axis_nan_policy(x, axis, rng, scale, nan_policy)
        values = retain_matrix_axis(x, axis, reduced)
    end function iqr_matrix_keepdims_nan_policy

    pure function median_abs_deviation_matrix_keepdims_nan_policy(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = median_abs_deviation_matrix_nan_policy(x, axis, scale, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function median_abs_deviation_matrix_keepdims_nan_policy

    pure function moment_matrix_keepdims_nan_policy(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = moment_matrix_nan_policy(x, order, axis, center, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function moment_matrix_keepdims_nan_policy

    pure function skew_matrix_keepdims_nan_policy(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = skew_matrix_nan_policy(x, axis, bias, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function skew_matrix_keepdims_nan_policy

    pure function kurtosis_matrix_keepdims_nan_policy(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 or 2; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp), allocatable :: reduced(:)

        reduced = kurtosis_matrix_nan_policy(x, axis, fisher, bias, nan_policy)
        values = retain_matrix_axis(x, default_matrix_axis(axis), reduced)
    end function kurtosis_matrix_keepdims_nan_policy


    pure function gmean_tensor(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = gmean_vector(x(:, j, k))
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = gmean_vector(x(i, :, k))
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = gmean_vector(x(i, j, :))
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function gmean_tensor

    pure function hmean_tensor(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = hmean_vector(x(:, j, k))
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = hmean_vector(x(i, :, k))
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = hmean_vector(x(i, j, :))
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function hmean_tensor

    pure function pmean_tensor(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = pmean_vector(x(:, j, k), p)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = pmean_vector(x(i, :, k), p)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = pmean_vector(x(i, j, :), p)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function pmean_tensor

    pure function gmean_tensor_weights1(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = gmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
    end function gmean_tensor_weights1

    pure function gmean_tensor_weights2(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = gmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
    end function gmean_tensor_weights2

    pure function gmean_tensor_weights3(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = gmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
    end function gmean_tensor_weights3

    pure function hmean_tensor_weights1(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = hmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
    end function hmean_tensor_weights1

    pure function hmean_tensor_weights2(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = hmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
    end function hmean_tensor_weights2

    pure function hmean_tensor_weights3(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = hmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
    end function hmean_tensor_weights3

    pure function pmean_tensor_weights1(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = pmean_tensor_weighted_core(broadcast_x, p, broadcast_weights, axis)
    end function pmean_tensor_weights1

    pure function pmean_tensor_weights2(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = pmean_tensor_weighted_core(broadcast_x, p, broadcast_weights, axis)
    end function pmean_tensor_weights2

    pure function pmean_tensor_weights3(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = pmean_tensor_weighted_core(broadcast_x, p, broadcast_weights, axis)
    end function pmean_tensor_weights3


    pure function sem_tensor(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = sem_vector(x(:, j, k), ddof)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = sem_vector(x(i, :, k), ddof)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = sem_vector(x(i, j, :), ddof)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function sem_tensor

    pure function variation_tensor(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = variation_vector(x(:, j, k), ddof)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = variation_vector(x(i, :, k), ddof)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = variation_vector(x(i, j, :), ddof)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function variation_tensor

    pure function median_abs_deviation_tensor(x, axis, scale) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = median_abs_deviation_vector(x(:, j, k), scale)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = median_abs_deviation_vector(x(i, :, k), scale)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = median_abs_deviation_vector(x(i, j, :), scale)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function median_abs_deviation_tensor

    pure function moment_tensor(x, order, axis, center) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = moment_vector(x(:, j, k), order, center)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = moment_vector(x(i, :, k), order, center)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = moment_vector(x(i, j, :), order, center)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function moment_tensor

    pure function skew_tensor(x, axis, bias) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = skew_vector(x(:, j, k), bias)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = skew_vector(x(i, :, k), bias)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = skew_vector(x(i, j, :), bias)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function skew_tensor

    pure function kurtosis_tensor(x, axis, fisher, bias) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = kurtosis_vector(x(:, j, k), fisher, bias)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = kurtosis_vector(x(i, :, k), fisher, bias)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = kurtosis_vector(x(i, j, :), fisher, bias)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function kurtosis_tensor

    pure function iqr_tensor_all(x, rng, scale) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations reduced over all axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        real(dp) :: value
        real(dp), allocatable :: flattened(:)

        flattened = reshape(x, [size(x)])
        value = iqr_vector(flattened, rng, scale)
    end function iqr_tensor_all

    pure function iqr_tensor_axis(x, axis, rng, scale) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 to 3
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        real(dp), allocatable :: values(:, :)

        integer :: i
        integer :: j
        integer :: k

        select case (axis)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = iqr_vector(x(:, j, k), rng, scale)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = iqr_vector(x(i, :, k), rng, scale)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = iqr_vector(x(i, j, :), rng, scale)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function iqr_tensor_axis

    pure function gmean_tensor_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = gmean_vector_nan_policy(x(:, j, k), nan_policy=nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = gmean_vector_nan_policy(x(i, :, k), nan_policy=nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = gmean_vector_nan_policy(x(i, j, :), nan_policy=nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function gmean_tensor_nan_policy

    pure function hmean_tensor_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = hmean_vector_nan_policy(x(:, j, k), nan_policy=nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = hmean_vector_nan_policy(x(i, :, k), nan_policy=nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = hmean_vector_nan_policy(x(i, j, :), nan_policy=nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function hmean_tensor_nan_policy

    pure function pmean_tensor_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = pmean_vector_nan_policy(x(:, j, k), p, nan_policy=nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = pmean_vector_nan_policy(x(i, :, k), p, nan_policy=nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = pmean_vector_nan_policy(x(i, j, :), p, nan_policy=nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function pmean_tensor_nan_policy

    pure function gmean_tensor_weights1_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = gmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
    end function gmean_tensor_weights1_nan_policy

    pure function gmean_tensor_weights2_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = gmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
    end function gmean_tensor_weights2_nan_policy

    pure function gmean_tensor_weights3_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = gmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
    end function gmean_tensor_weights3_nan_policy

    pure function hmean_tensor_weights1_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = hmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
    end function hmean_tensor_weights1_nan_policy

    pure function hmean_tensor_weights2_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = hmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
    end function hmean_tensor_weights2_nan_policy

    pure function hmean_tensor_weights3_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = hmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
    end function hmean_tensor_weights3_nan_policy

    pure function pmean_tensor_weights1_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = pmean_tensor_weighted_core( &
            broadcast_x, p, broadcast_weights, axis, nan_policy)
    end function pmean_tensor_weights1_nan_policy

    pure function pmean_tensor_weights2_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = pmean_tensor_weighted_core( &
            broadcast_x, p, broadcast_weights, axis, nan_policy)
    end function pmean_tensor_weights2_nan_policy

    pure function pmean_tensor_weights3_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)

        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            values = tensor_nan_reduced(shape(x), default_tensor_axis(axis))
            return
        end if
        values = pmean_tensor_weighted_core( &
            broadcast_x, p, broadcast_weights, axis, nan_policy)
    end function pmean_tensor_weights3_nan_policy


    pure function sem_tensor_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = sem_vector_nan_policy(x(:, j, k), ddof, nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = sem_vector_nan_policy(x(i, :, k), ddof, nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = sem_vector_nan_policy(x(i, j, :), ddof, nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function sem_tensor_nan_policy

    pure function variation_tensor_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = variation_vector_nan_policy(x(:, j, k), ddof, nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = variation_vector_nan_policy(x(i, :, k), ddof, nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = variation_vector_nan_policy(x(i, j, :), ddof, nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function variation_tensor_nan_policy

    pure function median_abs_deviation_tensor_nan_policy(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = median_abs_deviation_vector_nan_policy(x(:, j, k), scale, nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = median_abs_deviation_vector_nan_policy(x(i, :, k), scale, nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = median_abs_deviation_vector_nan_policy(x(i, j, :), scale, nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function median_abs_deviation_tensor_nan_policy

    pure function moment_tensor_nan_policy(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = moment_vector_nan_policy(x(:, j, k), order, center, nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = moment_vector_nan_policy(x(i, :, k), order, center, nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = moment_vector_nan_policy(x(i, j, :), order, center, nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function moment_tensor_nan_policy

    pure function skew_tensor_nan_policy(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = skew_vector_nan_policy(x(:, j, k), bias, nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = skew_vector_nan_policy(x(i, :, k), bias, nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = skew_vector_nan_policy(x(i, j, :), bias, nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function skew_tensor_nan_policy

    pure function kurtosis_tensor_nan_policy(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = kurtosis_vector_nan_policy(x(:, j, k), fisher, bias, nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = kurtosis_vector_nan_policy(x(i, :, k), fisher, bias, nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = kurtosis_vector_nan_policy(x(i, j, :), fisher, bias, nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function kurtosis_tensor_nan_policy

    pure function iqr_tensor_all_nan_policy(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations reduced over all axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp), allocatable :: flattened(:)

        flattened = reshape(x, [size(x)])
        value = iqr_vector_nan_policy(flattened, rng, scale, nan_policy)
    end function iqr_tensor_all_nan_policy

    pure function iqr_tensor_axis_nan_policy(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 to 3
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: i
        integer :: j
        integer :: k

        select case (axis)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(j, k) = iqr_vector_nan_policy(x(:, j, k), rng, scale, nan_policy)
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, k) = iqr_vector_nan_policy(x(i, :, k), rng, scale, nan_policy)
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j) = iqr_vector_nan_policy(x(i, j, :), rng, scale, nan_policy)
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function iqr_tensor_axis_nan_policy

    pure function gmean_tensor_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = gmean_tensor(x, axis)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function gmean_tensor_keepdims

    pure function hmean_tensor_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = hmean_tensor(x, axis)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function hmean_tensor_keepdims

    pure function pmean_tensor_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = pmean_tensor(x, p, axis)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function pmean_tensor_keepdims

    pure function gmean_tensor_keepdims_weights1(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = gmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function gmean_tensor_keepdims_weights1

    pure function gmean_tensor_keepdims_weights2(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = gmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function gmean_tensor_keepdims_weights2

    pure function gmean_tensor_keepdims_weights3(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = gmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function gmean_tensor_keepdims_weights3

    pure function hmean_tensor_keepdims_weights1(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = hmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function hmean_tensor_keepdims_weights1

    pure function hmean_tensor_keepdims_weights2(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = hmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function hmean_tensor_keepdims_weights2

    pure function hmean_tensor_keepdims_weights3(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = hmean_tensor_weighted_core(broadcast_x, broadcast_weights, axis)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function hmean_tensor_keepdims_weights3

    pure function pmean_tensor_keepdims_weights1(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = pmean_tensor_weighted_core(broadcast_x, p, broadcast_weights, axis)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function pmean_tensor_keepdims_weights1

    pure function pmean_tensor_keepdims_weights2(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = pmean_tensor_weighted_core(broadcast_x, p, broadcast_weights, axis)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function pmean_tensor_keepdims_weights2

    pure function pmean_tensor_keepdims_weights3(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = pmean_tensor_weighted_core(broadcast_x, p, broadcast_weights, axis)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function pmean_tensor_keepdims_weights3


    pure function sem_tensor_keepdims(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = sem_tensor(x, axis, ddof)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function sem_tensor_keepdims

    pure function variation_tensor_keepdims(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = variation_tensor(x, axis, ddof)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function variation_tensor_keepdims

    pure function median_abs_deviation_tensor_keepdims(x, axis, scale) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = median_abs_deviation_tensor(x, axis, scale)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function median_abs_deviation_tensor_keepdims

    pure function moment_tensor_keepdims(x, order, axis, center) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = moment_tensor(x, order, axis, center)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function moment_tensor_keepdims

    pure function skew_tensor_keepdims(x, axis, bias) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = skew_tensor(x, axis, bias)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function skew_tensor_keepdims

    pure function kurtosis_tensor_keepdims(x, axis, fisher, bias) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = kurtosis_tensor(x, axis, fisher, bias)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function kurtosis_tensor_keepdims

    pure function iqr_tensor_keepdims(x, axis, rng, scale) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimensions retained
        integer, intent(in), optional :: axis !! Fortran axis; absent reduces all dimensions
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        if (.not. present(axis)) then
            allocate(values(1, 1, 1))
            values(1, 1, 1) = iqr_tensor_all(x, rng, scale)
            return
        end if
        reduced = iqr_tensor_axis(x, axis, rng, scale)
        values = retain_tensor_axis(x, axis, reduced)
    end function iqr_tensor_keepdims

    pure function gmean_tensor_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = gmean_tensor_nan_policy(x, axis, nan_policy)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function gmean_tensor_keepdims_nan_policy

    pure function hmean_tensor_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = hmean_tensor_nan_policy(x, axis, nan_policy)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function hmean_tensor_keepdims_nan_policy

    pure function pmean_tensor_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = pmean_tensor_nan_policy(x, p, axis, nan_policy)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function pmean_tensor_keepdims_nan_policy

    pure function gmean_tensor_keepdims_weights1_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = gmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function gmean_tensor_keepdims_weights1_nan_policy

    pure function gmean_tensor_keepdims_weights2_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = gmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function gmean_tensor_keepdims_weights2_nan_policy

    pure function gmean_tensor_keepdims_weights3_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = gmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function gmean_tensor_keepdims_weights3_nan_policy

    pure function hmean_tensor_keepdims_weights1_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = hmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function hmean_tensor_keepdims_weights1_nan_policy

    pure function hmean_tensor_keepdims_weights2_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = hmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function hmean_tensor_keepdims_weights2_nan_policy

    pure function hmean_tensor_keepdims_weights3_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = hmean_tensor_weighted_core( &
            broadcast_x, broadcast_weights, axis, nan_policy)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function hmean_tensor_keepdims_weights3_nan_policy

    pure function pmean_tensor_keepdims_weights1_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = pmean_tensor_weighted_core( &
            broadcast_x, p, broadcast_weights, axis, nan_policy)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function pmean_tensor_keepdims_weights1_nan_policy

    pure function pmean_tensor_keepdims_weights2_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = pmean_tensor_weighted_core( &
            broadcast_x, p, broadcast_weights, axis, nan_policy)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function pmean_tensor_keepdims_weights2_nan_policy

    pure function pmean_tensor_keepdims_weights3_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; broadcast with weights before reduction
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in) :: weights(:, :, :) !! trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        logical :: valid
        real(dp), allocatable :: broadcast_x(:, :, :)
        real(dp), allocatable :: broadcast_weights(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        axis_value = default_tensor_axis(axis)
        call tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        if (.not. valid) then
            reduced = tensor_nan_reduced(shape(x), axis_value)
            values = retain_tensor_shape(shape(x), axis_value, reduced)
            return
        end if
        reduced = pmean_tensor_weighted_core( &
            broadcast_x, p, broadcast_weights, axis, nan_policy)
        values = retain_tensor_shape(shape(broadcast_x), axis_value, reduced)
    end function pmean_tensor_keepdims_weights3_nan_policy


    pure function sem_tensor_keepdims_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = sem_tensor_nan_policy(x, axis, ddof, nan_policy)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function sem_tensor_keepdims_nan_policy

    pure function variation_tensor_keepdims_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = variation_tensor_nan_policy(x, axis, ddof, nan_policy)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function variation_tensor_keepdims_nan_policy

    pure function median_abs_deviation_tensor_keepdims_nan_policy(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = median_abs_deviation_tensor_nan_policy(x, axis, scale, nan_policy)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function median_abs_deviation_tensor_keepdims_nan_policy

    pure function moment_tensor_keepdims_nan_policy(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = moment_tensor_nan_policy(x, order, axis, center, nan_policy)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function moment_tensor_keepdims_nan_policy

    pure function skew_tensor_keepdims_nan_policy(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = skew_tensor_nan_policy(x, axis, bias, nan_policy)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function skew_tensor_keepdims_nan_policy

    pure function kurtosis_tensor_keepdims_nan_policy(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        reduced = kurtosis_tensor_nan_policy(x, axis, fisher, bias, nan_policy)
        values = retain_tensor_axis(x, default_tensor_axis(axis), reduced)
    end function kurtosis_tensor_keepdims_nan_policy

    pure function iqr_tensor_keepdims_nan_policy(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; reduced with dimensions retained
        integer, intent(in), optional :: axis !! Fortran axis; absent reduces all dimensions
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp), allocatable :: reduced(:, :)

        if (.not. present(axis)) then
            allocate(values(1, 1, 1))
            values(1, 1, 1) = iqr_tensor_all_nan_policy(x, rng, scale, nan_policy)
            return
        end if
        reduced = iqr_tensor_axis_nan_policy(x, axis, rng, scale, nan_policy)
        values = retain_tensor_axis(x, axis, reduced)
    end function iqr_tensor_keepdims_nan_policy

    pure function gmean_rank4(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = gmean_vector(x(:, j, k, l))
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = gmean_vector(x(i, :, k, l))
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = gmean_vector(x(i, j, :, l))
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = gmean_vector(x(i, j, k, :))
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function gmean_rank4

    pure function hmean_rank4(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = hmean_vector(x(:, j, k, l))
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = hmean_vector(x(i, :, k, l))
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = hmean_vector(x(i, j, :, l))
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = hmean_vector(x(i, j, k, :))
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function hmean_rank4

    pure function pmean_rank4(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = pmean_vector(x(:, j, k, l), p)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = pmean_vector(x(i, :, k, l), p)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = pmean_vector(x(i, j, :, l), p)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = pmean_vector(x(i, j, k, :), p)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function pmean_rank4

    pure function sem_rank4(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = sem_vector(x(:, j, k, l), ddof)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = sem_vector(x(i, :, k, l), ddof)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = sem_vector(x(i, j, :, l), ddof)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = sem_vector(x(i, j, k, :), ddof)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function sem_rank4

    pure function variation_rank4(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = variation_vector(x(:, j, k, l), ddof)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = variation_vector(x(i, :, k, l), ddof)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = variation_vector(x(i, j, :, l), ddof)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = variation_vector(x(i, j, k, :), ddof)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function variation_rank4

    pure function median_abs_deviation_rank4(x, axis, scale) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = median_abs_deviation_vector(x(:, j, k, l), scale)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = median_abs_deviation_vector(x(i, :, k, l), scale)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = median_abs_deviation_vector(x(i, j, :, l), scale)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = median_abs_deviation_vector(x(i, j, k, :), scale)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function median_abs_deviation_rank4

    pure function moment_rank4(x, order, axis, center) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
                integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                real(dp), intent(in), optional :: center !! finite center; default slice mean
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = moment_vector(x(:, j, k, l), order, center)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = moment_vector(x(i, :, k, l), order, center)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = moment_vector(x(i, j, :, l), order, center)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = moment_vector(x(i, j, k, :), order, center)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function moment_rank4

    pure function skew_rank4(x, axis, bias) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = skew_vector(x(:, j, k, l), bias)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = skew_vector(x(i, :, k, l), bias)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = skew_vector(x(i, j, :, l), bias)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = skew_vector(x(i, j, k, :), bias)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function skew_rank4

    pure function kurtosis_rank4(x, axis, fisher, bias) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        logical, intent(in), optional :: fisher !! subtract 3, default true
                logical, intent(in), optional :: bias !! use biased sample moments, default true
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = kurtosis_vector(x(:, j, k, l), fisher, bias)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = kurtosis_vector(x(i, :, k, l), fisher, bias)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = kurtosis_vector(x(i, j, :, l), fisher, bias)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = kurtosis_vector(x(i, j, k, :), fisher, bias)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function kurtosis_rank4

    pure function gmean_rank4_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = gmean_vector_nan_policy(x(:, j, k, l), nan_policy=nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = gmean_vector_nan_policy(x(i, :, k, l), nan_policy=nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = gmean_vector_nan_policy(x(i, j, :, l), nan_policy=nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = gmean_vector_nan_policy(x(i, j, k, :), nan_policy=nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function gmean_rank4_nan_policy

    pure function hmean_rank4_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = hmean_vector_nan_policy(x(:, j, k, l), nan_policy=nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = hmean_vector_nan_policy(x(i, :, k, l), nan_policy=nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = hmean_vector_nan_policy(x(i, j, :, l), nan_policy=nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = hmean_vector_nan_policy(x(i, j, k, :), nan_policy=nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function hmean_rank4_nan_policy

    pure function pmean_rank4_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = pmean_vector_nan_policy(x(:, j, k, l), p, nan_policy=nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = pmean_vector_nan_policy(x(i, :, k, l), p, nan_policy=nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = pmean_vector_nan_policy(x(i, j, :, l), p, nan_policy=nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = pmean_vector_nan_policy(x(i, j, k, :), p, nan_policy=nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function pmean_rank4_nan_policy

    pure function sem_rank4_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = sem_vector_nan_policy(x(:, j, k, l), ddof, nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = sem_vector_nan_policy(x(i, :, k, l), ddof, nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = sem_vector_nan_policy(x(i, j, :, l), ddof, nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = sem_vector_nan_policy(x(i, j, k, :), ddof, nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function sem_rank4_nan_policy

    pure function variation_rank4_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = variation_vector_nan_policy(x(:, j, k, l), ddof, nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = variation_vector_nan_policy(x(i, :, k, l), ddof, nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = variation_vector_nan_policy(x(i, j, :, l), ddof, nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = variation_vector_nan_policy(x(i, j, k, :), ddof, nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function variation_rank4_nan_policy

    pure function median_abs_deviation_rank4_nan_policy(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = median_abs_deviation_vector_nan_policy(x(:, j, k, l), scale, nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = median_abs_deviation_vector_nan_policy(x(i, :, k, l), scale, nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = median_abs_deviation_vector_nan_policy(x(i, j, :, l), scale, nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = median_abs_deviation_vector_nan_policy(x(i, j, k, :), scale, nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function median_abs_deviation_rank4_nan_policy

    pure function moment_rank4_nan_policy(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
                integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = moment_vector_nan_policy(x(:, j, k, l), order, center, nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = moment_vector_nan_policy(x(i, :, k, l), order, center, nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = moment_vector_nan_policy(x(i, j, :, l), order, center, nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = moment_vector_nan_policy(x(i, j, k, :), order, center, nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function moment_rank4_nan_policy

    pure function skew_rank4_nan_policy(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = skew_vector_nan_policy(x(:, j, k, l), bias, nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = skew_vector_nan_policy(x(i, :, k, l), bias, nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = skew_vector_nan_policy(x(i, j, :, l), bias, nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = skew_vector_nan_policy(x(i, j, k, :), bias, nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function skew_rank4_nan_policy

    pure function kurtosis_rank4_nan_policy(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        logical, intent(in), optional :: fisher !! subtract 3, default true
                logical, intent(in), optional :: bias !! use biased sample moments, default true
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = default_rank4_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = kurtosis_vector_nan_policy(x(:, j, k, l), fisher, bias, nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = kurtosis_vector_nan_policy(x(i, :, k, l), fisher, bias, nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = kurtosis_vector_nan_policy(x(i, j, :, l), fisher, bias, nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = kurtosis_vector_nan_policy(x(i, j, k, :), fisher, bias, nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function kurtosis_rank4_nan_policy

    pure function iqr_rank4_all(x, rng, scale) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations reduced over all axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        real(dp) :: value
        real(dp), allocatable :: flattened(:)

        flattened = reshape(x, [size(x)])
        value = iqr_vector(flattened, rng, scale)
    end function iqr_rank4_all

    pure function iqr_rank4_axis(x, axis, rng, scale) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 to 4
                real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
                real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = axis
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = iqr_vector(x(:, j, k, l), rng, scale)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = iqr_vector(x(i, :, k, l), rng, scale)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = iqr_vector(x(i, j, :, l), rng, scale)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = iqr_vector(x(i, j, k, :), rng, scale)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function iqr_rank4_axis

    pure function iqr_rank4_all_nan_policy(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations reduced over all axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp), allocatable :: flattened(:)

        flattened = reshape(x, [size(x)])
        value = iqr_vector_nan_policy(flattened, rng, scale, nan_policy)
    end function iqr_rank4_all_nan_policy

    pure function iqr_rank4_axis_nan_policy(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 to 4
                real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
                real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        integer :: l

        axis_value = axis
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(j, k, l) = iqr_vector_nan_policy(x(:, j, k, l), rng, scale, nan_policy)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, k, l) = iqr_vector_nan_policy(x(i, :, k, l), rng, scale, nan_policy)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, l) = iqr_vector_nan_policy(x(i, j, :, l), rng, scale, nan_policy)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k) = iqr_vector_nan_policy(x(i, j, k, :), rng, scale, nan_policy)
                    end do
                end do
            end do
        case default
            allocate(values(0, 0, 0))
        end select
    end function iqr_rank4_axis_nan_policy

    pure function gmean_rank4_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = gmean_rank4(x, axis)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function gmean_rank4_keepdims

    pure function hmean_rank4_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = hmean_rank4(x, axis)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function hmean_rank4_keepdims

    pure function pmean_rank4_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        real(dp), intent(in) :: p !! finite power
                integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = pmean_rank4(x, p, axis)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function pmean_rank4_keepdims

    pure function sem_rank4_keepdims(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = sem_rank4(x, axis, ddof)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function sem_rank4_keepdims

    pure function variation_rank4_keepdims(x, axis, ddof) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = variation_rank4(x, axis, ddof)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function variation_rank4_keepdims

    pure function median_abs_deviation_rank4_keepdims(x, axis, scale) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = median_abs_deviation_rank4(x, axis, scale)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function median_abs_deviation_rank4_keepdims

    pure function moment_rank4_keepdims(x, order, axis, center) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
                integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                real(dp), intent(in), optional :: center !! finite center; default slice mean
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = moment_rank4(x, order, axis, center)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function moment_rank4_keepdims

    pure function skew_rank4_keepdims(x, axis, bias) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = skew_rank4(x, axis, bias)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function skew_rank4_keepdims

    pure function kurtosis_rank4_keepdims(x, axis, fisher, bias) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                logical, intent(in), optional :: fisher !! subtract 3, default true
                logical, intent(in), optional :: bias !! use biased sample moments, default true
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = kurtosis_rank4(x, axis, fisher, bias)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function kurtosis_rank4_keepdims

    pure function gmean_rank4_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = gmean_rank4_nan_policy(x, axis, nan_policy)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function gmean_rank4_keepdims_nan_policy

    pure function hmean_rank4_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = hmean_rank4_nan_policy(x, axis, nan_policy)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function hmean_rank4_keepdims_nan_policy

    pure function pmean_rank4_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        real(dp), intent(in) :: p !! finite power
                integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = pmean_rank4_nan_policy(x, p, axis, nan_policy)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function pmean_rank4_keepdims_nan_policy

    pure function sem_rank4_keepdims_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = sem_rank4_nan_policy(x, axis, ddof, nan_policy)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function sem_rank4_keepdims_nan_policy

    pure function variation_rank4_keepdims_nan_policy(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = variation_rank4_nan_policy(x, axis, ddof, nan_policy)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function variation_rank4_keepdims_nan_policy

    pure function median_abs_deviation_rank4_keepdims_nan_policy(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = median_abs_deviation_rank4_nan_policy(x, axis, scale, nan_policy)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function median_abs_deviation_rank4_keepdims_nan_policy

    pure function moment_rank4_keepdims_nan_policy(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
                integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = moment_rank4_nan_policy(x, order, axis, center, nan_policy)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function moment_rank4_keepdims_nan_policy

    pure function skew_rank4_keepdims_nan_policy(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = skew_rank4_nan_policy(x, axis, bias, nan_policy)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function skew_rank4_keepdims_nan_policy

    pure function kurtosis_rank4_keepdims_nan_policy(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
                logical, intent(in), optional :: fisher !! subtract 3, default true
                logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        reduced = kurtosis_rank4_nan_policy(x, axis, fisher, bias, nan_policy)
        values = retain_rank4_axis(x, default_rank4_axis(axis), reduced)
    end function kurtosis_rank4_keepdims_nan_policy

    pure function iqr_rank4_keepdims(x, axis, rng, scale) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimensions retained
        integer, intent(in), optional :: axis !! Fortran axis; absent reduces all dimensions
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        if (.not. present(axis)) then
            allocate(values(1, 1, 1, 1))
            values(1, 1, 1, 1) = iqr_rank4_all(x, rng, scale)
            return
        end if
        reduced = iqr_rank4_axis(x, axis, rng, scale)
        values = retain_rank4_axis(x, axis, reduced)
    end function iqr_rank4_keepdims

    pure function iqr_rank4_keepdims_nan_policy(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; reduced with dimensions retained
        integer, intent(in), optional :: axis !! Fortran axis; absent reduces all dimensions
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp), allocatable :: reduced(:, :, :)

        if (.not. present(axis)) then
            allocate(values(1, 1, 1, 1))
            values(1, 1, 1, 1) = iqr_rank4_all_nan_policy(x, rng, scale, nan_policy)
            return
        end if
        reduced = iqr_rank4_axis_nan_policy(x, axis, rng, scale, nan_policy)
        values = retain_rank4_axis(x, axis, reduced)
    end function iqr_rank4_keepdims_nan_policy

    pure function default_rank4_axis(axis) result(axis_value)
        integer, intent(in), optional :: axis !! optional Fortran axis
        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
    end function default_rank4_axis

    pure function retain_rank4_axis(x, axis, reduced) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! original array defining retained shape
        integer, intent(in) :: axis !! reduced Fortran axis, 1 to 4
        real(dp), intent(in) :: reduced(:, :, :) !! reduced values
        real(dp), allocatable :: values(:, :, :, :)

        select case (axis)
        case (1)
            if (all(shape(reduced) == [size(x, 2), size(x, 3), size(x, 4)])) then
                allocate(values(1, size(x, 2), size(x, 3), size(x, 4)))
                values(1, :, :, :) = reduced
            else
                allocate(values(0, 0, 0, 0))
            end if
        case (2)
            if (all(shape(reduced) == [size(x, 1), size(x, 3), size(x, 4)])) then
                allocate(values(size(x, 1), 1, size(x, 3), size(x, 4)))
                values(:, 1, :, :) = reduced
            else
                allocate(values(0, 0, 0, 0))
            end if
        case (3)
            if (all(shape(reduced) == [size(x, 1), size(x, 2), size(x, 4)])) then
                allocate(values(size(x, 1), size(x, 2), 1, size(x, 4)))
                values(:, :, 1, :) = reduced
            else
                allocate(values(0, 0, 0, 0))
            end if
        case (4)
            if (all(shape(reduced) == [size(x, 1), size(x, 2), size(x, 3)])) then
                allocate(values(size(x, 1), size(x, 2), size(x, 3), 1))
                values(:, :, :, 1) = reduced
            else
                allocate(values(0, 0, 0, 0))
            end if
        case default
            allocate(values(0, 0, 0, 0))
        end select
    end function retain_rank4_axis


    pure function covariance(x, y, ddof) result(value)
        real(dp), intent(in) :: x(:) !! first observations; all values must be finite
        real(dp), intent(in) :: y(:) !! second observations, same size as x
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp) :: value

        integer :: divisor_adjustment
        integer :: i
        real(dp) :: correction
        real(dp) :: factor
        real(dp) :: mean_x
        real(dp) :: mean_y
        real(dp) :: scale_x
        real(dp) :: scale_y
        real(dp) :: sum_cross
        real(dp) :: term
        real(dp) :: updated
        real(dp) :: xs
        real(dp) :: ys

        divisor_adjustment = 1
        if (present(ddof)) divisor_adjustment = ddof

        if (size(x) /= size(y) .or. size(x) == 0 .or. &
                divisor_adjustment < 0 .or. size(x) <= divisor_adjustment) then
            value = quiet_nan(0.0_dp)
            return
        end if

        scale_x = 0.0_dp
        scale_y = 0.0_dp
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i))) then
                value = quiet_nan(0.0_dp)
                return
            end if
            scale_x = max(scale_x, abs(x(i)))
            scale_y = max(scale_y, abs(y(i)))
        end do
        if (scale_x <= 0.0_dp .or. scale_y <= 0.0_dp) then
            value = 0.0_dp
            return
        end if

        mean_x = mean(x) / scale_x
        mean_y = mean(y) / scale_y
        sum_cross = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(x)
            xs = x(i) / scale_x - mean_x
            ys = y(i) / scale_y - mean_y
            term = xs * ys - correction
            updated = sum_cross + term
            correction = (updated - sum_cross) - term
            sum_cross = updated
        end do
        factor = sum_cross / real(size(x) - divisor_adjustment, dp)
        value = safe_product(safe_product(factor, min(scale_x, scale_y)), &
            max(scale_x, scale_y))
    end function covariance

    pure function pearson_correlation(x, y) result(value)
        real(dp), intent(in) :: x(:) !! first observations; all values must be finite
        real(dp), intent(in) :: y(:) !! second observations, same size as x
        real(dp) :: value

        integer :: i
        real(dp) :: correction_cross
        real(dp) :: correction_x
        real(dp) :: correction_y
        real(dp) :: mean_x
        real(dp) :: mean_y
        real(dp) :: scale_x
        real(dp) :: scale_y
        real(dp) :: sum_cross
        real(dp) :: sum_x2
        real(dp) :: sum_y2
        real(dp) :: term
        real(dp) :: updated
        real(dp) :: xs
        real(dp) :: ys

        if (size(x) /= size(y) .or. size(x) < 2) then
            value = quiet_nan(0.0_dp)
            return
        end if

        scale_x = 0.0_dp
        scale_y = 0.0_dp
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i))) then
                value = quiet_nan(0.0_dp)
                return
            end if
            scale_x = max(scale_x, abs(x(i)))
            scale_y = max(scale_y, abs(y(i)))
        end do
        if (scale_x <= 0.0_dp .or. scale_y <= 0.0_dp) then
            value = quiet_nan(0.0_dp)
            return
        end if

        mean_x = mean(x) / scale_x
        mean_y = mean(y) / scale_y
        sum_cross = 0.0_dp
        sum_x2 = 0.0_dp
        sum_y2 = 0.0_dp
        correction_cross = 0.0_dp
        correction_x = 0.0_dp
        correction_y = 0.0_dp
        do i = 1, size(x)
            xs = x(i) / scale_x - mean_x
            ys = y(i) / scale_y - mean_y

            term = xs * ys - correction_cross
            updated = sum_cross + term
            correction_cross = (updated - sum_cross) - term
            sum_cross = updated

            term = xs * xs - correction_x
            updated = sum_x2 + term
            correction_x = (updated - sum_x2) - term
            sum_x2 = updated

            term = ys * ys - correction_y
            updated = sum_y2 + term
            correction_y = (updated - sum_y2) - term
            sum_y2 = updated
        end do

        if (sum_x2 <= 0.0_dp .or. sum_y2 <= 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else
            value = sum_cross / sqrt(sum_x2 * sum_y2)
            value = max(-1.0_dp, min(1.0_dp, value))
        end if
    end function pearson_correlation

    pure function quantile_scalar(x, q) result(value)
        real(dp), intent(in) :: x(:) !! observations; NaN propagates
        real(dp), intent(in) :: q !! probability in [0, 1]
        real(dp) :: value

        integer :: i
        real(dp), allocatable :: work(:)

        if (size(x) == 0 .or. ieee_is_nan(q) .or. q < 0.0_dp .or. q > 1.0_dp) then
            value = quiet_nan(0.0_dp)
            return
        end if
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            end if
        end do

        allocate(work(size(x)))
        work = x
        call sort_values(work)
        value = sorted_quantile(work, q)
    end function quantile_scalar

    pure function quantile_array(x, q) result(values)
        real(dp), intent(in) :: x(:) !! observations; NaN propagates
        real(dp), intent(in) :: q(:) !! probabilities, each in [0, 1]
        real(dp), allocatable :: values(:)

        integer :: i
        real(dp), allocatable :: work(:)

        allocate(values(size(q)))
        if (size(x) == 0) then
            values = quiet_nan(0.0_dp)
            return
        end if
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) then
                values = quiet_nan(0.0_dp)
                return
            end if
        end do
        do i = 1, size(q)
            if (ieee_is_nan(q(i)) .or. q(i) < 0.0_dp .or. q(i) > 1.0_dp) then
                values = quiet_nan(0.0_dp)
                return
            end if
        end do

        allocate(work(size(x)))
        work = x
        call sort_values(work)
        do i = 1, size(q)
            values(i) = sorted_quantile(work, q(i))
        end do
    end function quantile_array

    pure function median(x) result(value)
        real(dp), intent(in) :: x(:) !! observations; NaN propagates
        real(dp) :: value

        value = quantile_scalar(x, 0.5_dp)
    end function median

    pure function rankdata(x, method) result(ranks)
        real(dp), intent(in) :: x(:) !! observations; NaN propagates to all ranks
        character(len=*), intent(in), optional :: method !! average, min, max, dense, or ordinal
        real(dp), allocatable :: ranks(:)

        character(len=7) :: selected_method
        integer :: dense_rank
        integer :: group_end
        integer :: group_start
        integer :: i
        integer, allocatable :: order(:)

        allocate(ranks(size(x)))
        if (size(x) == 0) return

        selected_method = 'average'
        if (present(method)) selected_method = method
        if (selected_method /= 'average' .and. selected_method /= 'min' .and. &
                selected_method /= 'max' .and. selected_method /= 'dense' .and. &
                selected_method /= 'ordinal') then
            ranks = quiet_nan(0.0_dp)
            return
        end if

        do i = 1, size(x)
            if (ieee_is_nan(x(i))) then
                ranks = quiet_nan(0.0_dp)
                return
            end if
        end do

        allocate(order(size(x)))
        do i = 1, size(x)
            order(i) = i
        end do
        call sort_indices_by_values(x, order)

        if (selected_method == 'ordinal') then
            do i = 1, size(x)
                ranks(order(i)) = real(i, dp)
            end do
            return
        end if

        group_start = 1
        dense_rank = 1
        do while (group_start <= size(x))
            group_end = group_start
            do while (group_end < size(x))
                if (.not. same_value(x(order(group_end)), x(order(group_end + 1)))) exit
                group_end = group_end + 1
            end do

            select case (selected_method)
            case ('average')
                do i = group_start, group_end
                    ranks(order(i)) = 0.5_dp * real(group_start + group_end, dp)
                end do
            case ('min')
                do i = group_start, group_end
                    ranks(order(i)) = real(group_start, dp)
                end do
            case ('max')
                do i = group_start, group_end
                    ranks(order(i)) = real(group_end, dp)
                end do
            case ('dense')
                do i = group_start, group_end
                    ranks(order(i)) = real(dense_rank, dp)
                end do
            end select

            dense_rank = dense_rank + 1
            group_start = group_end + 1
        end do
    end function rankdata

    pure function sorted_quantile(sorted, q) result(value)
        real(dp), intent(in) :: sorted(:) !! sorted observations, nonempty and without NaN
        real(dp), intent(in) :: q !! probability in [0, 1]
        real(dp) :: value

        integer :: lower_index
        real(dp) :: fraction
        real(dp) :: position
        real(dp) :: upper
        real(dp) :: lower

        if (size(sorted) == 1) then
            value = sorted(1)
            return
        end if

        position = q * real(size(sorted) - 1, dp)
        lower_index = int(floor(position)) + 1
        fraction = position - floor(position)
        if (fraction <= 0.0_dp .or. lower_index == size(sorted)) then
            value = sorted(lower_index)
            return
        end if

        lower = sorted(lower_index)
        upper = sorted(lower_index + 1)
        if (ieee_is_finite(lower) .and. ieee_is_finite(upper)) then
            value = (1.0_dp - fraction) * lower + fraction * upper
        else if (lower < 0.0_dp .and. upper < 0.0_dp) then
            value = negative_infinity(0.0_dp)
        else if (lower > 0.0_dp .and. upper > 0.0_dp) then
            value = positive_infinity(0.0_dp)
        else if (lower < 0.0_dp .and. ieee_is_finite(upper)) then
            value = negative_infinity(0.0_dp)
        else if (ieee_is_finite(lower) .and. upper > 0.0_dp) then
            value = positive_infinity(0.0_dp)
        else
            value = quiet_nan(0.0_dp)
        end if
    end function sorted_quantile

    pure subroutine sort_values(values)
        real(dp), intent(inout) :: values(:) !! values sorted ascending in place

        real(dp), allocatable :: temporary(:)

        if (size(values) <= 1) return
        allocate(temporary(size(values)))
        call merge_sort_values(values, temporary, 1, size(values))
    end subroutine sort_values

    recursive pure subroutine merge_sort_values(values, temporary, left, right)
        real(dp), intent(inout) :: values(:) !! values being sorted
        real(dp), intent(inout) :: temporary(:) !! merge workspace, same size as values
        integer, intent(in) :: left !! first index of the inclusive range
        integer, intent(in) :: right !! last index of the inclusive range

        integer :: i
        integer :: j
        integer :: k
        integer :: middle

        if (left >= right) return
        middle = left + (right - left) / 2
        call merge_sort_values(values, temporary, left, middle)
        call merge_sort_values(values, temporary, middle + 1, right)

        i = left
        j = middle + 1
        do k = left, right
            if (i > middle) then
                temporary(k) = values(j)
                j = j + 1
            else if (j > right) then
                temporary(k) = values(i)
                i = i + 1
            else if (values(i) <= values(j)) then
                temporary(k) = values(i)
                i = i + 1
            else
                temporary(k) = values(j)
                j = j + 1
            end if
        end do
        values(left:right) = temporary(left:right)
    end subroutine merge_sort_values

    pure subroutine sort_indices_by_values(values, order)
        real(dp), intent(in) :: values(:) !! source values used as sort keys
        integer, intent(inout) :: order(:) !! indices sorted by value, preserving tie order

        integer, allocatable :: temporary(:)

        if (size(order) <= 1) return
        allocate(temporary(size(order)))
        call merge_sort_indices(values, order, temporary, 1, size(order))
    end subroutine sort_indices_by_values

    recursive pure subroutine merge_sort_indices(values, order, temporary, left, right)
        real(dp), intent(in) :: values(:) !! source values used as sort keys
        integer, intent(inout) :: order(:) !! indices being sorted
        integer, intent(inout) :: temporary(:) !! merge workspace, same size as order
        integer, intent(in) :: left !! first index of the inclusive range
        integer, intent(in) :: right !! last index of the inclusive range

        integer :: i
        integer :: j
        integer :: k
        integer :: middle

        if (left >= right) return
        middle = left + (right - left) / 2
        call merge_sort_indices(values, order, temporary, left, middle)
        call merge_sort_indices(values, order, temporary, middle + 1, right)

        i = left
        j = middle + 1
        do k = left, right
            if (i > middle) then
                temporary(k) = order(j)
                j = j + 1
            else if (j > right) then
                temporary(k) = order(i)
                i = i + 1
            else if (values(order(i)) <= values(order(j))) then
                temporary(k) = order(i)
                i = i + 1
            else
                temporary(k) = order(j)
                j = j + 1
            end if
        end do
        order(left:right) = temporary(left:right)
    end subroutine merge_sort_indices

    pure logical function valid_nan_policy(nan_policy) result(valid)
        character(len=*), intent(in) :: nan_policy !! requested NaN handling policy

        select case (trim(adjustl(nan_policy)))
        case ('propagate', 'omit', 'raise')
            valid = .true.
        case default
            valid = .false.
        end select
    end function valid_nan_policy

    pure logical function nan_policy_raises(nan_policy) result(raise_nan)
        character(len=*), intent(in) :: nan_policy !! valid NaN handling policy

        raise_nan = trim(adjustl(nan_policy)) == 'raise'
    end function nan_policy_raises

    pure subroutine raise_if_nan_vector(x, weights)
        real(dp), intent(in) :: x(:) !! observations checked for NaNs
        real(dp), intent(in), optional :: weights(:) !! optional weights checked for NaNs

        if (any(ieee_is_nan(x))) then
            error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
        end if
        if (present(weights)) then
            if (any(ieee_is_nan(weights))) then
                error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
            end if
        end if
    end subroutine raise_if_nan_vector

    pure logical function nan_policy_omits(nan_policy) result(omit_nan)
        character(len=*), intent(in) :: nan_policy !! valid NaN handling policy

        omit_nan = trim(adjustl(nan_policy)) == 'omit'
    end function nan_policy_omits

    pure logical function valid_weight_shape(x, weights) result(valid)
        real(dp), intent(in) :: x(:) !! observations used to determine the required weight count
        real(dp), intent(in), optional :: weights(:) !! optional observation weights

        valid = size(x) > 0
        if (present(weights)) valid = valid .and. size(weights) == size(x)
    end function valid_weight_shape

    pure function observation_weight(weights, i) result(weight)
        real(dp), intent(in), optional :: weights(:) !! optional observation weights
        integer, intent(in) :: i !! one-based observation index
        real(dp) :: weight

        if (present(weights)) then
            weight = weights(i)
        else
            weight = 1.0_dp
        end if
    end function observation_weight

    pure subroutine weight_total(x, weights, total, sentinel)
        real(dp), intent(in) :: x(:) !! observations used only for their count
        real(dp), intent(in), optional :: weights(:) !! optional finite observation weights
        real(dp), intent(out) :: total !! sum of the observation weights
        real(dp), intent(out) :: sentinel !! zero on success and NaN on invalid weights

        integer :: i
        real(dp) :: correction
        real(dp) :: term
        real(dp) :: updated
        real(dp) :: weight

        total = 0.0_dp
        correction = 0.0_dp
        sentinel = 0.0_dp
        do i = 1, size(x)
            weight = observation_weight(weights, i)
            if (.not. ieee_is_finite(weight)) then
                sentinel = quiet_nan(0.0_dp)
                return
            end if
            term = weight - correction
            updated = total + term
            correction = (updated - total) - term
            total = updated
        end do
        if (total == 0.0_dp) sentinel = quiet_nan(0.0_dp)
    end subroutine weight_total

    pure function centered_moment(x, order, center) result(value)
        real(dp), intent(in) :: x(:) !! observations; all values must be finite
        integer, intent(in) :: order !! nonnegative moment order
        real(dp), intent(in) :: center !! finite center of the moment
        real(dp) :: value

        integer :: i
        integer :: k
        real(dp) :: correction
        real(dp) :: data_scale
        real(dp) :: scaled_center
        real(dp) :: sum_power
        real(dp) :: term
        real(dp) :: updated

        if (size(x) == 0 .or. order < 0 .or. .not. ieee_is_finite(center)) then
            value = quiet_nan(0.0_dp)
            return
        else if (order == 0) then
            value = 1.0_dp
            return
        end if

        data_scale = abs(center)
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            end if
            data_scale = max(data_scale, abs(x(i)))
        end do
        if (data_scale <= 0.0_dp) then
            value = 0.0_dp
            return
        end if

        scaled_center = center / data_scale
        sum_power = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(x)
            term = (x(i) / data_scale - scaled_center) ** order - correction
            updated = sum_power + term
            correction = (updated - sum_power) - term
            sum_power = updated
        end do
        value = sum_power / real(size(x), dp)
        do k = 1, order
            value = safe_product(value, data_scale)
            if (.not. ieee_is_finite(value)) exit
        end do
    end function centered_moment

    pure function divide_by_scale(numerator, denominator) result(value)
        real(dp), intent(in) :: numerator !! finite or infinite numerator
        real(dp), intent(in) :: denominator !! requested scale divisor
        real(dp) :: value

        if (ieee_is_nan(numerator) .or. ieee_is_nan(denominator)) then
            value = quiet_nan(0.0_dp)
        else if (denominator == 0.0_dp) then
            if (numerator == 0.0_dp) then
                value = quiet_nan(0.0_dp)
            else if (numerator > 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else
                value = negative_infinity(0.0_dp)
            end if
        else
            value = numerator / denominator
        end if
    end function divide_by_scale

    pure function nonnegative_power(x, exponent) result(value)
        real(dp), intent(in) :: x !! finite nonnegative base
        real(dp), intent(in) :: exponent !! finite nonzero exponent
        real(dp) :: value

        real(dp) :: log_value

        if (x == 0.0_dp) then
            if (exponent > 0.0_dp) then
                value = 0.0_dp
            else
                value = positive_infinity(0.0_dp)
            end if
            return
        end if

        log_value = safe_product(exponent, log(x))
        if (.not. ieee_is_finite(log_value)) then
            if (log_value > 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else
                value = 0.0_dp
            end if
        else if (log_value > log(huge(1.0_dp))) then
            value = positive_infinity(0.0_dp)
        else if (log_value < log(tiny(1.0_dp))) then
            value = 0.0_dp
        else
            value = exp(log_value)
        end if
    end function nonnegative_power


    pure function default_matrix_axis(axis) result(axis_value)
        integer, intent(in), optional :: axis !! optional Fortran axis
        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
    end function default_matrix_axis

    pure function retain_matrix_axis(x, axis, reduced) result(values)
        real(dp), intent(in) :: x(:, :) !! original matrix defining retained shape
        integer, intent(in) :: axis !! reduced Fortran axis, 1 or 2
        real(dp), intent(in) :: reduced(:) !! reduced values
        real(dp), allocatable :: values(:, :)

        if (axis == 1 .and. size(reduced) == size(x, 2)) then
            allocate(values(1, size(x, 2)))
            values(1, :) = reduced
        else if (axis == 2 .and. size(reduced) == size(x, 1)) then
            allocate(values(size(x, 1), 1))
            values(:, 1) = reduced
        else
            allocate(values(0, 0))
        end if
    end function retain_matrix_axis

    pure function matrix_weight_vector_valid(x, weights) result(valid)
        real(dp), intent(in) :: x(:, :) !! observations defining target shape
        real(dp), intent(in) :: weights(:) !! trailing one-dimensional weights
        logical :: valid

        valid = size(weights) == 1 .or. size(weights) == size(x, 2)
    end function matrix_weight_vector_valid

    pure function matrix_weight_matrix_valid(x, weights) result(valid)
        real(dp), intent(in) :: x(:, :) !! observations defining target shape
        real(dp), intent(in) :: weights(:, :) !! two-dimensional broadcast weights
        logical :: valid

        valid = (size(weights, 1) == 1 .or. size(weights, 1) == size(x, 1)) .and. &
            (size(weights, 2) == 1 .or. size(weights, 2) == size(x, 2))
    end function matrix_weight_matrix_valid

    pure subroutine matrix_weight_vector_slice(x, axis, slice_index, weights, work_weights)
        real(dp), intent(in) :: x(:, :) !! observations defining broadcast shape
        integer, intent(in) :: axis !! reduced Fortran axis, 1 or 2
        integer, intent(in) :: slice_index !! index along the retained axis
        real(dp), intent(in) :: weights(:) !! valid trailing one-dimensional weights
        real(dp), allocatable, intent(out) :: work_weights(:) !! expanded slice weights

        if (axis == 1) then
            allocate(work_weights(size(x, 1)))
            if (size(weights) == 1) then
                work_weights = weights(1)
            else
                work_weights = weights(slice_index)
            end if
        else
            allocate(work_weights(size(x, 2)))
            if (size(weights) == 1) then
                work_weights = weights(1)
            else
                work_weights = weights
            end if
        end if
    end subroutine matrix_weight_vector_slice

    pure subroutine matrix_weight_matrix_slice(x, axis, slice_index, weights, work_weights)
        real(dp), intent(in) :: x(:, :) !! observations defining broadcast shape
        integer, intent(in) :: axis !! reduced Fortran axis, 1 or 2
        integer, intent(in) :: slice_index !! index along the retained axis
        real(dp), intent(in) :: weights(:, :) !! valid two-dimensional broadcast weights
        real(dp), allocatable, intent(out) :: work_weights(:) !! expanded slice weights

        integer :: i
        integer :: row
        integer :: column

        if (axis == 1) then
            allocate(work_weights(size(x, 1)))
            column = merge(slice_index, 1, size(weights, 2) > 1)
            do i = 1, size(x, 1)
                row = merge(i, 1, size(weights, 1) > 1)
                work_weights(i) = weights(row, column)
            end do
        else
            allocate(work_weights(size(x, 2)))
            row = merge(slice_index, 1, size(weights, 1) > 1)
            do i = 1, size(x, 2)
                column = merge(i, 1, size(weights, 2) > 1)
                work_weights(i) = weights(row, column)
            end do
        end if
    end subroutine matrix_weight_matrix_slice


    pure subroutine broadcast_extent_pair(left, right, extent, valid)
        integer, intent(in) :: left !! first extent
        integer, intent(in) :: right !! second extent
        integer, intent(out) :: extent !! broadcast extent
        logical, intent(out) :: valid !! true when the extents are broadcast compatible

        if (left == right) then
            extent = left
            valid = .true.
        else if (left == 1) then
            extent = right
            valid = .true.
        else if (right == 1) then
            extent = left
            valid = .true.
        else
            extent = 0
            valid = .false.
        end if
    end subroutine broadcast_extent_pair

    pure function broadcast_source_index(index_value, source_extent) result(source_index)
        integer, intent(in) :: index_value !! index in the broadcast result
        integer, intent(in) :: source_extent !! corresponding source extent
        integer :: source_index

        if (source_extent == 1) then
            source_index = 1
        else
            source_index = index_value
        end if
    end function broadcast_source_index

    pure subroutine tensor_broadcast_weights1(x, weights, broadcast_x, broadcast_weights, valid)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: weights(:) !! trailing one-dimensional weights
        real(dp), allocatable, intent(out) :: broadcast_x(:, :, :) !! broadcast observations
        real(dp), allocatable, intent(out) :: broadcast_weights(:, :, :) !! broadcast weights
        logical, intent(out) :: valid !! true when x and weights are broadcast compatible

        integer :: i
        integer :: j
        integer :: k
        integer :: source_i
        integer :: source_j
        integer :: source_k
        integer :: target_shape(3)
        integer :: weight_k
        logical :: extent_valid

        target_shape = shape(x)
        call broadcast_extent_pair(size(x, 3), size(weights), target_shape(3), extent_valid)
        valid = extent_valid
        if (.not. valid) return

        allocate(broadcast_x(target_shape(1), target_shape(2), target_shape(3)))
        allocate(broadcast_weights(target_shape(1), target_shape(2), target_shape(3)))
        do k = 1, target_shape(3)
            source_k = broadcast_source_index(k, size(x, 3))
            weight_k = broadcast_source_index(k, size(weights))
            do j = 1, target_shape(2)
                source_j = broadcast_source_index(j, size(x, 2))
                do i = 1, target_shape(1)
                    source_i = broadcast_source_index(i, size(x, 1))
                    broadcast_x(i, j, k) = x(source_i, source_j, source_k)
                    broadcast_weights(i, j, k) = weights(weight_k)
                end do
            end do
        end do
    end subroutine tensor_broadcast_weights1

    pure subroutine tensor_broadcast_weights2(x, weights, broadcast_x, broadcast_weights, valid)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: weights(:, :) !! trailing two-dimensional weights
        real(dp), allocatable, intent(out) :: broadcast_x(:, :, :) !! broadcast observations
        real(dp), allocatable, intent(out) :: broadcast_weights(:, :, :) !! broadcast weights
        logical, intent(out) :: valid !! true when x and weights are broadcast compatible

        integer :: i
        integer :: j
        integer :: k
        integer :: source_i
        integer :: source_j
        integer :: source_k
        integer :: target_shape(3)
        integer :: weight_j
        integer :: weight_k
        logical :: valid_j
        logical :: valid_k

        target_shape = shape(x)
        call broadcast_extent_pair(size(x, 2), size(weights, 1), target_shape(2), valid_j)
        call broadcast_extent_pair(size(x, 3), size(weights, 2), target_shape(3), valid_k)
        valid = valid_j .and. valid_k
        if (.not. valid) return

        allocate(broadcast_x(target_shape(1), target_shape(2), target_shape(3)))
        allocate(broadcast_weights(target_shape(1), target_shape(2), target_shape(3)))
        do k = 1, target_shape(3)
            source_k = broadcast_source_index(k, size(x, 3))
            weight_k = broadcast_source_index(k, size(weights, 2))
            do j = 1, target_shape(2)
                source_j = broadcast_source_index(j, size(x, 2))
                weight_j = broadcast_source_index(j, size(weights, 1))
                do i = 1, target_shape(1)
                    source_i = broadcast_source_index(i, size(x, 1))
                    broadcast_x(i, j, k) = x(source_i, source_j, source_k)
                    broadcast_weights(i, j, k) = weights(weight_j, weight_k)
                end do
            end do
        end do
    end subroutine tensor_broadcast_weights2

    pure subroutine tensor_broadcast_weights3(x, weights, broadcast_x, broadcast_weights, valid)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: weights(:, :, :) !! three-dimensional weights
        real(dp), allocatable, intent(out) :: broadcast_x(:, :, :) !! broadcast observations
        real(dp), allocatable, intent(out) :: broadcast_weights(:, :, :) !! broadcast weights
        logical, intent(out) :: valid !! true when x and weights are broadcast compatible

        integer :: i
        integer :: j
        integer :: k
        integer :: source_i
        integer :: source_j
        integer :: source_k
        integer :: target_shape(3)
        integer :: weight_i
        integer :: weight_j
        integer :: weight_k
        logical :: valid_i
        logical :: valid_j
        logical :: valid_k

        call broadcast_extent_pair(size(x, 1), size(weights, 1), target_shape(1), valid_i)
        call broadcast_extent_pair(size(x, 2), size(weights, 2), target_shape(2), valid_j)
        call broadcast_extent_pair(size(x, 3), size(weights, 3), target_shape(3), valid_k)
        valid = valid_i .and. valid_j .and. valid_k
        if (.not. valid) return

        allocate(broadcast_x(target_shape(1), target_shape(2), target_shape(3)))
        allocate(broadcast_weights(target_shape(1), target_shape(2), target_shape(3)))
        do k = 1, target_shape(3)
            source_k = broadcast_source_index(k, size(x, 3))
            weight_k = broadcast_source_index(k, size(weights, 3))
            do j = 1, target_shape(2)
                source_j = broadcast_source_index(j, size(x, 2))
                weight_j = broadcast_source_index(j, size(weights, 2))
                do i = 1, target_shape(1)
                    source_i = broadcast_source_index(i, size(x, 1))
                    weight_i = broadcast_source_index(i, size(weights, 1))
                    broadcast_x(i, j, k) = x(source_i, source_j, source_k)
                    broadcast_weights(i, j, k) = weights(weight_i, weight_j, weight_k)
                end do
            end do
        end do
    end subroutine tensor_broadcast_weights3

    pure function tensor_nan_reduced(source_shape, axis) result(values)
        integer, intent(in) :: source_shape(3) !! shape before reduction
        integer, intent(in) :: axis !! reduced Fortran axis, 1 to 3
        real(dp), allocatable :: values(:, :)

        select case (axis)
        case (1)
            allocate(values(source_shape(2), source_shape(3)))
        case (2)
            allocate(values(source_shape(1), source_shape(3)))
        case (3)
            allocate(values(source_shape(1), source_shape(2)))
        case default
            allocate(values(0, 0))
            return
        end select
        values = quiet_nan(0.0_dp)
    end function tensor_nan_reduced

    pure function retain_tensor_shape(source_shape, axis, reduced) result(values)
        integer, intent(in) :: source_shape(3) !! broadcast source shape
        integer, intent(in) :: axis !! reduced Fortran axis, 1 to 3
        real(dp), intent(in) :: reduced(:, :) !! reduced values
        real(dp), allocatable :: values(:, :, :)

        select case (axis)
        case (1)
            if (all(shape(reduced) == [source_shape(2), source_shape(3)])) then
                allocate(values(1, source_shape(2), source_shape(3)))
                values(1, :, :) = reduced
            else
                allocate(values(0, 0, 0))
            end if
        case (2)
            if (all(shape(reduced) == [source_shape(1), source_shape(3)])) then
                allocate(values(source_shape(1), 1, source_shape(3)))
                values(:, 1, :) = reduced
            else
                allocate(values(0, 0, 0))
            end if
        case (3)
            if (all(shape(reduced) == [source_shape(1), source_shape(2)])) then
                allocate(values(source_shape(1), source_shape(2), 1))
                values(:, :, 1) = reduced
            else
                allocate(values(0, 0, 0))
            end if
        case default
            allocate(values(0, 0, 0))
        end select
    end function retain_tensor_shape

    pure function gmean_tensor_weighted_core(x, weights, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations after broadcasting
        real(dp), intent(in) :: weights(:, :, :) !! weights after broadcasting
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in), optional :: nan_policy !! optional propagate/omit/raise policy
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    if (present(nan_policy)) then
                        values(j, k) = gmean_vector_nan_policy( &
                            x(:, j, k), weights(:, j, k), nan_policy)
                    else
                        values(j, k) = gmean_vector(x(:, j, k), weights(:, j, k))
                    end if
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    if (present(nan_policy)) then
                        values(i, k) = gmean_vector_nan_policy( &
                            x(i, :, k), weights(i, :, k), nan_policy)
                    else
                        values(i, k) = gmean_vector(x(i, :, k), weights(i, :, k))
                    end if
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    if (present(nan_policy)) then
                        values(i, j) = gmean_vector_nan_policy( &
                            x(i, j, :), weights(i, j, :), nan_policy)
                    else
                        values(i, j) = gmean_vector(x(i, j, :), weights(i, j, :))
                    end if
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function gmean_tensor_weighted_core

    pure function hmean_tensor_weighted_core(x, weights, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations after broadcasting
        real(dp), intent(in) :: weights(:, :, :) !! weights after broadcasting
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in), optional :: nan_policy !! optional propagate/omit/raise policy
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    if (present(nan_policy)) then
                        values(j, k) = hmean_vector_nan_policy( &
                            x(:, j, k), weights(:, j, k), nan_policy)
                    else
                        values(j, k) = hmean_vector(x(:, j, k), weights(:, j, k))
                    end if
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    if (present(nan_policy)) then
                        values(i, k) = hmean_vector_nan_policy( &
                            x(i, :, k), weights(i, :, k), nan_policy)
                    else
                        values(i, k) = hmean_vector(x(i, :, k), weights(i, :, k))
                    end if
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    if (present(nan_policy)) then
                        values(i, j) = hmean_vector_nan_policy( &
                            x(i, j, :), weights(i, j, :), nan_policy)
                    else
                        values(i, j) = hmean_vector(x(i, j, :), weights(i, j, :))
                    end if
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function hmean_tensor_weighted_core

    pure function pmean_tensor_weighted_core(x, p, weights, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations after broadcasting
        real(dp), intent(in) :: p !! finite power
        real(dp), intent(in) :: weights(:, :, :) !! weights after broadcasting
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in), optional :: nan_policy !! optional propagate/omit/raise policy
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k

        axis_value = default_tensor_axis(axis)
        select case (axis_value)
        case (1)
            allocate(values(size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    if (present(nan_policy)) then
                        values(j, k) = pmean_vector_nan_policy( &
                            x(:, j, k), p, weights(:, j, k), nan_policy)
                    else
                        values(j, k) = pmean_vector(x(:, j, k), p, weights(:, j, k))
                    end if
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3)))
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    if (present(nan_policy)) then
                        values(i, k) = pmean_vector_nan_policy( &
                            x(i, :, k), p, weights(i, :, k), nan_policy)
                    else
                        values(i, k) = pmean_vector(x(i, :, k), p, weights(i, :, k))
                    end if
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    if (present(nan_policy)) then
                        values(i, j) = pmean_vector_nan_policy( &
                            x(i, j, :), p, weights(i, j, :), nan_policy)
                    else
                        values(i, j) = pmean_vector(x(i, j, :), p, weights(i, j, :))
                    end if
                end do
            end do
        case default
            allocate(values(0, 0))
        end select
    end function pmean_tensor_weighted_core

    pure function default_tensor_axis(axis) result(axis_value)
        integer, intent(in), optional :: axis !! optional Fortran axis
        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
    end function default_tensor_axis

    pure function retain_tensor_axis(x, axis, reduced) result(values)
        real(dp), intent(in) :: x(:, :, :) !! original tensor defining retained shape
        integer, intent(in) :: axis !! reduced Fortran axis, 1 to 3
        real(dp), intent(in) :: reduced(:, :) !! reduced values
        real(dp), allocatable :: values(:, :, :)

        select case (axis)
        case (1)
            if (all(shape(reduced) == [size(x, 2), size(x, 3)])) then
                allocate(values(1, size(x, 2), size(x, 3)))
                values(1, :, :) = reduced
            else
                allocate(values(0, 0, 0))
            end if
        case (2)
            if (all(shape(reduced) == [size(x, 1), size(x, 3)])) then
                allocate(values(size(x, 1), 1, size(x, 3)))
                values(:, 1, :) = reduced
            else
                allocate(values(0, 0, 0))
            end if
        case (3)
            if (all(shape(reduced) == [size(x, 1), size(x, 2)])) then
                allocate(values(size(x, 1), size(x, 2), 1))
                values(:, :, 1) = reduced
            else
                allocate(values(0, 0, 0))
            end if
        case default
            allocate(values(0, 0, 0))
        end select
    end function retain_tensor_axis

    pure function safe_product(x, y) result(value)
        real(dp), intent(in) :: x !! finite or infinite first factor
        real(dp), intent(in) :: y !! finite or infinite second factor
        real(dp) :: value

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = quiet_nan(0.0_dp)
        else if (abs(x) <= 0.0_dp .or. abs(y) <= 0.0_dp) then
            value = 0.0_dp
        else if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            if ((x > 0.0_dp) .eqv. (y > 0.0_dp)) then
                value = positive_infinity(0.0_dp)
            else
                value = negative_infinity(0.0_dp)
            end if
        else if (abs(x) <= 1.0_dp .or. abs(y) <= 1.0_dp) then
            value = x * y
        else if (abs(x) > huge(1.0_dp) / abs(y)) then
            if ((x > 0.0_dp) .eqv. (y > 0.0_dp)) then
                value = positive_infinity(0.0_dp)
            else
                value = negative_infinity(0.0_dp)
            end if
        else
            value = x * y
        end if
    end function safe_product

    pure logical function same_value(x, y) result(same)
        real(dp), intent(in) :: x !! first value, guaranteed not NaN
        real(dp), intent(in) :: y !! second value, guaranteed not NaN

        same = x <= y .and. y <= x
    end function same_value

end module scifort_descriptive_lowrank
