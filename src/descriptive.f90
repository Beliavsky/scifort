! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, negative_infinity, quiet_nan
    implicit none
    private

    interface quantile
        module procedure quantile_scalar
        module procedure quantile_array
    end interface quantile

    public :: central_moment
    public :: covariance
    public :: mean
    public :: median
    public :: pearson_correlation
    public :: quantile
    public :: rankdata
    public :: standard_deviation
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

end module scifort_descriptive
