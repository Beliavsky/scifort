! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_summary_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_descriptive_lowrank, only : central_moment, kurtosis, mean, skew
    implicit none
    private

    public :: describe
    public :: describe_all
    public :: mode
    public :: mode_all
    public :: mode_keepdims

contains

    pure subroutine mode(x, mode_value, count_value, axis, nan_policy)
        real(dp), intent(in) :: x(..) !! observations, rank 1 through 15
        real(dp), allocatable, intent(out) :: mode_value(..) !! modal values after reduction
        integer, allocatable, intent(out) :: count_value(..) !! counts of modal values
        integer, intent(in), optional :: axis !! one-based Fortran axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise

        integer :: axis_value
        real(dp), allocatable :: flat_modes(:)
        integer, allocatable :: flat_counts(:)
        integer, allocatable :: out_shape(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call mode_reduce_assumed_rank(x, axis_value, nan_policy, flat_modes, flat_counts, &
            out_shape, ok)
        if (.not. ok) error stop 'scifort: invalid mode reduction axis/rank', quiet=.true.
        call assign_real_output(mode_value, flat_modes, out_shape)
        call assign_integer_output(count_value, flat_counts, out_shape)
    end subroutine mode

    pure subroutine mode_keepdims(x, mode_value, count_value, axis, nan_policy)
        real(dp), intent(in) :: x(..) !! observations, rank 1 through 15
        real(dp), allocatable, intent(out) :: mode_value(..) !! modal values with reduced axis retained
        integer, allocatable, intent(out) :: count_value(..) !! counts with reduced axis retained
        integer, intent(in), optional :: axis !! one-based Fortran axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise

        integer :: axis_value
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        real(dp), allocatable :: flat_modes(:)
        integer, allocatable :: flat_counts(:)
        integer, allocatable :: keep_shape(:)
        integer, allocatable :: out_shape(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) error stop 'scifort: mode supports ranks 1-15', quiet=.true.
        if (axis_value < 1 .or. axis_value > size(dims)) then
            error stop 'scifort: invalid mode reduction axis', quiet=.true.
        end if
        call mode_reduce_assumed_rank(x, axis_value, nan_policy, flat_modes, flat_counts, &
            out_shape, ok)
        if (.not. ok) error stop 'scifort: invalid mode reduction axis/rank', quiet=.true.
        keep_shape = dims
        keep_shape(axis_value) = 1
        call assign_real_output(mode_value, flat_modes, keep_shape)
        call assign_integer_output(count_value, flat_counts, keep_shape)
    end subroutine mode_keepdims

    pure subroutine mode_all(x, mode_value, count_value, nan_policy)
        real(dp), intent(in) :: x(..) !! observations, rank 1 through 15
        real(dp), intent(out) :: mode_value !! modal value of the flattened input
        integer, intent(out) :: count_value !! count of the flattened-input mode
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) error stop 'scifort: mode supports ranks 1-15', quiet=.true.
        call mode_vector(flat, nan_policy, mode_value, count_value)
    end subroutine mode_all

    pure subroutine describe(x, nobs, minimum_value, maximum_value, mean_value, &
            variance_value, skewness_value, kurtosis_value, axis, ddof, bias, nan_policy)
        real(dp), intent(in) :: x(..) !! observations, rank 1 through 15
        integer, allocatable, intent(out) :: nobs(..) !! observations per reduced slice
        real(dp), allocatable, intent(out) :: minimum_value(..) !! slice minima
        real(dp), allocatable, intent(out) :: maximum_value(..) !! slice maxima
        real(dp), allocatable, intent(out) :: mean_value(..) !! slice means
        real(dp), allocatable, intent(out) :: variance_value(..) !! slice variances
        real(dp), allocatable, intent(out) :: skewness_value(..) !! slice skewness
        real(dp), allocatable, intent(out) :: kurtosis_value(..) !! Fisher kurtosis
        integer, intent(in), optional :: axis !! one-based Fortran axis, default 1
        integer, intent(in), optional :: ddof !! variance degrees-of-freedom adjustment, default 1
        logical, intent(in), optional :: bias !! use biased skew/kurtosis, default true
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise

        integer :: axis_value
        integer, allocatable :: flat_nobs(:)
        real(dp), allocatable :: flat_minimum(:)
        real(dp), allocatable :: flat_maximum(:)
        real(dp), allocatable :: flat_mean(:)
        real(dp), allocatable :: flat_variance(:)
        real(dp), allocatable :: flat_skewness(:)
        real(dp), allocatable :: flat_kurtosis(:)
        integer, allocatable :: out_shape(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call describe_reduce_assumed_rank(x, axis_value, ddof, bias, nan_policy, flat_nobs, &
            flat_minimum, flat_maximum, flat_mean, flat_variance, flat_skewness, &
            flat_kurtosis, out_shape, ok)
        if (.not. ok) error stop 'scifort: invalid describe reduction axis/rank', quiet=.true.
        call assign_integer_output(nobs, flat_nobs, out_shape)
        call assign_real_output(minimum_value, flat_minimum, out_shape)
        call assign_real_output(maximum_value, flat_maximum, out_shape)
        call assign_real_output(mean_value, flat_mean, out_shape)
        call assign_real_output(variance_value, flat_variance, out_shape)
        call assign_real_output(skewness_value, flat_skewness, out_shape)
        call assign_real_output(kurtosis_value, flat_kurtosis, out_shape)
    end subroutine describe

    pure subroutine describe_all(x, nobs, minimum_value, maximum_value, mean_value, &
            variance_value, skewness_value, kurtosis_value, ddof, bias, nan_policy)
        real(dp), intent(in) :: x(..) !! observations, rank 1 through 15
        integer, intent(out) :: nobs !! number of retained observations
        real(dp), intent(out) :: minimum_value !! minimum
        real(dp), intent(out) :: maximum_value !! maximum
        real(dp), intent(out) :: mean_value !! arithmetic mean
        real(dp), intent(out) :: variance_value !! variance
        real(dp), intent(out) :: skewness_value !! skewness
        real(dp), intent(out) :: kurtosis_value !! Fisher kurtosis
        integer, intent(in), optional :: ddof !! variance degrees-of-freedom adjustment, default 1
        logical, intent(in), optional :: bias !! use biased skew/kurtosis, default true
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        logical :: mask_invalid
        character(len=16) :: policy

        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) error stop 'scifort: describe supports ranks 1-15', quiet=.true.
        if (size(flat) == 0) error stop 'scifort: describe input must not be empty', quiet=.true.
        policy = 'propagate'
        if (present(nan_policy)) policy = lowercase(trim(adjustl(nan_policy)))
        if (.not. valid_nan_policy(policy)) error stop 'scifort: invalid nan_policy', quiet=.true.
        if (policy == 'raise' .and. any(ieee_is_nan(flat))) then
            error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
        end if
        mask_invalid = policy == 'omit' .and. any(ieee_is_nan(flat))
        call describe_vector(flat, ddof, bias, policy, mask_invalid, nobs, minimum_value, &
            maximum_value, mean_value, variance_value, skewness_value, kurtosis_value)
    end subroutine describe_all

    pure subroutine mode_reduce_assumed_rank(x, axis, nan_policy, modes, counts, out_shape, ok)
        real(dp), intent(in) :: x(..) !! observations
        integer, intent(in) :: axis !! one-based reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable, intent(out) :: modes(:) !! flattened modes
        integer, allocatable, intent(out) :: counts(:) !! flattened counts
        integer, allocatable, intent(out) :: out_shape(:) !! result extents
        logical, intent(out) :: ok !! valid rank and axis

        integer :: base_index, coordinate, d, i, linear, n_out, n_reduce, out_dim
        integer, allocatable :: dims(:), strides(:)
        real(dp), allocatable :: flat(:), slice(:)

        ok = .false.
        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(modes(0), counts(0), out_shape(0))
            return
        end if
        if (axis < 1 .or. axis > size(dims)) then
            allocate(modes(0), counts(0), out_shape(0))
            return
        end if
        allocate(out_shape(size(dims) - 1))
        out_dim = 0
        do d = 1, size(dims)
            if (d == axis) cycle
            out_dim = out_dim + 1
            out_shape(out_dim) = dims(d)
        end do
        allocate(strides(size(dims)))
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
        end do
        n_reduce = dims(axis)
        n_out = product(out_shape)
        allocate(modes(n_out), counts(n_out), slice(n_reduce))
        do i = 0, n_out - 1
            linear = i
            base_index = 1
            out_dim = 0
            do d = 1, size(dims)
                if (d == axis) cycle
                out_dim = out_dim + 1
                coordinate = mod(linear, out_shape(out_dim))
                linear = linear / out_shape(out_dim)
                base_index = base_index + coordinate * strides(d)
            end do
            do d = 1, n_reduce
                slice(d) = flat(base_index + (d - 1) * strides(axis))
            end do
            call mode_vector(slice, nan_policy, modes(i + 1), counts(i + 1))
        end do
        ok = .true.
    end subroutine mode_reduce_assumed_rank

    pure subroutine describe_reduce_assumed_rank(x, axis, ddof, bias, nan_policy, nobs, &
            minima, maxima, means, variances, skewness, kurtoses, out_shape, ok)
        real(dp), intent(in) :: x(..) !! observations
        integer, intent(in) :: axis !! one-based reduction axis
        integer, intent(in), optional :: ddof !! variance adjustment
        logical, intent(in), optional :: bias !! biased skew/kurtosis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, allocatable, intent(out) :: nobs(:) !! flattened counts
        real(dp), allocatable, intent(out) :: minima(:), maxima(:), means(:), variances(:)
        real(dp), allocatable, intent(out) :: skewness(:), kurtoses(:)
        integer, allocatable, intent(out) :: out_shape(:) !! result extents
        logical, intent(out) :: ok !! valid rank and axis

        integer :: base_index, coordinate, d, i, linear, n_out, n_reduce, out_dim
        integer, allocatable :: dims(:), strides(:)
        real(dp), allocatable :: flat(:), slice(:)
        logical :: mask_invalid
        character(len=16) :: policy

        ok = .false.
        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(nobs(0), minima(0), maxima(0), means(0), variances(0), skewness(0), &
                kurtoses(0), out_shape(0))
            return
        end if
        if (size(flat) == 0) error stop 'scifort: describe input must not be empty', quiet=.true.
        if (axis < 1 .or. axis > size(dims)) then
            allocate(nobs(0), minima(0), maxima(0), means(0), variances(0), skewness(0), &
                kurtoses(0), out_shape(0))
            return
        end if
        policy = 'propagate'
        if (present(nan_policy)) policy = lowercase(trim(adjustl(nan_policy)))
        if (.not. valid_nan_policy(policy)) error stop 'scifort: invalid nan_policy', quiet=.true.
        if (policy == 'raise' .and. any(ieee_is_nan(flat))) then
            error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
        end if
        mask_invalid = policy == 'omit' .and. any(ieee_is_nan(flat))
        allocate(out_shape(size(dims) - 1))
        out_dim = 0
        do d = 1, size(dims)
            if (d == axis) cycle
            out_dim = out_dim + 1
            out_shape(out_dim) = dims(d)
        end do
        allocate(strides(size(dims)))
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
        end do
        n_reduce = dims(axis)
        n_out = product(out_shape)
        allocate(nobs(n_out), minima(n_out), maxima(n_out), means(n_out), variances(n_out), &
            skewness(n_out), kurtoses(n_out), slice(n_reduce))
        do i = 0, n_out - 1
            linear = i
            base_index = 1
            out_dim = 0
            do d = 1, size(dims)
                if (d == axis) cycle
                out_dim = out_dim + 1
                coordinate = mod(linear, out_shape(out_dim))
                linear = linear / out_shape(out_dim)
                base_index = base_index + coordinate * strides(d)
            end do
            do d = 1, n_reduce
                slice(d) = flat(base_index + (d - 1) * strides(axis))
            end do
            call describe_vector(slice, ddof, bias, policy, mask_invalid, nobs(i + 1), &
                minima(i + 1), maxima(i + 1), means(i + 1), variances(i + 1), &
                skewness(i + 1), kurtoses(i + 1))
        end do
        ok = .true.
    end subroutine describe_reduce_assumed_rank

    pure subroutine mode_vector(x, nan_policy, mode_value, count_value)
        real(dp), intent(in) :: x(:) !! one reduction slice
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), intent(out) :: mode_value !! modal value
        integer, intent(out) :: count_value !! modal count

        real(dp), allocatable :: clean(:)
        integer :: current_count, i, run_start
        character(len=16) :: policy

        policy = 'propagate'
        if (present(nan_policy)) policy = lowercase(trim(adjustl(nan_policy)))
        if (.not. valid_nan_policy(policy)) error stop 'scifort: invalid nan_policy', quiet=.true.
        if (policy == 'raise') then
            if (any(ieee_is_nan(x))) then
                error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
            end if
            clean = x
        else if (policy == 'omit') then
            clean = pack(x, .not. ieee_is_nan(x))
        else
            clean = x
        end if
        if (size(clean) == 0) then
            mode_value = quiet_nan(0.0_dp)
            count_value = 0
            return
        end if
        call sort_values_nan_last(clean)
        mode_value = clean(1)
        count_value = 1
        run_start = 1
        do i = 2, size(clean) + 1
            if (i <= size(clean)) then
                if (same_mode_value(clean(i), clean(run_start))) cycle
            end if
            current_count = i - run_start
            if (current_count > count_value) then
                count_value = current_count
                mode_value = clean(run_start)
            end if
            run_start = i
        end do
    end subroutine mode_vector

    pure subroutine describe_vector(x, ddof, bias, policy, mask_invalid, nobs, minimum_value, &
            maximum_value, mean_value, variance_value, skewness_value, kurtosis_value)
        real(dp), intent(in) :: x(:) !! one reduction slice
        integer, intent(in), optional :: ddof !! variance adjustment
        logical, intent(in), optional :: bias !! biased skew/kurtosis
        character(len=*), intent(in) :: policy !! validated nan policy
        logical, intent(in) :: mask_invalid !! SciPy omit path masks NaN and infinities
        integer, intent(out) :: nobs !! retained observations
        real(dp), intent(out) :: minimum_value, maximum_value, mean_value, variance_value
        real(dp), intent(out) :: skewness_value, kurtosis_value

        real(dp), allocatable :: clean(:)
        logical :: use_bias

        use_bias = .true.
        if (present(bias)) use_bias = bias
        if (policy == 'propagate' .and. any(ieee_is_nan(x))) then
            nobs = size(x)
            minimum_value = quiet_nan(0.0_dp)
            maximum_value = quiet_nan(0.0_dp)
            mean_value = quiet_nan(0.0_dp)
            variance_value = quiet_nan(0.0_dp)
            skewness_value = quiet_nan(0.0_dp)
            kurtosis_value = quiet_nan(0.0_dp)
            return
        end if
        if (mask_invalid) then
            clean = pack(x, ieee_is_finite(x))
        else
            clean = x
        end if
        nobs = size(clean)
        if (nobs == 0) then
            minimum_value = quiet_nan(0.0_dp)
            maximum_value = quiet_nan(0.0_dp)
            mean_value = quiet_nan(0.0_dp)
            variance_value = quiet_nan(0.0_dp)
            skewness_value = quiet_nan(0.0_dp)
            kurtosis_value = quiet_nan(0.0_dp)
            return
        end if
        minimum_value = minval(clean)
        maximum_value = maxval(clean)
        mean_value = mean(clean)
        variance_value = describe_variance(clean, ddof)
        if (mask_invalid) then
            skewness_value = masked_skew(clean, use_bias)
            kurtosis_value = masked_kurtosis(clean, use_bias)
        else
            skewness_value = skew(clean, bias=use_bias)
            kurtosis_value = kurtosis(clean, fisher=.true., bias=use_bias)
        end if
    end subroutine describe_vector

    pure function describe_variance(x, ddof) result(value)
        real(dp), intent(in) :: x(:) !! nonempty observations
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        real(dp) :: value

        integer :: d
        real(dp) :: m2
        real(dp) :: factor

        d = 1
        if (present(ddof)) d = ddof
        m2 = central_moment(x, 2)
        if (ieee_is_nan(m2)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (size(x) == d) then
            if (m2 > 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else
                value = quiet_nan(0.0_dp)
            end if
            return
        end if
        factor = real(size(x), dp) / real(size(x) - d, dp)
        value = m2 * factor
    end function describe_variance

    pure function masked_skew(x, bias) result(value)
        real(dp), intent(in) :: x(:) !! finite retained observations
        logical, intent(in) :: bias !! use biased estimator
        real(dp) :: value

        real(dp) :: m2, m3, mu, n, resolution
        logical :: effectively_zero

        if (size(x) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        mu = mean(x)
        m2 = central_moment(x, 2)
        m3 = central_moment(x, 3)
        resolution = 10.0_dp ** (-precision(1.0_dp))
        effectively_zero = m2 <= (resolution * mu) ** 2
        if (effectively_zero) then
            value = 0.0_dp
            return
        end if
        value = m3 / (m2 * sqrt(m2))
        if (.not. bias .and. size(x) > 2) then
            n = real(size(x), dp)
            value = sqrt(n * (n - 1.0_dp)) / (n - 2.0_dp) * value
        end if
    end function masked_skew

    pure function masked_kurtosis(x, bias) result(value)
        real(dp), intent(in) :: x(:) !! finite retained observations
        logical, intent(in) :: bias !! use biased estimator
        real(dp) :: value

        real(dp) :: m2, m4, mu, n, pearson, resolution
        logical :: effectively_zero

        if (size(x) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        mu = mean(x)
        m2 = central_moment(x, 2)
        m4 = central_moment(x, 4)
        resolution = 10.0_dp ** (-precision(1.0_dp))
        effectively_zero = m2 <= (resolution * mu) ** 2
        if (effectively_zero) then
            value = -3.0_dp
            return
        end if
        pearson = m4 / (m2 * m2)
        if (.not. bias .and. size(x) > 3) then
            n = real(size(x), dp)
            pearson = ((n * n - 1.0_dp) * pearson - 3.0_dp * (n - 1.0_dp) ** 2) / &
                ((n - 2.0_dp) * (n - 3.0_dp)) + 3.0_dp
        end if
        value = pearson - 3.0_dp
    end function masked_kurtosis

    pure subroutine sort_values_nan_last(values)
        real(dp), intent(inout) :: values(:) !! ascending values with NaNs last
        real(dp), allocatable :: temporary(:)
        if (size(values) <= 1) return
        allocate(temporary(size(values)))
        call merge_sort_values_nan_last(values, temporary, 1, size(values))
    end subroutine sort_values_nan_last

    recursive pure subroutine merge_sort_values_nan_last(values, temporary, left, right)
        real(dp), intent(inout) :: values(:) !! values being sorted
        real(dp), intent(inout) :: temporary(:) !! merge workspace
        integer, intent(in) :: left !! first index
        integer, intent(in) :: right !! last index
        integer :: i, j, k, middle
        if (left >= right) return
        middle = left + (right - left) / 2
        call merge_sort_values_nan_last(values, temporary, left, middle)
        call merge_sort_values_nan_last(values, temporary, middle + 1, right)
        i = left
        j = middle + 1
        do k = left, right
            if (i > middle) then
                temporary(k) = values(j); j = j + 1
            else if (j > right) then
                temporary(k) = values(i); i = i + 1
            else if (comes_before(values(i), values(j))) then
                temporary(k) = values(i); i = i + 1
            else
                temporary(k) = values(j); j = j + 1
            end if
        end do
        values(left:right) = temporary(left:right)
    end subroutine merge_sort_values_nan_last

    pure logical function comes_before(left, right) result(before)
        real(dp), intent(in) :: left !! left comparison value
        real(dp), intent(in) :: right !! right comparison value
        if (ieee_is_nan(left)) then
            before = ieee_is_nan(right)
        else if (ieee_is_nan(right)) then
            before = .true.
        else
            before = left <= right
        end if
    end function comes_before

    pure logical function same_mode_value(left, right) result(same)
        real(dp), intent(in) :: left !! first value
        real(dp), intent(in) :: right !! second value
        if (ieee_is_nan(left)) then
            same = ieee_is_nan(right)
        else if (ieee_is_nan(right)) then
            same = .false.
        else
            same = left <= right .and. right <= left
        end if
    end function same_mode_value

    pure subroutine assign_real_output(output, flat, dims)
        real(dp), allocatable, intent(out) :: output(..) !! caller-provided result rank
        real(dp), intent(in) :: flat(:) !! flattened result
        integer, intent(in) :: dims(:) !! desired result shape
        select rank (output)
        rank (0)
            if (size(dims) /= 0 .or. size(flat) /= 1) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output)
            output = flat(1)
        rank (1)
            if (size(dims) /= 1) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1)))
            output = reshape(flat, shape(output))
        rank (2)
            if (size(dims) /= 2) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2)))
            output = reshape(flat, shape(output))
        rank (3)
            if (size(dims) /= 3) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3)))
            output = reshape(flat, shape(output))
        rank (4)
            if (size(dims) /= 4) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4)))
            output = reshape(flat, shape(output))
        rank (5)
            if (size(dims) /= 5) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5)))
            output = reshape(flat, shape(output))
        rank (6)
            if (size(dims) /= 6) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
            output = reshape(flat, shape(output))
        rank (7)
            if (size(dims) /= 7) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
            output = reshape(flat, shape(output))
        rank (8)
            if (size(dims) /= 8) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8)))
            output = reshape(flat, shape(output))
        rank (9)
            if (size(dims) /= 9) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9)))
            output = reshape(flat, shape(output))
        rank (10)
            if (size(dims) /= 10) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10)))
            output = reshape(flat, shape(output))
        rank (11)
            if (size(dims) /= 11) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), dims(11)))
            output = reshape(flat, shape(output))
        rank (12)
            if (size(dims) /= 12) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), &
                dims(11), dims(12)))
            output = reshape(flat, shape(output))
        rank (13)
            if (size(dims) /= 13) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), &
                dims(11), dims(12), dims(13)))
            output = reshape(flat, shape(output))
        rank (14)
            if (size(dims) /= 14) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), &
                dims(11), dims(12), dims(13), dims(14)))
            output = reshape(flat, shape(output))
        rank (15)
            if (size(dims) /= 15) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), &
                dims(11), dims(12), dims(13), dims(14), dims(15)))
            output = reshape(flat, shape(output))
        rank default
            error stop 'scifort: output rank exceeds 15', quiet=.true.
        end select
    end subroutine assign_real_output

    pure subroutine assign_integer_output(output, flat, dims)
        integer, allocatable, intent(out) :: output(..) !! caller-provided result rank
        integer, intent(in) :: flat(:) !! flattened result
        integer, intent(in) :: dims(:) !! desired result shape
        select rank (output)
        rank (0)
            if (size(dims) /= 0 .or. size(flat) /= 1) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output)
            output = flat(1)
        rank (1)
            if (size(dims) /= 1) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1)))
            output = reshape(flat, shape(output))
        rank (2)
            if (size(dims) /= 2) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2)))
            output = reshape(flat, shape(output))
        rank (3)
            if (size(dims) /= 3) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3)))
            output = reshape(flat, shape(output))
        rank (4)
            if (size(dims) /= 4) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4)))
            output = reshape(flat, shape(output))
        rank (5)
            if (size(dims) /= 5) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5)))
            output = reshape(flat, shape(output))
        rank (6)
            if (size(dims) /= 6) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
            output = reshape(flat, shape(output))
        rank (7)
            if (size(dims) /= 7) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
            output = reshape(flat, shape(output))
        rank (8)
            if (size(dims) /= 8) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8)))
            output = reshape(flat, shape(output))
        rank (9)
            if (size(dims) /= 9) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9)))
            output = reshape(flat, shape(output))
        rank (10)
            if (size(dims) /= 10) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10)))
            output = reshape(flat, shape(output))
        rank (11)
            if (size(dims) /= 11) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), dims(11)))
            output = reshape(flat, shape(output))
        rank (12)
            if (size(dims) /= 12) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), &
                dims(11), dims(12)))
            output = reshape(flat, shape(output))
        rank (13)
            if (size(dims) /= 13) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), &
                dims(11), dims(12), dims(13)))
            output = reshape(flat, shape(output))
        rank (14)
            if (size(dims) /= 14) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), &
                dims(11), dims(12), dims(13), dims(14)))
            output = reshape(flat, shape(output))
        rank (15)
            if (size(dims) /= 15) error stop 'scifort: output rank mismatch', quiet=.true.
            allocate(output(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7), dims(8), dims(9), dims(10), &
                dims(11), dims(12), dims(13), dims(14), dims(15)))
            output = reshape(flat, shape(output))
        rank default
            error stop 'scifort: output rank exceeds 15', quiet=.true.
        end select
    end subroutine assign_integer_output

    pure subroutine flatten_assumed_real(x, flat, dims)
        real(dp), intent(in) :: x(..) !! rank-1 through rank-15 real array
        real(dp), allocatable, intent(out) :: flat(:) !! Fortran-order flattened values
        integer, allocatable, intent(out) :: dims(:) !! input shape
        select rank (a => x)
        rank (1)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (2)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (3)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (4)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (5)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (6)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (7)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (8)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (9)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (10)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (11)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (12)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (13)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (14)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (15)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank default
            allocate(dims(0), flat(0))
        end select
    end subroutine flatten_assumed_real

    pure logical function valid_nan_policy(policy) result(valid)
        character(len=*), intent(in) :: policy !! requested policy
        valid = policy == 'propagate' .or. policy == 'omit' .or. policy == 'raise'
    end function valid_nan_policy

    pure function lowercase(text) result(lower)
        character(len=*), intent(in) :: text !! input text
        character(len=len(text)) :: lower
        integer :: code, i
        lower = text
        do i = 1, len(text)
            code = iachar(lower(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

end module scifort_descriptive_summary_engine
