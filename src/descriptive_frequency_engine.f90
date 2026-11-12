! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_frequency_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    type, public :: cumfreq_result
        real(dp), allocatable :: cumcount(:) !! cumulative histogram bin counts/weights
        real(dp) :: lowerlimit = 0.0_dp !! reported lower real limit
        real(dp) :: binsize = 0.0_dp !! histogram bin width
        integer :: extrapoints = 0 !! observations outside reported limits
    end type cumfreq_result

    type, public :: relfreq_result
        real(dp), allocatable :: frequency(:) !! histogram frequency divided by size(a,1)
        real(dp) :: lowerlimit = 0.0_dp !! reported lower real limit
        real(dp) :: binsize = 0.0_dp !! histogram bin width
        integer :: extrapoints = 0 !! observations outside reported limits
    end type relfreq_result

    public :: cumfreq
    public :: percentileofscore_scalar
    public :: percentileofscore_vector
    public :: relfreq
    public :: scoreatpercentile_all_assumed_rank
    public :: scoreatpercentile_reduce_assumed_rank
    public :: scoreatpercentile_vector
    public :: scoreatpercentile_vector_many

contains

    pure function percentileofscore_scalar(a, score, kind, nan_policy) result(value)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score !! score whose percentile rank is requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: clean(:)
        integer :: left_count
        integer :: matches
        integer :: n
        character(len=16) :: policy
        character(len=16) :: rank_kind
        integer :: right_count

        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))
        if (.not. valid_nan_policy(policy)) then
            error stop 'scifort: invalid nan_policy', quiet=.true.
        end if

        rank_kind = 'rank'
        if (present(kind)) rank_kind = trim(adjustl(kind))
        if (.not. valid_percentile_kind(rank_kind)) then
            error stop "scifort: kind must be 'rank', 'strict', 'weak', or 'mean'", quiet=.true.
        end if

        if (ieee_is_nan(score)) then
            value = quiet_nan(0.0_dp)
            return
        end if

        if (policy == 'raise') then
            if (any(ieee_is_nan(a))) then
                error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
            end if
            clean = a
        else if (policy == 'propagate') then
            if (any(ieee_is_nan(a))) then
                value = quiet_nan(0.0_dp)
                return
            end if
            clean = a
        else
            clean = pack(a, .not. ieee_is_nan(a))
        end if

        n = size(clean)
        if (n == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if

        left_count = count(clean < score)
        right_count = count(clean <= score)
        matches = right_count - left_count
        select case (rank_kind)
        case ('strict')
            value = 100.0_dp * real(left_count, dp) / real(n, dp)
        case ('weak')
            value = 100.0_dp * real(right_count, dp) / real(n, dp)
        case ('mean')
            value = 50.0_dp * real(left_count + right_count, dp) / real(n, dp)
        case default
            value = 50.0_dp * real(left_count + right_count + merge(1, 0, matches > 0), dp) / &
                real(n, dp)
        end select
    end function percentileofscore_scalar

    pure function percentileofscore_vector(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        real(dp), allocatable :: clean(:)
        integer :: i
        integer :: left_count
        integer :: matches
        integer :: n
        character(len=16) :: policy
        character(len=16) :: rank_kind
        integer :: right_count

        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))
        if (.not. valid_nan_policy(policy)) then
            error stop 'scifort: invalid nan_policy', quiet=.true.
        end if
        rank_kind = 'rank'
        if (present(kind)) rank_kind = trim(adjustl(kind))
        if (.not. valid_percentile_kind(rank_kind)) then
            error stop "scifort: kind must be 'rank', 'strict', 'weak', or 'mean'", quiet=.true.
        end if

        allocate(values(size(score)))
        if (policy == 'raise') then
            if (any(ieee_is_nan(a))) then
                error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
            end if
            clean = a
        else if (policy == 'propagate') then
            if (any(ieee_is_nan(a))) then
                values = quiet_nan(0.0_dp)
                return
            end if
            clean = a
        else
            clean = pack(a, .not. ieee_is_nan(a))
        end if

        n = size(clean)
        if (n == 0) then
            values = quiet_nan(0.0_dp)
            return
        end if
        do i = 1, size(score)
            if (ieee_is_nan(score(i))) then
                values(i) = quiet_nan(0.0_dp)
                cycle
            end if
            left_count = count(clean < score(i))
            right_count = count(clean <= score(i))
            matches = right_count - left_count
            select case (rank_kind)
            case ('strict')
                values(i) = 100.0_dp * real(left_count, dp) / real(n, dp)
            case ('weak')
                values(i) = 100.0_dp * real(right_count, dp) / real(n, dp)
            case ('mean')
                values(i) = 50.0_dp * real(left_count + right_count, dp) / real(n, dp)
            case default
                values(i) = 50.0_dp * &
                    real(left_count + right_count + merge(1, 0, matches > 0), dp) / real(n, dp)
            end select
        end do
    end function percentileofscore_vector

    pure function scoreatpercentile_vector(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:) !! observations in one reduction slice
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: selected(:)
        character(len=16) :: method
        integer :: n
        real(dp) :: q
        real(dp), allocatable :: sorted(:)

        if (ieee_is_nan(per) .or. per < 0.0_dp .or. per > 100.0_dp) then
            error stop 'scifort: percentile must be in the range [0, 100]', quiet=.true.
        end if
        method = 'fraction'
        if (present(interpolation_method)) method = trim(adjustl(interpolation_method))
        if (.not. valid_interpolation_method(method)) then
            error stop "scifort: interpolation_method must be 'fraction', 'lower', or 'higher'", quiet=.true.
        end if

        call select_score_values(x, selected, limit)
        n = size(selected)
        if (n == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        sorted = selected
        call sort_values_nan_last(sorted)
        q = per / 100.0_dp
        value = sorted_percentile(sorted, q, method)
    end function scoreatpercentile_vector

    pure function scoreatpercentile_vector_many(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:) !! observations in one reduction slice
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        real(dp), allocatable :: selected(:)
        character(len=16) :: method
        integer :: i
        real(dp), allocatable :: sorted(:)

        method = 'fraction'
        if (present(interpolation_method)) method = trim(adjustl(interpolation_method))
        if (.not. valid_interpolation_method(method)) then
            error stop "scifort: interpolation_method must be 'fraction', 'lower', or 'higher'", quiet=.true.
        end if
        do i = 1, size(per)
            if (ieee_is_nan(per(i)) .or. per(i) < 0.0_dp .or. per(i) > 100.0_dp) then
                error stop 'scifort: percentile must be in the range [0, 100]', quiet=.true.
            end if
        end do

        allocate(values(size(per)))
        call select_score_values(x, selected, limit)
        if (size(selected) == 0) then
            values = quiet_nan(0.0_dp)
            return
        end if
        sorted = selected
        call sort_values_nan_last(sorted)
        do i = 1, size(per)
            values(i) = sorted_percentile(sorted, per(i) / 100.0_dp, method)
        end do
    end function scoreatpercentile_vector_many

    pure function scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        if (rank(x) < 1 .or. rank(x) > 15) then
            allocate(values(size(per)))
            values = quiet_nan(0.0_dp)
            return
        end if
        call flatten_assumed_real(x, flat, dims)
        values = scoreatpercentile_vector_many(flat, per, limit, interpolation_method)
    end function scoreatpercentile_all_assumed_rank

    pure subroutine scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, &
            out_shape, ok, limit, interpolation_method)
        real(dp), intent(in) :: x(..) !! observations with rank 2 to 15
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), allocatable, intent(out) :: reduced(:) !! percentile-first flattened result
        integer, intent(out) :: out_shape(:) !! shape with reduction axis removed
        logical, intent(out) :: ok !! true when rank and axis are valid
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher

        integer :: base_index
        integer :: coordinate
        integer :: d
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        integer :: g
        integer :: j
        integer :: linear
        integer :: n_out
        integer :: n_reduce
        integer :: out_dim
        real(dp), allocatable :: slice(:)
        integer, allocatable :: strides(:)
        real(dp), allocatable :: slice_values(:)

        ok = .false.
        if (rank(x) < 2 .or. rank(x) > 15) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if
        if (size(out_shape) /= rank(x) - 1) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if
        if (axis < 1 .or. axis > rank(x)) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if

        call flatten_assumed_real(x, flat, dims)
        allocate(strides(size(dims)))
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
        end do

        out_dim = 0
        do d = 1, size(dims)
            if (d == axis) cycle
            out_dim = out_dim + 1
            out_shape(out_dim) = dims(d)
        end do
        n_reduce = dims(axis)
        n_out = product(out_shape)
        allocate(reduced(size(per) * n_out))
        allocate(slice(n_reduce))

        do g = 0, n_out - 1
            linear = g
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
            slice_values = scoreatpercentile_vector_many(slice, per, limit, interpolation_method)
            do j = 1, size(per)
                reduced(g * size(per) + j) = slice_values(j)
            end do
        end do
        ok = .true.
    end subroutine scoreatpercentile_reduce_assumed_rank

    pure function cumfreq(a, numbins, defaultreallimits, weights) result(result_value)
        real(dp), intent(in) :: a(..) !! input observations, flattened for histogramming
        integer, intent(in), optional :: numbins !! number of bins, default 10
        real(dp), intent(in), optional :: defaultreallimits(2) !! reported lower/upper real limits
        real(dp), intent(in), optional :: weights(..) !! weights with the same shape as a
        type(cumfreq_result) :: result_value

        real(dp), allocatable :: histogram(:)
        integer :: i
        real(dp) :: running
        real(dp) :: lower
        real(dp) :: binsize
        integer :: extrapoints

        call frequency_histogram(a, numbins, defaultreallimits, weights, histogram, &
            lower, binsize, extrapoints)
        allocate(result_value%cumcount(size(histogram)))
        running = 0.0_dp
        do i = 1, size(histogram)
            running = running + histogram(i)
            result_value%cumcount(i) = running
        end do
        result_value%lowerlimit = lower
        result_value%binsize = binsize
        result_value%extrapoints = extrapoints
    end function cumfreq

    pure function relfreq(a, numbins, defaultreallimits, weights) result(result_value)
        real(dp), intent(in) :: a(..) !! input observations, flattened for histogramming
        integer, intent(in), optional :: numbins !! number of bins, default 10
        real(dp), intent(in), optional :: defaultreallimits(2) !! reported lower/upper real limits
        real(dp), intent(in), optional :: weights(..) !! weights with the same shape as a
        type(relfreq_result) :: result_value

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        real(dp), allocatable :: histogram(:)
        real(dp) :: lower
        real(dp) :: binsize
        integer :: extrapoints
        real(dp) :: denominator

        call flatten_assumed_real(a, flat, dims)
        call frequency_histogram(a, numbins, defaultreallimits, weights, histogram, &
            lower, binsize, extrapoints)
        allocate(result_value%frequency(size(histogram)))
        if (size(dims) == 0) then
            denominator = 0.0_dp
        else
            denominator = real(dims(1), dp)
        end if
        if (denominator <= 0.0_dp) then
            result_value%frequency = quiet_nan(0.0_dp)
        else
            result_value%frequency = histogram / denominator
        end if
        result_value%lowerlimit = lower
        result_value%binsize = binsize
        result_value%extrapoints = extrapoints
    end function relfreq

    pure subroutine frequency_histogram(a, numbins, defaultreallimits, weights, histogram, &
            lowerlimit, binsize, extrapoints)
        real(dp), intent(in) :: a(..) !! observations to bin
        integer, intent(in), optional :: numbins !! number of bins, default 10
        real(dp), intent(in), optional :: defaultreallimits(2) !! requested real limits
        real(dp), intent(in), optional :: weights(..) !! weights with same shape as a
        real(dp), allocatable, intent(out) :: histogram(:) !! weighted bin counts
        real(dp), intent(out) :: lowerlimit !! reported lower real limit
        real(dp), intent(out) :: binsize !! actual histogram bin width
        integer, intent(out) :: extrapoints !! points outside reported limits

        real(dp) :: data_max
        real(dp) :: data_min
        integer, allocatable :: dims(:)
        integer, allocatable :: weight_dims(:)
        real(dp), allocatable :: flat(:)
        real(dp), allocatable :: flat_weights(:)
        real(dp) :: hist_lower
        real(dp) :: hist_upper
        integer :: i
        integer :: ibin
        real(dp) :: lower
        integer :: nbins
        real(dp) :: s
        real(dp) :: upper
        real(dp) :: w

        nbins = 10
        if (present(numbins)) nbins = numbins
        if (nbins <= 0) then
            error stop 'scifort: numbins must be positive', quiet=.true.
        end if

        call flatten_assumed_real(a, flat, dims)
        if (rank(a) < 1 .or. rank(a) > 15) then
            error stop 'scifort: frequency input rank must be 1 through 15', quiet=.true.
        end if
        if (present(weights)) then
            call flatten_assumed_real(weights, flat_weights, weight_dims)
            if (size(weight_dims) /= size(dims) .or. any(weight_dims /= dims)) then
                error stop 'scifort: weights must have the same shape as a', quiet=.true.
            end if
        else
            allocate(flat_weights(size(flat)))
            flat_weights = 1.0_dp
        end if

        if (present(defaultreallimits)) then
            lower = defaultreallimits(1)
            upper = defaultreallimits(2)
            if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper)) then
                error stop 'scifort: histogram limits must be finite', quiet=.true.
            end if
            if (upper < lower) then
                error stop 'scifort: histogram upper limit must not be below lower limit', quiet=.true.
            end if
        else if (size(flat) == 0) then
            lower = 0.0_dp
            upper = 1.0_dp
        else
            if (any(ieee_is_nan(flat))) then
                error stop 'scifort: automatic histogram limits require finite data', quiet=.true.
            end if
            data_min = minval(flat)
            data_max = maxval(flat)
            if (.not. ieee_is_finite(data_min) .or. .not. ieee_is_finite(data_max)) then
                error stop 'scifort: automatic histogram limits require finite data', quiet=.true.
            end if
            if (nbins == 1) then
                error stop 'scifort: automatic histogram limits are undefined for numbins=1', quiet=.true.
            end if
            s = (data_max - data_min) / (2.0_dp * real(nbins - 1, dp))
            lower = data_min - s
            upper = data_max + s
        end if

        lowerlimit = lower
        hist_lower = lower
        hist_upper = upper
        if (hist_upper <= hist_lower) then
            hist_lower = hist_lower - 0.5_dp
            hist_upper = hist_upper + 0.5_dp
        end if
        binsize = (hist_upper - hist_lower) / real(nbins, dp)
        allocate(histogram(nbins))
        histogram = 0.0_dp
        extrapoints = 0

        do i = 1, size(flat)
            if (ieee_is_nan(flat(i))) cycle
            if (flat(i) < lower .or. flat(i) > upper) extrapoints = extrapoints + 1
            if (flat(i) < hist_lower .or. flat(i) > hist_upper) cycle
            if (flat(i) >= hist_upper) then
                ibin = nbins
            else
                ibin = int((flat(i) - hist_lower) / binsize) + 1
                ibin = max(1, min(nbins, ibin))
            end if
            w = flat_weights(i)
            histogram(ibin) = histogram(ibin) + w
        end do
    end subroutine frequency_histogram

    pure subroutine select_score_values(x, selected, limit)
        real(dp), intent(in) :: x(:) !! observations before optional filtering
        real(dp), allocatable, intent(out) :: selected(:) !! values retained for percentile
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval

        logical, allocatable :: keep(:)
        integer :: i
        real(dp) :: lower
        real(dp) :: upper

        if (.not. present(limit)) then
            selected = x
            return
        end if
        lower = limit(1)
        upper = limit(2)
        if (ieee_is_nan(lower) .or. ieee_is_nan(upper) .or. upper < lower) then
            error stop 'scifort: invalid scoreatpercentile limit', quiet=.true.
        end if
        allocate(keep(size(x)))
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) then
                keep(i) = .false.
            else
                keep(i) = x(i) >= lower .and. x(i) <= upper
            end if
        end do
        selected = pack(x, keep)
    end subroutine select_score_values

    pure function sorted_percentile(sorted, q, interpolation_method) result(value)
        real(dp), intent(in) :: sorted(:) !! ascending observations with NaNs last
        real(dp), intent(in) :: q !! probability in [0, 1]
        character(len=*), intent(in) :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp) :: fraction
        real(dp) :: lower
        integer :: lower_index
        real(dp) :: position
        real(dp) :: upper
        integer :: upper_index

        if (size(sorted) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        position = q * real(size(sorted) - 1, dp)
        lower_index = int(floor(position)) + 1
        upper_index = int(ceiling(position)) + 1
        select case (trim(adjustl(interpolation_method)))
        case ('lower')
            value = sorted(lower_index)
            return
        case ('higher')
            value = sorted(upper_index)
            return
        end select

        fraction = position - floor(position)
        if (fraction <= 0.0_dp .or. lower_index == size(sorted)) then
            value = sorted(lower_index)
            return
        end if
        lower = sorted(lower_index)
        upper = sorted(lower_index + 1)
        value = safe_linear_interpolation(lower, upper, fraction)
    end function sorted_percentile

    pure function safe_linear_interpolation(lower, upper, fraction) result(value)
        real(dp), intent(in) :: lower !! lower sorted endpoint
        real(dp), intent(in) :: upper !! upper sorted endpoint
        real(dp), intent(in) :: fraction !! interpolation fraction in (0,1)
        real(dp) :: value

        if (ieee_is_nan(lower) .or. ieee_is_nan(upper)) then
            value = quiet_nan(0.0_dp)
        else if (ieee_is_finite(lower) .and. ieee_is_finite(upper)) then
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
    end function safe_linear_interpolation

    pure logical function valid_nan_policy(policy) result(valid)
        character(len=*), intent(in) :: policy !! requested NaN policy

        select case (trim(adjustl(policy)))
        case ('propagate', 'omit', 'raise')
            valid = .true.
        case default
            valid = .false.
        end select
    end function valid_nan_policy

    pure logical function valid_percentile_kind(kind) result(valid)
        character(len=*), intent(in) :: kind !! requested percentile rank convention

        select case (trim(adjustl(kind)))
        case ('rank', 'weak', 'strict', 'mean')
            valid = .true.
        case default
            valid = .false.
        end select
    end function valid_percentile_kind

    pure logical function valid_interpolation_method(method) result(valid)
        character(len=*), intent(in) :: method !! requested percentile interpolation convention

        select case (trim(adjustl(method)))
        case ('fraction', 'lower', 'higher')
            valid = .true.
        case default
            valid = .false.
        end select
    end function valid_interpolation_method

    pure subroutine sort_values_nan_last(values)
        real(dp), intent(inout) :: values(:) !! ascending values with NaNs placed last

        real(dp), allocatable :: temporary(:)

        if (size(values) <= 1) return
        allocate(temporary(size(values)))
        call merge_sort_values_nan_last(values, temporary, 1, size(values))
    end subroutine sort_values_nan_last

    recursive pure subroutine merge_sort_values_nan_last(values, temporary, left, right)
        real(dp), intent(inout) :: values(:) !! values being sorted
        real(dp), intent(inout) :: temporary(:) !! merge workspace
        integer, intent(in) :: left !! first index of inclusive range
        integer, intent(in) :: right !! last index of inclusive range

        integer :: i
        integer :: j
        integer :: k
        integer :: middle

        if (left >= right) return
        middle = left + (right - left) / 2
        call merge_sort_values_nan_last(values, temporary, left, middle)
        call merge_sort_values_nan_last(values, temporary, middle + 1, right)
        i = left
        j = middle + 1
        do k = left, right
            if (i > middle) then
                temporary(k) = values(j)
                j = j + 1
            else if (j > right) then
                temporary(k) = values(i)
                i = i + 1
            else if (comes_before(values(i), values(j))) then
                temporary(k) = values(i)
                i = i + 1
            else
                temporary(k) = values(j)
                j = j + 1
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

    pure subroutine flatten_assumed_real(x, flat, dims)
        real(dp), intent(in) :: x(..) !! array to flatten in Fortran element order
        real(dp), allocatable, intent(out) :: flat(:) !! flattened values
        integer, allocatable, intent(out) :: dims(:) !! source extents

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
            allocate(dims(0))
            allocate(flat(0))
        end select
    end subroutine flatten_assumed_real

end module scifort_descriptive_frequency_engine
