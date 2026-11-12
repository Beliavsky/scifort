! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_standardize_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: gstd_vector_core
    public :: gzscore_assumed_rank
    public :: gzscore_all_assumed_rank
    public :: zscore_assumed_rank
    public :: zscore_all_assumed_rank
    public :: zmap_assumed_rank
    public :: zmap_all_assumed_rank
    public :: gstd_reduce_assumed_rank
    public :: gstd_reduce_dims_keepdims_assumed_rank
    public :: gstd_all_assumed_rank

contains

    pure subroutine zscore_all_assumed_rank(x, ddof, nan_policy, values, dims, ok)
        real(dp), intent(in) :: x(..) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable, intent(out) :: values(:) !! flattened whole-array z-scores
        integer, allocatable, intent(out) :: dims(:) !! input extents
        logical, intent(out) :: ok !! true when the input rank is supported

        call standardize_all_assumed_rank(x, ddof, nan_policy, .false., values, dims, ok)
    end subroutine zscore_all_assumed_rank

    pure subroutine gzscore_all_assumed_rank(x, ddof, nan_policy, values, dims, ok)
        real(dp), intent(in) :: x(..) !! positive observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable, intent(out) :: values(:) !! flattened whole-array geometric z-scores
        integer, allocatable, intent(out) :: dims(:) !! input extents
        logical, intent(out) :: ok !! true when the input rank is supported

        call standardize_all_assumed_rank(x, ddof, nan_policy, .true., values, dims, ok)
    end subroutine gzscore_all_assumed_rank

    pure subroutine zmap_all_assumed_rank(scores, compare, ddof, nan_policy, values, &
            dims, ok)
        real(dp), intent(in) :: scores(..) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(..) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable, intent(out) :: values(:) !! flattened relative z-scores
        integer, allocatable, intent(out) :: dims(:) !! score extents
        logical, intent(out) :: ok !! true when ranks are supported and compatible

        integer :: divisor_adjustment
        real(dp) :: mean_scaled
        logical :: stats_ok
        character(len=9) :: policy
        real(dp) :: scale_value
        integer, allocatable :: compare_dims(:)
        real(dp), allocatable :: compare_flat(:)
        integer :: i
        real(dp), allocatable :: score_flat(:)
        real(dp) :: sd_scaled

        ok = .false.
        call flatten_assumed_real(scores, score_flat, dims)
        call flatten_assumed_real(compare, compare_flat, compare_dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(values(0))
            return
        end if
        if (size(compare_dims) /= size(dims)) then
            deallocate(dims)
            allocate(dims(0))
            allocate(values(0))
            return
        end if

        allocate(values(size(score_flat)))
        values = quiet_nan(0.0_dp)
        divisor_adjustment = 0
        if (present(ddof)) divisor_adjustment = ddof
        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))
        if (.not. valid_nan_policy(policy)) then
            ok = .true.
            return
        end if

        call comparison_mean_sd_core(compare_flat, divisor_adjustment, policy, &
            mean_scaled, sd_scaled, scale_value, stats_ok)
        if (.not. stats_ok) then
            ok = .true.
            return
        end if
        do i = 1, size(score_flat)
            values(i) = relative_standard_score(score_flat(i), mean_scaled, sd_scaled, &
                scale_value)
        end do
        ok = .true.
    end subroutine zmap_all_assumed_rank

    pure subroutine gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        real(dp), intent(in) :: x(..) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), intent(out) :: value !! geometric standard deviation over all elements
        logical, intent(out) :: ok !! true when the input rank is supported

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        ok = .false.
        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = gstd_vector_core(flat, ddof, nan_policy)
        ok = .true.
    end subroutine gstd_all_assumed_rank

    pure subroutine zscore_assumed_rank(x, axis, ddof, nan_policy, values, dims, ok)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: axis !! Fortran axis, 1 to rank(x)
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable, intent(out) :: values(:) !! flattened z-scores
        integer, allocatable, intent(out) :: dims(:) !! input extents
        logical, intent(out) :: ok !! true when rank and axis are valid

        call standardize_assumed_rank(x, axis, ddof, nan_policy, .false., values, dims, ok)
    end subroutine zscore_assumed_rank

    pure subroutine gzscore_assumed_rank(x, axis, ddof, nan_policy, values, dims, ok)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: axis !! Fortran axis, 1 to rank(x)
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable, intent(out) :: values(:) !! flattened geometric z-scores
        integer, allocatable, intent(out) :: dims(:) !! input extents
        logical, intent(out) :: ok !! true when rank and axis are valid

        call standardize_assumed_rank(x, axis, ddof, nan_policy, .true., values, dims, ok)
    end subroutine gzscore_assumed_rank

    pure subroutine zmap_assumed_rank(scores, compare, axis, ddof, nan_policy, &
            values, dims, ok)
        real(dp), intent(in) :: scores(..) !! values to standardize, rank 1 to 15
        real(dp), intent(in) :: compare(..) !! reference values supplying mean and variance
        integer, intent(in) :: axis !! Fortran axis reduced in compare, 1 to rank(scores)
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable, intent(out) :: values(:) !! flattened relative z-scores
        integer, allocatable, intent(out) :: dims(:) !! broadcast output extents
        logical, intent(out) :: ok !! true when ranks, axis, and shapes are compatible

        integer :: compare_base
        integer :: compare_coord
        integer, allocatable :: compare_dims(:)
        real(dp), allocatable :: compare_flat(:)
        integer, allocatable :: compare_strides(:)
        integer :: coordinate
        integer :: d
        integer :: divisor_adjustment
        integer :: i
        integer :: linear
        real(dp), allocatable :: mean_scaled(:)
        integer :: n_stats
        character(len=9) :: policy
        real(dp), allocatable :: scale_value(:)
        integer :: score_index
        integer, allocatable :: score_dims(:)
        real(dp), allocatable :: score_flat(:)
        integer, allocatable :: score_strides(:)
        real(dp), allocatable :: sd_scaled(:)
        real(dp), allocatable :: slice(:)
        integer :: stat_index
        integer, allocatable :: stat_strides(:)
        logical, allocatable :: stat_valid(:)

        ok = .false.
        call flatten_assumed_real(scores, score_flat, score_dims)
        call flatten_assumed_real(compare, compare_flat, compare_dims)
        if (size(score_dims) < 1 .or. size(score_dims) > 15) then
            allocate(values(0))
            allocate(dims(0))
            return
        end if
        if (size(compare_dims) /= size(score_dims)) then
            allocate(values(0))
            allocate(dims(0))
            return
        end if
        if (axis < 1 .or. axis > size(score_dims)) then
            allocate(values(0))
            allocate(dims(0))
            return
        end if

        divisor_adjustment = 0
        if (present(ddof)) divisor_adjustment = ddof
        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))

        allocate(dims(size(score_dims)))
        do d = 1, size(score_dims)
            if (d == axis) then
                dims(d) = score_dims(d)
            else if (score_dims(d) == compare_dims(d)) then
                dims(d) = score_dims(d)
            else if (score_dims(d) == 1) then
                dims(d) = compare_dims(d)
            else if (compare_dims(d) == 1) then
                dims(d) = score_dims(d)
            else
                deallocate(dims)
                allocate(dims(0))
                allocate(values(0))
                return
            end if
        end do

        allocate(values(product(dims)))
        if (.not. valid_nan_policy(policy)) then
            values = quiet_nan(0.0_dp)
            ok = .true.
            return
        end if
        if (size(values) == 0) then
            ok = .true.
            return
        end if

        allocate(score_strides(size(score_dims)))
        allocate(compare_strides(size(compare_dims)))
        allocate(stat_strides(size(compare_dims)))
        score_strides(1) = 1
        compare_strides(1) = 1
        do d = 2, size(score_dims)
            score_strides(d) = score_strides(d - 1) * score_dims(d - 1)
            compare_strides(d) = compare_strides(d - 1) * compare_dims(d - 1)
        end do
        stat_strides = 0
        n_stats = 1
        do d = 1, size(compare_dims)
            if (d == axis) cycle
            stat_strides(d) = n_stats
            n_stats = n_stats * compare_dims(d)
        end do

        allocate(mean_scaled(n_stats))
        allocate(sd_scaled(n_stats))
        allocate(scale_value(n_stats))
        allocate(stat_valid(n_stats))
        allocate(slice(compare_dims(axis)))

        do i = 0, n_stats - 1
            linear = i
            compare_base = 1
            do d = 1, size(compare_dims)
                if (d == axis) cycle
                compare_coord = mod(linear, compare_dims(d))
                linear = linear / compare_dims(d)
                compare_base = compare_base + compare_coord * compare_strides(d)
            end do
            do d = 1, compare_dims(axis)
                slice(d) = compare_flat(compare_base + (d - 1) * compare_strides(axis))
            end do
            call comparison_mean_sd_core(slice, divisor_adjustment, policy, &
                mean_scaled(i + 1), sd_scaled(i + 1), scale_value(i + 1), &
                stat_valid(i + 1))
        end do

        do i = 0, size(values) - 1
            linear = i
            score_index = 1
            stat_index = 1
            do d = 1, size(dims)
                coordinate = mod(linear, dims(d))
                linear = linear / dims(d)
                if (score_dims(d) /= 1) then
                    score_index = score_index + coordinate * score_strides(d)
                end if
                if (d /= axis .and. compare_dims(d) /= 1) then
                    stat_index = stat_index + coordinate * stat_strides(d)
                end if
            end do
            if (stat_valid(stat_index)) then
                values(i + 1) = relative_standard_score(score_flat(score_index), &
                    mean_scaled(stat_index), sd_scaled(stat_index), &
                    scale_value(stat_index))
            else
                values(i + 1) = quiet_nan(0.0_dp)
            end if
        end do
        ok = .true.
    end subroutine zmap_assumed_rank

    pure subroutine standardize_all_assumed_rank(x, ddof, nan_policy, geometric, &
            values, dims, ok)
        real(dp), intent(in) :: x(..) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        logical, intent(in) :: geometric !! use logarithms before standardization
        real(dp), allocatable, intent(out) :: values(:) !! flattened standardized data
        integer, allocatable, intent(out) :: dims(:) !! input extents
        logical, intent(out) :: ok !! true when the input rank is supported

        real(dp), allocatable :: flat(:)

        ok = .false.
        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(values(0))
            return
        end if
        if (geometric) then
            values = gzscore_vector_core(flat, ddof, nan_policy)
        else
            values = zscore_vector_core(flat, ddof, nan_policy)
        end if
        ok = .true.
    end subroutine standardize_all_assumed_rank

    pure subroutine standardize_assumed_rank(x, axis, ddof, nan_policy, geometric, &
            values, dims, ok)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: axis !! Fortran axis, 1 to rank(x)
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        logical, intent(in) :: geometric !! use logarithms before standardization
        real(dp), allocatable, intent(out) :: values(:) !! flattened standardized data
        integer, allocatable, intent(out) :: dims(:) !! input extents
        logical, intent(out) :: ok !! true when rank and axis are valid

        integer :: base_index
        integer :: coordinate
        integer :: d
        integer :: i
        integer :: linear
        integer :: n_out
        integer :: n_reduce
        integer :: out_dim
        integer, allocatable :: out_shape(:)
        real(dp), allocatable :: flat(:)
        real(dp), allocatable :: scores(:)
        real(dp), allocatable :: slice(:)
        integer, allocatable :: strides(:)

        ok = .false.
        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(values(0))
            return
        end if
        if (axis < 1 .or. axis > size(dims)) then
            allocate(values(0))
            return
        end if

        allocate(values(size(flat)))
        if (size(flat) == 0) then
            ok = .true.
            return
        end if

        allocate(strides(size(dims)))
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
        end do

        allocate(out_shape(size(dims) - 1))
        out_dim = 0
        do d = 1, size(dims)
            if (d == axis) cycle
            out_dim = out_dim + 1
            out_shape(out_dim) = dims(d)
        end do
        n_reduce = dims(axis)
        n_out = product(out_shape)
        allocate(slice(n_reduce))

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
            if (geometric) then
                scores = gzscore_vector_core(slice, ddof, nan_policy)
            else
                scores = zscore_vector_core(slice, ddof, nan_policy)
            end if
            do d = 1, n_reduce
                values(base_index + (d - 1) * strides(axis)) = scores(d)
            end do
        end do
        ok = .true.
    end subroutine standardize_assumed_rank

    pure subroutine gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, &
            reduced, out_shape, ok)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened keepdims results
        integer, intent(out) :: out_shape(:) !! input extents with reduced axes set to one
        logical, intent(out) :: ok !! true when rank, axes, and output rank are valid

        integer :: base_index
        integer :: coordinate
        integer :: d
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        integer :: i
        integer :: input_index
        logical, allocatable :: is_reduced(:)
        integer :: j
        integer :: linear
        integer :: n_out
        integer :: n_reduce
        real(dp), allocatable :: slice(:)
        integer, allocatable :: strides(:)

        ok = .false.
        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if
        if (size(out_shape) /= size(dims) .or. size(axes) == 0) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if

        allocate(is_reduced(size(dims)))
        is_reduced = .false.
        do i = 1, size(axes)
            if (axes(i) < 1 .or. axes(i) > size(dims)) then
                allocate(reduced(0))
                out_shape = 0
                return
            end if
            if (is_reduced(axes(i))) then
                allocate(reduced(0))
                out_shape = 0
                return
            end if
            is_reduced(axes(i)) = .true.
        end do

        out_shape = dims
        do d = 1, size(dims)
            if (is_reduced(d)) out_shape(d) = 1
        end do
        n_out = product(out_shape)
        n_reduce = 1
        do d = 1, size(dims)
            if (is_reduced(d)) n_reduce = n_reduce * dims(d)
        end do
        allocate(reduced(n_out))
        allocate(slice(n_reduce))
        if (n_out == 0) then
            ok = .true.
            return
        end if

        allocate(strides(size(dims)))
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
        end do

        do i = 0, n_out - 1
            linear = i
            base_index = 1
            do d = 1, size(dims)
                if (is_reduced(d)) cycle
                coordinate = mod(linear, dims(d))
                linear = linear / dims(d)
                base_index = base_index + coordinate * strides(d)
            end do
            do j = 0, n_reduce - 1
                linear = j
                input_index = base_index
                do d = 1, size(dims)
                    if (.not. is_reduced(d)) cycle
                    coordinate = mod(linear, dims(d))
                    linear = linear / dims(d)
                    input_index = input_index + coordinate * strides(d)
                end do
                slice(j + 1) = flat(input_index)
            end do
            reduced(i + 1) = gstd_vector_core(slice, ddof, nan_policy)
        end do
        ok = .true.
    end subroutine gstd_reduce_dims_keepdims_assumed_rank

    pure subroutine gstd_reduce_assumed_rank(x, axis, ddof, nan_policy, reduced, &
            out_shape, ok)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: axis !! Fortran axis, 1 to rank(x)
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened geometric standard deviations
        integer, intent(out) :: out_shape(:) !! extents after removing the reduced axis
        logical, intent(out) :: ok !! true when rank, axis, and output rank are valid

        integer :: base_index
        integer :: coordinate
        integer :: d
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        integer :: i
        integer :: linear
        integer :: n_out
        integer :: n_reduce
        integer :: out_dim
        real(dp), allocatable :: slice(:)
        integer, allocatable :: strides(:)

        ok = .false.
        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if
        if (size(out_shape) /= size(dims) - 1) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if
        if (axis < 1 .or. axis > size(dims)) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if

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
        allocate(reduced(n_out))
        allocate(slice(n_reduce))

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
            reduced(i + 1) = gstd_vector_core(slice, ddof, nan_policy)
        end do
        ok = .true.
    end subroutine gstd_reduce_assumed_rank

    pure function zscore_vector_core(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations to standardize
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: divisor_adjustment
        logical, allocatable :: valid(:)
        integer :: i
        real(dp) :: mean_scaled
        integer :: n_valid
        logical :: ok
        character(len=9) :: policy
        real(dp) :: scale_value
        real(dp) :: sd_scaled

        allocate(values(size(x)))
        values = quiet_nan(0.0_dp)
        divisor_adjustment = 0
        if (present(ddof)) divisor_adjustment = ddof
        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))
        if (.not. valid_nan_policy(policy)) return

        allocate(valid(size(x)))
        valid = .true.
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) then
                if (policy == 'omit') then
                    valid(i) = .false.
                else if (policy == 'raise') then
                    error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
                else
                    return
                end if
            else if (.not. ieee_is_finite(x(i))) then
                return
            end if
        end do
        n_valid = count(valid)
        if (n_valid == 0 .or. n_valid - divisor_adjustment <= 0) return

        call finite_scaled_mean_sd(x, valid, divisor_adjustment, mean_scaled, &
            sd_scaled, scale_value, ok)
        if (.not. ok .or. sd_scaled <= 0.0_dp) return
        do i = 1, size(x)
            if (valid(i)) then
                values(i) = (x(i) / scale_value - mean_scaled) / sd_scaled
            end if
        end do
    end function zscore_vector_core

    pure function gzscore_vector_core(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations transformed by the natural logarithm
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: i
        real(dp), allocatable :: logged(:)

        allocate(logged(size(x)))
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) then
                logged(i) = quiet_nan(0.0_dp)
            else if (x(i) < 0.0_dp) then
                logged(i) = quiet_nan(0.0_dp)
            else if (x(i) <= 0.0_dp) then
                logged(i) = negative_infinity(0.0_dp)
            else if (.not. ieee_is_finite(x(i))) then
                logged(i) = positive_infinity(0.0_dp)
            else
                logged(i) = log(x(i))
            end if
        end do
        values = zscore_vector_core(logged, ddof, nan_policy)
    end function gzscore_vector_core

    pure function gstd_vector_core(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations for geometric standard deviation
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        integer :: divisor_adjustment
        integer :: i
        real(dp), allocatable :: logged(:)
        logical, allocatable :: valid(:)
        real(dp) :: mean_scaled
        integer :: n_valid
        logical :: ok
        character(len=9) :: policy
        real(dp) :: scale_value
        real(dp) :: sd_scaled
        real(dp) :: sd_value

        divisor_adjustment = 1
        if (present(ddof)) divisor_adjustment = ddof
        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))
        if (.not. valid_nan_policy(policy)) then
            value = quiet_nan(0.0_dp)
            return
        end if

        allocate(logged(size(x)))
        allocate(valid(size(x)))
        valid = .true.
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) then
                if (policy == 'omit') then
                    valid(i) = .false.
                    logged(i) = 0.0_dp
                else if (policy == 'raise') then
                    error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
                else
                    value = quiet_nan(0.0_dp)
                    return
                end if
            else if (x(i) < 0.0_dp) then
                if (policy == 'omit') then
                    valid(i) = .false.
                    logged(i) = 0.0_dp
                else if (policy == 'raise') then
                    error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
                else
                    value = quiet_nan(0.0_dp)
                    return
                end if
            else if (x(i) <= 0.0_dp .or. .not. ieee_is_finite(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            else
                logged(i) = log(x(i))
            end if
        end do

        n_valid = count(valid)
        if (n_valid == 0 .or. n_valid - divisor_adjustment <= 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        call finite_scaled_mean_sd(logged, valid, divisor_adjustment, mean_scaled, &
            sd_scaled, scale_value, ok)
        if (.not. ok) then
            value = quiet_nan(0.0_dp)
            return
        end if
        sd_value = scale_value * sd_scaled
        if (sd_value > log(huge(1.0_dp))) then
            value = positive_infinity(0.0_dp)
        else
            value = exp(sd_value)
        end if
    end function gstd_vector_core

    pure subroutine comparison_mean_sd_core(x, ddof, policy, mean_scaled, &
            sd_scaled, scale_value, ok)
        real(dp), intent(in) :: x(:) !! comparison observations for one reduction slice
        integer, intent(in) :: ddof !! divisor adjustment
        character(len=*), intent(in) :: policy !! validated 'propagate', 'omit', or 'raise' policy
        real(dp), intent(out) :: mean_scaled !! mean divided by scale_value
        real(dp), intent(out) :: sd_scaled !! standard deviation divided by scale_value
        real(dp), intent(out) :: scale_value !! maximum absolute included comparison value
        logical, intent(out) :: ok !! true when comparison statistics are defined

        integer :: i
        integer :: n_valid
        logical, allocatable :: valid(:)

        ok = .false.
        mean_scaled = quiet_nan(0.0_dp)
        sd_scaled = quiet_nan(0.0_dp)
        scale_value = 0.0_dp
        allocate(valid(size(x)))
        valid = .true.
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) then
                if (policy == 'omit') then
                    valid(i) = .false.
                else if (policy == 'raise') then
                    error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
                else
                    return
                end if
            else if (.not. ieee_is_finite(x(i))) then
                return
            end if
        end do
        n_valid = count(valid)
        if (n_valid == 0 .or. n_valid - ddof <= 0) return
        call finite_scaled_mean_sd(x, valid, ddof, mean_scaled, sd_scaled, &
            scale_value, ok)
    end subroutine comparison_mean_sd_core

    pure function relative_standard_score(score, mean_scaled, sd_scaled, &
            scale_value) result(value)
        real(dp), intent(in) :: score !! value standardized against comparison statistics
        real(dp), intent(in) :: mean_scaled !! comparison mean divided by scale_value
        real(dp), intent(in) :: sd_scaled !! comparison standard deviation divided by scale_value
        real(dp), intent(in) :: scale_value !! comparison scaling magnitude
        real(dp) :: value

        real(dp) :: common_scale
        real(dp) :: denominator_scaled
        real(dp) :: numerator_scaled

        if (ieee_is_nan(score)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (.not. ieee_is_finite(score)) then
            if (score > 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else
                value = negative_infinity(0.0_dp)
            end if
            return
        end if

        common_scale = max(abs(score), scale_value)
        if (common_scale <= 0.0_dp) then
            value = quiet_nan(0.0_dp)
            return
        end if
        numerator_scaled = score / common_scale - &
            mean_scaled * (scale_value / common_scale)
        if (sd_scaled <= 0.0_dp) then
            if (numerator_scaled > 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else if (numerator_scaled < 0.0_dp) then
                value = negative_infinity(0.0_dp)
            else
                value = quiet_nan(0.0_dp)
            end if
            return
        end if

        denominator_scaled = sd_scaled * (scale_value / common_scale)
        if (denominator_scaled <= 0.0_dp) then
            if (numerator_scaled > 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else if (numerator_scaled < 0.0_dp) then
                value = negative_infinity(0.0_dp)
            else
                value = quiet_nan(0.0_dp)
            end if
        else if (denominator_scaled < 1.0_dp .and. &
                abs(numerator_scaled) > huge(1.0_dp) * denominator_scaled) then
            if (numerator_scaled > 0.0_dp) then
                value = positive_infinity(0.0_dp)
            else
                value = negative_infinity(0.0_dp)
            end if
        else
            value = numerator_scaled / denominator_scaled
        end if
    end function relative_standard_score

    pure logical function valid_nan_policy(policy) result(valid)
        character(len=*), intent(in) :: policy !! requested NaN handling policy

        select case (trim(adjustl(policy)))
        case ('propagate', 'omit', 'raise')
            valid = .true.
        case default
            valid = .false.
        end select
    end function valid_nan_policy

    pure subroutine finite_scaled_mean_sd(x, valid, ddof, mean_scaled, sd_scaled, &
            scale_value, ok)
        real(dp), intent(in) :: x(:) !! finite observations
        logical, intent(in) :: valid(:) !! observations included in the calculation
        integer, intent(in) :: ddof !! divisor adjustment
        real(dp), intent(out) :: mean_scaled !! mean divided by scale_value
        real(dp), intent(out) :: sd_scaled !! standard deviation divided by scale_value
        real(dp), intent(out) :: scale_value !! maximum absolute included observation
        logical, intent(out) :: ok !! true when the requested divisor is positive

        real(dp) :: correction
        integer :: i
        integer :: n_valid
        real(dp) :: sum_scaled
        real(dp) :: sum_squares
        real(dp) :: term
        real(dp) :: updated
        real(dp) :: centered

        ok = .false.
        mean_scaled = quiet_nan(0.0_dp)
        sd_scaled = quiet_nan(0.0_dp)
        scale_value = 0.0_dp
        if (size(x) /= size(valid)) return
        n_valid = count(valid)
        if (n_valid == 0 .or. n_valid - ddof <= 0) return

        do i = 1, size(x)
            if (valid(i)) then
                if (.not. ieee_is_finite(x(i))) return
                scale_value = max(scale_value, abs(x(i)))
            end if
        end do
        if (scale_value <= 0.0_dp) then
            mean_scaled = 0.0_dp
            sd_scaled = 0.0_dp
            ok = .true.
            return
        end if

        sum_scaled = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(x)
            if (.not. valid(i)) cycle
            term = x(i) / scale_value - correction
            updated = sum_scaled + term
            correction = (updated - sum_scaled) - term
            sum_scaled = updated
        end do
        mean_scaled = sum_scaled / real(n_valid, dp)

        sum_squares = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(x)
            if (.not. valid(i)) cycle
            centered = x(i) / scale_value - mean_scaled
            term = centered * centered - correction
            updated = sum_squares + term
            correction = (updated - sum_squares) - term
            sum_squares = updated
        end do
        sd_scaled = sqrt(max(0.0_dp, sum_squares / real(n_valid - ddof, dp)))
        ok = .true.
    end subroutine finite_scaled_mean_sd

    pure subroutine flatten_assumed_real(x, flat, dims)
        real(dp), intent(in) :: x(..) !! array to flatten in Fortran element order
        real(dp), allocatable, intent(out) :: flat(:) !! flattened observations
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

end module scifort_descriptive_standardize_engine
