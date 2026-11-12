! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_trimmed_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: trimmed_kind_mean = 1
    integer, parameter, public :: trimmed_kind_var = 2
    integer, parameter, public :: trimmed_kind_std = 3
    integer, parameter, public :: trimmed_kind_sem = 4
    integer, parameter, public :: trimmed_kind_min = 5
    integer, parameter, public :: trimmed_kind_max = 6

    public :: trimmed_reduce_assumed_rank
    public :: trimmed_stat_all_assumed_rank
    public :: trimmed_stat_vector

contains

    pure function trimmed_stat_vector(x, stat_kind, lower_limit, upper_limit, &
            inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations in one reduction slice
        integer, intent(in) :: stat_kind !! internal trimmed-statistic selector
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment for variance family
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: clean(:)
        integer :: ddof_value
        logical :: include_lower
        logical :: include_upper
        logical :: keep_value
        logical :: mask(size(x))
        real(dp) :: mean_value
        integer :: i
        integer :: n
        character(len=16) :: policy
        real(dp) :: variance_value

        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))
        if (.not. valid_nan_policy(policy)) then
            value = quiet_nan(0.0_dp)
            return
        end if

        if (policy == 'raise') then
            if (any(ieee_is_nan(x))) then
                error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
            end if
        else if (policy == 'propagate') then
            if (any(ieee_is_nan(x))) then
                value = quiet_nan(0.0_dp)
                return
            end if
        end if

        include_lower = .true.
        include_upper = .true.
        if (present(inclusive)) then
            include_lower = inclusive(1)
            include_upper = inclusive(2)
        end if

        mask = .false.
        do i = 1, size(x)
            if (policy == 'omit' .and. ieee_is_nan(x(i))) cycle
            keep_value = .true.
            if (present(lower_limit)) then
                if (include_lower) then
                    keep_value = x(i) >= lower_limit
                else
                    keep_value = x(i) > lower_limit
                end if
            end if
            if (keep_value .and. present(upper_limit)) then
                if (include_upper) then
                    keep_value = x(i) <= upper_limit
                else
                    keep_value = x(i) < upper_limit
                end if
            end if
            mask(i) = keep_value
        end do

        clean = pack(x, mask)
        n = size(clean)
        if (n == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if

        select case (stat_kind)
        case (trimmed_kind_mean)
            value = arithmetic_mean(clean)
        case (trimmed_kind_var, trimmed_kind_std, trimmed_kind_sem)
            ddof_value = 1
            if (present(ddof)) ddof_value = ddof
            if (n - ddof_value <= 0) then
                value = quiet_nan(0.0_dp)
                return
            end if
            mean_value = arithmetic_mean(clean)
            variance_value = sum((clean - mean_value) ** 2) / real(n - ddof_value, dp)
            if (variance_value < 0.0_dp .and. abs(variance_value) <= epsilon(1.0_dp) * &
                    max(1.0_dp, mean_value * mean_value)) then
                variance_value = 0.0_dp
            end if
            select case (stat_kind)
            case (trimmed_kind_var)
                value = variance_value
            case (trimmed_kind_std)
                value = sqrt(variance_value)
            case (trimmed_kind_sem)
                value = sqrt(variance_value) / sqrt(real(n, dp))
            end select
        case (trimmed_kind_min)
            value = minval(clean)
        case (trimmed_kind_max)
            value = maxval(clean)
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function trimmed_stat_vector

    pure subroutine trimmed_reduce_assumed_rank(x, axis, stat_kind, reduced, out_shape, &
            ok, lower_limit, upper_limit, inclusive, ddof, nan_policy)
        real(dp), intent(in) :: x(..) !! observations with rank 2 to 15
        integer, intent(in) :: axis !! one-based Fortran axis to reduce
        integer, intent(in) :: stat_kind !! internal trimmed-statistic selector
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened reduction result
        integer, intent(out) :: out_shape(:) !! shape with reduced axis removed
        logical, intent(out) :: ok !! true when rank and axis are valid
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment for variance family
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise

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
            reduced(i + 1) = trimmed_stat_vector(slice, stat_kind, lower_limit, &
                upper_limit, inclusive, ddof, nan_policy)
        end do
        ok = .true.
    end subroutine trimmed_reduce_assumed_rank

    pure function trimmed_stat_all_assumed_rank(x, stat_kind, lower_limit, upper_limit, &
            inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: stat_kind !! internal trimmed-statistic selector
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment for variance family
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        if (rank(x) < 1 .or. rank(x) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if
        call flatten_assumed_real(x, flat, dims)
        value = trimmed_stat_vector(flat, stat_kind, lower_limit, upper_limit, &
            inclusive, ddof, nan_policy)
    end function trimmed_stat_all_assumed_rank

    pure function arithmetic_mean(x) result(value)
        real(dp), intent(in) :: x(:) !! nonempty observations
        real(dp) :: value

        real(dp) :: correction
        integer :: i
        real(dp) :: sum_value
        real(dp) :: term
        real(dp) :: updated

        sum_value = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(x)
            term = x(i) - correction
            updated = sum_value + term
            correction = (updated - sum_value) - term
            sum_value = updated
        end do
        value = sum_value / real(size(x), dp)
    end function arithmetic_mean

    pure logical function valid_nan_policy(nan_policy) result(valid)
        character(len=*), intent(in) :: nan_policy !! requested NaN handling policy

        select case (trim(adjustl(nan_policy)))
        case ('propagate', 'omit', 'raise')
            valid = .true.
        case default
            valid = .false.
        end select
    end function valid_nan_policy

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

end module scifort_descriptive_trimmed_engine
