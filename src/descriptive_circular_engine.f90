! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_circular_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    integer, parameter, public :: circular_kind_mean = 1
    integer, parameter, public :: circular_kind_var = 2
    integer, parameter, public :: circular_kind_std = 3

    public :: circular_all_assumed_rank
    public :: circular_reduce_assumed_rank
    public :: circular_vector

contains

    pure function circular_vector(x, stat_kind, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:) !! angle observations in one reduction slice
        integer, intent(in) :: stat_kind !! internal selector for mean, variance, or std
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        real(dp), allocatable :: clean(:)
        real(dp) :: cos_mean
        real(dp) :: cos_sum
        real(dp) :: high_value
        integer :: i
        real(dp) :: low_value
        logical :: normalize_value
        real(dp) :: period
        character(len=16) :: policy
        real(dp) :: r
        real(dp) :: scale
        real(dp) :: sin_mean
        real(dp) :: sin_sum
        real(dp) :: theta

        if (stat_kind < circular_kind_mean .or. stat_kind > circular_kind_std) then
            value = quiet_nan(0.0_dp)
            return
        end if

        high_value = 2.0_dp * scifort_pi
        if (present(high)) high_value = high
        low_value = 0.0_dp
        if (present(low)) low_value = low
        period = high_value - low_value
        if (abs(period) < tiny(1.0_dp)) then
            error stop 'scifort: circular high and low must differ', quiet=.true.
        end if

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
            clean = x
        else if (policy == 'propagate') then
            if (any(ieee_is_nan(x))) then
                value = quiet_nan(0.0_dp)
                return
            end if
            clean = x
        else
            clean = pack(x, .not. ieee_is_nan(x))
        end if

        if (size(clean) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if

        scale = (2.0_dp * scifort_pi) / period
        sin_sum = 0.0_dp
        cos_sum = 0.0_dp
        do i = 1, size(clean)
            theta = clean(i) * scale
            sin_sum = sin_sum + sin(theta)
            cos_sum = cos_sum + cos(theta)
        end do

        if (stat_kind == circular_kind_mean) then
            theta = atan2(sin_sum, cos_sum)
            value = modulo(theta / scale - low_value, period) + low_value
            return
        end if

        sin_mean = sin_sum / real(size(clean), dp)
        cos_mean = cos_sum / real(size(clean), dp)
        r = sqrt(sin_mean * sin_mean + cos_mean * cos_mean)
        r = min(1.0_dp, r)
        if (stat_kind == circular_kind_var) then
            value = 1.0_dp - r
            return
        end if

        normalize_value = .false.
        if (present(normalize)) normalize_value = normalize
        if (r <= 0.0_dp) then
            value = positive_infinity(0.0_dp)
        else
            value = sqrt(-2.0_dp * log(r))
        end if
        if (.not. normalize_value) value = value / scale
    end function circular_vector

    pure subroutine circular_reduce_assumed_rank(x, axis, stat_kind, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        real(dp), intent(in) :: x(..) !! angle observations with rank 2 to 15
        integer, intent(in) :: axis !! one-based Fortran axis to reduce
        integer, intent(in) :: stat_kind !! internal selector for mean, variance, or std
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened reduction result
        integer, intent(out) :: out_shape(:) !! shape with reduced axis removed
        logical, intent(out) :: ok !! true when rank, axis, and selector are valid
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true

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
        if (stat_kind < circular_kind_mean .or. stat_kind > circular_kind_std) then
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
            reduced(i + 1) = circular_vector(slice, stat_kind, high, low, &
                nan_policy, normalize)
        end do
        ok = .true.
    end subroutine circular_reduce_assumed_rank

    pure function circular_all_assumed_rank(x, stat_kind, high, low, nan_policy, &
            normalize) result(value)
        real(dp), intent(in) :: x(..) !! angle observations with rank 1 to 15
        integer, intent(in) :: stat_kind !! internal selector for mean, variance, or std
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        if (rank(x) < 1 .or. rank(x) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if
        call flatten_assumed_real(x, flat, dims)
        value = circular_vector(flat, stat_kind, high, low, nan_policy, normalize)
    end function circular_all_assumed_rank

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
        real(dp), intent(in) :: x(..) !! input array of rank 1 to 15
        real(dp), allocatable, intent(out) :: flat(:) !! flattened values
        integer, allocatable, intent(out) :: dims(:) !! input shape

        allocate(dims(rank(x)))
        dims = shape(x)
        select rank (x)
        rank (1)
            flat = reshape(x, [size(x)])
        rank (2)
            flat = reshape(x, [size(x)])
        rank (3)
            flat = reshape(x, [size(x)])
        rank (4)
            flat = reshape(x, [size(x)])
        rank (5)
            flat = reshape(x, [size(x)])
        rank (6)
            flat = reshape(x, [size(x)])
        rank (7)
            flat = reshape(x, [size(x)])
        rank (8)
            flat = reshape(x, [size(x)])
        rank (9)
            flat = reshape(x, [size(x)])
        rank (10)
            flat = reshape(x, [size(x)])
        rank (11)
            flat = reshape(x, [size(x)])
        rank (12)
            flat = reshape(x, [size(x)])
        rank (13)
            flat = reshape(x, [size(x)])
        rank (14)
            flat = reshape(x, [size(x)])
        rank (15)
            flat = reshape(x, [size(x)])
        rank default
            allocate(flat(0))
        end select
    end subroutine flatten_assumed_real

end module scifort_descriptive_circular_engine
