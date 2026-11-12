! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_kstat_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: kstat_kind_stat = 1
    integer, parameter, public :: kstat_kind_var = 2

    public :: kstat_all_assumed_rank
    public :: kstat_reduce_assumed_rank
    public :: kstat_vector
    public :: kstatvar_vector

contains

    pure function kstat_vector(x, order, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations in one reduction slice
        integer, intent(in) :: order !! k-statistic order in the range 1 to 4
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: clean(:)
        real(dp), allocatable :: centered(:)
        real(dp) :: mean_value
        real(dp) :: m2
        real(dp) :: m3
        real(dp) :: m4
        integer :: n
        real(dp) :: nr
        character(len=16) :: policy

        if (order < 1 .or. order > 4) then
            error stop 'scifort: kstat order must be between 1 and 4', quiet=.true.
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

        n = size(clean)
        if (n < order) then
            value = quiet_nan(0.0_dp)
            return
        end if
        nr = real(n, dp)
        mean_value = arithmetic_mean(clean)
        if (order == 1) then
            value = mean_value
            return
        end if

        centered = clean - mean_value
        m2 = sum(centered ** 2) / nr
        if (order == 2) then
            value = nr * m2 / (nr - 1.0_dp)
            return
        end if

        m3 = sum(centered ** 3) / nr
        if (order == 3) then
            value = nr ** 2 * m3 / ((nr - 1.0_dp) * (nr - 2.0_dp))
            return
        end if

        m4 = sum(centered ** 4) / nr
        value = nr ** 2 * ((nr + 1.0_dp) * m4 - &
            3.0_dp * (nr - 1.0_dp) * m2 ** 2) / &
            ((nr - 1.0_dp) * (nr - 2.0_dp) * (nr - 3.0_dp))
    end function kstat_vector

    pure function kstatvar_vector(x, order, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations in one reduction slice
        integer, intent(in) :: order !! variance order, either 1 or 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: clean(:)
        real(dp) :: k2
        real(dp) :: k4
        integer :: n
        real(dp) :: nr
        character(len=16) :: policy

        if (order < 1 .or. order > 2) then
            error stop 'scifort: kstatvar order must be 1 or 2', quiet=.true.
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

        n = size(clean)
        if (order == 1) then
            if (n < 2) then
                value = quiet_nan(0.0_dp)
                return
            end if
            nr = real(n, dp)
            k2 = kstat_vector(clean, 2)
            value = k2 / nr
        else
            if (n < 4) then
                value = quiet_nan(0.0_dp)
                return
            end if
            nr = real(n, dp)
            k2 = kstat_vector(clean, 2)
            k4 = kstat_vector(clean, 4)
            value = (2.0_dp * nr * k2 ** 2 + (nr - 1.0_dp) * k4) / &
                (nr * (nr + 1.0_dp))
        end if
    end function kstatvar_vector

    pure subroutine kstat_reduce_assumed_rank(x, axis, order, stat_kind, reduced, &
            out_shape, ok, nan_policy)
        real(dp), intent(in) :: x(..) !! observations with rank 2 to 15
        integer, intent(in) :: axis !! one-based Fortran axis to reduce
        integer, intent(in) :: order !! requested k-statistic or variance order
        integer, intent(in) :: stat_kind !! internal selector for kstat or kstatvar
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened reduction result
        integer, intent(out) :: out_shape(:) !! shape with the reduced axis removed
        logical, intent(out) :: ok !! true when rank, axis, and selector are valid
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
        if (stat_kind /= kstat_kind_stat .and. stat_kind /= kstat_kind_var) then
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
            if (stat_kind == kstat_kind_stat) then
                reduced(i + 1) = kstat_vector(slice, order, nan_policy)
            else
                reduced(i + 1) = kstatvar_vector(slice, order, nan_policy)
            end if
        end do
        ok = .true.
    end subroutine kstat_reduce_assumed_rank

    pure function kstat_all_assumed_rank(x, order, stat_kind, nan_policy) result(value)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: order !! requested k-statistic or variance order
        integer, intent(in) :: stat_kind !! internal selector for kstat or kstatvar
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        if (rank(x) < 1 .or. rank(x) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if
        call flatten_assumed_real(x, flat, dims)
        if (stat_kind == kstat_kind_stat) then
            value = kstat_vector(flat, order, nan_policy)
        else if (stat_kind == kstat_kind_var) then
            value = kstatvar_vector(flat, order, nan_policy)
        else
            value = quiet_nan(0.0_dp)
        end if
    end function kstat_all_assumed_rank

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
            allocate(flat(0))
            allocate(dims(0))
        end select
    end subroutine flatten_assumed_real

end module scifort_descriptive_kstat_engine
