! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_proportion_trim_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: trim_mode_both = 1
    integer, parameter, public :: trim_mode_left = 2
    integer, parameter, public :: trim_mode_right = 3

    public :: proportion_trim_all_assumed_rank
    public :: proportion_trim_assumed_rank
    public :: trim_mean_all_assumed_rank
    public :: trim_mean_reduce_assumed_rank
    public :: trim_mean_vector
    public :: trim_mode_from_tail

contains

    pure function trim_mean_vector(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations in one reduction slice
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: clean(:)
        integer :: lowercut
        integer :: n
        character(len=16) :: policy
        real(dp), allocatable :: sorted(:)
        integer :: uppercut

        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))
        if (.not. valid_nan_policy(policy)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (.not. ieee_is_finite(proportiontocut) .or. proportiontocut < 0.0_dp) then
            error stop 'scifort: proportiontocut must be finite and nonnegative', quiet=.true.
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
        if (n == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        lowercut = int(proportiontocut * real(n, dp))
        uppercut = n - lowercut
        if (lowercut > uppercut) then
            error stop 'scifort: proportiontocut is too large', quiet=.true.
        end if
        if (lowercut == uppercut) then
            value = quiet_nan(0.0_dp)
            return
        end if

        sorted = clean
        call sort_values_nan_last(sorted)
        value = arithmetic_mean(sorted(lowercut + 1:uppercut))
    end function trim_mean_vector

    pure subroutine trim_mean_reduce_assumed_rank(x, axis, proportiontocut, reduced, &
            out_shape, ok, nan_policy)
        real(dp), intent(in) :: x(..) !! observations with rank 2 to 15
        integer, intent(in) :: axis !! one-based Fortran axis to reduce
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened reduction result
        integer, intent(out) :: out_shape(:) !! shape with reduced axis removed
        logical, intent(out) :: ok !! true when rank and axis are valid
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
            reduced(i + 1) = trim_mean_vector(slice, proportiontocut, nan_policy)
        end do
        ok = .true.
    end subroutine trim_mean_reduce_assumed_rank

    pure function trim_mean_all_assumed_rank(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        if (rank(x) < 1 .or. rank(x) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if
        call flatten_assumed_real(x, flat, dims)
        value = trim_mean_vector(flat, proportiontocut, nan_policy)
    end function trim_mean_all_assumed_rank

    pure subroutine proportion_trim_assumed_rank(x, axis, proportiontocut, trim_mode, &
            trimmed, out_shape, ok)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: axis !! one-based Fortran axis to trim
        real(dp), intent(in) :: proportiontocut !! fraction removed from selected tail(s)
        integer, intent(in) :: trim_mode !! trim_mode_both/left/right
        real(dp), allocatable, intent(out) :: trimmed(:) !! flattened trimmed result
        integer, intent(out) :: out_shape(:) !! result extents
        logical, intent(out) :: ok !! true when inputs are valid

        integer :: base_in
        integer :: base_out
        integer :: coordinate
        integer :: cut
        integer :: d
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        integer :: i
        integer :: j
        integer :: linear
        integer :: lowercut
        integer :: n_keep
        integer :: n_slices
        integer :: nobs
        integer :: out_dim
        integer, allocatable :: out_strides(:)
        real(dp), allocatable :: slice(:)
        integer, allocatable :: strides(:)
        integer :: uppercut

        ok = .false.
        if (rank(x) < 1 .or. rank(x) > 15) then
            allocate(trimmed(0))
            out_shape = 0
            return
        end if
        if (size(out_shape) /= rank(x)) then
            allocate(trimmed(0))
            out_shape = 0
            return
        end if
        if (axis < 1 .or. axis > rank(x)) then
            allocate(trimmed(0))
            out_shape = 0
            return
        end if
        if (.not. ieee_is_finite(proportiontocut) .or. proportiontocut < 0.0_dp) then
            error stop 'scifort: proportiontocut must be finite and nonnegative', quiet=.true.
        end if
        if (trim_mode < trim_mode_both .or. trim_mode > trim_mode_right) then
            allocate(trimmed(0))
            out_shape = 0
            return
        end if

        call flatten_assumed_real(x, flat, dims)
        if (size(flat) == 0) then
            out_shape = dims
            allocate(trimmed(0))
            ok = .true.
            return
        end if
        nobs = dims(axis)
        cut = int(proportiontocut * real(nobs, dp))
        select case (trim_mode)
        case (trim_mode_both)
            lowercut = cut
            uppercut = nobs - cut
            if (lowercut >= uppercut) then
                error stop 'scifort: proportiontocut is too large', quiet=.true.
            end if
        case (trim_mode_left)
            if (proportiontocut >= 1.0_dp) then
                lowercut = nobs
            else
                lowercut = cut
            end if
            uppercut = nobs
        case (trim_mode_right)
            lowercut = 0
            if (proportiontocut >= 1.0_dp) then
                uppercut = 0
            else
                uppercut = nobs - cut
            end if
        end select
        n_keep = max(0, uppercut - lowercut)
        out_shape = dims
        out_shape(axis) = n_keep
        allocate(trimmed(product(out_shape)))
        if (size(trimmed) == 0) then
            ok = .true.
            return
        end if

        allocate(strides(size(dims)))
        allocate(out_strides(size(dims)))
        strides(1) = 1
        out_strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
            out_strides(d) = out_strides(d - 1) * out_shape(d - 1)
        end do
        allocate(slice(nobs))
        if (nobs == 0) then
            ok = .true.
            return
        end if

        n_slices = 1
        do d = 1, size(dims)
            if (d /= axis) n_slices = n_slices * dims(d)
        end do
        do i = 0, n_slices - 1
            linear = i
            base_in = 1
            base_out = 1
            out_dim = 0
            do d = 1, size(dims)
                if (d == axis) cycle
                out_dim = out_dim + 1
                coordinate = mod(linear, dims(d))
                linear = linear / dims(d)
                base_in = base_in + coordinate * strides(d)
                base_out = base_out + coordinate * out_strides(d)
            end do
            do j = 1, nobs
                slice(j) = flat(base_in + (j - 1) * strides(axis))
            end do
            call sort_values_nan_last(slice)
            do j = 1, n_keep
                trimmed(base_out + (j - 1) * out_strides(axis)) = slice(lowercut + j)
            end do
        end do
        ok = .true.
    end subroutine proportion_trim_assumed_rank

    pure function proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode) result(trimmed)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        real(dp), intent(in) :: proportiontocut !! fraction removed from selected tail(s)
        integer, intent(in) :: trim_mode !! trim_mode_both/left/right
        real(dp), allocatable :: trimmed(:)

        integer :: cut
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        integer :: lowercut
        integer :: nobs
        integer :: uppercut

        if (rank(x) < 1 .or. rank(x) > 15) then
            allocate(trimmed(0))
            return
        end if
        if (.not. ieee_is_finite(proportiontocut) .or. proportiontocut < 0.0_dp) then
            error stop 'scifort: proportiontocut must be finite and nonnegative', quiet=.true.
        end if
        call flatten_assumed_real(x, flat, dims)
        nobs = size(flat)
        if (nobs == 0) then
            trimmed = flat
            return
        end if
        cut = int(proportiontocut * real(nobs, dp))
        select case (trim_mode)
        case (trim_mode_both)
            lowercut = cut
            uppercut = nobs - cut
            if (lowercut >= uppercut) then
                error stop 'scifort: proportiontocut is too large', quiet=.true.
            end if
        case (trim_mode_left)
            lowercut = merge(nobs, cut, proportiontocut >= 1.0_dp)
            uppercut = nobs
        case (trim_mode_right)
            lowercut = 0
            uppercut = merge(0, nobs - cut, proportiontocut >= 1.0_dp)
        case default
            allocate(trimmed(0))
            return
        end select
        call sort_values_nan_last(flat)
        trimmed = flat(lowercut + 1:uppercut)
    end function proportion_trim_all_assumed_rank

    pure function trim_mode_from_tail(tail) result(trim_mode)
        character(len=*), intent(in) :: tail !! left or right tail selector
        integer :: trim_mode

        character(len=len(tail)) :: normalized
        integer :: i
        integer :: code

        normalized = trim(adjustl(tail))
        do i = 1, len_trim(normalized)
            code = iachar(normalized(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                normalized(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
        select case (trim(normalized))
        case ('left')
            trim_mode = trim_mode_left
        case ('right')
            trim_mode = trim_mode_right
        case default
            error stop "scifort: tail must be 'left' or 'right'", quiet=.true.
        end select
    end function trim_mode_from_tail

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

end module scifort_descriptive_proportion_trim_engine
