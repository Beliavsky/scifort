! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_lmoment_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: lmoment_all_orders_assumed_rank
    public :: lmoment_all_scalar_assumed_rank
    public :: lmoment_reduce_orders_assumed_rank
    public :: lmoment_reduce_scalar_assumed_rank
    public :: lmoment_vector_orders
    public :: lmoment_vector_scalar

contains

    pure function lmoment_vector_scalar(x, order, sorted, standardize, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations in one reduction slice
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when x is already sorted ascending
        logical, intent(in), optional :: standardize !! divide orders 3+ by the second L-moment
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: values(:)

        call lmoment_vector_orders(x, [order], sorted, standardize, nan_policy, values)
        value = values(1)
    end function lmoment_vector_scalar

    pure subroutine lmoment_vector_orders(x, orders, sorted, standardize, nan_policy, values)
        real(dp), intent(in) :: x(:) !! observations in one reduction slice
        integer, intent(in) :: orders(:) !! positive requested L-moment orders
        logical, intent(in), optional :: sorted !! true when x is already sorted ascending
        logical, intent(in), optional :: standardize !! divide orders 3+ by the second L-moment
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable, intent(out) :: values(:) !! requested L-moments

        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: clean(:)
        real(dp) :: coefficient
        real(dp) :: next_coefficient
        integer :: i
        integer :: k
        integer :: max_order
        integer :: n
        integer :: r
        real(dp), allocatable :: raw(:)
        logical :: sorted_value
        logical :: standardize_value
        character(len=16) :: policy

        if (size(orders) == 0 .or. any(orders <= 0)) then
            error stop 'scifort: lmoment order must be a nonempty set of positive integers', quiet=.true.
        end if

        allocate(values(size(orders)))
        values = quiet_nan(0.0_dp)

        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))
        if (.not. valid_nan_policy(policy)) return

        if (policy == 'raise') then
            if (any(ieee_is_nan(x))) then
                error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
            end if
            clean = x
        else if (policy == 'propagate') then
            if (any(ieee_is_nan(x))) return
            clean = x
        else
            clean = pack(x, .not. ieee_is_nan(x))
        end if

        n = size(clean)
        if (n == 0) return

        sorted_value = .false.
        if (present(sorted)) sorted_value = sorted
        if (.not. sorted_value) call sort_values(clean)

        standardize_value = .true.
        if (present(standardize)) standardize_value = standardize

        max_order = min(maxval(orders), n)
        allocate(beta(max_order))
        allocate(raw(max_order))

        do r = 0, max_order - 1
            beta(r + 1) = probability_weighted_moment(clean, r)
        end do

        do r = 0, max_order - 1
            coefficient = merge(1.0_dp, -1.0_dp, mod(r, 2) == 0)
            raw(r + 1) = 0.0_dp
            do k = 0, r
                raw(r + 1) = safe_add(raw(r + 1), safe_multiply(coefficient, beta(k + 1)))
                if (k < r) then
                    next_coefficient = safe_multiply(coefficient, &
                        -real(r - k, dp) * real(r + k + 1, dp) / real(k + 1, dp) ** 2)
                    coefficient = next_coefficient
                end if
            end do
        end do

        do i = 1, size(orders)
            r = orders(i)
            if (r > n) cycle
            if (standardize_value .and. r >= 3) then
                values(i) = safe_divide(raw(r), raw(2))
            else
                values(i) = raw(r)
            end if
        end do
    end subroutine lmoment_vector_orders

    pure subroutine lmoment_reduce_orders_assumed_rank(x, axis, orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        integer, intent(in) :: orders(:) !! positive requested L-moment orders
        logical, intent(in), optional :: sorted !! true when slices are already sorted
        logical, intent(in), optional :: standardize !! standardize orders 3+
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable, intent(out) :: reduced(:) !! order-fast flattened results
        integer, intent(out) :: out_shape(:) !! shape with the reduced axis removed
        logical, intent(out) :: ok !! true when rank and axis are valid

        integer :: base_index
        integer :: coordinate
        integer :: d
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        integer :: i
        integer :: linear
        integer :: n_orders
        integer :: n_out
        integer :: n_reduce
        integer :: out_dim
        real(dp), allocatable :: slice(:)
        real(dp), allocatable :: slice_values(:)
        integer, allocatable :: strides(:)

        ok = .false.
        if (rank(x) < 1 .or. rank(x) > 15) then
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
        if (size(orders) == 0 .or. any(orders <= 0)) then
            error stop 'scifort: lmoment order must be a nonempty set of positive integers', quiet=.true.
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
        n_orders = size(orders)
        allocate(reduced(n_orders * n_out))
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
            call lmoment_vector_orders(slice, orders, sorted, standardize, nan_policy, slice_values)
            reduced(i * n_orders + 1:(i + 1) * n_orders) = slice_values
        end do
        ok = .true.
    end subroutine lmoment_reduce_orders_assumed_rank

    pure subroutine lmoment_reduce_scalar_assumed_rank(x, axis, order, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        integer, intent(in) :: order !! positive requested L-moment order
        logical, intent(in), optional :: sorted !! true when slices are already sorted
        logical, intent(in), optional :: standardize !! standardize orders 3+
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened result
        integer, intent(out) :: out_shape(:) !! shape with reduced axis removed
        logical, intent(out) :: ok !! true when rank and axis are valid

        call lmoment_reduce_orders_assumed_rank(x, axis, [order], sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
    end subroutine lmoment_reduce_scalar_assumed_rank

    pure subroutine lmoment_all_orders_assumed_rank(x, orders, sorted, standardize, &
            nan_policy, values)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15, flattened first
        integer, intent(in) :: orders(:) !! positive requested L-moment orders
        logical, intent(in), optional :: sorted !! true when flattened input is already sorted
        logical, intent(in), optional :: standardize !! standardize orders 3+
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable, intent(out) :: values(:) !! requested L-moments

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        if (rank(x) < 1 .or. rank(x) > 15) then
            allocate(values(size(orders)))
            values = quiet_nan(0.0_dp)
            return
        end if
        call flatten_assumed_real(x, flat, dims)
        call lmoment_vector_orders(flat, orders, sorted, standardize, nan_policy, values)
    end subroutine lmoment_all_orders_assumed_rank

    pure function lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, &
            nan_policy) result(value)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15, flattened first
        integer, intent(in) :: order !! positive requested L-moment order
        logical, intent(in), optional :: sorted !! true when flattened input is already sorted
        logical, intent(in), optional :: standardize !! standardize orders 3+
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, [order], sorted, standardize, nan_policy, values)
        value = values(1)
    end function lmoment_all_scalar_assumed_rank

    pure function probability_weighted_moment(x, r) result(value)
        real(dp), intent(in) :: x(:) !! sorted non-NaN observations
        integer, intent(in) :: r !! zero-based probability-weighted-moment order
        real(dp) :: value

        integer :: j
        integer :: n
        logical :: negative_inf_seen
        logical :: positive_inf_seen
        real(dp) :: sum_value
        real(dp) :: term
        real(dp) :: weight

        n = size(x)
        if (r < 0 .or. r >= n) then
            value = 0.0_dp
            return
        end if

        sum_value = 0.0_dp
        positive_inf_seen = .false.
        negative_inf_seen = .false.
        weight = 1.0_dp
        do j = n - 1, r, -1
            if (.not. ieee_is_finite(x(j + 1))) then
                if (x(j + 1) > 0.0_dp) then
                    positive_inf_seen = .true.
                else
                    negative_inf_seen = .true.
                end if
            else
                term = x(j + 1) * weight
                sum_value = safe_add(sum_value, term)
            end if
            if (j > r) weight = weight * real(j - r, dp) / real(j, dp)
        end do

        if (positive_inf_seen .and. negative_inf_seen) then
            value = quiet_nan(0.0_dp)
        else if (positive_inf_seen) then
            value = positive_infinity(0.0_dp)
        else if (negative_inf_seen) then
            value = negative_infinity(0.0_dp)
        else
            value = sum_value / real(n, dp)
        end if
    end function probability_weighted_moment

    pure function safe_add(a, b) result(value)
        real(dp), intent(in) :: a !! first addend
        real(dp), intent(in) :: b !! second addend
        real(dp) :: value

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            value = quiet_nan(0.0_dp)
        else if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            if (.not. ieee_is_finite(a) .and. .not. ieee_is_finite(b)) then
                if ((a > 0.0_dp) .eqv. (b > 0.0_dp)) then
                    value = a
                else
                    value = quiet_nan(0.0_dp)
                end if
            else if (.not. ieee_is_finite(a)) then
                value = a
            else
                value = b
            end if
        else if (a > 0.0_dp .and. b > 0.0_dp .and. a > huge(a) - b) then
            value = positive_infinity(0.0_dp)
        else if (a < 0.0_dp .and. b < 0.0_dp .and. a < -huge(a) - b) then
            value = negative_infinity(0.0_dp)
        else
            value = a + b
        end if
    end function safe_add

    pure function safe_multiply(a, b) result(value)
        real(dp), intent(in) :: a !! first factor
        real(dp), intent(in) :: b !! second factor
        real(dp) :: value

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            value = quiet_nan(0.0_dp)
        else if ((same_real(a, 0.0_dp) .and. .not. ieee_is_finite(b)) .or. &
                (same_real(b, 0.0_dp) .and. .not. ieee_is_finite(a))) then
            value = quiet_nan(0.0_dp)
        else if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            value = sign(positive_infinity(0.0_dp), a) * sign(1.0_dp, b)
        else if (same_real(a, 0.0_dp) .or. same_real(b, 0.0_dp)) then
            value = a * b
        else if (multiply_overflows(a, b)) then
            value = sign(positive_infinity(0.0_dp), a) * sign(1.0_dp, b)
        else
            value = a * b
        end if
    end function safe_multiply

    pure logical function multiply_overflows(a, b) result(overflows)
        real(dp), intent(in) :: a !! finite nonzero first factor
        real(dp), intent(in) :: b !! finite nonzero second factor

        if (abs(a) > 1.0_dp) then
            overflows = abs(b) > huge(a) / abs(a)
        else if (abs(b) > 1.0_dp) then
            overflows = abs(a) > huge(a) / abs(b)
        else
            overflows = .false.
        end if
    end function multiply_overflows

    pure function safe_divide(numerator, denominator) result(value)
        real(dp), intent(in) :: numerator !! numerator
        real(dp), intent(in) :: denominator !! denominator
        real(dp) :: value

        if (ieee_is_nan(numerator) .or. ieee_is_nan(denominator)) then
            value = quiet_nan(0.0_dp)
        else if (same_real(denominator, 0.0_dp)) then
            if (same_real(numerator, 0.0_dp)) then
                value = quiet_nan(0.0_dp)
            else
                value = sign(positive_infinity(0.0_dp), numerator) * &
                    sign(1.0_dp, denominator)
            end if
        else if (.not. ieee_is_finite(numerator) .and. .not. ieee_is_finite(denominator)) then
            value = quiet_nan(0.0_dp)
        else if (.not. ieee_is_finite(numerator)) then
            value = sign(positive_infinity(0.0_dp), numerator) * sign(1.0_dp, denominator)
        else if (.not. ieee_is_finite(denominator)) then
            value = sign(0.0_dp, numerator) * sign(1.0_dp, denominator)
        else
            value = numerator / denominator
        end if
    end function safe_divide

    pure logical function same_real(a, b) result(same)
        real(dp), intent(in) :: a !! first real value
        real(dp), intent(in) :: b !! second real value

        same = (a <= b) .and. (b <= a)
    end function same_real

    pure logical function valid_nan_policy(nan_policy) result(valid)
        character(len=*), intent(in) :: nan_policy !! requested NaN handling policy

        select case (trim(adjustl(nan_policy)))
        case ('propagate', 'omit', 'raise')
            valid = .true.
        case default
            valid = .false.
        end select
    end function valid_nan_policy

    pure subroutine sort_values(values)
        real(dp), intent(inout) :: values(:) !! values sorted in ascending order

        real(dp), allocatable :: temporary(:)

        if (size(values) <= 1) return
        allocate(temporary(size(values)))
        call merge_sort(values, temporary, 1, size(values))
    end subroutine sort_values

    recursive pure subroutine merge_sort(values, temporary, left, right)
        real(dp), intent(inout) :: values(:) !! values being sorted
        real(dp), intent(inout) :: temporary(:) !! merge workspace
        integer, intent(in) :: left !! first inclusive index
        integer, intent(in) :: right !! last inclusive index

        integer :: i
        integer :: j
        integer :: k
        integer :: middle

        if (left >= right) return
        middle = left + (right - left) / 2
        call merge_sort(values, temporary, left, middle)
        call merge_sort(values, temporary, middle + 1, right)
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
    end subroutine merge_sort

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

end module scifort_descriptive_lmoment_engine
