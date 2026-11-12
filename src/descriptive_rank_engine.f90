! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_rank_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    public :: rankdata_all_assumed_rank
    public :: rankdata_axis_assumed_rank
    public :: rankdata_vector_core
    public :: tiecorrect

contains

    pure function rankdata_vector_core(x, method, nan_policy) result(ranks)
        real(dp), intent(in) :: x(:) !! observations to rank
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: ranks(:)

        real(dp), allocatable :: clean(:)
        real(dp), allocatable :: clean_ranks(:)
        integer :: i
        integer :: k
        integer, allocatable :: locations(:)
        integer :: nclean
        character(len=9) :: policy

        policy = 'propagate'
        if (present(nan_policy)) policy = trim(adjustl(nan_policy))
        if (.not. valid_nan_policy(policy)) then
            error stop 'scifort: invalid nan_policy', quiet=.true.
        end if

        allocate(ranks(size(x)))
        if (size(x) == 0) return

        if (policy == 'raise') then
            if (any(ieee_is_nan(x))) then
                error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
            end if
            ranks = rank_no_nan(x, method)
            return
        end if
        if (policy == 'propagate') then
            if (any(ieee_is_nan(x))) then
                ranks = quiet_nan(0.0_dp)
                return
            end if
            ranks = rank_no_nan(x, method)
            return
        end if

        ranks = quiet_nan(0.0_dp)
        nclean = count(.not. ieee_is_nan(x))
        if (nclean == 0) return
        allocate(clean(nclean), locations(nclean))
        k = 0
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) cycle
            k = k + 1
            clean(k) = x(i)
            locations(k) = i
        end do
        clean_ranks = rank_no_nan(clean, method)
        do i = 1, nclean
            ranks(locations(i)) = clean_ranks(i)
        end do
    end function rankdata_vector_core

    pure function rankdata_all_assumed_rank(x, method, nan_policy) result(ranks)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15, flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: ranks(:)

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(ranks(0))
            return
        end if
        ranks = rankdata_vector_core(flat, method, nan_policy)
    end function rankdata_all_assumed_rank

    pure subroutine rankdata_axis_assumed_rank(x, axis, method, nan_policy, ranked, dims, ok)
        real(dp), intent(in) :: x(..) !! observations with rank 1 to 15
        integer, intent(in) :: axis !! one-based Fortran axis to rank independently
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable, intent(out) :: ranked(:) !! flattened ranked result
        integer, allocatable, intent(out) :: dims(:) !! input extents
        logical, intent(out) :: ok !! true when rank and axis are valid

        integer :: base_index
        integer :: coordinate
        integer :: d
        real(dp), allocatable :: flat(:)
        integer :: i
        integer :: linear
        integer :: n_slices
        integer :: n_reduce
        real(dp), allocatable :: slice(:)
        real(dp), allocatable :: slice_ranks(:)
        integer, allocatable :: strides(:)

        ok = .false.
        call flatten_assumed_real(x, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(ranked(0))
            return
        end if
        if (axis < 1 .or. axis > size(dims)) then
            allocate(ranked(0))
            return
        end if

        allocate(ranked(size(flat)))
        if (size(flat) == 0) then
            ok = .true.
            return
        end if
        allocate(strides(size(dims)))
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
        end do

        n_reduce = dims(axis)
        if (n_reduce == 0) then
            ok = .true.
            return
        end if
        n_slices = size(flat) / n_reduce
        allocate(slice(n_reduce))
        do i = 0, n_slices - 1
            linear = i
            base_index = 1
            do d = 1, size(dims)
                if (d == axis) cycle
                coordinate = mod(linear, dims(d))
                linear = linear / dims(d)
                base_index = base_index + coordinate * strides(d)
            end do
            do d = 1, n_reduce
                slice(d) = flat(base_index + (d - 1) * strides(axis))
            end do
            slice_ranks = rankdata_vector_core(slice, method, nan_policy)
            do d = 1, n_reduce
                ranked(base_index + (d - 1) * strides(axis)) = slice_ranks(d)
            end do
        end do
        ok = .true.
    end subroutine rankdata_axis_assumed_rank

    pure function tiecorrect(rankvals) result(factor)
        real(dp), intent(in) :: rankvals(:) !! one-dimensional rank values
        real(dp) :: factor

        integer :: group_end
        integer :: group_start
        integer :: n
        real(dp), allocatable :: sorted(:)
        real(dp) :: tie_sum
        integer :: tie_size

        n = size(rankvals)
        if (n < 2) then
            factor = 1.0_dp
            return
        end if
        sorted = rankvals
        call sort_nan_last(sorted)
        tie_sum = 0.0_dp
        group_start = 1
        do while (group_start <= n)
            if (ieee_is_nan(sorted(group_start))) exit
            group_end = group_start
            do while (group_end < n)
                if (ieee_is_nan(sorted(group_end + 1))) exit
                if (.not. same_rank_value(sorted(group_end + 1), sorted(group_start))) exit
                group_end = group_end + 1
            end do
            tie_size = group_end - group_start + 1
            if (tie_size > 1) then
                tie_sum = tie_sum + real(tie_size, dp)**3 - real(tie_size, dp)
            end if
            group_start = group_end + 1
        end do
        factor = 1.0_dp - tie_sum / (real(n, dp)**3 - real(n, dp))
    end function tiecorrect

    pure function rank_no_nan(x, method) result(ranks)
        real(dp), intent(in) :: x(:) !! NaN-free observations
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        real(dp), allocatable :: ranks(:)

        integer :: dense_rank
        integer :: group_end
        integer :: group_start
        integer :: i
        integer, allocatable :: order(:)
        character(len=7) :: selected_method

        allocate(ranks(size(x)))
        if (size(x) == 0) return
        selected_method = 'average'
        if (present(method)) selected_method = trim(adjustl(method))
        if (selected_method /= 'average' .and. selected_method /= 'min' .and. &
                selected_method /= 'max' .and. selected_method /= 'dense' .and. &
                selected_method /= 'ordinal') then
            ranks = quiet_nan(0.0_dp)
            return
        end if

        allocate(order(size(x)))
        do i = 1, size(x)
            order(i) = i
        end do
        call stable_sort_indices(x, order)
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
                if (.not. same_rank_value(x(order(group_end + 1)), x(order(group_start)))) exit
                group_end = group_end + 1
            end do
            select case (selected_method)
            case ('average')
                ranks(order(group_start:group_end)) = &
                    0.5_dp * real(group_start + group_end, dp)
            case ('min')
                ranks(order(group_start:group_end)) = real(group_start, dp)
            case ('max')
                ranks(order(group_start:group_end)) = real(group_end, dp)
            case ('dense')
                ranks(order(group_start:group_end)) = real(dense_rank, dp)
            end select
            dense_rank = dense_rank + 1
            group_start = group_end + 1
        end do
    end function rank_no_nan

    pure subroutine stable_sort_indices(x, order)
        real(dp), intent(in) :: x(:) !! values used as sort keys
        integer, intent(inout) :: order(:) !! permutation sorted by ascending x

        integer :: i
        integer :: j
        integer :: key

        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (x(order(j)) <= x(key)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
    end subroutine stable_sort_indices

    pure subroutine sort_nan_last(x)
        real(dp), intent(inout) :: x(:) !! values sorted ascending with NaNs last

        integer :: i
        integer :: j
        real(dp) :: key

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (ieee_is_nan(key)) exit
                if (ieee_is_nan(x(j))) then
                    x(j + 1) = x(j)
                    j = j - 1
                    cycle
                end if
                if (x(j) <= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_nan_last

    pure logical function same_rank_value(a, b)
        real(dp), intent(in) :: a !! first non-NaN rank value
        real(dp), intent(in) :: b !! second non-NaN rank value

        same_rank_value = a <= b .and. b <= a
    end function same_rank_value

    pure logical function valid_nan_policy(policy)
        character(len=*), intent(in) :: policy !! requested NaN policy

        valid_nan_policy = policy == 'propagate' .or. policy == 'omit' .or. policy == 'raise'
    end function valid_nan_policy

    pure subroutine flatten_assumed_real(x, flat, dims)
        real(dp), intent(in) :: x(..) !! rank-1 through rank-15 real array
        real(dp), allocatable, intent(out) :: flat(:) !! Fortran-order flattened values
        integer, allocatable, intent(out) :: dims(:) !! input shape

        select rank (x)
        rank (1)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (2)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (3)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (4)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (5)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (6)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (7)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (8)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (9)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (10)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (11)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (12)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (13)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (14)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank (15)
            dims = shape(x); flat = reshape(x, [size(x)])
        rank default
            allocate(dims(0), flat(0))
        end select
    end subroutine flatten_assumed_real

end module scifort_descriptive_rank_engine
