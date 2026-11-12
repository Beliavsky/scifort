! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_rank_1_8
    use scifort_descriptive_rank_engine, only : rankdata_all_assumed_rank, &
        rankdata_axis_assumed_rank, tiecorrect
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface rankdata
        module procedure rankdata_all_rank1
        module procedure rankdata_axis_rank1
        module procedure rankdata_all_rank2
        module procedure rankdata_axis_rank2
        module procedure rankdata_all_rank3
        module procedure rankdata_axis_rank3
        module procedure rankdata_all_rank4
        module procedure rankdata_axis_rank4
        module procedure rankdata_all_rank5
        module procedure rankdata_axis_rank5
        module procedure rankdata_all_rank6
        module procedure rankdata_axis_rank6
        module procedure rankdata_all_rank7
        module procedure rankdata_axis_rank7
        module procedure rankdata_all_rank8
        module procedure rankdata_axis_rank8
    end interface rankdata

    public :: rankdata
    public :: tiecorrect

contains

    pure function rankdata_all_rank1(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank1

    pure function rankdata_axis_rank1(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer, allocatable :: dims_out(:)
        logical :: ok
        real(dp), allocatable :: ranked(:)

        call rankdata_axis_assumed_rank(x, axis, method, nan_policy, ranked, dims_out, ok)
        if (.not. ok) then
            allocate(values, source=x)
            values = quiet_nan(0.0_dp)
            return
        end if
        values = reshape(ranked, shape(x))
    end function rankdata_axis_rank1

    pure function rankdata_all_rank2(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank2

    pure function rankdata_axis_rank2(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        integer, allocatable :: dims_out(:)
        logical :: ok
        real(dp), allocatable :: ranked(:)

        call rankdata_axis_assumed_rank(x, axis, method, nan_policy, ranked, dims_out, ok)
        if (.not. ok) then
            allocate(values, source=x)
            values = quiet_nan(0.0_dp)
            return
        end if
        values = reshape(ranked, shape(x))
    end function rankdata_axis_rank2

    pure function rankdata_all_rank3(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank3

    pure function rankdata_axis_rank3(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        integer, allocatable :: dims_out(:)
        logical :: ok
        real(dp), allocatable :: ranked(:)

        call rankdata_axis_assumed_rank(x, axis, method, nan_policy, ranked, dims_out, ok)
        if (.not. ok) then
            allocate(values, source=x)
            values = quiet_nan(0.0_dp)
            return
        end if
        values = reshape(ranked, shape(x))
    end function rankdata_axis_rank3

    pure function rankdata_all_rank4(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank4

    pure function rankdata_axis_rank4(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        integer, allocatable :: dims_out(:)
        logical :: ok
        real(dp), allocatable :: ranked(:)

        call rankdata_axis_assumed_rank(x, axis, method, nan_policy, ranked, dims_out, ok)
        if (.not. ok) then
            allocate(values, source=x)
            values = quiet_nan(0.0_dp)
            return
        end if
        values = reshape(ranked, shape(x))
    end function rankdata_axis_rank4

    pure function rankdata_all_rank5(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank5

    pure function rankdata_axis_rank5(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        integer, allocatable :: dims_out(:)
        logical :: ok
        real(dp), allocatable :: ranked(:)

        call rankdata_axis_assumed_rank(x, axis, method, nan_policy, ranked, dims_out, ok)
        if (.not. ok) then
            allocate(values, source=x)
            values = quiet_nan(0.0_dp)
            return
        end if
        values = reshape(ranked, shape(x))
    end function rankdata_axis_rank5

    pure function rankdata_all_rank6(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank6

    pure function rankdata_axis_rank6(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer, allocatable :: dims_out(:)
        logical :: ok
        real(dp), allocatable :: ranked(:)

        call rankdata_axis_assumed_rank(x, axis, method, nan_policy, ranked, dims_out, ok)
        if (.not. ok) then
            allocate(values, source=x)
            values = quiet_nan(0.0_dp)
            return
        end if
        values = reshape(ranked, shape(x))
    end function rankdata_axis_rank6

    pure function rankdata_all_rank7(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank7

    pure function rankdata_axis_rank7(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer, allocatable :: dims_out(:)
        logical :: ok
        real(dp), allocatable :: ranked(:)

        call rankdata_axis_assumed_rank(x, axis, method, nan_policy, ranked, dims_out, ok)
        if (.not. ok) then
            allocate(values, source=x)
            values = quiet_nan(0.0_dp)
            return
        end if
        values = reshape(ranked, shape(x))
    end function rankdata_axis_rank7

    pure function rankdata_all_rank8(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank8

    pure function rankdata_axis_rank8(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer, allocatable :: dims_out(:)
        logical :: ok
        real(dp), allocatable :: ranked(:)

        call rankdata_axis_assumed_rank(x, axis, method, nan_policy, ranked, dims_out, ok)
        if (.not. ok) then
            allocate(values, source=x)
            values = quiet_nan(0.0_dp)
            return
        end if
        values = reshape(ranked, shape(x))
    end function rankdata_axis_rank8

end module scifort_descriptive_rank_1_8
