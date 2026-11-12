! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_rank_9_15
    use scifort_descriptive_rank_engine, only : rankdata_all_assumed_rank, &
        rankdata_axis_assumed_rank, tiecorrect
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface rankdata
        module procedure rankdata_all_rank9
        module procedure rankdata_axis_rank9
        module procedure rankdata_all_rank10
        module procedure rankdata_axis_rank10
        module procedure rankdata_all_rank11
        module procedure rankdata_axis_rank11
        module procedure rankdata_all_rank12
        module procedure rankdata_axis_rank12
        module procedure rankdata_all_rank13
        module procedure rankdata_axis_rank13
        module procedure rankdata_all_rank14
        module procedure rankdata_axis_rank14
        module procedure rankdata_all_rank15
        module procedure rankdata_axis_rank15
    end interface rankdata

    public :: rankdata

contains

    pure function rankdata_all_rank9(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank9

    pure function rankdata_axis_rank9(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

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
    end function rankdata_axis_rank9

    pure function rankdata_all_rank10(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank10

    pure function rankdata_axis_rank10(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

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
    end function rankdata_axis_rank10

    pure function rankdata_all_rank11(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank11

    pure function rankdata_axis_rank11(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

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
    end function rankdata_axis_rank11

    pure function rankdata_all_rank12(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank12

    pure function rankdata_axis_rank12(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

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
    end function rankdata_axis_rank12

    pure function rankdata_all_rank13(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank13

    pure function rankdata_axis_rank13(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

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
    end function rankdata_axis_rank13

    pure function rankdata_all_rank14(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank14

    pure function rankdata_axis_rank14(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

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
    end function rankdata_axis_rank14

    pure function rankdata_all_rank15(x, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations flattened before ranking
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        values = rankdata_all_assumed_rank(x, method, nan_policy)
    end function rankdata_all_rank15

    pure function rankdata_axis_rank15(x, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations ranked independently by slice
        integer, intent(in) :: axis !! one-based Fortran ranking axis
        character(len=*), intent(in), optional :: method !! average/min/max/dense/ordinal
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

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
    end function rankdata_axis_rank15

end module scifort_descriptive_rank_9_15
