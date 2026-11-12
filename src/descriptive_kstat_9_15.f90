! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_kstat_9_15
    use scifort_descriptive_kstat_engine, only : kstat_all_assumed_rank, &
        kstat_kind_stat, kstat_kind_var, kstat_reduce_assumed_rank, &
        kstat_vector, kstatvar_vector
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface kstat
        module procedure kstat_all_rank9
        module procedure kstat_axis_rank9
        module procedure kstat_all_rank10
        module procedure kstat_axis_rank10
        module procedure kstat_all_rank11
        module procedure kstat_axis_rank11
        module procedure kstat_all_rank12
        module procedure kstat_axis_rank12
        module procedure kstat_all_rank13
        module procedure kstat_axis_rank13
        module procedure kstat_all_rank14
        module procedure kstat_axis_rank14
        module procedure kstat_all_rank15
        module procedure kstat_axis_rank15
    end interface kstat

    interface kstat_all
        module procedure kstat_all_rank9
        module procedure kstat_all_rank10
        module procedure kstat_all_rank11
        module procedure kstat_all_rank12
        module procedure kstat_all_rank13
        module procedure kstat_all_rank14
        module procedure kstat_all_rank15
    end interface kstat_all

    interface kstat_keepdims
        module procedure kstat_keepdims_rank9
        module procedure kstat_keepdims_rank10
        module procedure kstat_keepdims_rank11
        module procedure kstat_keepdims_rank12
        module procedure kstat_keepdims_rank13
        module procedure kstat_keepdims_rank14
        module procedure kstat_keepdims_rank15
    end interface kstat_keepdims

    interface kstatvar
        module procedure kstatvar_all_rank9
        module procedure kstatvar_axis_rank9
        module procedure kstatvar_all_rank10
        module procedure kstatvar_axis_rank10
        module procedure kstatvar_all_rank11
        module procedure kstatvar_axis_rank11
        module procedure kstatvar_all_rank12
        module procedure kstatvar_axis_rank12
        module procedure kstatvar_all_rank13
        module procedure kstatvar_axis_rank13
        module procedure kstatvar_all_rank14
        module procedure kstatvar_axis_rank14
        module procedure kstatvar_all_rank15
        module procedure kstatvar_axis_rank15
    end interface kstatvar

    interface kstatvar_all
        module procedure kstatvar_all_rank9
        module procedure kstatvar_all_rank10
        module procedure kstatvar_all_rank11
        module procedure kstatvar_all_rank12
        module procedure kstatvar_all_rank13
        module procedure kstatvar_all_rank14
        module procedure kstatvar_all_rank15
    end interface kstatvar_all

    interface kstatvar_keepdims
        module procedure kstatvar_keepdims_rank9
        module procedure kstatvar_keepdims_rank10
        module procedure kstatvar_keepdims_rank11
        module procedure kstatvar_keepdims_rank12
        module procedure kstatvar_keepdims_rank13
        module procedure kstatvar_keepdims_rank14
        module procedure kstatvar_keepdims_rank15
    end interface kstatvar_keepdims

    public :: kstat
    public :: kstat_all
    public :: kstat_keepdims
    public :: kstatvar
    public :: kstatvar_all
    public :: kstatvar_keepdims

contains

    pure function kstat_all_rank9(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank9

    pure function kstat_axis_rank9(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank9

    pure function kstat_keepdims_rank9(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: order_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank9

    pure function kstatvar_all_rank9(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank9

    pure function kstatvar_axis_rank9(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank9

    pure function kstatvar_keepdims_rank9(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: order_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank9

    pure function kstat_all_rank10(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank10

    pure function kstat_axis_rank10(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank10

    pure function kstat_keepdims_rank10(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: order_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank10

    pure function kstatvar_all_rank10(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank10

    pure function kstatvar_axis_rank10(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank10

    pure function kstatvar_keepdims_rank10(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: order_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank10

    pure function kstat_all_rank11(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank11

    pure function kstat_axis_rank11(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank11

    pure function kstat_keepdims_rank11(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: order_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank11

    pure function kstatvar_all_rank11(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank11

    pure function kstatvar_axis_rank11(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank11

    pure function kstatvar_keepdims_rank11(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: order_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank11

    pure function kstat_all_rank12(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank12

    pure function kstat_axis_rank12(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank12

    pure function kstat_keepdims_rank12(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: order_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank12

    pure function kstatvar_all_rank12(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank12

    pure function kstatvar_axis_rank12(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank12

    pure function kstatvar_keepdims_rank12(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: order_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank12

    pure function kstat_all_rank13(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank13

    pure function kstat_axis_rank13(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank13

    pure function kstat_keepdims_rank13(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: order_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank13

    pure function kstatvar_all_rank13(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank13

    pure function kstatvar_axis_rank13(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank13

    pure function kstatvar_keepdims_rank13(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: order_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank13

    pure function kstat_all_rank14(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank14

    pure function kstat_axis_rank14(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank14

    pure function kstat_keepdims_rank14(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: order_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank14

    pure function kstatvar_all_rank14(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank14

    pure function kstatvar_axis_rank14(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank14

    pure function kstatvar_keepdims_rank14(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: order_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank14

    pure function kstat_all_rank15(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank15

    pure function kstat_axis_rank15(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank15

    pure function kstat_keepdims_rank15(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: order_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank15

    pure function kstatvar_all_rank15(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank15

    pure function kstatvar_axis_rank15(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank15

    pure function kstatvar_keepdims_rank15(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: order_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank15

end module scifort_descriptive_kstat_9_15
